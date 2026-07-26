# ============================================================
# 9. VALIDATION AGENT: GRAPH MODEL VALIDATION AND BENCHMARKING
# ============================================================

safe_divide <- function(num, den) {
  ifelse(is.finite(den) & den != 0, num / den, NA_real_)
}

confusion_metrics <- function(y, probability, threshold = 0.50) {
  ok <- is.finite(y) & is.finite(probability)
  y <- as.integer(y[ok])
  probability <- as.numeric(probability[ok])

  if (length(y) == 0) {
    return(tibble(
      n = 0L, TP = 0L, TN = 0L, FP = 0L, FN = 0L,
      Accuracy = NA_real_, Sensitivity = NA_real_, Specificity = NA_real_,
      Precision = NA_real_, F1 = NA_real_, MCC = NA_real_, ROC_AUC = NA_real_
    ))
  }

  prediction <- as.integer(probability >= threshold)
  tp <- sum(prediction == 1L & y == 1L)
  tn <- sum(prediction == 0L & y == 0L)
  fp <- sum(prediction == 1L & y == 0L)
  fn <- sum(prediction == 0L & y == 1L)

  sensitivity <- safe_divide(tp, tp + fn)
  specificity <- safe_divide(tn, tn + fp)
  precision <- safe_divide(tp, tp + fp)
  f1 <- ifelse(
    is.finite(precision) && is.finite(sensitivity) && (precision + sensitivity) > 0,
    2 * precision * sensitivity / (precision + sensitivity),
    NA_real_
  )
  mcc_den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))

  tibble(
    n = length(y),
    TP = tp,
    TN = tn,
    FP = fp,
    FN = fn,
    Accuracy = safe_divide(tp + tn, length(y)),
    Sensitivity = sensitivity,
    Specificity = specificity,
    Precision = precision,
    F1 = f1,
    MCC = safe_divide(tp * tn - fp * fn, mcc_den),
    ROC_AUC = binary_auc(y, probability)
  )
}

make_stratified_folds <- function(y, requested_folds = 5L, seed = 42L) {
  y <- as.integer(y)
  class_counts <- table(y)
  if (length(class_counts) < 2 || min(class_counts) < 2) return(NULL)

  k <- max(2L, min(as.integer(requested_folds), as.integer(min(class_counts))))
  set.seed(seed)
  fold_id <- integer(length(y))

  for (class_value in sort(unique(y))) {
    idx <- which(y == class_value)
    idx <- sample(idx, length(idx), replace = FALSE)
    fold_id[idx] <- rep(seq_len(k), length.out = length(idx))
  }

  list(fold_id = fold_id, k = k)
}

validation_feature_signal <- function(X) {
  required <- c(
    "MeanSynergy", "Delta_vs_Reference", "MeanObserved",
    "CombinationObservedGain_vs_Monotherapy", "MeanMonotherapyResponse",
    "MonotherapyAUC", "MaxCRISPR", "MeanExpressionPercentile",
    "n_targets", "AnyTargetMutation"
  )
  missing <- setdiff(required, colnames(X))
  if (length(missing) > 0) stop("Validation feature matrix is missing: ", paste(missing, collapse = ", "))

  0.25 * safe_graph_rescale(X[, "MeanSynergy"]) +
    0.25 * safe_graph_rescale(X[, "Delta_vs_Reference"]) +
    0.10 * safe_graph_rescale(X[, "MeanObserved"]) +
    0.10 * safe_graph_rescale(X[, "CombinationObservedGain_vs_Monotherapy"]) +
    0.10 * safe_graph_rescale(X[, "MeanMonotherapyResponse"]) +
    0.05 * safe_graph_rescale(X[, "MonotherapyAUC"]) +
    0.05 * safe_graph_rescale(X[, "MaxCRISPR"]) +
    0.04 * safe_graph_rescale(X[, "MeanExpressionPercentile"]) +
    0.04 * safe_graph_rescale(X[, "n_targets"]) +
    0.02 * safe_graph_rescale(X[, "AnyTargetMutation"])
}

supervised_graph_propagation_predict <- function(
    A_norm,
    X,
    y,
    train_idx,
    alpha = 0.80,
    iterations = 50L
) {
  feature_signal <- validation_feature_signal(X)
  prior <- mean(y[train_idx], na.rm = TRUE)
  if (!is.finite(prior)) prior <- 0.5

  probability <- 0.55 * feature_signal + 0.45 * prior
  seed_signal <- probability
  seed_signal[train_idx] <- 0.90 * y[train_idx] + 0.10 * feature_signal[train_idx]

  for (iteration in seq_len(as.integer(iterations))) {
    probability <- alpha * as.numeric(A_norm %*% probability) + (1 - alpha) * feature_signal
    probability[train_idx] <- seed_signal[train_idx]
  }

  pmin(pmax(probability, 0.001), 0.999)
}

train_gcn_validation_fold <- function(
    A_norm,
    X,
    y,
    train_idx,
    epochs = 250,
    hidden_dim = 16,
    seed = 42
) {
  if (!has_torch) stop("The optional R package 'torch' is not installed.")
  if (length(unique(y[train_idx])) < 2) stop("The training fold contains only one outcome class.")

  torch <- asNamespace("torch")
  set.seed(seed)
  torch$torch_manual_seed(seed)

  x_t <- torch$torch_tensor(X, dtype = torch$torch_float())
  a_t <- torch$torch_tensor(A_norm, dtype = torch$torch_float())
  y_t <- torch$torch_tensor(as.numeric(y), dtype = torch$torch_float())
  train_index_t <- torch$torch_tensor(as.integer(train_idx), dtype = torch$torch_long())

  ValidationGCN <- torch$nn_module(
    "ValidationDrugGCN",
    initialize = function(n_features, hidden) {
      self$lin1 <- torch$nn_linear(n_features, hidden)
      self$lin2 <- torch$nn_linear(hidden, 1)
    },
    forward = function(x, adjacency) {
      h <- torch$torch_matmul(adjacency, x)
      h <- self$lin1(h)
      h <- torch$nnf_relu(h)
      h <- torch$nnf_dropout(h, p = 0.20, training = self$training)
      h <- torch$torch_matmul(adjacency, h)
      self$lin2(h)$squeeze(2)
    }
  )

  model <- ValidationGCN(ncol(X), as.integer(hidden_dim))
  optimizer <- torch$optim_adam(model$parameters, lr = 0.01, weight_decay = 5e-4)

  for (epoch in seq_len(as.integer(epochs))) {
    model$train()
    optimizer$zero_grad()
    logits <- model(x_t, a_t)
    train_logits <- logits$index_select(1, train_index_t)
    train_labels <- y_t$index_select(1, train_index_t)
    loss <- torch$nnf_binary_cross_entropy_with_logits(train_logits, train_labels)
    loss$backward()
    optimizer$step()
  }

  model$eval()
  as.numeric(torch$torch_sigmoid(model(x_t, a_t))$detach()$cpu())
}

fit_pharmacology_benchmark <- function(X, y, train_idx, test_idx) {
  benchmark_features <- intersect(
    c("MeanSynergy", "Delta_vs_Reference", "MeanObserved", "MeanMonotherapyResponse", "MonotherapyAUC"),
    colnames(X)
  )

  train_dat <- as.data.frame(X[train_idx, benchmark_features, drop = FALSE])
  train_dat$Label <- y[train_idx]
  test_dat <- as.data.frame(X[test_idx, benchmark_features, drop = FALSE])

  fit <- suppressWarnings(stats::glm(Label ~ ., data = train_dat, family = stats::binomial()))
  prediction <- suppressWarnings(stats::predict(fit, newdata = test_dat, type = "response"))
  prediction[!is.finite(prediction)] <- mean(y[train_idx], na.rm = TRUE)
  pmin(pmax(as.numeric(prediction), 0.001), 0.999)
}

roc_curve_points <- function(y, probability, model_label, model_id, cell_line) {
  thresholds <- sort(unique(c(Inf, probability, -Inf)), decreasing = TRUE)
  map_dfr(thresholds, function(threshold) {
    prediction <- as.integer(probability >= threshold)
    tp <- sum(prediction == 1L & y == 1L)
    tn <- sum(prediction == 0L & y == 0L)
    fp <- sum(prediction == 1L & y == 0L)
    fn <- sum(prediction == 0L & y == 1L)
    tibble(
      ModelID = model_id,
      CellLineName = cell_line,
      Benchmark = model_label,
      Threshold = threshold,
      FalsePositiveRate = safe_divide(fp, fp + tn),
      TruePositiveRate = safe_divide(tp, tp + fn)
    )
  }) %>%
    mutate(
      FalsePositiveRate = replace_na(FalsePositiveRate, 0),
      TruePositiveRate = replace_na(TruePositiveRate, 0)
    ) %>%
    arrange(FalsePositiveRate, TruePositiveRate)
}

decision_curve_points <- function(y, probability, model_label, model_id, cell_line) {
  n <- length(y)
  prevalence <- mean(y == 1L)
  thresholds <- seq(0.05, 0.95, by = 0.05)

  map_dfr(thresholds, function(threshold) {
    prediction <- as.integer(probability >= threshold)
    tp <- sum(prediction == 1L & y == 1L)
    fp <- sum(prediction == 1L & y == 0L)
    odds <- threshold / (1 - threshold)

    tibble(
      ModelID = model_id,
      CellLineName = cell_line,
      Benchmark = model_label,
      Threshold = threshold,
      NetBenefit = tp / n - fp / n * odds,
      TreatAllNetBenefit = prevalence - (1 - prevalence) * odds,
      TreatNoneNetBenefit = 0
    )
  })
}

combination2_count <- function(x) {
  x <- as.numeric(x)
  x * (x - 1) / 2
}

rand_statistics <- function(reference, cluster) {
  tab <- table(reference, cluster)
  n <- sum(tab)
  if (n < 2) return(c(RandIndex = NA_real_, AdjustedRandIndex = NA_real_))

  a <- sum(combination2_count(tab))
  row_pairs <- sum(combination2_count(rowSums(tab)))
  col_pairs <- sum(combination2_count(colSums(tab)))
  total_pairs <- combination2_count(n)
  agreements_same <- a
  agreements_different <- total_pairs - row_pairs - col_pairs + a
  rand_index <- (agreements_same + agreements_different) / total_pairs

  expected <- row_pairs * col_pairs / total_pairs
  max_index <- 0.5 * (row_pairs + col_pairs)
  adjusted_rand <- ifelse(max_index == expected, NA_real_, (a - expected) / (max_index - expected))

  c(RandIndex = rand_index, AdjustedRandIndex = adjusted_rand)
}

mean_silhouette_score <- function(distance_matrix, clusters) {
  clusters <- as.integer(as.factor(clusters))
  n <- length(clusters)
  if (n < 3 || length(unique(clusters)) < 2) return(NA_real_)

  scores <- vapply(seq_len(n), function(i) {
    same <- which(clusters == clusters[i] & seq_len(n) != i)
    other_clusters <- setdiff(unique(clusters), clusters[i])
    if (length(same) == 0 || length(other_clusters) == 0) return(0)

    a <- mean(distance_matrix[i, same])
    b <- min(vapply(other_clusters, function(k) {
      mean(distance_matrix[i, clusters == k])
    }, numeric(1)))

    if (!is.finite(a) || !is.finite(b) || max(a, b) == 0) return(0)
    (b - a) / max(a, b)
  }, numeric(1))

  mean(scores, na.rm = TRUE)
}

graph_cluster_validation <- function(A_norm, X_scaled, y, seed = 42) {
  embedding <- cbind(A_norm %*% X_scaled, A_norm)
  embedding <- safe_standardise_matrix(embedding)

  set.seed(seed)
  clustering <- tryCatch(
    stats::kmeans(embedding, centers = 2, nstart = 25)$cluster,
    error = function(e) rep(1L, length(y))
  )

  distance_matrix <- as.matrix(stats::dist(embedding))
  rand <- rand_statistics(y, clustering)

  tibble(
    SilhouetteScore = mean_silhouette_score(distance_matrix, clustering),
    RandIndex = unname(rand[["RandIndex"]]),
    AdjustedRandIndex = unname(rand[["AdjustedRandIndex"]])
  )
}

run_validation_agent <- function(
    graph_features,
    synergy_method,
    edge_threshold = 0.20,
    requested_folds = 5L,
    decision_threshold = 0.50,
    epochs = 250,
    hidden_dim = 16,
    engine = "auto"
) {
  if (is.null(graph_features) || nrow(graph_features) == 0) {
    return(list(
      predictions = tibble(), fold_metrics = tibble(), summary = tibble(),
      confusion = tibble(), roc = tibble(), decision_curve = tibble(),
      clustering = tibble()
    ))
  }

  feature_names <- c(
    "MeanSynergy", "MaxSynergy", "MinSynergy", "SDSynergy",
    "ReferenceMeanSynergy", "Delta_vs_Reference", "MeanExpected", "MeanObserved",
    "CombinationObservedGain_vs_Monotherapy",
    "CombinationExpectedGain_vs_Monotherapy",
    "MeanMonotherapyResponse", "MaxMonotherapyResponse", "MonotherapyAUC",
    "log10EC50", "MeanCRISPR", "MaxCRISPR", "MeanExpression",
    "MaxExpression", "MeanExpressionPercentile", "n_targets",
    "AnyTargetMutation", "MutationCountTotal"
  )

  model_outputs <- map(unique(graph_features$ModelID), function(model_id) {
    dat <- graph_features %>%
      filter(ModelID == model_id) %>%
      arrange(ThirdDrugCompoundID)

    y <- as.integer(dat$Label)
    cell_line <- dat$CellLineName[[1]]
    fold_spec <- make_stratified_folds(
      y,
      requested_folds = requested_folds,
      seed = make_seed(model_id, paste0(synergy_method, "_validation"))
    )

    if (is.null(fold_spec)) {
      return(list(
        predictions = tibble(), fold_metrics = tibble(), roc = tibble(),
        decision_curve = tibble(), clustering = tibble(
          ModelID = model_id, CellLineName = cell_line,
          SynergyMethod = synergy_method_label(synergy_method),
          n_nodes = nrow(dat), SilhouetteScore = NA_real_,
          RandIndex = NA_real_, AdjustedRandIndex = NA_real_,
          ValidationNote = "Cross-validation unavailable because both classes require at least two candidates."
        )
      ))
    }

    graph <- build_drug_graph(
      model_tbl = dat,
      edge_threshold = edge_threshold,
      synergy_method = synergy_method
    )
    A_norm <- normalise_adjacency(graph$adjacency)
    X_raw <- dat %>% select(all_of(feature_names)) %>% as.data.frame()
    X_scaled <- safe_standardise_matrix(X_raw)
    colnames(X_scaled) <- feature_names

    selected_engine <- engine
    if (selected_engine == "auto") selected_engine <- if (has_torch) "gcn" else "propagation"

    fold_predictions <- map_dfr(seq_len(fold_spec$k), function(fold_number) {
      test_idx <- which(fold_spec$fold_id == fold_number)
      train_idx <- which(fold_spec$fold_id != fold_number)

      graph_probability <- tryCatch(
        {
          if (selected_engine == "gcn") {
            train_gcn_validation_fold(
              A_norm = A_norm,
              X = X_scaled,
              y = y,
              train_idx = train_idx,
              epochs = epochs,
              hidden_dim = hidden_dim,
              seed = make_seed(model_id, paste0(synergy_method, "_fold_", fold_number))
            )
          } else {
            supervised_graph_propagation_predict(A_norm, X_raw, y, train_idx)
          }
        },
        error = function(e) supervised_graph_propagation_predict(A_norm, X_raw, y, train_idx)
      )

      pharmacology_probability <- tryCatch(
        fit_pharmacology_benchmark(X_scaled, y, train_idx, test_idx),
        error = function(e) rep(mean(y[train_idx]), length(test_idx))
      )
      prevalence_probability <- rep(mean(y[train_idx]), length(test_idx))

      bind_rows(
        tibble(
          ModelID = model_id,
          CellLineName = cell_line,
          ThirdDrugCompoundID = dat$ThirdDrugCompoundID[test_idx],
          ThirdDrug = dat$ThirdDrug[test_idx],
          SynergyMethod = synergy_method_label(synergy_method),
          Fold = fold_number,
          Benchmark = ifelse(selected_engine == "gcn" && has_torch, "Graph convolutional network", "Graph propagation"),
          ObservedLabel = y[test_idx],
          Probability = graph_probability[test_idx]
        ),
        tibble(
          ModelID = model_id,
          CellLineName = cell_line,
          ThirdDrugCompoundID = dat$ThirdDrugCompoundID[test_idx],
          ThirdDrug = dat$ThirdDrug[test_idx],
          SynergyMethod = synergy_method_label(synergy_method),
          Fold = fold_number,
          Benchmark = "Pharmacology-only logistic regression",
          ObservedLabel = y[test_idx],
          Probability = pharmacology_probability
        ),
        tibble(
          ModelID = model_id,
          CellLineName = cell_line,
          ThirdDrugCompoundID = dat$ThirdDrugCompoundID[test_idx],
          ThirdDrug = dat$ThirdDrug[test_idx],
          SynergyMethod = synergy_method_label(synergy_method),
          Fold = fold_number,
          Benchmark = "Prevalence baseline",
          ObservedLabel = y[test_idx],
          Probability = prevalence_probability
        )
      )
    })

    fold_metrics <- fold_predictions %>%
      group_by(ModelID, CellLineName, SynergyMethod, Benchmark, Fold) %>%
      group_modify(~ confusion_metrics(.x$ObservedLabel, .x$Probability, decision_threshold)) %>%
      ungroup()

    roc <- split(fold_predictions, fold_predictions$Benchmark) %>%
      imap_dfr(~ roc_curve_points(
        .x$ObservedLabel, .x$Probability, .y,
        model_id = model_id, cell_line = cell_line
      ))

    decision_curve <- split(fold_predictions, fold_predictions$Benchmark) %>%
      imap_dfr(~ decision_curve_points(
        .x$ObservedLabel, .x$Probability, .y,
        model_id = model_id, cell_line = cell_line
      ))

    clustering <- graph_cluster_validation(
      A_norm = A_norm,
      X_scaled = X_scaled,
      y = y,
      seed = make_seed(model_id, paste0(synergy_method, "_cluster"))
    ) %>%
      mutate(
        ModelID = model_id,
        CellLineName = cell_line,
        SynergyMethod = synergy_method_label(synergy_method),
        n_nodes = nrow(dat),
        ValidationNote = NA_character_,
        .before = 1
      )

    list(
      predictions = fold_predictions,
      fold_metrics = fold_metrics,
      roc = roc,
      decision_curve = decision_curve,
      clustering = clustering
    )
  })

  predictions <- map_dfr(model_outputs, "predictions")
  fold_metrics <- map_dfr(model_outputs, "fold_metrics")
  roc <- map_dfr(model_outputs, "roc")
  decision_curve <- map_dfr(model_outputs, "decision_curve")
  clustering <- map_dfr(model_outputs, "clustering")

  summary <- predictions %>%
    group_by(SynergyMethod, Benchmark) %>%
    group_modify(~ confusion_metrics(.x$ObservedLabel, .x$Probability, decision_threshold)) %>%
    ungroup() %>%
    left_join(
      fold_metrics %>%
        group_by(SynergyMethod, Benchmark) %>%
        summarise(
          CV_Folds = n_distinct(paste(ModelID, Fold, sep = "_")),
          MeanFoldAccuracy = safe_mean(Accuracy),
          SDFoldAccuracy = stats::sd(Accuracy, na.rm = TRUE),
          MeanFoldROC_AUC = safe_mean(ROC_AUC),
          SDFoldROC_AUC = stats::sd(ROC_AUC, na.rm = TRUE),
          MeanFoldF1 = safe_mean(F1),
          MeanFoldMCC = safe_mean(MCC),
          .groups = "drop"
        ),
      by = c("SynergyMethod", "Benchmark")
    ) %>%
    mutate(DecisionThreshold = decision_threshold)

  confusion <- predictions %>%
    mutate(
      PredictedLabel = as.integer(Probability >= decision_threshold),
      Observed = ifelse(ObservedLabel == 1, "Stronger than reference", "Not stronger"),
      Predicted = ifelse(PredictedLabel == 1, "Stronger than reference", "Not stronger")
    ) %>%
    count(SynergyMethod, Benchmark, Observed, Predicted, name = "Count")

  list(
    predictions = predictions,
    fold_metrics = fold_metrics,
    summary = summary,
    confusion = confusion,
    roc = roc,
    decision_curve = decision_curve,
    clustering = clustering
  )
}
