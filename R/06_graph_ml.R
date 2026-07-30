# ============================================================
# 7. GRAPH NEURAL NETWORK PRIORITISATION MODULE
# ============================================================

# This module constructs one drug graph per cervical cancer model. Each node is
# a candidate third drug. Edges encode shared targets, shared mechanism/pathway
# annotations, and similarity in CRISPR dependency, expression and mutation
# context. A two-layer graph convolutional network (GCN) predicts whether a
# candidate will exceed the model-specific clinical reference interaction.
#
# IMPORTANT: labels are generated from the current pharmacological screen
# (Interaction == "Stronger_than_reference"). Therefore, this is an internal
# graph-learning prioritisation layer, not an externally validated clinical
# response predictor. External combination datasets are required for a fully
# independent predictive model.

tokenise_graph_text <- function(x) {
  value <- as.character(x %||% "")
  value[is.na(value)] <- ""
  value <- stringr::str_to_lower(value)
  value <- stringr::str_replace_all(value, "[^a-z0-9]+", " ")
  tokens <- unlist(stringr::str_split(value, "\\s+"), use.names = FALSE)
  tokens <- tokens[!is.na(tokens) & nzchar(tokens) & nchar(tokens) >= 3]
  unique(tokens)
}

jaccard_similarity <- function(a, b) {
  a <- unique(a); b <- unique(b)
  if (length(a) == 0 && length(b) == 0) return(0)
  union_n <- length(union(a, b))
  if (union_n == 0) return(0)
  length(intersect(a, b)) / union_n
}

scaled_numeric_similarity <- function(a, b, scale_value = 1) {
  if (!is.finite(a) || !is.finite(b)) return(0)
  1 / (1 + abs(a - b) / max(scale_value, 1e-8))
}

empty_graph_edges <- function() {
  tibble(
    from = character(),
    to = character(),
    weight = double(),
    SynergyMethod = character(),
    shared_target = double(),
    shared_mechanism = double(),
    selected_synergy_similarity = double(),
    delta_reference_similarity = double(),
    observed_response_similarity = double(),
    expected_response_similarity = double(),
    monotherapy_similarity = double(),
    dependency_similarity = double(),
    expression_similarity = double(),
    mutation_similarity = double()
  )
}

build_graph_feature_table <- function(third_ranking, synergy_method = NULL) {
  if (is.null(third_ranking) || nrow(third_ranking) == 0) return(tibble())

  # The pharmacology module normally writes SynergyMethod into third_ranking.
  # The explicit argument from the server is retained as a defensive fallback.
  if (!"SynergyMethod" %in% names(third_ranking)) {
    third_ranking$SynergyMethod <- synergy_method %||% "Selected synergy model"
  } else if (!is.null(synergy_method)) {
    third_ranking$SynergyMethod <- as.character(synergy_method)
  }

  target_tbl <- add_target_annotation(third_ranking)
  if (nrow(target_tbl) == 0) return(tibble())

  evidence_tbl <- target_tbl %>%
    mutate(
      Functional_Profile = map2(ModelID, TargetGene, get_crispr_dependency),
      Molecular_Profile = map2(ModelID, TargetGene, get_target_expression),
      Biomarker_Profile = map2(ModelID, TargetGene, get_target_mutation)
    ) %>%
    unnest(Functional_Profile) %>%
    unnest(Molecular_Profile) %>%
    unnest(Biomarker_Profile)

  evidence_summary <- evidence_tbl %>%
    group_by(ModelID, CellLineName, ThirdDrugCompoundID, ThirdDrug) %>%
    summarise(
      TargetGenes = paste(sort(unique(na.omit(TargetGene))), collapse = ";"),
      n_targets = n_distinct(TargetGene[!is.na(TargetGene)]),
      MeanCRISPR = safe_mean(CRISPR_Dependency),
      MaxCRISPR = safe_max(CRISPR_Dependency),
      MeanExpression = safe_mean(Expression_log2TPM),
      MaxExpression = safe_max(Expression_log2TPM),
      MeanExpressionPercentile = safe_mean(ExpressionPercentile),
      AnyTargetMutation = as.integer(any(MutationStatus == "Mutated", na.rm = TRUE)),
      MutationCountTotal = sum(MutationCount, na.rm = TRUE),
      .groups = "drop"
    )

  third_ranking %>%
    left_join(
      evidence_summary,
      by = c("ModelID", "CellLineName", "ThirdDrugCompoundID", "ThirdDrug")
    ) %>%
    left_join(
      PortalCompounds %>%
        select(CompoundID, GeneSymbolOfTargets, TargetOrMechanism) %>%
        distinct(CompoundID, .keep_all = TRUE),
      by = c("ThirdDrugCompoundID" = "CompoundID")
    ) %>%
    mutate(
      TargetGenes = coalesce(TargetGenes, as.character(GeneSymbolOfTargets), ""),
      MechanismText = coalesce(as.character(TargetOrMechanism), ""),
      Label = as.integer(Interaction == "Stronger_than_reference"),
      log10EC50 = ifelse(is.finite(EC50) & EC50 > 0, log10(EC50), NA_real_),
      across(
        c(
          MeanSynergy, MaxSynergy, MinSynergy, SDSynergy,
          ReferenceMeanSynergy, Delta_vs_Reference, MeanExpected, MeanObserved,
          CombinationObservedGain_vs_Monotherapy,
          CombinationExpectedGain_vs_Monotherapy,
          MeanMonotherapyResponse, MaxMonotherapyResponse, MonotherapyAUC,
          log10EC50, MeanCRISPR, MaxCRISPR, MeanExpression, MaxExpression,
          MeanExpressionPercentile, n_targets, AnyTargetMutation,
          MutationCountTotal
        ),
        ~ replace_na(as.numeric(.x), 0)
      )
    )
}

robust_graph_scale <- function(x) {
  x <- as.numeric(x)
  value <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(value) || value == 0) value <- stats::IQR(x, na.rm = TRUE)
  if (!is.finite(value) || value == 0) value <- 1
  value
}

build_drug_graph <- function(model_tbl, edge_threshold = 0.20, synergy_method = NULL) {
  if (nrow(model_tbl) < 2) {
    return(list(
      adjacency = matrix(1, nrow(model_tbl), nrow(model_tbl)),
      edges = empty_graph_edges()
    ))
  }

  if (!"SynergyMethod" %in% names(model_tbl)) {
    model_tbl$SynergyMethod <- synergy_method %||% "Selected synergy model"
  }

  target_tokens <- map(model_tbl$TargetGenes, tokenise_graph_text)
  mechanism_tokens <- map(model_tbl$MechanismText, tokenise_graph_text)

  # These scales are estimated within each cervical model. MeanSynergy,
  # Delta_vs_Reference, MeanObserved and MeanExpected are generated by the
  # selected interaction model; therefore they make the adjacency matrix
  # explicitly Bliss-, HSA-, Loewe- or ZIP-specific.
  synergy_scale <- robust_graph_scale(model_tbl$MeanSynergy)
  delta_scale <- robust_graph_scale(model_tbl$Delta_vs_Reference)
  observed_scale <- robust_graph_scale(model_tbl$MeanObserved)
  expected_scale <- robust_graph_scale(model_tbl$MeanExpected)
  mono_scale <- robust_graph_scale(model_tbl$MeanMonotherapyResponse)
  dep_scale <- robust_graph_scale(model_tbl$MeanCRISPR)
  expr_scale <- robust_graph_scale(model_tbl$MeanExpression)

  n <- nrow(model_tbl)
  A <- diag(1, n)
  edge_rows <- vector("list", n * (n - 1) / 2)
  k <- 0L

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      target_sim <- jaccard_similarity(target_tokens[[i]], target_tokens[[j]])
      mechanism_sim <- jaccard_similarity(mechanism_tokens[[i]], mechanism_tokens[[j]])

      synergy_sim <- scaled_numeric_similarity(
        model_tbl$MeanSynergy[[i]], model_tbl$MeanSynergy[[j]], synergy_scale
      )
      delta_sim <- scaled_numeric_similarity(
        model_tbl$Delta_vs_Reference[[i]], model_tbl$Delta_vs_Reference[[j]], delta_scale
      )
      observed_sim <- scaled_numeric_similarity(
        model_tbl$MeanObserved[[i]], model_tbl$MeanObserved[[j]], observed_scale
      )
      expected_sim <- scaled_numeric_similarity(
        model_tbl$MeanExpected[[i]], model_tbl$MeanExpected[[j]], expected_scale
      )
      mono_sim <- scaled_numeric_similarity(
        model_tbl$MeanMonotherapyResponse[[i]],
        model_tbl$MeanMonotherapyResponse[[j]],
        mono_scale
      )
      dep_sim <- scaled_numeric_similarity(
        model_tbl$MeanCRISPR[[i]], model_tbl$MeanCRISPR[[j]], dep_scale
      )
      expr_sim <- scaled_numeric_similarity(
        model_tbl$MeanExpression[[i]], model_tbl$MeanExpression[[j]], expr_scale
      )
      mutation_sim <- as.numeric(
        model_tbl$AnyTargetMutation[[i]] == model_tbl$AnyTargetMutation[[j]]
      )

      # Selected-synergy evidence contributes 40% directly through MeanSynergy
      # and Delta_vs_Reference. Expected and observed response contribute an
      # additional 10%, so changing the selected interaction method can change
      # edge retention, weights, centrality, communities and GCN predictions.
      weight <-
        0.20 * target_sim +
        0.12 * mechanism_sim +
        0.20 * synergy_sim +
        0.20 * delta_sim +
        0.05 * observed_sim +
        0.05 * expected_sim +
        0.06 * mono_sim +
        0.06 * dep_sim +
        0.04 * expr_sim +
        0.02 * mutation_sim

      if (is.finite(weight) && weight >= edge_threshold) {
        A[i, j] <- weight
        A[j, i] <- weight
        k <- k + 1L
        edge_rows[[k]] <- tibble(
          from = model_tbl$ThirdDrug[[i]],
          to = model_tbl$ThirdDrug[[j]],
          weight = weight,
          SynergyMethod = as.character(model_tbl$SynergyMethod[[i]]),
          shared_target = target_sim,
          shared_mechanism = mechanism_sim,
          selected_synergy_similarity = synergy_sim,
          delta_reference_similarity = delta_sim,
          observed_response_similarity = observed_sim,
          expected_response_similarity = expected_sim,
          monotherapy_similarity = mono_sim,
          dependency_similarity = dep_sim,
          expression_similarity = expr_sim,
          mutation_similarity = mutation_sim
        )
      }
    }
  }

  edges <- if (k == 0L) {
    empty_graph_edges()
  } else {
    bind_rows(edge_rows[seq_len(k)])
  }

  list(adjacency = A, edges = edges)
}

normalise_adjacency <- function(A) {
  d <- rowSums(A)
  d_inv_sqrt <- ifelse(d > 0, 1 / sqrt(d), 0)
  diag(d_inv_sqrt) %*% A %*% diag(d_inv_sqrt)
}

safe_standardise_matrix <- function(X) {
  X <- as.matrix(X)
  X[!is.finite(X)] <- 0
  means <- colMeans(X)
  sds <- apply(X, 2, sd)
  sds[!is.finite(sds) | sds == 0] <- 1
  sweep(sweep(X, 2, means, "-"), 2, sds, "/")
}

safe_graph_rescale <- function(x, to = c(0, 1)) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- 0
  if (length(x) == 0) return(numeric())
  if (length(unique(x)) <= 1) return(rep(mean(to), length(x)))
  scales::rescale(x, to = to)
}

binary_auc <- function(y, p) {
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]
  p <- p[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  pos <- p[y == 1]
  neg <- p[y == 0]
  if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
  mean(outer(pos, neg, ">")) + 0.5 * mean(outer(pos, neg, "=="))
}

# Transparent graph-propagation fallback used when torch is unavailable or the
# model contains too few labelled examples for stable GCN fitting.
graph_propagation_predict <- function(A_norm, X, y) {
  required <- c(
    "MeanSynergy", "Delta_vs_Reference", "MeanObserved",
    "CombinationObservedGain_vs_Monotherapy", "MeanMonotherapyResponse",
    "MonotherapyAUC", "MaxCRISPR", "MeanExpressionPercentile",
    "n_targets", "AnyTargetMutation"
  )
  missing <- setdiff(required, colnames(X))
  if (length(missing) > 0) {
    stop("Graph feature matrix is missing: ", paste(missing, collapse = ", "))
  }

  base_signal <-
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

  propagated <- as.numeric(A_norm %*% base_signal)
  safe_graph_rescale(propagated, to = c(0.01, 0.99))
}

train_gcn_model <- function(A_norm, X, y, epochs = 250, hidden_dim = 16, seed = 42) {
  if (!has_torch) stop("The optional R package 'torch' is not installed.")
  if (length(unique(y)) < 2 || length(y) < 10) {
    stop("Insufficient class variation for GCN training.")
  }

  torch <- asNamespace("torch")
  set.seed(seed)
  torch$torch_manual_seed(seed)

  x_t <- torch$torch_tensor(X, dtype = torch$torch_float())
  a_t <- torch$torch_tensor(A_norm, dtype = torch$torch_float())
  y_t <- torch$torch_tensor(as.numeric(y), dtype = torch$torch_float())

  GCN <- torch$nn_module(
    "DrugGCN",
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

  model <- GCN(ncol(X), as.integer(hidden_dim))
  optimizer <- torch$optim_adam(model$parameters, lr = 0.01, weight_decay = 5e-4)

  for (epoch in seq_len(as.integer(epochs))) {
    model$train()
    optimizer$zero_grad()
    logits <- model(x_t, a_t)
    loss <- torch$nnf_binary_cross_entropy_with_logits(logits, y_t)
    loss$backward()
    optimizer$step()
  }

  model$eval()
  probs <- as.numeric(torch$torch_sigmoid(model(x_t, a_t))$detach()$cpu())
  list(probability = probs, train_idx = seq_along(y), test_idx = integer(0))
}

run_graph_learning <- function(
    third_ranking,
    synergy_method = NULL,
    edge_threshold = 0.20,
    epochs = 250,
    hidden_dim = 16,
    engine = "auto"
) {
  features <- build_graph_feature_table(
    third_ranking = third_ranking,
    synergy_method = synergy_method
  )

  if (nrow(features) == 0) {
    return(list(
      ranking = tibble(),
      diagnostics = tibble(),
      edges = tibble(),
      features = tibble()
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

  model_results <- map(unique(features$ModelID), function(model_id) {
    dat <- features %>%
      filter(ModelID == model_id) %>%
      arrange(ThirdDrugCompoundID)

    graph <- build_drug_graph(
      model_tbl = dat,
      edge_threshold = edge_threshold,
      synergy_method = synergy_method
    )
    A_norm <- normalise_adjacency(graph$adjacency)
    X_raw <- dat %>% select(all_of(feature_names)) %>% as.data.frame()
    X_scaled <- safe_standardise_matrix(X_raw)
    colnames(X_scaled) <- feature_names
    y <- dat$Label

    selected_engine <- engine
    if (selected_engine == "auto") {
      selected_engine <- if (has_torch) "gcn" else "propagation"
    }

    fit <- tryCatch(
      {
        if (selected_engine == "gcn") {
          train_gcn_model(
            A_norm, X_scaled, y,
            epochs = epochs,
            hidden_dim = hidden_dim,
            seed = make_seed(model_id, synergy_method %||% "GCN")
          )
        } else {
          list(
            probability = graph_propagation_predict(A_norm, X_raw, y),
            train_idx = seq_len(nrow(dat)),
            test_idx = integer(0)
          )
        }
      },
      error = function(e) {
        list(
          probability = graph_propagation_predict(A_norm, X_raw, y),
          train_idx = seq_len(nrow(dat)),
          test_idx = integer(0),
          fallback_reason = conditionMessage(e)
        )
      }
    )

    engine_used <- if (!is.null(fit$fallback_reason) || selected_engine == "propagation") {
      "Graph propagation"
    } else {
      "Two-layer GCN"
    }

    ranked <- dat %>%
      mutate(
        GraphProbability = pmin(pmax(fit$probability, 0), 1),
        GraphPredictedClass = ifelse(
          GraphProbability >= 0.5,
          "Predicted stronger than reference",
          "Predicted not stronger"
        ),
        GraphRank = rank(-GraphProbability, ties.method = "first"),
        GraphEngine = engine_used
      ) %>%
      arrange(GraphRank)

    test_idx <- fit$test_idx
    eval_y <- if (length(test_idx) > 0) y[test_idx] else y
    eval_p <- if (length(test_idx) > 0) fit$probability[test_idx] else fit$probability

    diagnostics <- tibble(
      ModelID = model_id,
      CellLineName = dat$CellLineName[[1]],
      SynergyMethod = as.character(dat$SynergyMethod[[1]]),
      n_nodes = nrow(dat),
      n_edges = nrow(graph$edges),
      GraphDensity = if (nrow(dat) > 1) {
        2 * nrow(graph$edges) / (nrow(dat) * (nrow(dat) - 1))
      } else {
        0
      },
      MeanEdgeWeight = if (nrow(graph$edges) > 0) mean(graph$edges$weight) else 0,
      positive_labels = sum(y == 1),
      negative_labels = sum(y == 0),
      Engine = engine_used,
      EvaluationSet = ifelse(length(test_idx) > 0, "Held-out 20%", "Training/internal"),
      Accuracy = mean((eval_p >= 0.5) == eval_y),
      AUC = binary_auc(eval_y, eval_p),
      FallbackReason = fit$fallback_reason %||% NA_character_
    )

    edges <- if (nrow(graph$edges) == 0) {
      graph$edges %>%
        mutate(
          ModelID = character(),
          CellLineName = character(),
          .before = 1
        )
    } else {
      graph$edges %>%
        mutate(
          ModelID = model_id,
          CellLineName = dat$CellLineName[[1]],
          .before = 1
        )
    }

    list(ranking = ranked, diagnostics = diagnostics, edges = edges)
  })

  list(
    ranking = map_dfr(model_results, "ranking"),
    diagnostics = map_dfr(model_results, "diagnostics"),
    edges = map_dfr(model_results, "edges"),
    features = features
  )
}
