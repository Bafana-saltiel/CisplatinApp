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

build_graph_feature_table <- function(third_ranking) {
  if (is.null(third_ranking) || nrow(third_ranking) == 0) return(tibble())

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
    left_join(evidence_summary,
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
          MeanMonotherapyResponse, MaxMonotherapyResponse, MonotherapyAUC,
          log10EC50, MeanCRISPR, MaxCRISPR, MeanExpression, MaxExpression,
          MeanExpressionPercentile, n_targets, AnyTargetMutation,
          MutationCountTotal
        ),
        ~ replace_na(as.numeric(.x), 0)
      )
    )
}

build_drug_graph <- function(model_tbl, edge_threshold = 0.20) {
  if (nrow(model_tbl) < 2) {
    return(list(adjacency = matrix(1, nrow(model_tbl), nrow(model_tbl)), edges = tibble()))
  }

  target_tokens <- map(model_tbl$TargetGenes, tokenise_graph_text)
  mechanism_tokens <- map(model_tbl$MechanismText, tokenise_graph_text)

  dep_scale <- stats::sd(model_tbl$MeanCRISPR, na.rm = TRUE)
  expr_scale <- stats::sd(model_tbl$MeanExpression, na.rm = TRUE)
  if (!is.finite(dep_scale) || dep_scale == 0) dep_scale <- 1
  if (!is.finite(expr_scale) || expr_scale == 0) expr_scale <- 1

  n <- nrow(model_tbl)
  A <- diag(1, n)
  edge_rows <- vector("list", n * (n - 1) / 2)
  k <- 0L

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      target_sim <- jaccard_similarity(target_tokens[[i]], target_tokens[[j]])
      mechanism_sim <- jaccard_similarity(mechanism_tokens[[i]], mechanism_tokens[[j]])
      dep_sim <- scaled_numeric_similarity(model_tbl$MeanCRISPR[[i]], model_tbl$MeanCRISPR[[j]], dep_scale)
      expr_sim <- scaled_numeric_similarity(model_tbl$MeanExpression[[i]], model_tbl$MeanExpression[[j]], expr_scale)
      mutation_sim <- as.numeric(model_tbl$AnyTargetMutation[[i]] == model_tbl$AnyTargetMutation[[j]])

      weight <-
        0.45 * target_sim +
        0.25 * mechanism_sim +
        0.15 * dep_sim +
        0.10 * expr_sim +
        0.05 * mutation_sim

      if (is.finite(weight) && weight >= edge_threshold) {
        A[i, j] <- weight
        A[j, i] <- weight
        k <- k + 1L
        edge_rows[[k]] <- tibble(
          from = model_tbl$ThirdDrug[[i]],
          to = model_tbl$ThirdDrug[[j]],
          weight = weight,
          shared_target = target_sim,
          shared_mechanism = mechanism_sim,
          dependency_similarity = dep_sim,
          expression_similarity = expr_sim,
          mutation_similarity = mutation_sim
        )
      }
    }
  }

  edges <- if (k == 0L) tibble() else bind_rows(edge_rows[seq_len(k)])
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

binary_auc <- function(y, p) {
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]; p <- p[ok]
  if (length(unique(y)) < 2) return(NA_real_)
  pos <- p[y == 1]; neg <- p[y == 0]
  if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
  mean(outer(pos, neg, ">")) + 0.5 * mean(outer(pos, neg, "=="))
}

# Transparent graph-propagation fallback used when torch is unavailable or the
# model contains too few labelled examples for stable GCN fitting.
graph_propagation_predict <- function(A_norm, X, y) {
  base_signal <-
    0.30 * scales::rescale(X[, "MeanMonotherapyResponse"], to = c(0, 1)) +
    0.20 * scales::rescale(X[, "MonotherapyAUC"], to = c(0, 1)) +
    0.20 * scales::rescale(X[, "MaxCRISPR"], to = c(0, 1)) +
    0.15 * scales::rescale(X[, "MeanExpressionPercentile"], to = c(0, 1)) +
    0.10 * scales::rescale(X[, "n_targets"], to = c(0, 1)) +
    0.05 * scales::rescale(X[, "AnyTargetMutation"], to = c(0, 1))
  propagated <- as.numeric(A_norm %*% base_signal)
  scales::rescale(propagated, to = c(0.01, 0.99))
}

train_gcn_model <- function(A_norm, X, y, epochs = 250, hidden_dim = 16, seed = 42) {
  if (!has_torch) stop("The optional R package 'torch' is not installed.")
  if (length(unique(y)) < 2 || length(y) < 10) stop("Insufficient class variation for GCN training.")

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
    edge_threshold = 0.20,
    epochs = 250,
    hidden_dim = 16,
    engine = "auto"
) {
  features <- build_graph_feature_table(third_ranking)
  if (nrow(features) == 0) {
    return(list(ranking = tibble(), diagnostics = tibble(), edges = tibble(), features = tibble()))
  }

  feature_names <- c(
    "MeanMonotherapyResponse", "MaxMonotherapyResponse", "MonotherapyAUC",
    "log10EC50", "MeanCRISPR", "MaxCRISPR", "MeanExpression",
    "MaxExpression", "MeanExpressionPercentile", "n_targets",
    "AnyTargetMutation", "MutationCountTotal"
  )

  model_results <- map(unique(features$ModelID), function(model_id) {
    dat <- features %>% filter(ModelID == model_id) %>% arrange(ThirdDrugCompoundID)
    graph <- build_drug_graph(dat, edge_threshold = edge_threshold)
    A_norm <- normalise_adjacency(graph$adjacency)
    X_raw <- dat %>% select(all_of(feature_names)) %>% as.data.frame()
    X_scaled <- safe_standardise_matrix(X_raw)
    colnames(X_scaled) <- feature_names
    y <- dat$Label

    selected_engine <- engine
    if (selected_engine == "auto") selected_engine <- if (has_torch) "gcn" else "propagation"

    fit <- tryCatch(
      {
        if (selected_engine == "gcn") {
          train_gcn_model(A_norm, X_scaled, y, epochs, hidden_dim, seed = make_seed(model_id, "GCN"))
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
        GraphPredictedClass = ifelse(GraphProbability >= 0.5, "Predicted stronger than reference", "Predicted not stronger"),
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
      n_nodes = nrow(dat),
      n_edges = nrow(graph$edges),
      positive_labels = sum(y == 1),
      negative_labels = sum(y == 0),
      Engine = engine_used,
      EvaluationSet = ifelse(length(test_idx) > 0, "Held-out 20%", "Training/internal"),
      Accuracy = mean((eval_p >= 0.5) == eval_y),
      AUC = binary_auc(eval_y, eval_p),
      FallbackReason = fit$fallback_reason %||% NA_character_
    )

    edges <- graph$edges %>%
      mutate(ModelID = model_id, CellLineName = dat$CellLineName[[1]], .before = 1)

    list(ranking = ranked, diagnostics = diagnostics, edges = edges)
  })

  list(
    ranking = map_dfr(model_results, "ranking"),
    diagnostics = map_dfr(model_results, "diagnostics"),
    edges = map_dfr(model_results, "edges"),
    features = features
  )
}

