# Clear all environmental variables
rm(list = ls())

start_time           <- Sys.time()
#=======================================================================================================================#
#                                                   Libraries:
#=======================================================================================================================#
library(dplyr)
library(cluster)
library(tidyr)
library(tibble)
library(here)
library(ggplot2)
library(cowplot)
library(progress)

#=======================================================================================================================#
#                                               Import & Clean Data:
#=======================================================================================================================#
folderLoc <- here("AdvancedClustering")
dfraw <- readRDS(paste0(folderLoc, "\\Testing\\Valarray.rds")) %>% as.data.frame()
dfMaster <- readRDS(paste0(folderLoc, "\\Testing\\InputOutputTable.rds")) %>% as.data.frame()

lenRow <- dim(dfraw)[2]

# Extract day and sample numbers
col_info <- do.call(rbind, strsplit(colnames(dfraw), "\\.sample_"))
col_info <- data.frame(
  day = as.integer(col_info[,1]),
  sample = as.integer(col_info[,2])
)

C  <- 12
T  <- 366
N  <- lenRow / T

df <- array(
  as.numeric(as.matrix(dfraw)),
  dim      = c(C, T, N),
  dimnames = list(
    compartment = c(
      "Susceptible_Gentoo",    "Infected_Gentoo",    "Recovered_Gentoo",    "HPAI_Killed_Gentoo",
      "Susceptible_Predators", "Infected_Predators", "Recovered_Predators", "HPAI_Killed_Predators",
      "Susceptible_Albatross", "Infected_Albatross", "Recovered_Albatross", "HPAI_Killed_Albatross"),
    time        = paste0("Day_", 1:T),
    sample      = paste0("Sample_", 1:N)
  )
)

#=======================================================================================================================#
#                                                      Cluster:
#=======================================================================================================================#
# Utility: null coalesce
`%||%` <- function(a, b) if (!is.null(a)) a else b

# 1) Build features: per-compartment histograms with shared breaks
build_histogram_features <- function(
    valarray,                      # 3D: [compartment, time, sample]
    bins = 52,
    popvar = NULL,                 # list/df with total_gentoo/total_predators/total_albatross
    relativize_to_species_totals = FALSE,
    compartments_use = c(
      "Infected_Gentoo", "Infected_Predators", "Infected_Albatross",
      "HPAI_Killed_Gentoo", "HPAI_Killed_Predators", "HPAI_Killed_Albatross"
    )
){
  dn <- dimnames(valarray)
  if (is.null(dn) || is.null(dn$compartment)) {
    stop("valarray must have dimnames with a 'compartment' entry.")
  }
  if (!all(compartments_use %in% dn$compartment)) {
    missing <- setdiff(compartments_use, dn$compartment)
    stop("Missing compartments in valarray: ", paste(missing, collapse = ", "))
  }
  
  C <- length(compartments_use)
  T <- dim(valarray)[2]
  N <- dim(valarray)[3]
  
  # Optional per-capita relativization
  if (isTRUE(relativize_to_species_totals)) {
    if (is.null(popvar)) stop("relativize_to_species_totals=TRUE requires popvar.")
    read_total <- function(pv, name) {
      if (is.list(pv) && !is.data.frame(pv))      as.numeric(pv[[name]])
      else if (is.data.frame(pv) && name %in% names(pv)) as.numeric(pv[[name]][1])
      else stop("Could not locate '", name, "' in popvar.")
    }
    total_gentoo    <- read_total(popvar, "total_gentoo")
    total_predators <- read_total(popvar, "total_predators")
    total_albatross <- read_total(popvar, "total_albatross")
    if (any(!is.finite(c(total_gentoo, total_predators, total_albatross))) ||
        any(c(total_gentoo, total_predators, total_albatross) <= 0)) {
      stop("Species totals in popvar must be positive, finite numbers.")
    }
    species_for_comp <- c("Gentoo", "Predators", "Albatross",
                          "Gentoo", "Predators", "Albatross")
    totals <- list(Gentoo=total_gentoo, Predators=total_predators, Albatross=total_albatross)
    for (j in seq_along(compartments_use)) {
      comp <- compartments_use[j]
      denom <- totals[[ species_for_comp[j] ]]
      valarray[comp, , ] <- valarray[comp, , ] / denom
    }
  }
  
  # Shared breaks per compartment
  comp_breaks <- vector("list", C)
  names(comp_breaks) <- compartments_use
  for (j in seq_along(compartments_use)) {
    comp <- compartments_use[j]
    x_all <- as.numeric(valarray[comp, , ])
    x_all <- x_all[is.finite(x_all)]
    if (!length(x_all)) stop("All values non-finite for compartment: ", comp)
    rng <- range(x_all, na.rm = TRUE)
    if (rng[1] == rng[2]) {
      eps <- if (rng[1] == 0) 1e-8 else abs(rng[1]) * 1e-6
      rng <- rng + c(-eps, +eps)
    }
    comp_breaks[[j]] <- seq(rng[1], rng[2], length.out = bins + 1)
  }
  
  # Build feature matrix (density histograms)
  feature_mat <- matrix(NA_real_, nrow = dim(valarray)[3], ncol = 0)
  rn <- dimnames(valarray)$sample %||% paste0("Sample_", seq_len(dim(valarray)[3]))
  rownames(feature_mat) <- rn
  
  for (j in seq_along(compartments_use)) {
    comp <- compartments_use[j]
    brks <- comp_breaks[[j]]
    block <- matrix(NA_real_, nrow = length(rn), ncol = bins)
    colnames(block) <- paste0(comp, "_bin", seq_len(bins))
    for (i in seq_along(rn)) {
      ts_ij <- as.numeric(valarray[comp, , i])
      h <- hist(ts_ij, breaks = brks, plot = FALSE, right = TRUE, include.lowest = TRUE)
      dens <- h$density
      if (any(!is.finite(dens))) {
        cnt <- h$counts; s <- sum(cnt); dens <- if (s > 0) cnt / s else rep(0, bins)
      }
      block[i, ] <- dens
    }
    feature_mat <- cbind(feature_mat, block)
  }
  
  # Drop zero-variance features (across runs)
  keep_cols <- which(apply(feature_mat, 2, function(v) stats::sd(v, na.rm=TRUE) > 0))
  feature_mat <- feature_mat[, keep_cols, drop = FALSE]
  
  list(
    features = feature_mat,
    breaks   = comp_breaks,
    compartments = compartments_use
  )
}

# 2) Evaluate k over a range using multiple criteria (no external deps)
evaluate_kmeans_k <- function(
    feature_mat,
    k_min = 2, k_max = 20,
    nstart = 50, iter.max = 100,
    scale_features = TRUE,
    compute_gap = TRUE,  # set FALSE if speed is critical
    gap_B = 20,          # number of reference datasets for Gap statistic
    seed = 1L
){
  stopifnot(nrow(feature_mat) >= k_max)
  set.seed(seed)
  
  X <- feature_mat
  if (isTRUE(scale_features)) {
    X <- scale(X)
  }
  
  N <- nrow(X)
  p <- ncol(X)
  # Precompute pairwise distances once for silhouette; we recompute within each k for cluster labels
  D <- as.matrix(dist(X, method = "euclidean"))
  
  # total sum of squares (about the grand mean)
  grand_mean <- colMeans(X)
  TSS <- sum(rowSums((X - matrix(grand_mean, nrow=N, ncol=p, byrow=TRUE))^2))
  
  ks <- k_min:k_max
  out <- data.frame(
    k = ks,
    WSS = NA_real_,                # within-cluster sum of squares
    CH  = NA_real_,                # Calinski–Harabasz (higher is better)
    DB  = NA_real_,                # Davies–Bouldin (lower is better)
    Silhouette = NA_real_,         # mean silhouette width (higher is better)
    Gap = NA_real_,                # Tibshirani gap (higher is better)
    Gap_se = NA_real_
  )
  
  # helper: compute cluster stats used for CH and DB
  cluster_stats <- function(X, cl, centers) {
    k <- length(unique(cl))
    # WSS
    WSS <- 0
    S_i <- numeric(k)     # average within-cluster distance to center (for DB)
    for (j in seq_len(k)) {
      idx <- which(cl == j)
      if (length(idx) == 0) next
      Cj <- centers[j, , drop = FALSE]
      diffs <- X[idx, , drop = FALSE] - matrix(Cj, nrow=length(idx), ncol=ncol(X), byrow=TRUE)
      sq <- rowSums(diffs^2)
      WSS <- WSS + sum(sq)
      # S_i: mean Euclidean distance to center
      S_i[j] <- mean(sqrt(sq))
    }
    list(WSS = WSS, S_i = S_i)
  }
  
  # helper: Davies–Bouldin
  davies_bouldin <- function(centers, S_i) {
    k <- nrow(centers)
    if (k == 1) return(NA_real_)
    # inter-center distances
    Cdist <- as.matrix(dist(centers))
    # for each i, compute max_j ( (S_i + S_j) / d(c_i,c_j) )
    R_i <- numeric(k)
    for (i in seq_len(k)) {
      vals <- (S_i[i] + S_i[-i]) / Cdist[i, -i]
      R_i[i] <- max(vals[is.finite(vals)])
    }
    mean(R_i)
  }
  
  mean_silhouette <- function(D, cl) {
    # D: full NxN distance matrix (Euclidean)
    # cl: integer cluster labels (length N)
    labs <- sort(unique(cl))
    k <- length(labs)
    if (k < 2) return(NA_real_)  # silhouette undefined for k=1
    
    N <- nrow(D)
    a <- numeric(N)  # mean intra-cluster distance (exclude self)
    b <- numeric(N)  # mean distance to the nearest *other* cluster
    
    for (lab in labs) {
      idx <- which(cl == lab)
      # a_i: mean distance to points in own cluster (exclude self)
      if (length(idx) > 1) {
        Di <- D[idx, idx, drop = FALSE]
        diag(Di) <- NA_real_                   # exclude self
        a[idx] <- rowMeans(Di, na.rm = TRUE)
      } else {
        a[idx] <- 0                            # singleton cluster: define a_i = 0
      }
      
      # b_i: for each other cluster, the mean distance to that cluster; take min
      other_labs <- setdiff(labs, lab)
      if (length(other_labs) == 0) {
        b[idx] <- 0
      } else {
        other_means <- lapply(other_labs, function(l2) {
          jdx <- which(cl == l2)
          rowMeans(D[idx, jdx, drop = FALSE], na.rm = TRUE)
        })
        Bmat <- do.call(cbind, other_means)
        if (is.null(dim(Bmat))) {
          # Only one other cluster -> Bmat is a vector
          b[idx] <- as.numeric(Bmat)
        } else {
          b[idx] <- apply(Bmat, 1, min)
        }
      }
    }
    
    s <- (b - a) / pmax(a, b)
    mean(s[is.finite(s)])
  }
  
  # helper: Gap statistic (Tibshirani et al. 2001)
  compute_gap_stat <- function(X, k, B = 20) {
    N <- nrow(X); p <- ncol(X)
    # fit kmeans on X to get W_k
    km <- kmeans_kl_fast(X, k = k, nstart=25)
    W  <- as.numeric(km$tot.withinss)   # total within-cluster SS (equivalent to sum WSS over clusters)
    
    # uniform reference in bounding box of X
    mins <- apply(X, 2, min); maxs <- apply(X, 2, max)
    lW <- numeric(B)
    for (b in seq_len(B)) {
      U <- matrix(runif(N * p, mins, maxs), nrow=N, ncol=p)
      kmb <- kmeans(U, centers = k, nstart=10)
      lW[b] <- log(as.numeric(kmb$tot.withinss))
    }
    gap    <- mean(lW) - log(W)
    gap_se <- sqrt(1 + 1/B) * sd(lW)
    c(gap=gap, se=gap_se)
  }
  
  pb <- progress_bar$new(
    format = "Cluster Evaluation [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
    total = length(ks), clear = FALSE, width = 90
  )
  
  for (idx in seq_along(ks)) {
    k <- ks[idx]
    km <- kmeans_kl_fast(X, k = k, nstart = nstart)
    stats <- cluster_stats(X, km$cluster, km$centers)
    WSS <- stats$WSS
    # Calinski–Harabasz
    CH <- ((TSS - WSS) / (k - 1)) / (WSS / (N - k))
    # Davies–Bouldin
    DB <- davies_bouldin(km$centers, stats$S_i)
    # Silhouette (mean)
    Sil <- mean_silhouette(D, km$cluster)
    # Gap (optional)
    Gap <- NA_real_; Gap_se <- NA_real_
    if (isTRUE(compute_gap)) {
      gs <- compute_gap_stat(X, k, B = gap_B)
      Gap <- gs["gap"]; Gap_se <- gs["se"]
    }
    out[idx, c("WSS","CH","DB","Silhouette","Gap","Gap_se")] <- c(WSS, CH, DB, Sil, Gap, Gap_se)
    pb$tick()
  }
  
  # Simple recommendations by criterion:
  # - CH: maximize
  # - Silhouette: maximize
  # - DB: minimize
  # - Gap: choose smallest k s.t. Gap(k) >= Gap(k+1) - se(k+1) (Tibshirani rule)
  rec <- list()
  rec$k_CH  <- out$k[ which.max(out$CH) ]
  rec$k_Sil <- out$k[ which.max(out$Silhouette) ]
  rec$k_DB  <- out$k[ which.min(out$DB) ]
  if (isTRUE(compute_gap)) {
    kstar <- out$k[1]
    for (i in seq_len(nrow(out)-1)) {
      if (!is.na(out$Gap[i]) && !is.na(out$Gap[i+1]) && !is.na(out$Gap_se[i+1])) {
        if (out$Gap[i] >= out$Gap[i+1] - out$Gap_se[i+1]) { kstar <- out$k[i]; break }
      }
    }
    rec$k_Gap <- kstar
  } else {
    rec$k_Gap <- NA_integer_
  }
  
  list(metrics = out, recommended_k = rec, scaled_features = X)
}

# 3) (Optional) Fit kmeans at a chosen k after scanning
fit_kmeans_at_k <- function(X_scaled, k, nstart=100, iter.max=200, seed=1L) {
  set.seed(seed)
  kmeans_kl_fast(X_scaled, k = k, nstart = nstart)
}

#=======================================================================================================================#
#                                                 Generate Graphs:
#=======================================================================================================================#
vardistplots                  <- function(dfMaster, cnum){
  inputdf <- dfMaster %>%
    transmute(
      prev_infection_preds = ifelse(prev_infection_preds != 0, prev_infection_preds, NA),
      prev_infection_albs  = ifelse(prev_infection_albs  != 0, prev_infection_albs,  NA),
      midseason_pred_day   = ifelse(midseason_pred_day   != -1, midseason_pred_day,   NA),
      midseason_albs_day   = ifelse(midseason_albs_day   != -1, midseason_albs_day,   NA),
      
      g_InfDur  = 1 / gamma_gentoo,                    
      g_InfMort = (1 / gamma_gentoo) * muI_gentoo,
      
      p_InfDur  = 1 / gamma_predators,
      p_InfMort = (1 / gamma_predators) * muI_predators,
      
      a_InfDur  = 1 / gamma_albatross,
      a_InfMort = (1 / gamma_albatross) * muI_albatross,
      
      B_GentGent = Beta_GentGent,
      B_GentPred = Beta_GentPred,
      B_PredGent = Beta_PredGent,
      B_PredPred = Beta_PredPred,
      B_PredAlbs = Beta_PredAlbs,
      B_AlbsPred = Beta_AlbsPred,
      B_AlbsAlbs = Beta_AlbsAlbs,
      
      R0 = R0,
      clusternum = Cluster
    )
  
  df_long <- inputdf %>%
    mutate(row_id = row_number()) %>%
    pivot_longer(
      cols = -c(row_id, clusternum),
      names_to = "variable",
      values_to = "value"
    ) %>%
    mutate(value = suppressWarnings(as.numeric(value))) %>%
    filter(is.finite(value))
  
  df_red <- df_long %>% filter(clusternum == cnum)
  
  # ---- plot: constant x so each facet's box fills the width ----
  ggplot(df_long, aes(x = 0, y = value)) +
    geom_violin(
      width = 0.9, trim = FALSE,
      fill = "lightblue", color = "darkblue"
    ) +
    geom_jitter(
      data  = df_red,
      aes(x = 0, y = value),
      width = 0.25, height = 0,
      color = "red", size = 1.5, alpha = 0.8
    ) +
    facet_wrap(~ variable, scales = "free_y", ncol = 4) +
    scale_x_continuous(breaks = NULL, labels = NULL, expand = expansion(mult = c(0, 0))) +
    labs(
      title = paste0("Violin Plot of Variable Space — Cluster ", cnum, " in red"),
      x = NULL, y = "Value"
    ) +
    theme_minimal() +
    theme(
      strip.text  = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x= element_blank()
    )
}
outdistgraph                  <- function(valarray, clusternum, cnum, compartments, popvar) {
  dn <- dimnames(valarray)
  time_labels <- dn$time
  
  to_int_time <- function(x) {
    x2 <- gsub("[^0-9\\-]+", "", x)         # strip non-digits
    suppressWarnings(as.integer(x2))
  }
  
  speccomp     <- length(compartments)/3
  yvals        <- c(popvar$total_gentoo * 1.5, popvar$total_preds * 1.5, popvar$total_albatross * 1.5)
  colorsr      <- c("red", "blue", "orange", "darkgreen")
  plots        <- vector("list", length(compartments))
  names(plots) <- compartments
  i = 1
  for (comp in compartments) {
    # slice: [time x sample] for this compartment
    mat_ts <- drop(valarray[comp, , , drop = FALSE])    # [T x N]
    if (nrow(mat_ts) != length(time_labels)) {
      stop("Row count of mat_ts does not match length of time_labels for compartment ", comp)
    }
    
    # keep only the chosen cluster
    keep <- which(clusternum == cnum)
    if (length(keep) == 0) {
      warning("No samples in cluster ", cnum, " for compartment ", comp, ". Skipping.")
      next
    }
    mat_ts_cluster <- mat_ts[, keep, drop = FALSE]      # [T x K]
    
    # build df with correct axes: rows=time, cols=runs
    df <- as.data.frame(mat_ts_cluster, check.names = FALSE)
    rownames(df) <- time_labels
    
    # long format: time column from rownames; runs come from column names
    df_all <- df %>%
      rownames_to_column(var = "time") %>%
      pivot_longer(cols = -time, names_to = "run", values_to = "value") %>%
      mutate(time = to_int_time(time)) %>%
      filter(is.finite(time), !is.na(value))
    
    # average over runs for each time
    df_avg <- df_all %>%
      group_by(time) %>%
      summarise(avg_value = mean(value, na.rm = TRUE), .groups = "drop")
    
    # plot
    p <- ggplot(df_all, aes(x = time, y = value, group = run)) +
      geom_line(color = "gray60", linetype = "dotted") +
      geom_line(data = df_avg, aes(x = time, y = avg_value),
                inherit.aes = FALSE, color = colorsr[(i %% speccomp) + 1], linewidth = 1.2) +
      scale_y_continuous(
        limits = c(-10, yvals[((i-1) %/% speccomp) + 1])
      ) +
      labs(
        title = paste0(comp, " (Cluster ", cnum, ")"),
        x = "Time (days)",
        y = "Number of Individuals"
      ) +
      theme_minimal()
    
    plots[[i]] <- p 
    i = i + 1
  }
  return(plots)
}


validate_clustering           <- function(scan, km_final, chosen_k, features_scaled = NULL, 
                                          file = "k_scan_log.txt"){
  # ---- basic validations ----
  stopifnot(is.list(scan),
            !is.null(scan$metrics),
            !is.null(scan$recommended_k),
            is.list(km_final),
            is.numeric(chosen_k), length(chosen_k) == 1L, chosen_k >= 1L)
  
  if (!is.null(km_final$centers)) {
    stopifnot(nrow(km_final$centers) == as.integer(chosen_k))
  }
  
  # prefer provided features; else use what's inside scan (if present)
  if (is.null(features_scaled) && !is.null(scan$scaled_features)) {
    features_scaled <- scan$scaled_features
  }
  
  N <- if (!is.null(features_scaled)) nrow(features_scaled) else length(km_final$cluster)
  p <- if (!is.null(features_scaled)) ncol(features_scaled) else if (!is.null(km_final$centers)) ncol(km_final$centers) else NA_integer_
  
  # extract metrics at chosen_k (if present)
  met <- scan$metrics
  has_k <- "k" %in% names(met)
  row_k <- if (has_k) which(met$k == as.integer(chosen_k)) else integer(0)
  
  metric_at_k <- if (length(row_k) == 1L) {
    as.list(met[row_k, , drop = FALSE])
  } else list()
  
  # cluster sizes from the provided fit
  cl_sizes <- sort(table(km_final$cluster))
  cluster_df <- data.frame(
    cluster = as.integer(names(cl_sizes)),
    size    = as.integer(cl_sizes)
  )
  
  # open connection and write
  con <- file(file, open = "wt"); on.exit(close(con), add = TRUE)
  w <- function(...) writeLines(paste0(...), con)
  
  w("==== K-SCAN AND FINAL CLUSTERING LOG (from existing fit) ====")
  w(sprintf("Timestamp: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
  if (!is.null(N)) w(sprintf("N (rows): %s", N))
  if (!is.null(p) && !is.na(p)) w(sprintf("p (features): %s", p))
  w("")
  
  # Recommendations
  w("== Recommended k values (by criterion) ==")
  for (nm in names(scan$recommended_k)) {
    w(sprintf("  %-12s : %s", nm, as.character(scan$recommended_k[[nm]])))
  }
  w("")
  
  # Chosen k and summary from the provided fit
  w("== Final choice and provided kmeans fit ==")
  w(sprintf("  chosen_k        : %d", chosen_k))
  # km_final diagnostics (present in base kmeans objects)
  if (!is.null(km_final$tot.withinss)) w(sprintf("  tot.withinss    : %.6f", km_final$tot.withinss))
  if (!is.null(km_final$betweenss))    w(sprintf("  betweenss       : %.6f", km_final$betweenss))
  if (!is.null(km_final$totss))        w(sprintf("  totss           : %.6f", km_final$totss))
  if (!is.null(km_final$iter))         w(sprintf("  iterations used : %s", km_final$iter))
  w("")
  
  # ---- Nicely aligned metrics table ----
  w("== Metrics across k ==")
  tmp <- scan$metrics
  
  # Format all columns to the same width
  col_widths <- vapply(tmp, function(col) {
    max(nchar(c(names(tmp)[1], as.character(col))), na.rm = TRUE)
  }, numeric(1))
  
  # Create a formatted text matrix
  formatted_lines <- apply(tmp, 1, function(row) {
    paste(
      vapply(seq_along(row), function(i)
        format(row[i], width = col_widths[i], justify = "right"), ""),
      collapse = "  "
    )
  })
  
  # Write header and rows
  header <- paste(
    vapply(names(tmp), function(nm)
      format(nm, width = col_widths[nm], justify = "right"), ""),
    collapse = "  "
  )
  
  writeLines(header, con)
  writeLines(formatted_lines, con)
  w("")
  
  # Metrics at chosen_k (if available)
  if (length(metric_at_k)) {
    w("== Metrics at chosen_k ==")
    for (nm in names(metric_at_k)) {
      val <- metric_at_k[[nm]]
      if (is.numeric(val)) {
        w(sprintf("  %-12s : %s", nm, sprintf("%.6f", val)))
      } else {
        w(sprintf("  %-12s : %s", nm, as.character(val)))
      }
    }
    w("")
  }
  
  # Cluster sizes from the provided fit
  w("== Final cluster size distribution (chosen_k) ==")
  utils::write.table(cluster_df, con, sep = "\t", row.names = FALSE, quote = FALSE)
  w("")
  
  invisible(list(file = normalizePath(file), cluster_sizes = cluster_df))
}
#=======================================================================================================================#
#                                                   Main:
#=======================================================================================================================#
compartments <- c(
  "Susceptible_Gentoo",    "Infected_Gentoo",    "Recovered_Gentoo",    "HPAI_Killed_Gentoo",
  "Susceptible_Predators", "Infected_Predators", "Recovered_Predators", "HPAI_Killed_Predators",
  "Susceptible_Albatross", "Infected_Albatross", "Recovered_Albatross", "HPAI_Killed_Albatross"
)

# X: N x B matrix; each row nonnegative and will be row-normalized to sum=1
# k: number of clusters
# eps: small smoothing to avoid log(0)
kmeans_kl      <- function(X, k, max_iter = 100, nstart = 10, eps = 1e-12, seed = 1L) {
  N <- nrow(X); B <- ncol(X)
  # row-normalize to probabilities
  row_norm <- function(M) { s <- rowSums(M); M / pmax(s, eps) }
  Xp <- pmax(X, 0); Xp <- row_norm(Xp)
  
  # KL(P||Q) for many Q: returns length-k vector for one P
  kl_vec <- function(P, Qlist) sapply(Qlist, function(Q) sum(P * (log(pmax(P,eps)) - log(pmax(Q,eps)))))
  
  best <- NULL; set.seed(seed)
  
  pb1  <- progress_bar$new(
    format = "Iterations [:bar] :percent | Elapsed: :elapsed | ETA: :eta",
    total = nstart, clear = FALSE, width = 90
  )
  
  for (start in seq_len(nstart)) {
    # init with k random rows as centroids
    Cidx <- sample.int(N, k)
    Q <- lapply(Cidx, function(i) Xp[i, ])
    
    cl_old <- rep(0L, N)
    for (it in seq_len(max_iter)) {
      # assign
      cl <- apply(Xp, 1, function(P) which.min(kl_vec(P, Q)))
      if (identical(cl, cl_old)) break
      cl_old <- cl
      # update centroids: mean distribution of members (then renormalize)
      Q <- lapply(seq_len(k), function(r) {
        idx <- which(cl == r)
        if (length(idx) == 0) { # empty cluster: re-seed
          Xp[sample.int(N, 1), ]
        } else {
          q <- colMeans(Xp[idx, , drop = FALSE])
          q / sum(q)
        }
      })
    }
    # objective
    obj <- sum(mapply(function(i, r) sum(Xp[i, ] * (log(pmax(Xp[i, ],eps)) - log(pmax(Q[[r]],eps)))), seq_len(N), cl))
    fit <- list(cluster = cl, centers = do.call(rbind, Q), iter = it, objective = obj)
    if (is.null(best) || fit$objective < best$objective) best <- fit
    pb1$tick()
  }
  best
}

kmeans_kl_fast <- function(X, k, max_iter = 100, nstart = 10, batch_size = NULL,
                           eps = 1e-12, tol = 1e-6, seed = 1L) {
  set.seed(seed)
  N <- nrow(X); B <- ncol(X)
  
  # Row-normalize to probabilities
  X <- pmax(X, 0); rs <- rowSums(X); X <- X / pmax(rs, eps)
  
  # Precompute row constants sum(P*logP) (not needed for argmin, but useful for objective)
  row_const <- rowSums(X * log(pmax(X, eps)))
  
  best <- NULL
  for (s in seq_len(nstart)) {
    # --- simple init: pick k distinct rows
    centers_idx <- sample.int(N, k)
    Q <- X[centers_idx, , drop = FALSE]
    Q <- sweep(Q, 1, rowSums(Q), "/")
    
    cl_old <- integer(N)
    it <- 0L
    
    repeat {
      it <- it + 1L
      # ---- ASSIGN (vectorized): choose r maximizing X %*% t(log Q)
      logQ <- log(pmax(Q, eps))                # k x B
      score <- X %*% t(logQ)                   # N x k  (fast BLAS)
      cl <- max.col(score, ties.method = "first")
      
      # ---- check convergence
      if (identical(cl, cl_old) || it > max_iter) break
      cl_old <- cl
      
      # ---- UPDATE centers
      if (is.null(batch_size)) {
        # full batch
        for (r in seq_len(k)) {
          idx <- which(cl == r)
          if (length(idx) == 0L) {             # reseed if empty
            Q[r, ] <- X[sample.int(N, 1L), ]
          } else {
            q <- colMeans(X[idx, , drop = FALSE])
            Q[r, ] <- q / sum(q)
          }
        }
      } else {
        # mini-batch: sample a subset of rows to update
        mb <- sample.int(N, min(batch_size, N))
        for (r in seq_len(k)) {
          idx <- mb[cl[mb] == r]
          if (length(idx) > 0L) {
            q <- colMeans(X[idx, , drop = FALSE])
            Q[r, ] <- q / sum(q)
          }
        }
      }
      
      # ---- early stop on centroid change
      # (cheap proxy: recompute score and see if assignments stable next step)
      if (it >= 2L && !is.null(tol) && tol > 0) {
        logQ2 <- log(pmax(Q, eps))
        score2 <- X %*% t(logQ2)
        cl2 <- max.col(score2, ties.method = "first")
        if (mean(cl2 != cl) < tol) { cl <- cl2; break }
      }
    }
    
    # objective (lower is better)
    logQ <- log(pmax(Q, eps))
    obj <- sum(row_const - score[cbind(seq_len(N), cl)])
    
    fit <- list(cluster = cl, centers = Q, iter = it, objective = obj)
    if (is.null(best) || fit$objective < best$objective) best <- fit
  }
  
  best
}

# 1) Build features from your 12x365x1500 array
feat <- build_histogram_features(
  df)

# 2) Scan k = 2..12 and compute metrics
scan <- evaluate_kmeans_k(
  feature_mat     = feat$features,
  k_min           = 2,
  k_max           = 30,
  nstart          = 100,
  iter.max        = 500,
  scale_features  = TRUE,
  compute_gap     = TRUE,   # set FALSE if you need a faster pass
  gap_B           = 20)

# See the full metrics table
scan$metrics
# See recommended k values per criterion
scan$recommended_k

# 3) Choose a k (e.g., by majority across CH / Silhouette / Gap)
chosen_k <- stats::na.omit(unlist(scan$recommended_k)) |>
  as.integer() |>
  (\(v) as.integer(names(sort(table(v), decreasing = TRUE))[1]))()

# 4) Fit the final k-means at your chosen k
km_final <- fit_kmeans_at_k(scan$scaled_features, k = chosen_k, nstart = 200, iter.max = 300)

# Cluster labels per run:
clusters <- km_final$cluster
#table(clusters)

stopifnot(!is.null(dimnames(df)$compartment),
          !is.null(dimnames(df)$time),
          !is.null(dimnames(df)$sample))

dfMaster$Cluster <- clusters

outputLoc <- paste0(folderLoc, "\\Output\\", format(Sys.time(), "O_%Y-%m-%d_%H-%M-%S"))
dir.create(outputLoc)

validate_clustering(scan, km_final, chosen_k, file = paste0(outputLoc, "\\k_scan_log.txt"))

# 5) Generate graphs
# 5a) Violin plots of input variables, highlighting one cluster at a time
for (cnum in 1:chosen_k) {
  p <- vardistplots(dfMaster, cnum)
  ggsave(paste0(outputLoc, "\\InputDistrib_Cluster", cnum, ".png"), p, width = 14, height = 8)}

# 5b) Time series plots of compartment trajectories, highlighting one cluster at a time
# Define population sizes (for y-axis scaling)
popvar <- list(
  total_gentoo    =  2000,
  total_preds     =   750,
  total_albatross = 50000)

for (cnum in 1:chosen_k) {
  plots <- outdistgraph(df, clusters, cnum, compartments, popvar)
  p     <- plot_grid(plotlist = plots, ncol = length(compartments)/3, align = "hv")
  ggsave(paste0(outputLoc, "\\CompartmentGraphs_Cluster_", cnum,".png"), p, width = 15, height = 10)
}

saveRDS(dfMaster, file = paste0(outputLoc, "\\InputOutputFile.rds"))

#=======================================================================================================================#
#                                                   Clean Up:
#=======================================================================================================================#
# Measure time of execution
end_time <- Sys.time()
execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
hs <- as.integer(execution_time %/% 3600)
ms <- as.integer((execution_time %% 3600) %/% 60)
ss <- execution_time %% 60
print(cat(sprintf("Execution Time: %02d:%02d:%05.2f\n", hs, ms, ss)))

# Clear all environmental variables
rm(list = ls())
