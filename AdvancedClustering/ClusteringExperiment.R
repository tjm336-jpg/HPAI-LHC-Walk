start_time           <- Sys.time()
#=======================================================================================================================#
#                                               Libraries:
#=======================================================================================================================#
library(dplyr)
library(cluster)
library(tidyr)
library(tibble)
library(here)

#=======================================================================================================================#
#                                               Import Data:
#=======================================================================================================================#
folderLoc <- here("AdvancedClustering")
dfraw <- readRDS(paste0(folderLoc, "\\Testing\\InputOutputTable.rds")) %>% as.data.frame()

#=======================================================================================================================#
#                                             Data Cleaning:
#=======================================================================================================================#
drop_patterns <- c(
  "^gamma_", "^sigma_",
  "^muI_", "^muH_",
  "^Beta_",
  "^prev_infection", "midseason_",
  "^t_", "^mig_day_", "^breed_start", "^hatch_start", 
  "^R0" #, "run_id" #, "Cluster"
)
print(names(dfraw))
drop_regex <- paste0("(", paste(drop_patterns, collapse = "|") , ")")
df <- dfraw %>% select(-matches(drop_regex))

df$terI_Gentoo    <- round(df$terI_Gentoo)
df$terI_Predators <- round(df$terI_Predators)
df$terI_Gentoo    <- ifelse(df$terI_Gentoo    < 0, 0, df$terI_Gentoo)
df$terI_Predators <- ifelse(df$terI_Predators < 0, 0, df$terI_Predators)


#=======================================================================================================================#
#                                             Histogram Variables:
#=======================================================================================================================#
# Derive a global time range across all chosen peak lists
peak_cols <- c("I_G_Peak", "I_P_Peak", "I_A_Peak")
all_peaks <- unlist(df[peak_cols], recursive = TRUE, use.names = FALSE)
T_max <- suppressWarnings(max(as.numeric(all_peaks), na.rm = TRUE))
if (!is.finite(T_max)) T_max <- 365  # fallback horizon if no peaks

# Weekly bins 
bin_width <- 7
breaks <- seq(0, ceiling(365/bin_width)*bin_width, by = bin_width)

#=======================================================================================================================#
#                                             Feature Extraction:
#=======================================================================================================================#
peak_hist <- function(times, amps, breaks, normalize = TRUE, prefix = "h") {
  times <- unlist(times)
  amps  <- unlist(amps)
  
  # Count peaks per bin
  counts <- hist(times, breaks = breaks, plot = FALSE)$counts
  
  # Weighted sum of amplitudes per bin
  sums   <- hist(times, breaks = breaks, weights = amps, plot = FALSE)$counts
  
  if (normalize && sum(counts) > 0) {
    counts <- counts / sum(counts)  # probability distribution over bins
  }
  if (normalize && sum(sums) > 0) {
    sums <- sums / sum(sums)        # normalized amplitude distribution
  }
  
  names(counts) <- paste0(prefix, "_count_bin", seq_along(counts))
  names(sums)   <- paste0(prefix, "_amp_bin",   seq_along(sums))
  
  c(counts, sums)
}
peak_summaries <- function(times, amps, prefix = "s") {
  times <- sort(unlist(times))
  amps  <- unlist(amps)
  
  if (length(times) == 0) {
    out <- c(count = 0, first = NA, last = NA, mean_time = NA,
             sd_time = NA, iqr_time = NA, median_iei = NA,
             mean_amp = NA, sd_amp = NA, max_amp = NA)
  } else {
    iei <- diff(times)
    out <- c(
      count     = length(times),
      first     = times[1],
      last      = times[length(times)],
      mean_time = mean(times),
      sd_time   = stats::sd(times),
      iqr_time  = IQR(times),
      median_iei= if (length(iei)) stats::median(iei) else NA_real_,
      mean_amp  = mean(amps, na.rm = TRUE),
      sd_amp    = stats::sd(amps, na.rm = TRUE),
      max_amp   = max(amps, na.rm = TRUE)
    )
  }
  
  names(out) <- paste0(prefix, "_", names(out))
  out
}

# Extract features for each run
features <- df %>%
  select(run_id, I_G_Peak, I_G_PeakValue,
         I_P_Peak, I_P_PeakValue,
         I_A_Peak, I_A_PeakValue) %>%
  rowwise() %>%
  mutate(
    hist_I_G = list(peak_hist(I_G_Peak, I_G_PeakValue, breaks, normalize=TRUE, prefix="IG")),
    hist_I_P = list(peak_hist(I_P_Peak, I_P_PeakValue, breaks, normalize=TRUE, prefix="IP")),
    hist_I_A = list(peak_hist(I_A_Peak, I_A_PeakValue, breaks, normalize=TRUE, prefix="IA")),
    
    sum_I_G  = list(peak_summaries(I_G_Peak, I_G_PeakValue, prefix="IGs")),
    sum_I_P  = list(peak_summaries(I_P_Peak, I_P_PeakValue, prefix="IPs")),
    sum_I_A  = list(peak_summaries(I_A_Peak, I_A_PeakValue, prefix="IAs"))
  ) %>%
  ungroup() %>%
  unnest_wider(hist_I_G) %>%
  unnest_wider(hist_I_P) %>%
  unnest_wider(hist_I_A) %>%
  unnest_wider(sum_I_G)  %>%
  unnest_wider(sum_I_P)  %>%
  unnest_wider(sum_I_A)

# Find clusters on the features
cluster_features <- function(features,
                             id_col = "run_id",
                             k_min = 2, k_max = 30,
                             method = c("pam", "kmeans"),
                             scale_features = TRUE,
                             impute = c("median", "drop")) {
  method <- match.arg(method)
  impute <- match.arg(impute)
  
  stopifnot(id_col %in% names(features))
  # keep ID separate
  ids <- features[[id_col]]
  
  # Keep only numeric feature columns
  X <- features %>% select(where(is.numeric))
  # If run_id is numeric, drop it from the feature set
  if (id_col %in% names(X)) X[[id_col]] <- NULL
  
  # Drop all-NA and ~zero-variance columns
  all_na <- vapply(X, function(x) all(is.na(x)), logical(1))
  nzvar  <- vapply(X, function(x) {
    x2 <- x[is.finite(x)]
    if (length(x2) < 2) return(FALSE)
    var(x2) > 1e-12
  }, logical(1))
  X <- X[, !(all_na | !nzvar), drop = FALSE]
  if (ncol(X) < 2) stop("Not enough usable numeric features to cluster.")
  
  # Simple NA handling
  if (impute == "drop") {
    ok <- stats::complete.cases(X)
    X  <- X[ok, , drop = FALSE]
    ids_used <- ids[ok]
  } else {  # median imputation
    for (j in seq_along(X)) {
      x <- X[[j]]
      if (anyNA(x)) {
        med <- stats::median(x, na.rm = TRUE)
        x[is.na(x)] <- med
        X[[j]] <- x
      }
    }
    ids_used <- ids
  }
  
  # Scale features (recommended for Euclidean / PAM)
  Xs <- if (scale_features) scale(X) else as.matrix(X)
  
  # Choose k by average silhouette
  ks <- k_min:k_max
  sil_scores <- numeric(length(ks))
  models <- vector("list", length(ks))
  
  if (method == "pam") {
    # PAM uses a dissimilarity; Euclidean on scaled features is fine here.
    D <- dist(Xs)  # numeric-only, Euclidean
    for (i in seq_along(ks)) {
      k <- ks[i]
      fit <- pam(D, k = k, diss = TRUE)
      si  <- silhouette(fit$clustering, D)
      sil_scores[i] <- mean(si[, "sil_width"])
      models[[i]] <- list(fit = fit, sil = si)
    }
  } else { # kmeans
    for (i in seq_along(ks)) {
      k <- ks[i]
      fit <- kmeans(Xs, centers = k, nstart = 50, iter.max = 1000)
      D   <- dist(Xs)
      si  <- silhouette(fit$cluster, D)
      sil_scores[i] <- mean(si[, "sil_width"])
      models[[i]] <- list(fit = fit, sil = si)
    }
  }
  
  best_idx <- which.max(sil_scores)
  best_k   <- ks[best_idx]
  best_fit <- models[[best_idx]]$fit
  best_sil <- models[[best_idx]]$sil
  
  # Assemble outputs
  if (method == "pam") {
    clusters <- factor(best_fit$clustering)
    centers  <- best_fit$medoids
    centers_note <- "medoids are row indices in the distance matrix"
  } else {
    clusters <- factor(best_fit$cluster)
    # unscale kmeans centers back to original units
    centers <- sweep(best_fit$centers, 2, attr(Xs, "scaled:scale"), `*`)
    centers <- sweep(centers,        2, attr(Xs, "scaled:center"), `+`)
    centers_note <- "centers are in original feature units"
  }
  
  result <- list(
    ids_used         = ids_used,
    feature_names    = colnames(X),
    X_scaled         = Xs,
    k_grid           = ks,
    silhouette_table = data.frame(k = ks, avg_silhouette = sil_scores),
    best_k           = best_k,
    clusters         = clusters,
    model            = best_fit,
    silhouette       = best_sil,
    centers          = centers,
    centers_note     = centers_note,
    method           = method
  )
  
  # tidy data frame with cluster labels for your runs
  clustered_df <- tibble(!!id_col := ids_used, cluster = clusters, row.names = NULL)
  
  list(result = result, clustered_df = clustered_df)
}

out <- cluster_features(features,
                        id_col = "run_id",
                        k_min = 80, k_max = 180,
                        method = "pam",        # or "kmeans"
                        scale_features = TRUE,
                        impute = "median")     # or "drop"

# Inspect chosen k and silhouette curve
print(out$result$best_k)
print(out$result$silhouette_table)

# Cluster assignments per run_id
head(out$clustered_df)

# Medoids (PAM) or centers (kmeans)
out$result$centers_note
out$result$centers

# Save results if you like
saveRDS(out, file = paste0(folderLoc, "\\peak_features_clustering_results_80-1.rds"))

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

