start_time           <- Sys.time()

library(here)
library(ggplot2)
library(GGally)
library(dplyr)
library(rlang)
library(ragg)
library(progress)

inputLoc   <- here("AdvancedClustering/Testing")
dfMaster   <- readRDS(paste0(inputLoc, "//InputOutputTable.rds"))
clustervar <- c("maxI_Gentoo",	        "time_maxI_Gentoo",	
              "maxI_Predators",	      "time_maxI_Predators",
              "maxI_Albatross",	      "time_maxI_Albatross",	
              "infperiod_Gentoo",     "infperiod_Predators",  "infperiod_Albatross",
              "deadGentoo",	          "deadPredators",        "deadAlbatross",	
              "gentpred_transevents",	"predgent_transevents",	"predalbs_transevents", "albspred_transevents")

missing_vars <- setdiff(clustervar, names(dfMaster))
if(length(missing_vars) > 0){
  stop(
    paste(
      "Missing the following variables:
      ",
      paste(missing_vars, collapse = ","), "
       ", length(missing_vars), "/", length(clustervar)
    )
  )
}

print("All Variables Present")

df_use <- dfMaster %>%
  select(all_of(clustervar)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(across(everything(), ~ ifelse(is.finite(.), ., NA_real_)))

evaluate_kmeans_k <- function(feature_mat,
                              k_min = 2, k_max = 20,
                              nstart = 5000, iter.max = 10000,
                              compute_gap = TRUE,  # set FALSE if speed is critical
                              gap_B = 20,          # number of reference datasets for Gap statistic
                              seed = 1L){
  stopifnot(nrow(feature_mat) >= k_max)
  set.seed(seed)
  
  X <- feature_mat
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
    km <- kmeans(X, centers=k, nstart=25, iter.max=100)
    W <- km$tot.withinss   # total within-cluster SS (equivalent to sum WSS over clusters)
    
    # uniform reference in bounding box of X
    mins <- apply(X, 2, min); maxs <- apply(X, 2, max)
    lW <- numeric(B)
    for (b in seq_len(B)) {
      U <- matrix(runif(N * p, mins, maxs), nrow=N, ncol=p)
      kmb <- kmeans(U, centers=k, nstart=10, iter.max=100)
      lW[b] <- log(kmb$tot.withinss)
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
    km <- kmeans(X, centers = k, nstart = nstart, iter.max = iter.max)
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
  
  list(metrics = out, recommended_k = rec, features = X)
}

scan <- evaluate_kmeans_k(df_use)
print(scan$recommended_k)

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