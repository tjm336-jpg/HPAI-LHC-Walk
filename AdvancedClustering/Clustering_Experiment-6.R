rm(list = ls())
start_time <- Sys.time()

#=======================================================================================================================#
#                                                   Libraries:
#=======================================================================================================================#
library(dirichletprocess)
library(dplyr)
library(ggplot2)
library(here)

#=======================================================================================================================#
#                                                Generate Variables:
#=======================================================================================================================#
generate_popvar          <- function(){
  popvar                       <- data.frame(total_gentoo = 2000)                         # Total Gentoo population [Source: Dr. Gamble]
  popvar$ResPer_gentoo         <- 1                                                       # Percent of Gentoo that are resident to the island [Source: 1]
  popvar$Res_gentoo            <- popvar$total_gentoo * popvar$ResPer_gentoo              # Calculate the number of resident Gentoo
  popvar$gentoo_migrate_in     <- popvar$total_gentoo * (1 - popvar$ResPer_gentoo)        # Calculate the number of migratory Gentoo
  popvar$gentoo_chicksurvival  <- 0.51                                                    # Breeding success rate [Source: 7]
  
  popvar$total_preds           <- 750                                                     # Total predator population [Source: Dr. Gamble]
  popvar$ResPer_preds          <- 0.1                                                     # Percent of predator that are resident to the island [Source: 1]
  popvar$Res_preds             <- popvar$total_preds * popvar$ResPer_preds                # Calculate the number of resident predators 
  popvar$preds_migrate_in      <- popvar$total_preds * (1 - popvar$ResPer_preds)          # Calculate the number of migratory predators
  popvar$preds_chicksurvival   <- 0.34                                                    # Breeding success rate [Source: 10]                                                 
  
  popvar$total_albatross       <- 50000                                                   # Total Albatross population [Source: Dr. Gamble]
  popvar$ResPer_albatross      <- 0                                                       # Percent of Albatross that are resident to the island [Source: 1]
  popvar$Res_albatross         <- popvar$total_albatross * popvar$ResPer_albatross        # Resident Albatross 
  popvar$albatross_migrate_in  <- popvar$total_albatross * (1 - popvar$ResPer_albatross)  # 100% of albatross are migratory
  popvar$albs_chicksurvival    <- 0.47                                                    # Breeding success rate [Source: 9]
  
  return(popvar)
}

#=======================================================================================================================#
#                                                   Cluster Functions:
#=======================================================================================================================#
build_histogram_features <- function(
    valarray,                      # 3D: [compartment, time, sample]
    bins = 52,
    popvar = NULL,                 # list/df with total_gentoo/total_predators/total_albatross
    relativize_to_species_totals = TRUE,
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
    total_predators <- read_total(popvar, "total_preds")
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
    features     = feature_mat,
    breaks       = comp_breaks,
    compartments = compartments_use
  )
}


#=======================================================================================================================#
#                                                   Import Data:
#=======================================================================================================================#
folderLoc <- here("AdvancedClustering")
dfraw     <- readRDS(paste0(folderLoc, "\\Testing\\Valarray.rds")) %>% as.data.frame()
dfMaster  <- readRDS(paste0(folderLoc, "\\Testing\\InputOutputTable.rds")) %>% as.data.frame()

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
#                                                     Main:
#=======================================================================================================================#
popvar   <- generate_popvar()
feat     <- build_histogram_features(df, popvar = popvar)
features <- feat$features 
stopifnot(is.matrix(features))

N        <- nrow(features)
keep     <- apply(features, 2, function(v) sd(v, na.rm = TRUE) > 0)
X        <- scale(features[, keep, drop = FALSE])

pca      <- prcomp(X, center = TRUE, scale. = FALSE)
var_exp  <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
k_pc     <- min(which(var_exp >= 0.95)[1],15)
scores   <- pca$x[, seq_len(k_pc), drop = FALSE]

dp <- DirichletProcessMvnormal(scores)
dp <- Fit(dp, its = 10, progressBar = TRUE)
labels <- if (!is.null(dp$clusterLabels)) dp$clusterLabels else dp$labels
labels <- as.integer(labels)
cat("Estimated number of clusters:", length(unique(labels)), "\n")
print(table(labels))
df_plot <- data.frame(PC1 = scores[,1], PC2 = scores[,2], cluster = factor(labels))
ggplot(df_plot, aes(PC1, PC2, color = cluster)) +
  geom_point(alpha = 0.7) + theme_minimal() +
  ggtitle("Dirichlet Process MVN Clustering on PCA Scores")

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
