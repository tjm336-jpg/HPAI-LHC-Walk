start_time           <- Sys.time()

library(here)
library(ggplot2)
library(GGally)
library(dplyr)
library(rlang)
library(ragg)

inputLoc <- here("AdvancedClustering/Testing")

dfMaster <- readRDS(paste0(inputLoc, "//InputOutputTable.rds"))

checkvar <- c("maxI_Gentoo",	        "time_maxI_Gentoo",	
              "maxI_Predators",	      "time_maxI_Predators",
              "maxI_Albatross",	      "time_maxI_Albatross",	
              "outbreakdur",	        
              "infperiod_Gentoo",     "infperiod_Predators",  "infperiod_Albatross",
              "deadGentoo",	          "deadPredators",        "deadAlbatross",	
              "gentpred_transevents",	"predgent_transevents",	"predalbs_transevents", "albspred_transevents",	
              "cuminf_Gentoo",	      "cuminf_Predators",	    "cuminf_Albatross")

missing_vars <- setdiff(checkvar, names(dfMaster))
if(length(missing_vars) > 0){
  stop(
    paste(
      "Missing the following variables:
      ",
      paste(missing_vars, collapse = ","), "
       ", length(missing_vars), "/", length(checkvar)
    )
  )
}

print("All Variables Present")

df_use <- dfMaster %>%
  select(all_of(checkvar)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(across(everything(), ~ ifelse(is.finite(.), ., NA_real_)))

save_pairwise_from_df <- function(df_use,
                                  checkvar,
                                  out_dir = "pairwise_plots",
                                  width = 6, height = 5, dpi = 300,
                                  prefix = "") {
  # ---- validate variables ----
  missing_vars <- setdiff(checkvar, names(df_use))
  if (length(missing_vars) > 0) {
    stop("The following variable(s) are not in df_use: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }
  
  # ---- prep output dir ----
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ---- helper: sanitize filenames ----
  sanitize <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)
  
  # ---- generate all unique pairs (i < j) ----
  pairs <- t(combn(checkvar, 2))
  
  # ---- iterate and save ----
  apply(pairs, 1, function(v) {
    x <- v[1]; y <- v[2]
    
    # subset & coerce to numeric safely
    dat <- df_use[, c(x, y)]
    dat[[x]] <- suppressWarnings(as.numeric(dat[[x]]))
    dat[[y]] <- suppressWarnings(as.numeric(dat[[y]]))
    dat <- dat[is.finite(dat[[x]]) & is.finite(dat[[y]]), , drop = FALSE]
    
    # compute R² (needs at least 3 rows)
    r2 <- NA_real_
    if (nrow(dat) >= 3) {
      r2 <- summary(lm(reformulate(x, y), data = dat))$r.squared
    }
    
    # build plot
    p <- ggplot2::ggplot(dat, ggplot2::aes_string(x = x, y = y)) +
      ggplot2::geom_point(alpha = 0.6, size = 1.2) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
      ggplot2::annotate("text", x = Inf, y = Inf,
                        label = sprintf("R² = %s",
                                        ifelse(is.finite(r2), sprintf("%.3f", r2), "NA")),
                        hjust = 1.05, vjust = 1.3, size = 4) +
      ggplot2::labs(title = paste(y, "vs", x), x = x, y = y) +
      ggplot2::theme_minimal()
    
    # filename
    fname <- file.path(out_dir, paste0(
      if (nzchar(prefix)) paste0(sanitize(prefix), "_") else "",
      sanitize(y), "_vs_", sanitize(x), ".png"
    ))
    
    # save (use ragg if available to avoid rendering quirks)
    if (requireNamespace("ragg", quietly = TRUE)) {
      ggplot2::ggsave(filename = fname, plot = p,
                      width = width, height = height, dpi = dpi,
                      device = ragg::agg_png, bg = "white")
    } else {
      ggplot2::ggsave(filename = fname, plot = p,
                      width = width, height = height, dpi = dpi,
                      bg = "white")
    }
    
    invisible(fname)
  })
}

fileLoc = here("AdvancedClustering/Output")
folderLoc = paste0(fileLoc, "/", format(Sys.time(), "R2_%Y-%m-%d_%H-%M-%S"))
dir.create(folderLoc)

save_pairwise_from_df(df_use,
                      checkvar,
                      out_dir = folderLoc)

# Measure time of execution
end_time <- Sys.time()
execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
hs <- as.integer(execution_time %/% 3600)
ms <- as.integer((execution_time %% 3600) %/% 60)
ss <- execution_time %% 60
print(cat(sprintf("Execution Time: %02d:%02d:%05.2f\n", hs, ms, ss)))

# Clear all environmental variables
rm(list = ls())