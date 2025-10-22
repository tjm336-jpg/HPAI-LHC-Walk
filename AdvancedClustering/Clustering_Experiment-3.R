library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(cowplot)
library(progress)
library(ggh4x)
library(here)
library(ragg)

inputLoc <- here("AdvancedClustering/Testing")

dfMaster <- readRDS(paste0(inputLoc, "//InputOutputTable.rds"))
valarray <- readRDS(paste0(inputLoc, "//Valarray.rds"))
popvar   <- data.frame(total_gentoo = 2000)
popvar$total_preds <- 750
popvar$total_albatross <- 50000

vardistplots <- function(dfMaster, varout, thresholdvalue, ncol = 4, direction, checkvar) {
  inputdf <- dfMaster %>%
    dplyr::transmute(
      prev_infection_preds = dplyr::na_if(prev_infection_preds, 0),
      prev_infection_albs  = dplyr::na_if(prev_infection_albs,  0),
      midseason_pred_day   = dplyr::na_if(midseason_pred_day,  -1),
      midseason_albs_day   = dplyr::na_if(midseason_albs_day,  -1),
      
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
      
      R0 = R0
    )
  
  df_long <- inputdf %>%
    dplyr::mutate(row_id = dplyr::row_number()) %>%
    tidyr::pivot_longer(
      cols = -row_id,
      names_to = "variable",
      values_to = "value"
    ) %>%
    dplyr::mutate(value = suppressWarnings(as.numeric(value))) %>%
    dplyr::filter(is.finite(value))
  
  df_highlight <- df_long %>% dplyr::filter(row_id %in% varout)
  
  ggplot2::ggplot(df_long, ggplot2::aes(x = 0, y = value)) +
    ggplot2::geom_violin(width = 0.9, trim = FALSE,
                         fill = "lightblue", color = "darkblue") +
    ggplot2::geom_jitter(
      data  = df_highlight,
      ggplot2::aes(x = 0, y = value),
      width = 0.25, height = 0,
      color = "red", size = 1.5, alpha = 0.8
    ) +
    ggplot2::facet_wrap(~ variable, scales = "free_y", ncol = ncol) +
    ggplot2::scale_x_continuous(breaks = NULL, labels = NULL,
                                expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(
      title = paste0("Violin Plot — samples with outlier ", checkvar, " values (", ifelse(direction == "above", ">", 
                                                                                          ifelse(direction== "below", "<",
                                                                                                 stop("Direction Input Invalid.")
                                                                                                 )
                                                                                          ), thresholdvalue, ") in red"),
      x = NULL, y = "Value"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      strip.text  = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x= ggplot2::element_blank()
    )
}
outdistgraph <- function(valarray, varout, checkvar, compartments, popvar) {
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
    if (length(varout) == 0) {
      warning("No samples in variable ", checkvar, " for compartment ", comp, ". Skipping.")
      next
    }
    mat_ts_cluster <- mat_ts[, varout, drop = FALSE]      # [T x K]
    
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
        title = comp,
        x = "Time (days)",
        y = "Number of Individuals"
      ) +
      theme_minimal()
    
    plots[[i]] <- p 
    i = i + 1
  }
  return(plots)
}

variable_dist_check <- function(checkvar, threshold = c("95", "66.8", "99.7"), folderLoc){
  threshold <- match.arg(as.character(threshold), choices = c("95","66.8","99.7"))
  if(!checkvar %in% names(dfMaster)){
    warning(paste(checkvar, "not in the InputOutputFile provided. Check Spelling"))
    stop
  }
  if(threshold == 66.8){
    stdevnum <- 1
  } else if (threshold == 95){
    stdevnum <- 2
  } else if (threshold == 99.7){
    stdevnum <- 3
  } else{
    stop("Invalid Threshold set")
  }
  
  varmean <- mean(dfMaster[[checkvar]])
  varstdv <- sd(dfMaster[[checkvar]])
  hvarout <- which(dfMaster[[checkvar]] > varmean + (stdevnum * varstdv))
  lvarout <- which(dfMaster[[checkvar]] < varmean - (stdevnum * varstdv))
  
  if(length(hvarout != 0)){
    print(paste("Number of outliers above mean in", checkvar, "-", length(hvarout)))
    p <- vardistplots(dfMaster, hvarout , round(varmean + (stdevnum * varstdv)), 
                      direction = "above", checkvar = checkvar)
    ggsave(paste0(folderLoc, "\\InputDistrib_Outlier_",stdevnum,"above_", checkvar, ".png"), p, width = 14, height = 8)
    plots <- outdistgraph(valarray, hvarout, checkvar, compartments, popvar) 
    pgrid <- plot_grid(plotlist = plots, ncol = length(compartments)/3, align = "hv")
    labelMain = paste0(checkvar, " values ", stdevnum, " standard deviations above mean (>", round(varmean + (stdevnum * varstdv)),"),  n=", length(hvarout))
    pfinal <- ggdraw() +
      draw_plot(pgrid, x = 0, y = 0, height = 0.94) +
      draw_label(labelMain, fontface = "bold", size = 16,
                 x = 0.5, y = 0.97, hjust = 0.5)
    ggsave(paste0(folderLoc, "\\CompartmentGraphs_Outlier_",stdevnum,"above_", checkvar,".png"), 
           pfinal, width = 15, height = 10, dpi = 300,
           device = ragg::agg_png, bg = "white")
  }

  if(length(lvarout != 0)){
    print(paste("Number of outliers below mean in", checkvar, "-", length(lvarout)))
    p <- vardistplots(dfMaster, lvarout , round(varmean - (stdevnum * varstdv)), 
                      direction = "below", checkvar = checkvar)
    ggsave(paste0(folderLoc, "\\InputDistrib_Outlier_",stdevnum,"below_", checkvar, ".png"), p, width = 14, height = 8)
    plots <- outdistgraph(valarray, lvarout, checkvar, compartments, popvar) 
    pgrid <- plot_grid(plotlist = plots, ncol = length(compartments)/3, align = "hv")
    labelMain = paste0(checkvar, " values ", stdevnum, " standard deviations below mean (<", round(varmean - (stdevnum * varstdv)),"),  n=", length(lvarout))
    pfinal <- ggdraw() +
      draw_plot(pgrid, x = 0, y = 0, height = 0.94) +
      draw_label(labelMain, fontface = "bold", size = 16,
                 x = 0.5, y = 0.97, hjust = 0.5)
    ggsave(paste0(folderLoc, "\\CompartmentGraphs_Outlier_",stdevnum,"below_", checkvar,".png"), 
           pfinal, width = 15, height = 10, dpi = 300,
           device = ragg::agg_png, bg = "white")
  }
}

checkcompartments <- c("deadGentoo", "deadPredators", "deadAlbatross", 
                       "maxI_Gentoo", "maxI_Predators", "maxI_Albatross",
                       "time_maxI_Gentoo", "time_maxI_Predators", "time_maxI_Albatross",
                       "outbreakdur", "cuminf_Total")
compartments      <- c(
  "Susceptible_Gentoo",    "Infected_Gentoo",    "Recovered_Gentoo",    "HPAI_Killed_Gentoo",
  "Susceptible_Predators", "Infected_Predators", "Recovered_Predators", "HPAI_Killed_Predators",
  "Susceptible_Albatross", "Infected_Albatross", "Recovered_Albatross", "HPAI_Killed_Albatross"
)

fileLoc = here("AdvancedClustering/Output")
folderLoc = paste0(fileLoc, "/", format(Sys.time(), "O_%Y-%m-%d_%H-%M-%S"))
print(folderLoc)
dir.create(folderLoc)


for (i in seq_along(checkcompartments)){
  variable_dist_check(checkcompartments[i], folderLoc = folderLoc)
}