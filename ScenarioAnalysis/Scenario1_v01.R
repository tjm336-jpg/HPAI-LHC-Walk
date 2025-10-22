start_time           <- Sys.time()
#=======================================================================================================================#
#                                               Control Variables:
#=======================================================================================================================#
# Control Variable
scenario <- 4

#=======================================================================================================================#
#                                                   Libraries:
#=======================================================================================================================#
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(cowplot)
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
#                                                 Generate Graphs:
#=======================================================================================================================#
vardistplots_kl               <- function(dfMaster, cnum, bins = 30, eps = 1e-8, ncol = 4) {
  inputdf <- dfMaster %>%
    transmute(
      prev_infection_preds = na_if(prev_infection_preds, 0),
      prev_infection_albs  = na_if(prev_infection_albs,  0),
      midseason_pred_day   = na_if(midseason_pred_day,   -1),
      midseason_albs_day   = na_if(midseason_albs_day,   -1),
      
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
    mutate(row_id = dplyr::row_number()) %>%
    pivot_longer(
      cols = -c(row_id, clusternum),
      names_to = "variable", values_to = "value"
    ) %>%
    mutate(value = suppressWarnings(as.numeric(value))) %>%
    filter(is.finite(value))
  
  df_red <- df_long %>% filter(clusternum == cnum)
  
  # ---- KL helpers (no .data pronoun used) ----
  kl_div <- function(p, q, eps = 1e-8) {
    p <- p / pmax(sum(p), eps); q <- q / pmax(sum(q), eps)
    p <- pmax(p, eps); q <- pmax(q, eps)
    sum(p * log(p / q))
  }
  make_breaks <- function(x, bins) {
    r <- range(x, na.rm = TRUE)
    if (!all(is.finite(r)) || r[1] == r[2]) r <- r + c(-0.5, 0.5)
    seq(r[1], r[2], length.out = bins + 1)
  }
  
  vars <- unique(df_long$variable)
  kl_list <- lapply(vars, function(vn) {
    v_all <- df_long$value[df_long$variable == vn]
    v_clu <- df_red$value[df_red$variable == vn]
    br    <- make_breaks(v_all, bins)
    P     <- hist(v_clu, breaks = br, plot = FALSE)$counts
    Q     <- hist(v_all, breaks = br, plot = FALSE)$counts
    KL    <- kl_div(P, Q, eps)
    tibble(variable = vn, KL = KL,
           KL_label = sprintf("%s  (KL=%.3f)", vn, KL))
  })
  kl_tbl <- bind_rows(kl_list)
  lab_map <- setNames(kl_tbl$KL_label, kl_tbl$variable)
  
  p <- ggplot(df_long, aes(x = 0, y = value)) +
    geom_violin(width = 0.9, trim = FALSE, fill = "lightblue", color = "darkblue") +
    geom_jitter(data = df_red, aes(x = 0, y = value),
                width = 0.25, height = 0, color = "red", size = 1.5, alpha = 0.8) +
    facet_wrap(~ variable, scales = "free_y", ncol = ncol,
               labeller = as_labeller(lab_map)) +
    scale_x_continuous(breaks = NULL, labels = NULL, expand = expansion(mult = c(0, 0))) +
    labs(
      title = paste0("Inputs vs. Cluster ", cnum, " — KL(P_cluster || Q_all) per variable"),
      x = NULL, y = "Value"
    ) +
    theme_minimal() +
    theme(
      strip.text  = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x= element_blank()
    )
  
  list(plot = p, kl_table = select(kl_tbl, variable, KL))
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

#=======================================================================================================================#
#                                                 Model Grouping:
#=======================================================================================================================#
#=======================================================================================================================#
#                                                     Main:
#=======================================================================================================================#
if (scenario == 1){
  inputLoc <- here("ScenarioAnalysis/Testing/Prev_Preds")
} else if (scenario == 2){
  inputLoc <- here("ScenarioAnalysis/Testing/Prev_Albs")
} else if (scenario == 3){
  inputLoc <- here("ScenarioAnalysis/Testing/Mids_Preds")
} else if (scenario == 4){
  inputLoc <- here("ScenarioAnalysis/Testing/Mids_Albs")
} 

Valarray   <- readRDS(paste0(inputLoc, "/Valarray.rds"))
dfMaster   <- readRDS(paste0(inputLoc, "/InputOutputTable.rds"))

compartments <- c(
  "Susceptible_Gentoo",    "Infected_Gentoo",    "Recovered_Gentoo",    "HPAI_Killed_Gentoo",
  "Susceptible_Predators", "Infected_Predators", "Recovered_Predators", "HPAI_Killed_Predators",
  "Susceptible_Albatross", "Infected_Albatross", "Recovered_Albatross", "HPAI_Killed_Albatross"
)
popvar <- generate_popvar()

dfMaster$Cluster <- ifelse(dfMaster$overwintering == 1, 1, 0)

outputLoc <- paste0(here("ScenarioAnalysis/Output"),"/",scenario, "-", format(Sys.time(), "O_%Y-%m-%d_%H-%M-%S"))
dir.create(outputLoc)

res   <- vardistplots_kl(dfMaster, 1)
p     <- res$plot
ggsave(paste0(outputLoc, "/InputDistrib_Cluster_Overwinter.png"), p, width = 14, height = 8)
plots <- outdistgraph(Valarray, dfMaster$Cluster, 1, compartments, popvar)
p     <- plot_grid(plotlist = plots, ncol = length(compartments)/3, align = "hv")
ggsave(paste0(outputLoc, "/CompartmentGraphs_Overwinter.png"), p, width = 15, height = 10)  

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
