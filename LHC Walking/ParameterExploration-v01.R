#=======================================================================================================================#
#                                               Compatability Notes:
#=======================================================================================================================#
# Compatible with LHC-Walk-v02.R

#=======================================================================================================================#
#                                                   Libraries:
#=======================================================================================================================#
library(ggplot2)
library(tidyverse)

#=======================================================================================================================#
#                                                   Load data:
#=======================================================================================================================#
fileLoc   = "C:\\Users\\tjmcg\\OneDrive - Cornell University\\01 Lab Work\\01 Rotations\\01 Bento-Gamble\\04-TriSpecies-HPAI\\LHC Walking"
inputFile = "Foundation.csv"
rawdf     <- as.data.frame(read.csv(paste0(fileLoc, "\\", inputFile)))

#=======================================================================================================================#
#                                                   Isolate Inputs:
#=======================================================================================================================#
inputdf                     <- data.frame(prev_infection_preds = rawdf$prev_infection_preds)
inputdf$prev_infection_albs <- rawdf$prev_infection_albs
inputdf$gentoo_inf_dur      <- 1 / rawdf$gamma_genetoo
inputdf$gentoo_inf_mort     <- inputdf$gentoo_inf_dur * rawdf$muI_gentoo
inputdf$preds_inf_dur       <- 1 / rawdf$gamma_predators
inputdf$preds_inf_mort      <- inputdf$preds_inf_dur * rawdf$muI_predators
inputdf$albs_inf_dur        <- 1 / rawdf$gamma_albatross
inputdf$albs_inf_mort       <- inputdf$albs_inf_dur * rawdf$muI_albatross
inputdf$Beta_gentgent       <- rawdf$Beta_GentGent
inputdf$Beta_gentpred       <- rawdf$Beta_GentPred
inputdf$Beta_predgent       <- rawdf$Beta_PredGent
inputdf$Beta_predpred       <- rawdf$Beta_PredPred
inputdf$Beta_predalbs       <- rawdf$Beta_PredAlbs
inputdf$Beta_albspred       <- rawdf$Beta_AlbsPred
inputdf$Beta_albsalbs       <- rawdf$Beta_AlbsAlbs

# Reshape inputdf to long format
df_long <- inputdf %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value")

# Plot: one boxplot per variable, faceted with free Y-axis
ggplot(df_long, aes(x = variable, y = value)) +
  geom_boxplot(fill = "lightblue", color = "darkblue", outlier.color = "red", outlier.shape = 16) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +
  labs(title = "Box and Whisker Plots of Input Variables",
       x = NULL, y = "Value") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )