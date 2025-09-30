#=======================================================================================================================#
#                                               Compatability Notes:
#=======================================================================================================================#
# Compatible with LHC-Walk-v02.R

#=======================================================================================================================#
#                                                   Libraries:
#=======================================================================================================================#
library(ggplot2)
library(patchwork)
library(tidyverse)

#=======================================================================================================================#
#                                                   Load data:
#=======================================================================================================================#
fileLoc   = "C:\\Users\\tjm336\\OneDrive - Cornell University\\01 Lab Work\\01 Rotations\\01 Bento-Gamble\\04-TriSpecies-HPAI\\LHC Walking"
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

overwintering               <- rawdf$overwintering

# 2. Pivot to long + attach the overwintering flag
df_long <- inputdf %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(-row_id, names_to = "variable", values_to = "value") %>%
  left_join(
    rawdf %>% mutate(row_id = row_number()) %>% select(row_id, overwintering),
    by = "row_id"
  )

# 3. Plot boxplots per variable, with overwinterers as red dots
ggplot(df_long, aes(x = "", y = value)) +
  geom_boxplot(fill = "lightblue",
               color = "darkblue",
               outlier.shape = NA) +             # hide default outliers
  geom_jitter(
    data     = filter(df_long, overwintering == 1),
    aes(x = "", y = value),
    color    = "red",
    size     = 1.5,
    width    = 0.2,
    alpha    = 0.8
  ) +
  facet_wrap(~ variable, scales = "free_y", ncol = 3) +
  labs(
    title = "Box & Whisker Plots of Variable Space\nRed points = overwintering events",
    x     = NULL,
    y     = "Value"
  ) +
  theme_minimal() +
  theme(
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank()
  )