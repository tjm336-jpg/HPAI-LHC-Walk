timevar                   <- data.frame(t_arrival = 81) # Day at which migratory birds begin arriving                  [Source: 9]
timevar$mig_day_arrival   <- 10                         # Number of days between arrival and 24 days before egg-laying [Source: 9]
timevar$breed_start       <- 115                        # When egg-laying/incubation starts [Source: 7/8]
timevar$hatch_start       <- 153                        # When hatching starts              [Source: 7/8]
timevar$hatch_duration    <- 7                          # Hatching duration                 [Source: 8]
timevar$t_departure       <- 300                        # 30 days post-fledging in April    [Source: 8]
timevar$mig_day_departure <- 31
timevar$horizon           <- 365

popvar                       <- data.frame(total_gentoo = 2000)                         # Total Gentoo population [Source: Dr. Gamble]
popvar$ResPer_gentoo         <- 1                                                       # Percent of Gentoo that are resident to the island [Source: 1]
popvar$Res_gentoo            <- popvar$total_gentoo * popvar$ResPer_gentoo              # Calculate the number of resident Gentoo
popvar$gentoo_migrate_in     <- popvar$total_gentoo * (1 - popvar$ResPer_gentoo)        # Calculate the number of migratory Gentoo
popvar$gentoo_chicksurvival  <- 0.51                                                    # Breeding success rate [Source: 7]

popvar$total_preds           <- 750                                                     # Total predator population [Source: Dr. Gamble]
popvar$ResPer_preds          <- 0.1                                                     # Percent of predator that are resident to the island [Source: 1]
popvar$Res_preds             <- popvar$total_preds * popvar$ResPer_preds                # Calculate the number of resident predators 
popvar$preds_migrate_in      <- popvar$total_preds * (1 - popvar$ResPer_preds)          # Calculate the number of migratory predators
popvar$prev_infection_preds  <- 0
popvar$preds_chicksurvival   <- 0.34                                                    # Breeding success rate [Source: 10]                                                 

popvar$total_albatross       <- 50000                                                   # Total Albatross population [Source: Dr. Gamble]
popvar$ResPer_albatross      <- 0                                                       # Percent of Albatross that are resident to the island [Source: 1]
popvar$Res_albatross         <- popvar$total_albatross * popvar$ResPer_albatross        # Resident Albatross 
popvar$albatross_migrate_in  <- popvar$total_albatross * (1 - popvar$ResPer_albatross)  # 100% of albatross are migratory
popvar$prev_infection_albs   <- 0
popvar$albs_chicksurvival    <- 0.47                                                    # Breeding success rate [Source: 9]

Initial_values <- c(
  S_Gentoo     = popvar$Res_gentoo,     # Susceptible

  S_Predators  = popvar$Res_preds,       # Susceptible

  S_Albatross  = popvar$Res_albatross,   # Susceptible

  timer        = 0)


min_Life_Gentoo       <- 10                                         # [Source: 2]
max_Life_Gentoo       <- 15                                         # [Source: 2]
avg_Life_Gentoo       <- (max_Life_Gentoo + min_Life_Gentoo) / 2    # Calculate average lifespan
avg_Life_Gentoo_Days  <- avg_Life_Gentoo * 365                      # Calculate average lifespan in days
muH_Gentoo            <- 1 / avg_Life_Gentoo_Days                   # Calculate MuH for the Gentoo

min_Life_Preds        <- 0                                          # [Source: NA]
max_Life_Preds        <- 28.8                                       # [Source: 4]
avg_Life_Preds        <- (max_Life_Preds + min_Life_Preds) / 2      # Calculate average lifespan
avg_Life_Preds_Days   <- avg_Life_Preds * 365                       # Calculate average lifespan in days
muH_Preds             <- 1 / avg_Life_Preds_Days                    # Calculate MuH for the Predators

min_Life_Albs         <- 0                                          # [Source: NA]
max_Life_Albs         <- 47                                         # [Source: 3]
avg_Life_Albs         <- (max_Life_Albs + min_Life_Albs) / 2        # Calculate average lifespan
avg_Life_Albs_Days    <- avg_Life_Albs * 365                        # Calculate average lifespan in days
muH_Albs              <- 1 / avg_Life_Albs_Days                     # Calculate MuH for the Albatross

library(deSolve)
library(ggplot2)

validate_births <- function(values){
  sus_cols <- c("Susceptible_Gentoo",
                "Susceptible_Predators",
                "Susceptible_Albatross")
  
  df_long <- as.data.frame(values) %>%
    select(Time, all_of(sus_cols)) %>%
    pivot_longer(-Time, names_to = "compartment", values_to = "value") %>%
    mutate(
      Species = sub("^Susceptible_", "", compartment),
      Species = factor(Species, levels = c("Gentoo", "Predators", "Albatross"))
    )
  
  p <- ggplot(df_long, aes(x = Time, y = value)) +
    geom_line(color = "steelblue") +
    facet_wrap(~ Species, ncol = 1, scales = "free_y") +
    scale_y_continuous(
      limits = c(-1, NA),                 # force min=0, let max be computed per facet
      expand = expansion(mult = c(0, 0.05)) # small headroom above max
    ) +
    labs(
      title = "Migration, Birth, and Death Dynamics without Infection",
      x = "Time (days)",
      y = "Susceptible"
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold"),
      panel.spacing.y = unit(0.8, "lines")
    )
  
  return(p)
}
birthrate_dyn   <- function(time, state, parameters){
  # Extract state
  S_Gentoo    <- round(state["S_Gentoo"]   , 4)
  S_Predators <- round(state["S_Predators"], 4)
  S_Albatross <- round(state["S_Albatross"], 4)
  
  timer       <- state["timer"]
  
  # Extract parameters
  muH_Gentoo  <- parameters["muH_Gentoo"]
  muH_Preds   <- parameters["muH_Preds"]
  muH_Albs    <- parameters["muH_Albs"]
  
  t_departure <- parameters["t_departure"]
  
  dS_Gentoo = 
    - muH_Gentoo * S_Gentoo
  dS_Predators = 
    - muH_Preds * S_Predators
  dS_Albatross = 
    ifelse(time < t_departure, - muH_Albs * S_Albatross, 0)
  dTimer = 1
  
  list(c(dS_Gentoo, dS_Predators, dS_Albatross, dTimer))
}

create_migin_eventdata  <- function(timevar, popvar){
  # Duration (days) for each species' arrival window
  L <- timevar$mig_day_arrival
  if (L <= 0) stop("timevar$mig_day_arrival must be > 0")
  
  # --- Arrival windows (sequential) ---
  # Albatross arrive first
  t_alb_start <- timevar$t_arrival
  t_alb_end   <- t_alb_start + L - 1
  alb_times   <- seq(from = t_alb_start, to = t_alb_end, length.out = L)
  
  # Predators arrive after albatrosses are done, for the same duration
  t_pred_start <- t_alb_end + 1
  t_pred_end   <- t_pred_start + L - 1
  pred_times   <- seq(from = t_pred_start, to = t_pred_end, length.out = L)
  
  # --- Per-day arrivals (constant per day over each window) ---
  pred_S_arrival <- (1 - popvar$prev_infection_preds) * popvar$preds_migrate_in     / L
#  pred_I_arrival <- (    popvar$prev_infection_preds) * popvar$preds_migrate_in     / L
  alb_S_arrival  <- (1 - popvar$prev_infection_albs)  * popvar$albatross_migrate_in / L
#  alb_I_arrival  <- (    popvar$prev_infection_albs)  * popvar$albatross_migrate_in / L
  
  # --- Build eventdata1 (sequential blocks) ---
  df_alb <- data.frame(
    var    = rep(c("S_Albatross"), each = length(alb_times)),
    time   = rep(alb_times, times = 1),
    value  = rep(c(alb_S_arrival), each = length(alb_times)),
    method = "add"
  )
  
  df_pred <- data.frame(
    var    = rep(c("S_Predators"), each = length(pred_times)),
    time   = rep(pred_times, times = 2),
    value  = rep(c(pred_S_arrival), each = length(pred_times)),
    method = "add"
  )
  
  migineventdata      <- rbind(df_alb, df_pred)
  migineventdata$time <- round(migineventdata$time, 3)
  migineventdata      <- migineventdata[order(migineventdata$time), ]
  
  return(migineventdata)
}
create_hatch_eventdata  <- function(timevar, popvar, N_Gent, N_Pred, N_Albs){
  # Duration (days) for each species' arrival window
  L <- timevar$hatch_duration
  if (L <= 0) stop("timevar$hatch_duration must be > 0")
  
  # --- Hatching window ---
  t_hatch_start <- timevar$hatch_start
  t_hatch_end   <- t_hatch_start + L - 1
  hatch_times   <- seq(from = t_hatch_start, to = t_hatch_end, length.out = L)

  # --- Per-day hatches (constant per day over each window) ---
  gent_S_hatch <- (N_Gent / 2) * popvar$gentoo_chicksurvival / L
  pred_S_hatch <- (N_Pred / 2) * popvar$preds_chicksurvival  / L
  albs_S_hatch <- (N_Albs / 2) * popvar$albs_chicksurvival   / L
  
  # --- Build hatcheventdata ---
  hatcheventdata <- data.frame(
    var    = rep(c("S_Gentoo","S_Predators", "S_Albatross"),  each = length(hatch_times)),
    time   = rep(hatch_times, times = 3),
    value  = rep(c(gent_S_hatch, pred_S_hatch, albs_S_hatch), each = length(hatch_times)),
    method = "add"
  )
  
  hatcheventdata$time <- round(hatcheventdata$time, 3)
  hatcheventdata      <- hatcheventdata[order(hatcheventdata$time), ]
  
  return(hatcheventdata)
}
create_migout_eventdata <- function(timevar, popvar, statedeparturevector){
  # Duration (days) for each species' arrival window
  L  <- timevar$mig_day_departure
  if (L <= 0) stop("timevar$t_departure must be > 0")
  
  # --- Hatching window ---
  t_migout_start <- timevar$t_departure
  t_migout_end   <- t_migout_start + L - 1
  migout_times   <- seq(from = t_migout_start, to = t_migout_end, length.out = L)
  
  S_Gentoo    <- statedeparturevector["S_Gentoo"]
  S_Predators <- statedeparturevector["S_Predators"]
  S_Albatross <- statedeparturevector["S_Albatross"]

  S_Gentoo_Migout    <- S_Gentoo    * (1 - popvar$ResPer_gentoo)    / L
  S_Predators_Migout <- S_Predators * (1 - popvar$ResPer_preds)     / L
  S_Albatross_Migout <- S_Albatross * (1 - popvar$ResPer_albatross) / L

  migout_eventdata <- data.frame(
    var    = rep(c("S_Gentoo","S_Predators", "S_Albatross"),  each = length(migout_times)),
    time   = rep(migout_times, times = 3),
    value  = rep(c(- S_Gentoo_Migout, - S_Predators_Migout, - S_Albatross_Migout), each = length(migout_times)),
    method = "add"
  )
  
  migout_eventdata$time <- round(migout_eventdata$time, 3)
  migout_eventdata      <- migout_eventdata[order(migout_eventdata$time), ]
  
  return(migout_eventdata)
}

run_Model     <- function(timevar, popvar, Initial_values, parameters){
  time1 <- seq(from = 0, to = timevar$hatch_start, by = 1)
  time2 <- seq(from = timevar$hatch_start, to = timevar$t_departure, by = 1)
  time3 <- seq(from = timevar$t_departure, to = timevar$horizon, by = 1)
  
  migin_eventdata <- create_migin_eventdata(timevar, popvar)
  print(Initial_values)
  print(parameters)
  print(time1)
  print(migin_eventdata)
  
  migin_eventdata <- create_migin_eventdata(timevar, popvar)
  solved1 <- ode(
    y      = Initial_values,
    func   = birthrate_dyn,
    parms  = parameters,
    times  = time1,
    events = list(data = migin_eventdata),
    method = "lsoda"
  )
  
  solved1 <- as.data.frame(solved1)
  statedeparturevector<- as.numeric(solved1[ nrow(solved1), -1])
  names(statedeparturevector) <- names(Initial_values)
  N_Gent <- statedeparturevector["S_Gentoo"]
  N_Pred <- statedeparturevector["S_Predators"]
  N_Albs <- statedeparturevector["S_Albatross"]
  
  hatcheventdata <- create_hatch_eventdata(timevar, popvar, N_Gent, N_Pred, N_Albs)

  solved2 <- ode(
    y      = statedeparturevector,
    func   = birthrate_dyn,
    parms  = parameters,
    times  = time2,
    events = list(data = hatcheventdata), 
    method = "lsoda"
  )
  solved2 <- as.data.frame(solved2)
  
  statedeparturevector<- as.numeric(solved2[ nrow(solved2), -1])
  names(statedeparturevector) <- names(Initial_values)
  N_Gent <- statedeparturevector["S_Gentoo"]
  N_Pred <- statedeparturevector["S_Predators"]
  N_Albs <- statedeparturevector["S_Albatross"]
  
  migout_eventdata <- create_migout_eventdata(timevar, popvar, statedeparturevector)

  solved3 <- ode(
    y      = statedeparturevector,
    func   = birthrate_dyn,
    parms  = parameters,
    times  = time3,
    events = list(data = migout_eventdata),
    method = "lsoda"
  )
  solved3 <- as.data.frame(solved3)
  
  out <- NULL
  out <- rbind(solved1, solved2, solved3)
  
  names(out) <- c(
    "Time",
    "Susceptible_Gentoo",
    "Susceptible_Predators",
    "Susceptible_Albatross",
    "Timer"
  )
  
  out <- out[!duplicated(out$Time),]
  
  return(out)
}

parameters <- c(
  muH_Gentoo     = muH_Gentoo,
  muH_Preds      = muH_Preds,
  muH_Albs       = muH_Albs,
  
  t_departure    = timevar$t_departure
)

values <- run_Model(timevar, popvar, Initial_values, parameters)
p <- validate_births(values)
print(p)

#=======================================================================================================================#
#                                                         Sources:
#=======================================================================================================================#
# [10] Catry, Paulo, Ana Almeida, Miguel Lecoq, José Pedro Granadeiro, and Rafael Matias. “Low Breeding Success and Sharp
#         Population Decline at the Largest Known Falkland Skua Colony.” Polar Biology 34, no. 8 (2011): 1239–41. 
#         https://doi.org/10.1007/s00300-011-0978-1.
