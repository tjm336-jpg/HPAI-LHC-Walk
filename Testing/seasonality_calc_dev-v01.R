timevar                   <- data.frame(t_arrival = 81) # Day at which migratory birds begin arriving                  [Source: 9]
timevar$mig_day_arrival   <- 10                         # Number of days between arrival and 24 days before egg-laying [Source: 9]
timevar$breed_start       <- 115                        # When egg-laying/incubation starts [Source: 7/8]
timevar$hatch_start       <- 153                        # When hatching starts              [Source: 7/8]
timevar$hatch_duration    <- 7                          # Hatching duration                 [Source: 8]
timevar$t_departure       <- 300                        # 30 days post-fledging in April    [Source: 8]
timevar$mig_day_departure <- 31
timevar$horizon           <- 365

time <- seq(from = 0, to = timevar$horizon, by = 1)

seasonality_transmission_calc <- function(time, timevar){
  if(time < (timevar$t_arrival + (timevar$mig_day_arrival / 2))){
    1
  } else{
    if(time < timevar$breed_start){
      2
    } else{
      if(time < timevar$hatch_start){
        1
      } else{
        if(time < timevar$t_departure){
          2
        } else{
          1
        }
      }
    }
  }
}

validate_Seasonality          <- function(func, time, timevar){
  # Calculate values
  values <- vapply(time,
                   function(x) func(x, timevar),
                   FUN.VALUE = numeric(1))
  
  df <- data.frame(
    time = time,
    value = values
  )
  
  # Generate Plot for Values
  p <- ggplot(df, aes(x = time, y = value)) +
    geom_line(color = "steelblue", size = 1) +
    labs(
      title = "Seasonality Transmission Over Time",
      x     = "Time Interval (days)",
      y     = "Computed Value"
    ) +
    theme_minimal()
  
  # Return the values calculated
  list(data = df, plot = p)
}                  # This function validates the "seasonality_transmission_calc"

seas_Val <- validate_Seasonality(seasonality_transmission_calc, time, timevar)
print(seas_Val$p)

#=======================================================================================================================#
#                                                         Sources:
#=======================================================================================================================#
# [7] Williams, T. D. “Annual Variation in Breeding Biology of Gentoo Penguins, Pygoscelis Papua, at Bird Island, South 
#         Georgia.” Journal of Zoology 222, no. 2 (1990): 247–58. https://doi.org/10.1111/j.1469-7998.1990.tb05675.x.
# [8] Catry, Paulo, Jaume Forcada, and Ana Almeida. “Demographic Parameters of Black-Browed Albatrosses Thalassarche 
#         Melanophris from the Falkland Islands.” Polar Biology 34, no. 8 (2011): 1221–29. 
#         https://doi.org/10.1007/s00300-011-0984-3.
# [9] Tickell, W. L. N., and R. Pinder. “BREEDING BIOLOGY OF THE BLACK‐BROWED ALBATROSS DIOMEDEA MELANOPHRIS AND 
#         GREY‐HEADED ALBATROSS D. CHRYSOSTOMA AT BIRD ISLAND, SOUTH GEORGIA.” Ibis 117, no. 4 (1975): 433–51. 
#         https://doi.org/10.1111/j.1474-919X.1975.tb04237.x.
