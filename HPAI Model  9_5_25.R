####PACKAGES####
#Install packages with install.packages("name of the package")

library(deSolve) #Used to solve ordinary differential equations (ODEs). More info at https://cran.r-project.org/web/packages/deSolve/deSolve.pdf
library(ggplot2) #Used to plot results from the model. More info at https://cran.r-project.org/web/packages/ggplot2/ggplot2.pdf
#MODEL OVERVIEW#
#You can write a model overview here.

# Time points
horizon<-365
time=seq(from=0,to=horizon,by=1)

###SWITCH for seasonality in transmission: ON/OFF
switch<- 1
#Set the starting values (these values indicate the initial population values; i.e., the number of individuals that would be
#in each compartment at the beginning of the simulation (time=0))

#Populations
Res_gulls <- (74000*0.1)/12 #10% of gulls are resident
Res_gentoo<-200000/20 #Resident gentoo
albatross <- 0 #Resident Albatross 

gulls_migrate_in <-(84000*0.9)/20  #90% of gulls are migratory

albatross_migrate_in <-(951000)/20  #100% of albatross are migratory

Initial_values <-c(
  ###Compartments
  S_gulls= Res_gulls,    #Susceptible
  I_gulls= 0,    #Infected                   
  R_gulls= 0,    #Recovered 
  D_gulls = 0, 
  
  S_res_gentoo= Res_gentoo,    #Susceptible
  I_res_gentoo= 0,            #Infected                   
  R_res_gentoo= 0,            #Recovered
  D_res_gentoo = 0, 
  
  S_albatross= albatross,    #Susceptible
  I_albatross= 0,            #Infected                   
  R_albatross= 0,            #Recovered 
  D_albatross = 0, 
  



  
  timer=0
  
)

# Parameters for the first sigmoid curve
k1 = 0.1
x01 = 31
a1 = 2
b1 = 1

# Parameters for the second sigmoid curve
k2 = 0.1
x02 = 193
a2 = 2
b2 = 1

k3 = 0.1
x03 = 357.5  # Midpoint of the third sigmoid
a3 = 1.05    # Upper bound of the third sigmoid (previously `top`)
b3 = 0.98      # Lower bound of the third sigmoid (previously `bottom`)

#b4=0

#Gull Parameters   
prev_infection_gulls <- 0.01

mig_day_arrival<- 31

#scalingfactor <- 0.175

beta_gulls <- 0.000003

gentoogull_trans <- 0.000003

mu_gulls <- 0.0001

immunity_duration_gulls <- 182.5

sigma_gulls <- 1/immunity_duration_gulls #How is SIGMA being estimated?? This should be based on the period in R, right?

infectious_period_gulls <- 8

gamma_gulls <- 1/infectious_period_gulls

death_prob_gulls <- 0

death_rate_gulls <- death_prob_gulls/infectious_period_gulls

##Gentoo Parameters
mu_gentoo <- 0.0001 #Natural mortality rate for gentoo (per bird per day) 

beta_gentoo_density <- 0.000003

gullgentoo_trans <- 0.000003

immunity_duration_gentoo <- 182.5

sigma_gentoo <- 1/immunity_duration_gentoo #This parameter can be set to 0 to go from an "SIRS" to an "SIR" model

infectious_period_gentoo<- 8 #Duration of the infectious period (days)

gamma_gentoo <- 1/infectious_period_gentoo #Recovery rate from HPAI infection in gentoos (per bird per day)

death_prob_gentoo <- 0.5 #Probability of death from HPAI infection in gentoos ##[need to be adjusted]##

death_rate_gentoo <- death_prob_gentoo/infectious_period_gentoo  #Death rate from HPAI infection in gentoos (per bird per day)

##Albatross Parameters
prev_infection_albatross = 0

mu_albatross <- 0.0000626 #Natural mortality rate for albatross (per bird per day) 

beta_albatross <- 0.000003 #HPAI transmission rate between albatross (per bird per day)

gullalbatross_trans <- 0.000003

albatrossgull_trans <- 0.000003

immunity_duration_albatross <- 182.5

sigma_albatross <- 1/immunity_duration_albatross #This parameter can be set to 0 to go from an "SIRS" to an "SIR" model

infectious_period_albatross <- 8 #Duration of the infectious period (days)

gamma_albatross <- 1/infectious_period_albatross #Recovery rate from HPAI infection in albatross (per bird per day)

death_prob_albatross <- 0.5 #Probability of death from HPAI infection in albatross ##[need to be adjusted]##

death_rate_albatross <- death_prob_albatross/infectious_period_albatross  #Death rate from HPAI infection in albatross (per bird per day)

parameters<- c(prev_infection_gulls,
               mu_gulls, 
               sigma_gulls,
               beta_gulls, 
               gentoogull_trans, 
               gamma_gulls,
               death_rate_gulls,
               
               mu_gentoo,
               sigma_gentoo, 
               beta_gentoo_density,
               gullgentoo_trans,
               gamma_gentoo,
               death_rate_gentoo, 
               
               mu_albatross,
               sigma_albatross, 
               beta_albatross,
               gullalbatross_trans,
               albatrossgull_trans,
               gamma_albatross,
               death_rate_albatross
)


#This function can be used to establish migration season
PiecewiseFlux<-function(time, t1, t2, b) {ifelse(t1<=floor(time %% 365) && floor(time %% 365) <=t2,b,0)} #This is to incorporate gulls' migration

PiecewiseFlux_albatross <- function(time, t1, t2, b) {
  ifelse(t1 <= floor(time %% 365) && floor(time %% 365) <= t2, b, 0)
}

t_arrival=62

post_migration_reduction_gulls <- 0.1 #90% reduction
post_migration_reduction_albatross <- 0 #

# Define the event function
event_function <- function(time, state, parms) {
  with(as.list(state), {
    if (time == 300) {
      # Reduce the populations by 90%
      S_gulls <- S_gulls * post_migration_reduction_gulls
      I_gulls <- I_gulls * post_migration_reduction_gulls
      R_gulls <- R_gulls * post_migration_reduction_gulls
      S_albatross <- S_albatross * post_migration_reduction_albatross
      I_albatross <- I_albatross * post_migration_reduction_albatross
      R_albatross <- R_albatross * post_migration_reduction_albatross
    }
    return(c(S_gulls, I_gulls, R_gulls, D_gulls,
             S_res_gentoo, I_res_gentoo, R_res_gentoo, D_res_gentoo, 
             S_albatross, I_albatross, R_albatross, D_albatross,
             timer)) # Return updated state
  })
}

#class(time)
# print(time)

####MODEL EQUATIONS#####
HPAI_dyn <- function(time,state,parameters){   
  if(!is.numeric(time)) {
    stop("Time should be numeric")
  }
  with(as.list(c(state,parameters)),{
    
    if(switch==0){seasonality_transmission = rep(1, length(time))
    }else{seasonality_transmission = ifelse(time <= 100,
                                            b1 + (a1 - b1) / (1 + exp(-k1 * (time - x01))),
                                            ifelse(time <= 349,
                                                   a2 - (a2 - b2) * (1 / (1 + exp(-k2 * (time - x02)))),
                                                   b3 + (a3 - b3) / (1 + exp(-k3 * (time - x03)))))}
    
    #=== Gulls ===#
    dS_gulls = (1-prev_infection_gulls) * PiecewiseFlux(time, t_arrival, t_arrival+mig_day_arrival, gulls_migrate_in/mig_day_arrival) + #Arrival of susceptible gulls through migration
      ifelse(time >= 93 & time <= 225, 2.5 * mu_gulls * (S_gulls + I_gulls + R_gulls), 0) + # Breeding term for gulls
      sigma_gulls*R_gulls - #Loss of immunity in gulls
      
      seasonality_transmission * ifelse(S_gulls + I_gulls + R_gulls == 0, 0,  (beta_gulls * I_gulls * S_gulls)) - #Between-gulls HPAI transmission 
      seasonality_transmission * ifelse(S_res_gentoo + I_res_gentoo + R_res_gentoo == 0, 0,  (gentoogull_trans * I_res_gentoo * S_gulls)) - #gentoo-to-gull HPAI transmission (frequency-dependent)
      seasonality_transmission * ifelse(S_albatross + I_albatross + R_albatross == 0, 0,  (albatrossgull_trans * I_albatross * S_gulls)) -
      
      mu_gulls*S_gulls  #Natural mortality in gulls

    
    dI_gulls = prev_infection_gulls * PiecewiseFlux(time, t_arrival, t_arrival+ mig_day_arrival, gulls_migrate_in/mig_day_arrival) + #Arrival of infectious gulls through migration
      
      seasonality_transmission * ifelse(S_gulls + I_gulls + R_gulls  == 0, 0,  (beta_gulls * I_gulls * S_gulls)) + #Between-gulls HPAI transmission
      seasonality_transmission * ifelse(S_res_gentoo + I_res_gentoo + R_res_gentoo == 0, 0, (gentoogull_trans * I_res_gentoo * S_gulls)) + #gentoo-to-gull HPAI transmission (frequency-dependent)
      seasonality_transmission * ifelse(S_albatross + I_albatross + R_albatross == 0, 0, (albatrossgull_trans * I_albatross * S_gulls)) -
      
      (gamma_gulls + mu_gulls + death_rate_gulls) * I_gulls #Recovery, natural mortality, and mortality due to HPAI in infectious gulls
    
    dR_gulls = gamma_gulls * I_gulls - #Recovery from infection in gulls
      mu_gulls * R_gulls - #Natural mortality in gulls
      sigma_gulls * R_gulls #Loss of immunity in gulls
    
    dD_gulls =  death_rate_gulls * I_gulls
    
    #=== Resident gentoos ===#
    dS_res_gentoo = ifelse(time >= 31 & time <= 212, 2.0 * mu_gentoo * (S_res_gentoo+I_res_gentoo+R_res_gentoo), 0) +
      sigma_gentoo*R_res_gentoo - #Loss of immunity in gentoos- #Breeding season in gentoos. #2.985 try to change this number to reach equilibrium
      
      seasonality_transmission * beta_gentoo_density*I_res_gentoo*S_res_gentoo - 
      seasonality_transmission * ifelse(S_gulls + I_gulls + R_gulls == 0, 0,  (gullgentoo_trans * I_gulls * S_res_gentoo)) - #Gull-to-gentoo HPAI transmission (frequency-dependent)
     
      mu_gentoo*S_res_gentoo
    
    dI_res_gentoo = seasonality_transmission * beta_gentoo_density * I_res_gentoo * S_res_gentoo + #Between-gentoos HPAI transmission 
      seasonality_transmission * ifelse(S_gulls + I_gulls + R_gulls == 0, 0,  (gullgentoo_trans * I_gulls * S_res_gentoo)) - #Gull-to-gentoo HPAI transmission (frequency-dependent)
      (gamma_gentoo + mu_gentoo + death_rate_gentoo) * I_res_gentoo
    
    dR_res_gentoo = gamma_gentoo * I_res_gentoo - #Recovery from infection in gentoos
      mu_gentoo * R_res_gentoo - #Natural mortality in gentoos
      sigma_gentoo * R_res_gentoo #Loss of immunity in gentoos

    dD_res_gentoo = death_rate_gentoo * I_res_gentoo
    
    #=== Albatross ===#
    dS_albatross = (1 - prev_infection_albatross) * PiecewiseFlux_albatross(time, t_arrival, t_arrival + mig_day_arrival, albatross_migrate_in/mig_day_arrival) +
      ifelse(time >= 31 & time <= 245, 1.7 * mu_albatross * (S_albatross+I_albatross+R_albatross), 0) +
      sigma_albatross * R_albatross -
      (seasonality_transmission * beta_albatross * I_albatross * S_albatross) - 
      seasonality_transmission * ifelse(S_gulls + I_gulls + R_gulls == 0, 0,   (gullalbatross_trans * I_gulls * S_albatross)) - 
      mu_albatross * S_albatross 

    
    dI_albatross = (seasonality_transmission * beta_albatross * I_albatross * S_albatross) + 
      seasonality_transmission * ifelse(S_gulls + I_gulls + R_gulls == 0, 0,   (gullalbatross_trans * I_gulls * S_albatross)) -
      (gamma_albatross + mu_albatross + death_rate_albatross) * I_albatross
    
    dR_albatross =    gamma_albatross * I_albatross - #Recovery from infection in albatross
      mu_albatross * R_albatross - #Natural mortality in albatross
      sigma_albatross * R_albatross #Loss of immunity in albatross
    
    dD_albatross = death_rate_albatross * I_albatross

    dtimer = 1
    
    return(list(c(dS_gulls, dI_gulls, dR_gulls, dD_gulls, 
                  dS_res_gentoo, dI_res_gentoo, dR_res_gentoo, dD_res_gentoo,
                  dS_albatross, dI_albatross, dR_albatross, dD_albatross,

                  dtimer)))
  }
  )
}

output <- as.data.frame(lsode(
  y = Initial_values,
  func = HPAI_dyn,
  parms = parameters,
  time = time,
  events = list(func = event_function, time = 300)  # Add event at day 300
))

colnames(output)<- c("time", "Susceptible_Gulls", "Infectious_Gulls", "Recovered_Gulls","Dead_Gulls",
                     "Susceptible_Resident_gentoos","Infectious_Resident_gentoos","Recovered_Resident_gentoos",  "Dead_Gentoo",
                     "Susceptible_Albatross", "Infectious_Albatross", "Recovered_Albatross", "Dead_Albatross",
                       "Timer")


head(output) #Shows you the first 6 rows of the dataset


if(switch == 0){
  seasonality_transmission = rep(1, length(time))
} else{seasonality_transmission = ifelse(time <= 100,
                                         b1 + (a1 - b1) / (1 + exp(-k1 * (time - x01))),
                                         ifelse(time <= 349,
                                                a2 - (a2 - b2) * (1 / (1 + exp(-k2 * (time - x02)))),
                                                b3 + (a3 - b3) / (1 + exp(-k3 * (time - x03)))))}

final_values = tail(output, 1)

death_counts_table <- data.frame(
  Species = c("Dead Gulls", "Dead Gentoo", "Dead Albatross"),
  Deaths = c(final_values$Dead_Gulls, final_values$Dead_Gentoo, final_values$Dead_Albatross)
)

# Print the table
print(death_counts_table)

# library(reshape2)
# library(ggplot2)
# library(dplyr)
#


library(ggplot2)
library(dplyr)
library(tidyr)

#== PLOTTING EVERYONE == #
sir_data <- output %>%
  select(time,
         Susceptible_Gulls, Infectious_Gulls, Recovered_Gulls,
         Susceptible_Resident_gentoos, Infectious_Resident_gentoos, Recovered_Resident_gentoos,
         Susceptible_Albatross, Infectious_Albatross, Recovered_Albatross) %>%
  pivot_longer(-time, names_to = "Compartment", values_to = "Population")

# Add species column
sir_data <- sir_data %>%
  mutate(
    Species = case_when(
      grepl("Gulls", Compartment) ~ "Gulls (Apex)",
      grepl("gentoo", Compartment) ~ "Gentoo Penguins",
      grepl("Albatross", Compartment) ~ "Albatross (Meso)"
    ),
    State = case_when(
      grepl("Susceptible", Compartment) ~ "Susceptible",
      grepl("Infectious", Compartment) ~ "Infectious",
      grepl("Recovered", Compartment) ~ "Recovered"
    )
  )

# Month labels
month_day <- c(1, 32, 62, 93, 123, 154, 184, 214, 245, 275, 305, 335, 365)
month_labels <- c("Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug")

# Plot
ggplot(sir_data, aes(x = time, y = Population, color = State, linetype = Species)) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = c(
      "Susceptible" = "#4169E1",  # Blue
      "Infectious" = "#DC143C",   # Red
      "Recovered" = "#228B22"     # Green
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Gentoo Penguins" = "solid",
      "Gulls (Apex)" = "dashed",
      "Albatross (Meso)" = "dotted"
    )
  ) +
  scale_x_continuous(
    breaks = month_day,
    labels = month_labels,
    name = "Month of Year"
  ) +
  ylab("Number of individuals") +
  labs(
    color = "Compartment",
    linetype = "Species",
    title = "All species dynamics"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5)
  )

# == PLOTTING GULLS AND PENGUINS ONLY == # 
out_long <- reshape2::melt(output[, c("time", 
                                      "Susceptible_Gulls", "Infectious_Gulls", "Recovered_Gulls",
                                      "Susceptible_Resident_gentoos", "Infectious_Resident_gentoos", "Recovered_Resident_gentoos")],
                           id = "time")

# Assign species based on variable names
out_long$species <- ifelse(grepl("Gulls", out_long$variable), "Apex predators", "Penguins")

# Assign compartment (state) based on variable names
out_long$state <- ifelse(grepl("Susceptible", out_long$variable), "Susceptible",
                         ifelse(grepl("Infectious", out_long$variable), "Infectious", "Recovered"))

# Factor order for species and states
out_long$species <- factor(out_long$species, levels = c("Penguins", "Apex predators"))
out_long$state <- factor(out_long$state, levels = c("Susceptible", "Infectious", "Recovered"))

# Define month breaks and labels for the x-axis
month_day <- c(1, 32, 62, 93, 123, 154, 184, 214, 245, 275, 305, 335, 365)
month_labels <- c("Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug")

# Plot
ggplot(out_long, aes(x = time, y = value, color = state)) +
  geom_line(aes(linetype = species), linewidth = 1) +
  scale_color_manual(
    values = c(
      "Susceptible" = "#4169E1",  # Blue
      "Infectious" = "#DC143C",   # Red
      "Recovered" = "#228B22"     # Green
    )
  ) +
  scale_linetype_manual(
    values = c(
      "Penguins" = "solid",
      "Apex predators" = "dashed"
    )
  ) +
  scale_x_continuous(
    breaks = month_day,
    labels = month_labels,
    name = "Month of Year"
  ) +
  ylab("Number of individuals") +
  labs(color = "Compartment", linetype = "Bird group") +
  ggtitle("Penguins and Gulls") +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 12, hjust = 0.5)
  )

# == Mesopredators only no Susceptibles == #
out_long <- reshape2::melt(output[, c("time", "Infectious_Albatross", "Recovered_Albatross")],
                           id = "time")

# Assign species and state based on variable names
out_long$species <- "Mesopredators"  # All are Albatross
out_long$state <- ifelse(grepl("Infectious", out_long$variable), "Infectious", "Recovered")

# Set factor levels
out_long$species <- factor(out_long$species, levels = c("Mesopredators"))
out_long$state <- factor(out_long$state, levels = c("Infectious", "Recovered"))

# Define month breaks and labels
month_day <- c(1, 32, 62, 93, 123, 154, 184, 214, 245, 275, 305, 335, 365)
month_labels <- c("Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug")

# Plot
ggplot(out_long, aes(x = time, y = value, color = state)) +
  geom_line(aes(linetype = species), linewidth = 1) +
  scale_color_manual(
    values = c(
      "Infectious" = "#DC143C",  # Red
      "Recovered" = "#228B22"    # Green
    )
  ) +
  scale_linetype_manual(
    values = c("Mesopredators" = "dotted")
  ) +
  scale_x_continuous(
    breaks = month_day,
    labels = month_labels,
    name = "Month of Year"
  ) +
  ylab("Number of individuals") +
  labs(color = "Compartment", linetype = "Bird group") +
  ggtitle("Infectious and Recovered Albatross") +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 12, hjust = 0.5)
  )

# == Gulls Only No Susceptible == #
gulls_data <- output %>%
  select(time, Infectious_Gulls, Recovered_Gulls)

# Reshape the data for ggplot (long format)
gulls_long <- gulls_data %>%
  gather(key = "Compartment", value = "Population", -time)

# Plot the Infectious and Recovered Gulls populations
ggplot(gulls_long, aes(x = time, y = Population, color = Compartment)) +
  geom_line(aes(linetype = "Gulls"), size = 1) +  # Dashed line for Gulls
  scale_color_manual(
    values = c(
      "Infectious_Gulls" = "#DC143C",  # Red for Infectious
      "Recovered_Gulls" = "#228B22"    # Green for Recovered
    )
  ) +
  scale_linetype_manual(
    values = c("Gulls" = "dashed")  # Dashed line for Gulls
  ) +
  scale_x_continuous(
    breaks = month_day,
    labels = month_labels,
    name = "Month of Year"
  ) +
  ylab("Number of individuals") +
  labs(color = "Compartment", linetype = "Bird group") +
  ggtitle("Infectious and Recovered Gulls") + 
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5)
  )


# == Penguins Only No Susceptible == # 
gentoo_data <- output %>%
  select(time, Infectious_Resident_gentoos, Recovered_Resident_gentoos)

# Reshape the data for ggplot (long format)
gentoo_long <- gentoo_data %>%
  gather(key = "Compartment", value = "Population", -time)

# Plot the Infectious and Recovered Gentoo populations
ggplot(gentoo_long, aes(x = time, y = Population, color = Compartment)) +
  geom_line(aes(linetype = "Gentoo Penguins"), size = 1) +  # Solid line for Penguins
  scale_color_manual(
    values = c(
      "Infectious_Resident_gentoos" = "#DC143C",  # Red for Infectious
      "Recovered_Resident_gentoos" = "#228B22"    # Green for Recovered
    ),
    labels = c("Infectious_Resident_gentoos" = "Infectious Gentoo", 
               "Recovered_Resident_gentoos" = "Recovered Gentoo")  # Shorten legend labels
  ) +
  scale_linetype_manual(
    values = c("Gentoo Penguins" = "solid")  # Solid line for Gentoo Penguins
  ) +
  scale_x_continuous(
    breaks = month_day,
    labels = month_labels,
    name = "Month of Year"
  ) +
  ylab("Number of individuals") +
  labs(
    color = "Compartment", 
    linetype = "Bird group", 
    title = "Infectious and Recovered Gentoo Penguins"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 14, hjust = 0.5)
  )

# == Max Cumulative Dead on a Given Day?? == #
# output$Daily_Dead_Gulls <- c(0, diff(output$Dead_Gulls))
# output$Daily_Dead_Gentoo <- c(0, diff(output$Dead_Gentoo))
# output$Daily_Dead_Albatross <- c(0, diff(output$Dead_Albatross))
# 
# peak_day_albatross <- output[which.max(output$Daily_Dead_Albatross), ]
# peak_day_gulls <- output[which.max(output$Daily_Dead_Gulls), ]
# peak_day_gentoo <- output[which.max(output$Daily_Dead_Gentoo), ]
# 
# # Print the peak daily death day for each species
# cat("Peak daily death for Albatross on day", peak_day_albatross$time, "with", peak_day_albatross$Daily_Dead_Albatross, "deaths\n")
# cat("Peak daily death for Gulls on day", peak_day_gulls$time, "with", peak_day_gulls$Daily_Dead_Gulls, "deaths\n")
# cat("Peak daily death for Gentoo on day", peak_day_gentoo$time, "with", peak_day_gentoo$Daily_Dead_Gentoo, "deaths\n")


# Now, print the first few rows to check
#head(output)

print(death_counts_table)

