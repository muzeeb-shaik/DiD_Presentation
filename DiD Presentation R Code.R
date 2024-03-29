
# 0. Setup ---- 

# Clear the workspace to start fresh
rm(list=ls())

# Check if the 'librarian' package is installed, install it if not, and load necessary packages
if (!require("librarian")) install.packages("librarian")
librarian::shelf(here, ggiplot,gsynth,synthdid, ggthemes, ggplot2, dplyr, 
                 haven, epiDisplay, broom, corrr, knitr, kableExtra, readr, tidyr, 
                 readxl, stargazer, fixest, fst, tibble, panelView, tidyverse)

# Set the working directory to the current project directory
setwd(here())

# 1. TWFE DiD Implementation ----

# Load the simulated dataset for difference-in-differences analysis
did_simulation_data <- read.csv("https://raw.githubusercontent.com/muzeeb-shaik/DiD_Presentation/main/did_simulation_data.csv")

# Visualization of treat_post effect over time using the panelview function (package may be missing or misspelled)
panelview(outcome ~ treat_post, data = did_simulation_data,  index = c("individual","time"), pre.post = TRUE) 


# Two-Way Fixed Effects (TWFE) Model Section
# Estimate the TWFE model without time-varying effects, controlling for individual and time fixed effects
twfe_model <- feols(outcome ~ treat_post | individual + time,cluster = ~individual, data = did_simulation_data)

# Display the table of estimates for the TWFE model
etable(twfe_model)



# Parallel Trends Plot Section
# Calculate the average outcomes by time and treat_post status
avg_outcomes <- did_simulation_data %>%
  group_by(time, treat) %>%
  summarize(average_outcome = mean(outcome), .groups = 'drop')

# Plot the parallel trends with treat_post status differentiated by color
ggplot(avg_outcomes, aes(x = time, y = average_outcome, color = as.factor(treat), group = treat)) +
  geom_line() + 
  geom_point() +
  scale_color_manual(values = c("0" = "blue", "1" = "red"), 
                     labels = c("0" = "Untreated", "1" = "Treated"), 
                     name = "Treatment Status") +
  labs(title = "Parallel Trends Plot", x = "Time Period", y = "Average Outcome") +
  theme_minimal() +
  geom_vline(xintercept = 5.5, linetype="dashed", color = "grey") +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 12)) +
  scale_x_continuous(breaks = 1:10) # Ensure all time periods are shown


# Extended TWFE Model with Time to treat_post Interaction

# Estimate a TWFE model incorporating interaction between time to treatment and treatment status, including individual and time fixed effects, and clustering standard errors at the individual level
mod_twfe = feols(outcome ~ i(time_to_treat, treat, ref = -1) | individual + time, cluster = ~individual, data = did_simulation_data)

# Interactive plot for the event study analysis showing the effect of treat_post over time
iplot(mod_twfe, xlab = 'Time to Treatment', main = 'Event study: Two Way Fixed Effects Model')


# 2. Staggered DiD Implementation ----


## Staggered DiD panelview

data(gsynth)

staggered_data <- turnout

panelview(turnout ~ policy_edr, data = staggered_data, 
          index = c("abb","year"), pre.post = TRUE, 
          by.timing = TRUE) 


## Staggered DiD TWFE Event Study

staggered_did_simulation_data <- get(data(base_stagg))

staggered_did_simulation_data <- staggered_did_simulation_data %>%   mutate(treat_post = if_else(treated==1 & time_to_treatment>=0, 1, 0))

panelview(y  ~ treat_post , data = staggered_did_simulation_data,  index = c("id","year"),by.timing = TRUE
          
          , pre.post = TRUE) 

# "Naive" TWFE DiD (note that the time to treatment for the never treated is -1000)
# (by using ref = c(-1, -1000) we exclude the period just before the treatment and 
# the never treated)
res_twfe = feols(y ~ x1 + i(time_to_treatment, ref = c(-1, -1000)) | id + year, staggered_did_simulation_data)

iplot(res_twfe)


# Sun & Abraham (2020) correction

# To implement the Sun and Abraham (2020) method,
# we use the sunab(cohort, period) function
res_sa20 = feols(y ~ x1 + sunab(year_treated, year) | id + year, staggered_did_simulation_data)

summary(res_sa20, agg = "att")

summary(res_sa20, agg = "cohort")


res_twfe_n = feols(y ~ x1 + treat_post | id + year, staggered_did_simulation_data)

summary(res_twfe_n)

# Plot the two TWFE results
iplot(list(res_twfe, res_sa20), sep = 0.5)

###

data(gsynth)

staggered_data <- turnout

panelview(turnout ~ policy_edr, data = staggered_data, 
          index = c("abb","year"), pre.post = TRUE, 
          by.timing = TRUE) 


# 3. Synthetic Control Implementation ----

data(gsynth)

synth_data <- simdata
synth_data <- synth_data %>%
  group_by(id) %>%
  mutate(treated = max(D))

# 5 treated units, 45 control units, and 30 time periods. 
#The treatment kicks at Period 21 for all treated units

panelview(Y  ~ D , data = synth_data,  index = c("id","time"), pre.post = TRUE) 

avg_outcomes_synth <- synth_data %>%
  group_by(time, treated) %>%
  summarize(average_Y = mean(Y), .groups = 'drop')


# Plot the parallel trends with treatment status differentiated by color
ggplot(avg_outcomes_synth, aes(x = time, y = average_Y, color = as.factor(treated), group = treated)) +
  geom_line() + 
  geom_point() +
  scale_color_manual(values = c("0" = "blue", "1" = "red"), 
                     labels = c("0" = "Control", "1" = "Treated"), 
                     name = "Treatment Status") +
  labs(title = "Parallel Trends Plot", x = "Time Period", y = "Average Outcome") +
  theme_minimal() +
  geom_vline(xintercept = 20.5, linetype="dashed", color = "blue") +
  theme(legend.title = element_text(size = 12), legend.text = element_text(size = 12)) 



## Generalized Synthetic Control

gsynth_out <- gsynth(Y ~ D + X1 + X2, data = synth_data, 
              index = c("id","time"), force = "two-way", 
              CV = TRUE, r = c(0, 5), se = TRUE, 
              inference = "parametric", nboots = 1000, 
              parallel = TRUE, seed = 12345)

gsynth_out_1 <- gsynth(Y ~ D , data = synth_data, 
                     index = c("id","time"), force = "unit", 
                     CV = TRUE, r = c(0, 5), se = TRUE, 
                     inference = "parametric", nboots = 1000, 
                     parallel = TRUE, seed = 12345)


gsynth_out$Ntr

gsynth_out$Nco

gsynth_out$est.avg

gsynth_out$est.att

plot(gsynth_out, type = "gap", xlab = "Time to Treatment", 
     main = "ATT Plot")
plot(gsynth_out, type="counterfactual")

plot(gsynth_out, type="counterfactual", id = "102")
plot(gsynth_out, type = "gap", xlab = "Time to Treatment", id = "102",
     main = "ATT Plot- Unit 101")

implied_weights <- gsynth_out$wgt.implied

# 4. SynthDiD Implementation ----


# Converts the 'synth_data' dataframe into a list of matrices suitable for synthetic control and difference-in-differences analyses.

setup <- panel.matrices(as.data.frame(synth_data))

# Estimates the treatment effect using the Synthetic Difference-in-Differences (SynthDiD) method.
#'tau.hat' represents the estimated treatment effect. 
tau.hat = synthdid_estimate(setup$Y, setup$N0, setup$T0)


print(summary(tau.hat))

# Calculates the standard error of the estimated treatment effect using a placebo method
# which is a way to assess the robustness of the treatment effect by simulating a scenario 
# where no treatment is administered.
se = sqrt(vcov(tau.hat, method='placebo'))

# Prints the point estimate of the treatment effect and its 95% confidence interval

sprintf('point estimate: %1.2f', tau.hat)
sprintf('95%% CI (%1.2f, %1.2f)', tau.hat - 1.96 * se, tau.hat + 1.96 * se)


# Generates a plot to visualize the estimated treatment effect (tau.hat) and its uncertainty (standard error)

plot(tau.hat, se.method='placebo')

# Estimates the treatment effect using 
# the Synthetic Control (SC) method and the traditional Difference-in-Differences (DiD) approach for comparison purposes.
tau.sc   = sc_estimate(setup$Y, setup$N0, setup$T0)
tau.did  = did_estimate(setup$Y, setup$N0, setup$T0)

# Stores the estimates from different methods in a list for easy comparison.

estimates = list(tau.did, tau.sc, tau.hat)
names(estimates) = c('Diff-in-Diff', 'Synthetic Control', 'Synthetic Diff-in-Diff')


# Estimates a Two-Way Fixed Effects (TWFE) model using the 'feols' function from the 'fixest' package. 

twfe_model_new <- feols(Y ~ D | id + time, cluster = ~id, data = synth_data)

# Displays a table of the model's estimates

etable(twfe_model_new)
