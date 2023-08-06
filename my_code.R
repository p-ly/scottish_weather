# Peter Ly (p-ly)

# Load function definitions
source("my_functions.R")
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(StatCompLab))

# General data
station_info <- ghcnd_stations
saveRDS(station_info, file = "data/station_info.rds")

ghcnd <- left_join(ghcnd_values, ghcnd_stations, by = "ID")
saveRDS(ghcnd, file = "data/ghcnd.rds")

## Monthly trend
# Data
seasonal_TMIN <- ghcnd %>%
  pivot_wider(names_from = Element, values_from = Value) %>%
  filter(!is.na(TMIN)) %>%
  group_by(ID, Name, Month, Year) %>%
  summarise(TMINBAR = mean(TMIN), .groups = "drop")
saveRDS(seasonal_TMIN, file = "data/seasonal_TMIN.rds")

# Monthly trend estimation ALL
trend_est <- NULL
ID_vect <- ghcnd_stations$ID
for(j in ID_vect){
  for(i in 1:12){
    trend_est <- rbind.data.frame(trend_est, trend(j, i, seasonal_TMIN))
  }
}
trend_est <- left_join(ghcnd_stations[,1:2], trend_est, by = "ID")
saveRDS(trend_est, file = "data/trend_est.rds")

#January prediction data
Jan_pred <-NULL
for(j in ID_vect){
  Jan_pred <- rbind.data.frame(Jan_pred, Month_pred(1, j, seasonal_TMIN))
}

Jan_pred <- left_join(ghcnd_stations[,1:2], Jan_pred, by = "ID")
saveRDS(Jan_pred, file = "data/Jan_pred.rds")

#October predicition data
Oct_pred <-NULL
for(j in ID_vect){
  Oct_pred <- rbind.data.frame(Oct_pred, Month_pred(10, j, seasonal_TMIN))
}

Oct_pred <- left_join(ghcnd_stations[,1:2], Oct_pred, by = "ID")
saveRDS(Oct_pred, file = "data/Oct_pred.rds")

#November predicition data
Nov_pred <-NULL
for(j in ID_vect){
  Nov_pred <- rbind.data.frame(Nov_pred, Month_pred(11, j, seasonal_TMIN))
}

Nov_pred <- left_join(ghcnd_stations[,1:2], Nov_pred, by = "ID")
saveRDS(Nov_pred, file = "data/Nov_pred.rds")

#December predicition data
Dec_pred <-NULL
for(j in ID_vect){
  Dec_pred <- rbind.data.frame(Dec_pred, Month_pred(12, j, seasonal_TMIN))
}

Dec_pred <- left_join(ghcnd_stations[,1:2], Dec_pred, by = "ID")
saveRDS(Dec_pred, file = "data/Dec_pred.rds")

## Seasonally varying variability
#Balmoral residual data
Bal_residual_data <- ghcnd %>%
  pivot_wider(names_from = Element, values_from = Value) %>%
  filter(!is.na(TMIN), ID == "UKE00105875") %>%
  group_by(ID, Name, Month, Year)%>%
  select(ID, Name, Year, Month, Day, DecYear, TMIN) %>%
  mutate(TMINBAR = mean(TMIN)) %>%
  mutate(residual = (TMIN - TMINBAR))

# Observed test statistic
H0 <- var(Bal_residual_data$residual)
T_Month <- NULL
for(mon in 1:12){
  T_Month[mon] <- var((Bal_residual_data %>%
                         filter(Month == mon))$residual)
}
T_Month <- cbind.data.frame(Month = 1:12, T_Month, H0)
T_Month <- T_Month %>%
  mutate(T_stat = T_Month - H0)
saveRDS(T_Month, file = "data/T_Month.rds")

# Permutated test statistic
Test_stats <- Test_stat_sampling(2500, Bal_residual_data)

# pvalues for randomisation test
pvalue <- NULL
for(mon in 1:12){
  pvalue[mon] <- (1- mean(T_Month[mon,4] >= Test_stats[mon,]))
}
pvalues <- cbind.data.frame(Month = 1:12, pvalue = pvalue)
saveRDS(pvalues, file = "data/pvalues.rds")

## Estimation
# Estimation Data
data_est <- ghcnd %>%
  pivot_wider(names_from = Element, values_from = Value) %>%
  group_by(ID, Month, Year) %>%
  select(ID, Year, Month, Day, DecYear, TMAX) %>%
  pivot_wider(names_from = ID, values_from = TMAX) %>%
  filter(!is.na(UKE00105874), !is.na(UKE00105875), !is.na(UKE00105884),
         !is.na(UKE00105885), !is.na(UKE00105886), !is.na(UKE00105887),
         !is.na(UKE00105888), !is.na(UKE00105930))

# Model estimation for Edinburgh weather station
Edin_model  <- estim_model(data_est)
saveRDS(Edin_model, file = "data/Edin_model.rds")

# 12 Month prediction data for Edinburgh weather station
Edin_pred <- NULL
for(mon in 1:12){
  Edin_pred <- rbind.data.frame(Edin_pred, pred_model(mon,data_est))
}
saveRDS(Edin_pred, file = "data/Edin_pred.rds")

# November prediction data for Edinburgh weather station
Edin_Nov_pred <- pred_model(11, data_est)
saveRDS(Edin_Nov_pred, file = "data/Edin_Nov_pred.rds")

# May prediction data for Edinburgh weather station
Edin_May_pred <- pred_model(5, data_est)
saveRDS(Edin_May_pred, file = "data/Edin_May_pred.rds")

## Assessment
set.seed(123)

# Proper scores for Edinburgh weather station monthly models
score_summary <- NULL
for(mon in 1:12){
  score_summary <- rbind.data.frame(score_summary, cross_validation(mon, data_est))
}
saveRDS(score_summary, file = "data/score_summary.rds")

# Proper scores plot data
score_summary_confint <-  score_summary %>%
  mutate(ae_lwr = mean_ae + sd_ae *qt(.025,df=9),
         ae_upr = mean_ae + sd_ae *qt(.975,df=9),
         ds_lwr = mean_ds + sd_ds *qt(.025,df=9),
         ds_upr = mean_ds + sd_ds *qt(.975,df=9))
saveRDS(score_summary_confint, file = "data/score_summary_confint.rds")
