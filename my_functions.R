# Peter Ly (p-ly)

# Place your function definitions that may be needed in my_code.R and report.Rmd
# in this file, including documentation.
# You can also include any needed library() calls here

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(StatCompLab))

## Monthly Trends
#' Trend estimation function
#'
#' Estimates the trend by modeling the average TMIN of the month against the
#' covariate year for a weather station ID.
#'
#' @param id The identification number of the weather station to be modelled for
#' @param mon The month which the model should estimate the trend for
#' @param data The data frame which is being used that contains the required
#'   data in order to estimate the trend.
#'
#' @return A single row with 6 columns. The output is the ID, Month, trend,
#'   standard deviation of the trend, the lower and upper 95% confidence interval
#'   of the trend, using the t-distribution.

trend <- function(id, mon, data){
  trend_est<-NULL
  ID <- id
  Month <- mon
  fit <- lm(TMINBAR ~ 1 + Year, data = data %>%
              filter(Month == mon, ID == id))
  trend <- coef(fit)["Year"]
  sd <- sqrt(vcov(fit)["Year", "Year"])
  CI <-t(confint(fit)[2,])
  result <- cbind.data.frame(ID, Month, trend, sd, CI)
  row.names(result) <- NULL
  trend_est<- rbind.data.frame(trend_est, result)
  colnames(trend_est)[5] <- "lwr"
  colnames(trend_est)[6] <- "upr"
  as.data.frame(trend_est)
}

#' Monthly prediction function over years for TMINBAR
#'
#' Estimates the trend by modeling the average TMIN of the month against the
#' covariate year for a weather station ID and uses the model to predict the
#' observed values.
#'
#' @param id The identification number of the weather station to be modelled for
#' @param mon The month which the model should estimate the trend for
#' @param data The data frame which is being used that contains the required
#'   data in order to estimate the trend.
#'
#' @return Estimates a linear model for a the average TMIN month value and
#'   weather station. Using the linear model, it will predict the observed
#'   values in the data and its 95% prediciton interval.

Month_pred <- function(mon, id, data){
  Mon_pred <-NULL
  ID <- id
  Month <- mon
  Mon_fit <-lm(TMINBAR ~ 1 + Year, data = data %>%
                 filter(Month == mon, ID == id))
  newdata <- data.frame(Year = 1960:2018)
  Mon_pred <- rbind(Mon_pred, cbind(ID, Month,newdata,
                                    data.frame(pred = predict(Mon_fit,
                                                              newdata = newdata,
                                                              interval = "prediction"))))
}

## Seasonally varying variability
#' Generator for random test statics for one station
#'
#' The function permutes the data into random months and calculates the
#' test statistic for the variance of the monthly residual against the variance
#' of the set of residuals
#'
#' @param samples Number of times to permute/randomise the data.
#' @param data The data frame which is being used that contains the required
#'   data in order to permute.
#'
#' @return The output is the 12 rows x number of samples columns.
#'   Generates test statistics for all months in a year.

Test_stat_sampling <- function(samples, data){
  H0 <- var(data$residual)
  T_Month_sample <- matrix(NA, 12, samples)
  for(i in 1:samples){
    sampling <- sample(data$Month,
                       size = nrow(data),
                       replace = FALSE)
    newdata <- data
    newdata$Month <- sampling
    for(mon in 1:12){
      var_resid_month <- var((newdata %>%
                                filter(Month == mon))$residual)
      H0 <- H0
      T_Month_sample[mon, i] <- var_resid_month - H0
    }
  }
  T_Month_sample
}

## Estimation
#' Linear model estimator for the Edinburgh:ROYAL BOTANIC GARDE
#'
#' The function estimates a 12 linear models, for each month, using all other stations
#' and year as a covariates for the TMAX.
#'
#' @param data The data frame which is being used that contains the required
#'   estimate a linear model
#'
#' @return The output is the 12 x 12 data frame.
#'   This data frame contains the estimates for each covariate and intercept
#'   for each month.

estim_model <- function(data){
  model <- NULL
  for(i in 1:12){
    Month <- i
    ID <- "UKE00105888"
    Name <- "EDINBURGH: ROYAL BOTANIC GARDE"
    est <- lm(UKE00105888 ~ 1 + UKE00105874 + UKE00105875 + UKE00105884 + UKE00105885 + UKE00105886 + UKE00105887 + UKE00105930 + Year,
              data = data %>%
                filter(Month == i))
    coeff <- t(coef(est))
    result <- cbind.data.frame(ID, Name, Month, coeff)
    model <- rbind(model, result)
  }
  model
}

#' Predicts the daily TMAX for the Edinburgh:ROYAL BOTANIC GARDE by month
#'
#' The function estimates a linear model by month and predicts the TMAX value
#' using all other stations and year as a covariates.
#'
#' @param mon The month which the function should estimate and predict TMAX
#'   values for
#' @param data The data frame which is being used that contains the required
#'   estimate a linear model
#'
#' @return The day, month and year with the observed and predicted values using
#' the estimated model.

pred_model <- function(mon, data){
  data <- data %>%
    filter(Month == mon)
  Month <- mon
  ID <- "UKE00105888"
  Name <- "EDINBURGH: ROYAL BOTANIC GARDE"
  est <- lm(UKE00105888 ~ 1 + UKE00105874 + UKE00105875 + UKE00105884 + UKE00105885 + UKE00105886 + UKE00105887 + UKE00105930 + Year,
            data = data)
  Mon_pred <- cbind.data.frame(data,
                               data.frame(pred = predict(est,
                                                         data)))
  Mon_pred <- Mon_pred %>%
    select(Year, Month, Day, DecYear, True_value = UKE00105888, pred)
  Mon_pred
}


## Assessment
#' Cross validation function
#'
#' The function is used to 10-fold cross validate the month of the model.
#' It splits the data into 10 random parts, 9 are used to train the model.
#' The other subset of data is used to validate and produce a absolute-error
#' score and Dawid-Sebastiani score.
#'
#' @param mon The month for which the model should estimate and cross validate.
#' @param data The data frame which is being used that contains the required
#'   estimate a linear model, cross validate and produce scores for.
#'
#' @return The mean and the standard deviation of the absolute error score
#'   and Dawid-Sebastiani score for the monthly models.

cross_validation <- function(mon, data){
  sampled_data <- data %>%
    filter(Month == mon)
  index <- sample(rep_len(1:10, length.out = nrow(sampled_data)))
  indexed_data <-  cbind(sampled_data, indx = index)
  score_summary <- NULL
  for(k in 1:10){
    training_data <- indexed_data %>%
      filter(!indx == k)
    valid_data <- indexed_data %>%
      filter(indx == k)
    estim <- lm(UKE00105888 ~ 1 + UKE00105874 + UKE00105875 + UKE00105884 + UKE00105885 + UKE00105886 + UKE00105887 + UKE00105930 + Year,
                data = training_data)
    pred <- cbind.data.frame(valid_data, predict(estim,
                                                 newdata = valid_data,
                                                 interval = "prediction",
                                                 se.fit = TRUE))
    pred <- pred %>%
      mutate(sd = sqrt(se.fit^2 + residual.scale^2)) %>%
      select(Year, Month, Day, True_value = UKE00105888, fit = fit.fit, lwr = fit.lwr, upr = fit.upr, sd) %>%
      mutate(
        ae = proper_score("ae", True_value, median = fit),
        ds = proper_score("ds", True_value, mean = fit, sd = sd))
    score_summary <- rbind.data.frame(score_summary,
                                      pred %>%
                                        summarise(ae = mean(ae),
                                                  ds = mean(ds)))
  }
  Month <- mon
  score_summary <- score_summary %>%
    summarise(mean_ae = mean(ae),
              sd_ae = sqrt(var(ae)/10),
              mean_ds = mean(ds),
              sd_ds = sqrt(var(ds)/10))
  score_summary <- cbind.data.frame(Month, score_summary)
  score_summary
}
