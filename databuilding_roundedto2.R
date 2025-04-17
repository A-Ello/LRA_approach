# Packages
library(survival)
library(ggplot2)
library(tictoc)
library(patchwork)
library(ggpubr)
library(pander)
library(systemfonts)
library(dplyr)
library(stringr)
library(gridExtra)
library(car)


# Fixed parameters
sample_size_values = seq(10, 1000, 10) # sample size (number of at risk patients for ex.)
hazard_surv = 1

# To obtain the values of the survival
vec_sample_size <- NULL
for (i in 1:length(sample_size_values)) {
  vec_sample_size <- c(vec_sample_size, rep(sample_size_values[i], sample_size_values[i] - 1))
}
sample_size_and_survival <- matrix(0, ncol = 2, nrow=length(vec_sample_size), dimnames = list(NULL, c("sample_size", "survival")))
sample_size_and_survival[,"sample_size"] <- vec_sample_size

ndeaths = 1 # number of deaths
for (i in 1:nrow(sample_size_and_survival)) {
  sample_size_and_survival[i, "survival"] = (sample_size_and_survival[i,"sample_size"] - ndeaths)/sample_size_and_survival[i,"sample_size"]
  
  if (ndeaths == sample_size_and_survival[i,"sample_size"]-1) {
    ndeaths = 1}  else {ndeaths = ndeaths+1}
}

t0 = -log(sample_size_and_survival[, "survival"])
# sum(is.na(t0)) # if 0 >> all good


# Find hazard of censoring for a given censoring rate
# inflation_coeff = (hazard_surv/(hazard_surv+hazard_cens))*(exp((hazard_surv)*t0)-1)/(exp(t0*He)-1) # from Nemes et al., Statistica, 2020
# censoring = (hazard_surv/(hazard_surv+hazard_cens))(1-exp(-(hazard_surv+hazard_cens)*t0))


find_hazard_cens <- function(t0, cens_value) {
  # Define the function to find the root of
  func <- function(hazard_cens) {
    (hazard_cens / (hazard_surv + hazard_cens)) * (1 - exp(- (hazard_surv + hazard_cens) * t0)) - cens_value
  }
  
  # Use uniroot to solve for hazard_cens in the range [lower, upper]
  hazard_cens_solution <- uniroot(func, lower = - cens_value, upper = 1000, tol = 1e-06)$root
  return(hazard_cens_solution)
}

hazard_cens_10 <- sapply(t0, function(t) find_hazard_cens(t, 0.10))
hazard_cens_25 <- sapply(t0, function(t) find_hazard_cens(t, 0.25))
hazard_cens_50 <- sapply(t0, function(t) find_hazard_cens(t, 0.50))

# Inflation coefficients
inflation_coeff0 = rep(1, length(t0)) # ones because the standard error is not inflated when there is no censoring so it is multiplied by 1
# We also find  1 1 1 1 1 when we use the code above as we did for the other censoring levels

inflation_coeff10 = (hazard_surv/(hazard_surv+hazard_cens_10))*(exp((hazard_surv+hazard_cens_10)*t0) - 1)/(exp(t0*hazard_surv) - 1) 
inflation_coeff25 = (hazard_surv/(hazard_surv+hazard_cens_25))*(exp((hazard_surv+hazard_cens_25)*t0) - 1)/(exp(t0*hazard_surv) - 1) 
inflation_coeff50 = (hazard_surv/(hazard_surv+hazard_cens_50))*(exp((hazard_surv+hazard_cens_50)*t0) - 1)/(exp(t0*hazard_surv) - 1) 


# Checkpoint

max((hazard_cens_10 / (hazard_surv + hazard_cens_10)) * (1 - exp(- (hazard_surv + hazard_cens_10) * t0)) - 0.10)
max((hazard_cens_25 / (hazard_surv + hazard_cens_25)) * (1 - exp(- (hazard_surv + hazard_cens_25) * t0)) - 0.25)
max((hazard_cens_50 / (hazard_surv + hazard_cens_50)) * (1 - exp(- (hazard_surv + hazard_cens_50) * t0)) - 0.50)


min((hazard_cens_10 / (hazard_surv + hazard_cens_10)) * (1 - exp(- (hazard_surv + hazard_cens_10) * t0)) - 0.10)
min((hazard_cens_25 / (hazard_surv + hazard_cens_25)) * (1 - exp(- (hazard_surv + hazard_cens_25) * t0)) - 0.25)
min((hazard_cens_50 / (hazard_surv + hazard_cens_50)) * (1 - exp(- (hazard_surv + hazard_cens_50) * t0)) - 0.50)


var_names_data_cens <- c("sample_size", "survival_initial", "survival_rounded", 
                         "t0", "inflation_coeff", "se_id_initial", "se_id_infl",
                         
                         "id_low", "id_up", "truncated_id", "AA_id", "RA_id",
                         "log_low", "log_up", "truncated_log", "AA_log", "RA_log",
                         "cloglog_low", "cloglog_up", "truncated_cloglog", "AA_cloglog", "RA_cloglog",
                         "logit_low", "logit_up", "truncated_logit", "AA_logit", "RA_logit",
                         "arcsine_low", "arcsine_up", "truncated_arcsine", "AA_arcsine", "RA_arcsine",
                         
                         "err_log_instead_of_id", "err_cloglog_instead_of_id", "err_logit_instead_of_id", "err_arcsine_instead_of_id", "true_id",
                         "err_id_instead_of_log", "err_cloglog_instead_of_log", "err_logit_instead_of_log", "err_arcsine_instead_of_log", "true_log",
                         "err_id_instead_of_cloglog", "err_log_instead_of_cloglog","err_logit_instead_of_cloglog", "err_arcsine_instead_of_cloglog", "true_cloglog",
                         "err_id_instead_of_logit", "err_log_instead_of_logit","err_cloglog_instead_of_logit", "err_arcsine_instead_of_logit", "true_logit",
                         "err_id_instead_of_arcsine", "err_log_instead_of_arcsine","err_cloglog_instead_of_arcsine", "err_logit_instead_of_arcsine", "true_arcsine")

#err_log_instead_of_id = extraction using log transformation instead of identity 
#true_id = extraction and computation function are both identity
#AA = Absolute asymmetry : Upper bound + lower bound - 2survival = (UB - survival) - (survival - LB)
#RA = Relative asymmetry : Upper bound / lower bound

data_cens0_r2  <- as.data.frame(matrix(0, nrow = nrow(sample_size_and_survival), ncol = length(var_names_data_cens)))
data_cens10_r2 <- as.data.frame(matrix(0, nrow = nrow(sample_size_and_survival), ncol = length(var_names_data_cens)))
data_cens25_r2 <- as.data.frame(matrix(0, nrow = nrow(sample_size_and_survival), ncol = length(var_names_data_cens)))
data_cens50_r2 <- as.data.frame(matrix(0, nrow = nrow(sample_size_and_survival), ncol = length(var_names_data_cens)))

colnames(data_cens0_r2) <- colnames(data_cens10_r2) <- colnames(data_cens25_r2) <- colnames(data_cens50_r2) <- var_names_data_cens

data_cens0_r2$sample_size <- data_cens10_r2$sample_size <- data_cens25_r2$sample_size <- data_cens50_r2$sample_size <- sample_size_and_survival[,"sample_size"]
data_cens0_r2$survival_initial <- data_cens10_r2$survival_initial <- data_cens25_r2$survival_initial <- data_cens50_r2$survival_initial <- sample_size_and_survival[,"survival"]
data_cens0_r2$survival_rounded <- data_cens10_r2$survival_rounded <- data_cens25_r2$survival_rounded <- data_cens50_r2$survival_rounded <- round(sample_size_and_survival[,"survival"], 2)
data_cens0_r2$t0 <- data_cens10_r2$t0 <- data_cens25_r2$t0 <- data_cens50_r2$t0 <- t0

data_cens0_r2$inflation_coeff  <- inflation_coeff0
data_cens10_r2$inflation_coeff <- inflation_coeff10
data_cens25_r2$inflation_coeff <- inflation_coeff25
data_cens50_r2$inflation_coeff <- inflation_coeff50



# No censoring ----
for (i in 1:nrow(data_cens0_r2)) {
  
  print(i)
  ## standard error on different scales
  data_cens0_r2$se_id_initial[i] = sqrt((data_cens0_r2$survival_initial[i]*(1-data_cens0_r2$survival_initial[i]))/data_cens0_r2$sample_size[i])
  data_cens0_r2$se_id_infl[i] = data_cens0_r2$se_id_initial[i]* data_cens0_r2$inflation_coeff[i]
  
  se_id = data_cens0_r2$se_id_infl[i]
  se_log = se_id/data_cens0_r2$survival_initial[i]
  se_cloglog = - se_id / (data_cens0_r2$survival_initial[i]*log(data_cens0_r2$survival_initial[i]))
  se_logit = se_id / (data_cens0_r2$survival_initial[i]*(1-data_cens0_r2$survival_initial[i]))
  se_arcsine =  se_id / (2*sqrt((data_cens0_r2$survival_initial[i])*(1-data_cens0_r2$survival_initial[i])))
  
  
  ## 95% confidence interval of survival computed using the different methods
  # If the lower bound is equal to 0 or the upper bound equal to 1 whether due to rounding or truncation, we consider that the interval is truncated/trimmed
  #Thus, we don't need the bounds not rounded for the following steps, we only need the rounded bounds
  
  data_cens0_r2$id_low[i] = round(data_cens0_r2$survival_initial[i] - qnorm(0.975) * se_id, 2)
  data_cens0_r2$id_up[i]  = round(data_cens0_r2$survival_initial[i] + qnorm(0.975) * se_id, 2)
  
  data_cens0_r2$log_low[i] = round(data_cens0_r2$survival_initial[i] * exp(-qnorm(0.975)*se_log), 2)
  data_cens0_r2$log_up[i]  = round(data_cens0_r2$survival_initial[i] * exp(qnorm(0.975)*se_log), 2)
  
  data_cens0_r2$cloglog_low[i] = round(data_cens0_r2$survival_initial[i] ^ (exp(qnorm(0.975)*se_cloglog)), 2)
  data_cens0_r2$cloglog_up[i]  = round(data_cens0_r2$survival_initial[i] ^ (exp(-qnorm(0.975)*se_cloglog)), 2)
  
  
  # exp(x)/(1+exp(x)) = 1/(exp(-x)+1)
  data_cens0_r2$logit_low[i] = round(1/(exp(-(logit(data_cens0_r2$survival_initial[i], adjust = F) - qnorm(0.975)*se_logit)) + 1), 2)
  data_cens0_r2$logit_up[i]  = round(1/(exp(-(logit(data_cens0_r2$survival_initial[i], adjust = F) + qnorm(0.975)*se_logit)) + 1), 2)
  
  data_cens0_r2$arcsine_low[i] = round(((sin(max(c(0,asin(sqrt(data_cens0_r2$survival_initial[i])) - qnorm(0.975)*se_arcsine))))^2), 2)
  data_cens0_r2$arcsine_up[i]  = round(((sin(min(c(pi/2,asin(sqrt(data_cens0_r2$survival_initial[i])) + qnorm(0.975)*se_arcsine))))^2), 2)
  
  
  #truncation
  if(data_cens0_r2$id_low[i] <= 0) {
    data_cens0_r2$id_low[i] = 0
    data_cens0_r2$truncated_id[i] = 1} 
  
  if(data_cens0_r2$log_low[i] <= 0)  {
    data_cens0_r2$log_low[i] = 0
    data_cens0_r2$truncated_log[i] = 1} 
  
  if(data_cens0_r2$cloglog_low[i] <= 0)  {
    data_cens0_r2$cloglog_low[i] = 0
    data_cens0_r2$truncated_cloglog[i] = 1}
  
  if(data_cens0_r2$logit_low[i] <= 0)  {
    data_cens0_r2$logit_low[i] = 0
    data_cens0_r2$truncated_logit[i] = 1}
  
  if(data_cens0_r2$arcsine_low[i] <= 0)  {
    data_cens0_r2$arcsine_low[i] = 0
    data_cens0_r2$truncated_arcsine[i] = 1}
  
  
  if(data_cens0_r2$id_up[i] >= 1) {
    data_cens0_r2$id_up[i] = 1
    data_cens0_r2$truncated_id[i] = 1} 
  
  if(data_cens0_r2$log_up[i] >= 1)  {
    data_cens0_r2$log_up[i] = 1
    data_cens0_r2$truncated_log[i] = 1} 
  
  if(data_cens0_r2$cloglog_up[i] >= 1)  {
    data_cens0_r2$cloglog_up[i] = 1
    data_cens0_r2$truncated_cloglog[i] = 1} 
  
  if(data_cens0_r2$logit_up[i] >= 1)  {
    data_cens0_r2$logit_up[i] = 1
    data_cens0_r2$truncated_logit[i] = 1} 
  
  if(data_cens0_r2$arcsine_up[i] >= 1)  {
    data_cens0_r2$arcsine_up[i] = 1
    data_cens0_r2$truncated_arcsine[i] = 1} 
  
  # Absolute asymmetry 
  data_cens0_r2$AA_id[i]      = data_cens0_r2$id_up[i] + data_cens0_r2$id_low[i]- 2*data_cens0_r2$survival_rounded[i]
  data_cens0_r2$AA_log[i]     = data_cens0_r2$log_up[i] + data_cens0_r2$log_low[i]- 2*data_cens0_r2$survival_rounded[i]
  data_cens0_r2$AA_cloglog[i] = data_cens0_r2$cloglog_up[i] + data_cens0_r2$cloglog_low[i]- 2*data_cens0_r2$survival_rounded[i]
  data_cens0_r2$AA_logit[i]   = data_cens0_r2$logit_up[i] + data_cens0_r2$logit_low[i]- 2*data_cens0_r2$survival_rounded[i]
  data_cens0_r2$AA_arcsine[i] = data_cens0_r2$arcsine_up[i] + data_cens0_r2$arcsine_low[i]- 2*data_cens0_r2$survival_rounded[i]
  
  # Relative asymmetry 
  data_cens0_r2$RA_id[i]      = (data_cens0_r2$id_up[i] - data_cens0_r2$survival_rounded[i]) / (data_cens0_r2$survival_rounded[i] - data_cens0_r2$id_low[i])
  data_cens0_r2$RA_log[i]     = (data_cens0_r2$log_up[i] - data_cens0_r2$survival_rounded[i]) / (data_cens0_r2$survival_rounded[i] - data_cens0_r2$log_low[i])
  data_cens0_r2$RA_cloglog[i] = (data_cens0_r2$cloglog_up[i] - data_cens0_r2$survival_rounded[i]) / (data_cens0_r2$survival_rounded[i] - data_cens0_r2$cloglog_low[i])
  data_cens0_r2$RA_logit[i]   = (data_cens0_r2$logit_up[i] - data_cens0_r2$survival_rounded[i]) / (data_cens0_r2$survival_rounded[i] - data_cens0_r2$logit_low[i])
  data_cens0_r2$RA_arcsine[i] = (data_cens0_r2$arcsine_up[i] - data_cens0_r2$survival_rounded[i]) / (data_cens0_r2$survival_rounded[i] - data_cens0_r2$arcsine_low[i])
  
  
  # Extraction transformation is identity
  
  if(data_cens0_r2$id_low[i] == 0){
    se_true_id = (data_cens0_r2$id_up[i] - data_cens0_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens0_r2$id_up[i] == 1) {
    se_true_id = (data_cens0_r2$survival_rounded[i] - data_cens0_r2$id_low[i])/qnorm(0.975)
    
  } else {se_true_id = (data_cens0_r2$id_up[i] - data_cens0_r2$id_low[i])/(2*qnorm(0.975))}
  data_cens0_r2$true_id[i] = (se_true_id - se_id)/se_id
  
  
  if(data_cens0_r2$log_low[i] == 0){
    se_id_instead_of_log = (data_cens0_r2$log_up[i] - data_cens0_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens0_r2$log_up[i] == 1){
    se_id_instead_of_log = (data_cens0_r2$survival_rounded[i] - data_cens0_r2$log_low[i])/qnorm(0.975)
    
  } else {se_id_instead_of_log = (data_cens0_r2$log_up[i] - data_cens0_r2$log_low[i])/(2*qnorm(0.975))}
  data_cens0_r2$err_id_instead_of_log[i] = (se_id_instead_of_log - se_id)/se_id
  
  
  if(data_cens0_r2$cloglog_low[i] == 0){
    se_id_instead_of_cloglog = (data_cens0_r2$cloglog_up[i] - data_cens0_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens0_r2$cloglog_up[i] == 1){
    se_id_instead_of_cloglog = (data_cens0_r2$survival_rounded[i] - data_cens0_r2$cloglog_low[i])/qnorm(0.975)
    
  } else {se_id_instead_of_cloglog = (data_cens0_r2$cloglog_up[i] - data_cens0_r2$cloglog_low[i])/(2*qnorm(0.975))}
  data_cens0_r2$err_id_instead_of_cloglog[i] = (se_id_instead_of_cloglog - se_id)/se_id
  
  
  if(data_cens0_r2$logit_low[i] == 0){
    se_id_instead_of_logit = (data_cens0_r2$logit_up[i] - data_cens0_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens0_r2$logit_up[i] == 1){
    se_id_instead_of_logit = (data_cens0_r2$survival_rounded[i] - data_cens0_r2$logit_low[i])/qnorm(0.975)
    
  } else {se_id_instead_of_logit = (data_cens0_r2$logit_up[i] - data_cens0_r2$logit_low[i])/(2*qnorm(0.975))}
  data_cens0_r2$err_id_instead_of_logit[i] = (se_id_instead_of_logit - se_id)/se_id
  
  
  if(data_cens0_r2$arcsine_low[i] == 0){
    se_id_instead_of_arcsine = (data_cens0_r2$arcsine_up[i] - data_cens0_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens0_r2$arcsine_up[i] == 1){
    se_id_instead_of_arcsine = (data_cens0_r2$survival_rounded[i] - data_cens0_r2$arcsine_low[i])/qnorm(0.975)
    
  } else {se_id_instead_of_arcsine = (data_cens0_r2$arcsine_up[i] - data_cens0_r2$arcsine_low[i])/(2*qnorm(0.975))}
  data_cens0_r2$err_id_instead_of_arcsine[i] = (se_id_instead_of_arcsine - se_id)/se_id
  
  
  
  # Extraction transformation is logarithm
  if (data_cens0_r2$log_low[i] == 0) {
    se_true_log = data_cens0_r2$survival_rounded[i]*((log(data_cens0_r2$log_up[i]) - log(data_cens0_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens0_r2$log_up[i] == 1) {
    se_true_log = data_cens0_r2$survival_rounded[i]*((log(data_cens0_r2$survival_rounded[i]) - log(data_cens0_r2$log_low[i]))/(qnorm(0.975)))
    
  } else {se_true_log = data_cens0_r2$survival_rounded[i]*((log(data_cens0_r2$log_up[i]) - log(data_cens0_r2$log_low[i]))/(2*qnorm(0.975)))}
  data_cens0_r2$true_log[i] = (se_true_log - se_id)/se_id
  
  
  if (data_cens0_r2$id_low[i] == 0) {
    se_log_instead_of_id = data_cens0_r2$survival_rounded[i]*((log(data_cens0_r2$id_up[i]) - log(data_cens0_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens0_r2$id_up[i] == 1) {
    se_log_instead_of_id = data_cens0_r2$survival_rounded[i]*((log(data_cens0_r2$survival_rounded[i]) - log(data_cens0_r2$id_low[i]))/(qnorm(0.975)))
    
  } else {se_log_instead_of_id = data_cens0_r2$survival_rounded[i]*((log(data_cens0_r2$id_up[i]) - log(data_cens0_r2$id_low[i]))/(2*qnorm(0.975)))}
  data_cens0_r2$err_log_instead_of_id[i] = (se_log_instead_of_id - se_id)/se_id
  
  
  if (data_cens0_r2$cloglog_low[i] == 0) {
    se_log_instead_of_cloglog = data_cens0_r2$survival_rounded[i]*((log(data_cens0_r2$cloglog_up[i]) - log(data_cens0_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens0_r2$cloglog_up[i] == 1) {
    se_log_instead_of_cloglog = data_cens0_r2$survival_rounded[i]*((log(data_cens0_r2$survival_rounded[i]) - log(data_cens0_r2$cloglog_low[i]))/qnorm(0.975))
    
  } else {se_log_instead_of_cloglog = data_cens0_r2$survival_rounded[i]*((log(data_cens0_r2$cloglog_up[i]) - log(data_cens0_r2$cloglog_low[i]))/(2*qnorm(0.975)))}
  data_cens0_r2$err_log_instead_of_cloglog[i] = (se_log_instead_of_cloglog - se_id)/se_id
  
  
  if (data_cens0_r2$logit_low[i] == 0) {
    se_log_instead_of_logit = data_cens0_r2$survival_rounded[i]*((log(data_cens0_r2$logit_up[i]) - log(data_cens0_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens0_r2$logit_up[i] == 1) {
    se_log_instead_of_logit = data_cens0_r2$survival_rounded[i]*((log(data_cens0_r2$survival_rounded[i]) - log(data_cens0_r2$logit_low[i]))/qnorm(0.975))
    
  } else {se_log_instead_of_logit = data_cens0_r2$survival_rounded[i]*((log(data_cens0_r2$logit_up[i]) - log(data_cens0_r2$logit_low[i]))/(2*qnorm(0.975)))}
  data_cens0_r2$err_log_instead_of_logit[i] = (se_log_instead_of_logit - se_id)/se_id
  
  
  if (data_cens0_r2$arcsine_low[i] == 0) {
    se_log_instead_of_arcsine = data_cens0_r2$survival_rounded[i]*((log(data_cens0_r2$arcsine_up[i]) - log(data_cens0_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens0_r2$arcsine_up[i] == 1) {
    se_log_instead_of_arcsine = data_cens0_r2$survival_rounded[i]*((log(data_cens0_r2$survival_rounded[i]) - log(data_cens0_r2$arcsine_low[i]))/qnorm(0.975))
    
  } else {se_log_instead_of_arcsine = data_cens0_r2$survival_rounded[i]*((log(data_cens0_r2$arcsine_up[i]) - log(data_cens0_r2$arcsine_low[i]))/(2*qnorm(0.975)))}
  data_cens0_r2$err_log_instead_of_arcsine[i] = (se_log_instead_of_arcsine - se_id)/se_id
  
  
  
  # Extraction transformation is complementary log-log
  if (data_cens0_r2$cloglog_low[i] == 0) {
    se_true_cloglog = (- data_cens0_r2$survival_rounded[i]*log(data_cens0_r2$survival_rounded[i]))*((log(-log(data_cens0_r2$survival_rounded[i])) - log(-log(data_cens0_r2$cloglog_up[i])))/qnorm(0.975))
    
  } else if (data_cens0_r2$cloglog_up[i] ==1) {
    se_true_cloglog = (- data_cens0_r2$survival_rounded[i]*log(data_cens0_r2$survival_rounded[i]))*((log(-log(data_cens0_r2$cloglog_low[i])) - log(-log(data_cens0_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else {se_true_cloglog = (- data_cens0_r2$survival_rounded[i]*log(data_cens0_r2$survival_rounded[i]))*((log(-log(data_cens0_r2$cloglog_low[i])) - log(-log(data_cens0_r2$cloglog_up[i])))/(2*qnorm(0.975)))}
  data_cens0_r2$true_cloglog[i] = (se_true_cloglog - se_id)/se_id
  
  
  if (data_cens0_r2$id_low[i] == 0) {
    se_cloglog_instead_of_id = (- data_cens0_r2$survival_rounded[i]*log(data_cens0_r2$survival_rounded[i]))*((log(-log(data_cens0_r2$survival_rounded[i])) - log(-log(data_cens0_r2$id_up[i])))/(qnorm(0.975)))
    
  } else if (data_cens0_r2$id_up[i] == 1) {
    se_cloglog_instead_of_id = (- data_cens0_r2$survival_rounded[i]*log(data_cens0_r2$survival_rounded[i]))*((log(-log(data_cens0_r2$id_low[i])) - log(-log(data_cens0_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else{se_cloglog_instead_of_id = (- data_cens0_r2$survival_rounded[i]*log(data_cens0_r2$survival_rounded[i]))*((log(-log(data_cens0_r2$id_low[i])) - log(-log(data_cens0_r2$id_up[i])))/(2*qnorm(0.975)))}
  data_cens0_r2$err_cloglog_instead_of_id[i] = (se_cloglog_instead_of_id - se_id)/se_id
  
  
  if (data_cens0_r2$log_low[i] == 0) {
    se_cloglog_instead_of_log = (- data_cens0_r2$survival_rounded[i]*log(data_cens0_r2$survival_rounded[i]))*((log(-log(data_cens0_r2$survival_rounded[i])) - log(-log(data_cens0_r2$log_up[i])))/qnorm(0.975))
    
  } else if (data_cens0_r2$log_up[i] ==1) {
    se_cloglog_instead_of_log = (- data_cens0_r2$survival_rounded[i]*log(data_cens0_r2$survival_rounded[i]))*((log(-log(data_cens0_r2$log_low[i])) - log(-log(data_cens0_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else{se_cloglog_instead_of_log = (- data_cens0_r2$survival_rounded[i]*log(data_cens0_r2$survival_rounded[i]))*((log(-log(data_cens0_r2$log_low[i])) - log(-log(data_cens0_r2$log_up[i])))/(2*qnorm(0.975)))}
  data_cens0_r2$err_cloglog_instead_of_log[i] = (se_cloglog_instead_of_log - se_id)/se_id
  
  if (data_cens0_r2$logit_low[i] == 0) {
    se_cloglog_instead_of_logit = (- data_cens0_r2$survival_rounded[i]*log(data_cens0_r2$survival_rounded[i]))*((log(-log(data_cens0_r2$survival_rounded[i])) - log(-log(data_cens0_r2$logit_up[i])))/qnorm(0.975))
    
  } else if (data_cens0_r2$logit_up[i] ==1) {
    se_cloglog_instead_of_logit = (- data_cens0_r2$survival_rounded[i]*log(data_cens0_r2$survival_rounded[i]))*((log(-log(data_cens0_r2$logit_low[i])) - log(-log(data_cens0_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else{se_cloglog_instead_of_logit = (- data_cens0_r2$survival_rounded[i]*log(data_cens0_r2$survival_rounded[i]))*((log(-log(data_cens0_r2$logit_low[i])) - log(-log(data_cens0_r2$logit_up[i])))/(2*qnorm(0.975)))}
  data_cens0_r2$err_cloglog_instead_of_logit[i] = (se_cloglog_instead_of_logit - se_id)/se_id
  
  
  if (data_cens0_r2$arcsine_low[i] == 0) {
    se_cloglog_instead_of_arcsine = (- data_cens0_r2$survival_rounded[i]*log(data_cens0_r2$survival_rounded[i]))*((log(-log(data_cens0_r2$survival_rounded[i])) - log(-log(data_cens0_r2$arcsine_up[i])))/qnorm(0.975))
    
  } else if (data_cens0_r2$arcsine_up[i] ==1) {
    se_cloglog_instead_of_arcsine = (- data_cens0_r2$survival_rounded[i]*log(data_cens0_r2$survival_rounded[i]))*((log(-log(data_cens0_r2$arcsine_low[i])) - log(-log(data_cens0_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else{se_cloglog_instead_of_arcsine = (- data_cens0_r2$survival_rounded[i]*log(data_cens0_r2$survival_rounded[i]))*((log(-log(data_cens0_r2$arcsine_low[i])) - log(-log(data_cens0_r2$arcsine_up[i])))/(2*qnorm(0.975)))}
  data_cens0_r2$err_cloglog_instead_of_arcsine[i] = (se_cloglog_instead_of_arcsine - se_id)/se_id
  
  
  
  # Extraction transformation is logit
  if (data_cens0_r2$logit_low[i] == 0) {
    se_true_logit = ((logit(data_cens0_r2$logit_up[i], adjust = F) - logit(data_cens0_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i])
    
  } else if (data_cens0_r2$logit_up[i] ==1) {
    se_true_logit = ((logit(data_cens0_r2$survival_rounded[i], adjust = F) - logit(data_cens0_r2$logit_low[i], adjust = F))/qnorm(0.975))*data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i])
    
  } else {se_true_logit = ((logit(data_cens0_r2$logit_up[i], adjust = F) - logit(data_cens0_r2$logit_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i])}
  data_cens0_r2$true_logit[i] = (se_true_logit - se_id)/se_id
  
  
  if (data_cens0_r2$id_low[i] == 0) {
    se_logit_instead_of_id = ((logit(data_cens0_r2$id_up[i], adjust = F) - logit(data_cens0_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i])
    
  } else if (data_cens0_r2$id_up[i] == 1) {
    se_logit_instead_of_id = ((logit(data_cens0_r2$survival_rounded[i], adjust = F) - logit(data_cens0_r2$id_low[i], adjust = F))/qnorm(0.975))*data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i])
    
  } else{se_logit_instead_of_id = ((logit(data_cens0_r2$id_up[i], adjust = F) - logit(data_cens0_r2$id_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i])}
  data_cens0_r2$err_logit_instead_of_id[i] = (se_logit_instead_of_id - se_id)/se_id
  
  
  if (data_cens0_r2$log_low[i] == 0) {
    se_logit_instead_of_log = ((logit(data_cens0_r2$log_up[i], adjust = F) - logit(data_cens0_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i])
    
  } else if (data_cens0_r2$log_up[i] ==1) {
    se_logit_instead_of_log = ((logit(data_cens0_r2$survival_rounded[i], adjust = F) - logit(data_cens0_r2$log_low[i], adjust = F))/qnorm(0.975))*data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i])
    
  } else{se_logit_instead_of_log = ((logit(data_cens0_r2$log_up[i], adjust = F) - logit(data_cens0_r2$log_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i])}
  data_cens0_r2$err_logit_instead_of_log[i] = (se_logit_instead_of_log - se_id)/se_id
  
  
  if (data_cens0_r2$cloglog_low[i] == 0) {
    se_logit_instead_of_cloglog = ((logit(data_cens0_r2$cloglog_up[i], adjust = F) - logit(data_cens0_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i])
    
  } else if (data_cens0_r2$cloglog_up[i] ==1) {
    se_logit_instead_of_cloglog = ((logit(data_cens0_r2$survival_rounded[i], adjust = F) - logit(data_cens0_r2$cloglog_low[i], adjust = F))/qnorm(0.975))*data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i])
    
  } else{se_logit_instead_of_cloglog = ((logit(data_cens0_r2$cloglog_up[i], adjust = F) - logit(data_cens0_r2$cloglog_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i])}
  data_cens0_r2$err_logit_instead_of_cloglog[i] = (se_logit_instead_of_cloglog - se_id)/se_id
  
  
  if (data_cens0_r2$arcsine_low[i] == 0) {
    se_logit_instead_of_arcsine = ((logit(data_cens0_r2$arcsine_up[i], adjust = F) - logit(data_cens0_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i])
    
  } else if (data_cens0_r2$arcsine_up[i] ==1) {
    se_logit_instead_of_arcsine = ((logit(data_cens0_r2$survival_rounded[i], adjust = F) - logit(data_cens0_r2$arcsine_low[i], adjust = F))/qnorm(0.975))*data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i])
    
  } else{se_logit_instead_of_arcsine = ((logit(data_cens0_r2$arcsine_up[i], adjust = F) - logit(data_cens0_r2$arcsine_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i])}
  data_cens0_r2$err_logit_instead_of_arcsine[i] = (se_logit_instead_of_arcsine - se_id)/se_id
  
  
  # Extraction transformation is arcsine
  if (data_cens0_r2$arcsine_low[i] == 0) {
    se_true_arcsine = ((asin(sqrt(data_cens0_r2$arcsine_up[i])) - asin(sqrt(data_cens0_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i]))
    
  } else if (data_cens0_r2$arcsine_up[i] == 1) {
    se_true_arcsine = ((asin(sqrt(data_cens0_r2$survival_rounded[i])) - asin(sqrt(data_cens0_r2$arcsine_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i]))
    
  } else {se_true_arcsine = ((asin(sqrt(data_cens0_r2$arcsine_up[i])) - asin(sqrt(data_cens0_r2$arcsine_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i]))}
  data_cens0_r2$true_arcsine[i] = (se_true_arcsine - se_id)/se_id
  
  
  if (data_cens0_r2$id_low[i] == 0) {
    se_arcsine_instead_of_id = ((asin(sqrt(data_cens0_r2$id_up[i])) - asin(sqrt(data_cens0_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i]))
    
  } else if (data_cens0_r2$id_up[i] == 1) {
    se_arcsine_instead_of_id = ((asin(sqrt(data_cens0_r2$survival_rounded[i])) - asin(sqrt(data_cens0_r2$id_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i]))
    
  } else{se_arcsine_instead_of_id = ((asin(sqrt(data_cens0_r2$id_up[i])) - asin(sqrt(data_cens0_r2$id_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i]))}
  data_cens0_r2$err_arcsine_instead_of_id[i] = (se_arcsine_instead_of_id - se_id)/se_id
  
  
  if (data_cens0_r2$log_low[i] == 0) {
    se_arcsine_instead_of_log = ((asin(sqrt(data_cens0_r2$log_up[i])) - asin(sqrt(data_cens0_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i]))
    
  } else if (data_cens0_r2$log_up[i] ==1) {
    se_arcsine_instead_of_log = ((asin(sqrt(data_cens0_r2$survival_rounded[i])) - asin(sqrt(data_cens0_r2$log_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i]))
    
  } else{se_arcsine_instead_of_log = ((asin(sqrt(data_cens0_r2$log_up[i])) - asin(sqrt(data_cens0_r2$log_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i]))}
  data_cens0_r2$err_arcsine_instead_of_log[i] = (se_arcsine_instead_of_log - se_id)/se_id
  
  
  if (data_cens0_r2$cloglog_low[i] == 0) {
    se_arcsine_instead_of_cloglog = ((asin(sqrt(data_cens0_r2$cloglog_up[i])) - asin(sqrt(data_cens0_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i]))
    
  } else if (data_cens0_r2$arcsine_up[i] ==1) {
    se_arcsine_instead_of_cloglog = ((asin(sqrt(data_cens0_r2$survival_rounded[i])) - asin(sqrt(data_cens0_r2$cloglog_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i]))
    
  } else{se_arcsine_instead_of_cloglog = ((asin(sqrt(data_cens0_r2$cloglog_up[i])) - asin(sqrt(data_cens0_r2$cloglog_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i]))}
  data_cens0_r2$err_arcsine_instead_of_cloglog[i] = (se_arcsine_instead_of_cloglog - se_id)/se_id
  
  
  if (data_cens0_r2$logit_low[i] == 0) {
    se_arcsine_instead_of_logit = ((asin(sqrt(data_cens0_r2$logit_up[i])) - asin(sqrt(data_cens0_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i]))
    
  } else if (data_cens0_r2$logit_up[i] ==1) {
    se_arcsine_instead_of_logit = ((asin(sqrt(data_cens0_r2$survival_rounded[i])) - asin(sqrt(data_cens0_r2$logit_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i]))
    
  } else{se_arcsine_instead_of_logit = ((asin(sqrt(data_cens0_r2$logit_up[i])) - asin(sqrt(data_cens0_r2$logit_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens0_r2$survival_rounded[i]*(1-data_cens0_r2$survival_rounded[i]))}
  data_cens0_r2$err_arcsine_instead_of_logit[i] = (se_arcsine_instead_of_logit - se_id)/se_id
  
  
}
# Transform extraction errors into percentages by multiplying them by 100 
for(i in names(data_cens0_r2)) {
  if (startsWith(i, "err_") | startsWith(i, "true_")) {
    data_cens0_r2[[i]] <- data_cens0_r2[[i]] * 100
  }
}
save(data_cens0_r2,  file = "your_directory/datacens0_extr_erorr_r2.Rdata")


# 10% censoring  ----
for (i in 1:nrow(data_cens10_r2)) {
  
  print(i)
  ## standard error on different scales
  data_cens10_r2$se_id_initial[i] = sqrt((data_cens10_r2$survival_initial[i]*(1-data_cens10_r2$survival_initial[i]))/data_cens10_r2$sample_size[i])
  data_cens10_r2$se_id_infl[i] = data_cens10_r2$se_id_initial[i]* data_cens10_r2$inflation_coeff[i]
  
  se_id = data_cens10_r2$se_id_infl[i]
  se_log = se_id/data_cens10_r2$survival_initial[i]
  se_cloglog = - se_id / (data_cens10_r2$survival_initial[i]*log(data_cens10_r2$survival_initial[i]))
  se_logit = se_id / (data_cens10_r2$survival_initial[i]*(1-data_cens10_r2$survival_initial[i]))
  se_arcsine =  se_id / (2*sqrt((data_cens10_r2$survival_initial[i])*(1-data_cens10_r2$survival_initial[i])))
  
  
  ## 95% confidence interval of survival computed using the different methods
  # If the lower bound is equal to 0 or the upper bound equal to 1 whether due to rounding or truncation, we consider that the interval is truncated/trimmed
  #Thus, we don't need the bounds not rounded for the following steps, we only need the rounded bounds
  
  data_cens10_r2$id_low[i] = round(data_cens10_r2$survival_initial[i] - qnorm(0.975) * se_id, 2)
  data_cens10_r2$id_up[i]  = round(data_cens10_r2$survival_initial[i] + qnorm(0.975) * se_id, 2)
  
  data_cens10_r2$log_low[i] = round(data_cens10_r2$survival_initial[i] * exp(-qnorm(0.975)*se_log), 2)
  data_cens10_r2$log_up[i]  = round(data_cens10_r2$survival_initial[i] * exp(qnorm(0.975)*se_log), 2)
  
  data_cens10_r2$cloglog_low[i] = round(data_cens10_r2$survival_initial[i] ^ (exp(qnorm(0.975)*se_cloglog)), 2)
  data_cens10_r2$cloglog_up[i]  = round(data_cens10_r2$survival_initial[i] ^ (exp(-qnorm(0.975)*se_cloglog)), 2)
  
  data_cens10_r2$logit_low[i] = round(1/(exp(-(logit(data_cens10_r2$survival_initial[i], adjust = F) - qnorm(0.975)*se_logit)) + 1), 2)
  data_cens10_r2$logit_up[i]  = round(1/(exp(-(logit(data_cens10_r2$survival_initial[i], adjust = F) + qnorm(0.975)*se_logit)) + 1), 2)
  
  data_cens10_r2$arcsine_low[i] = round(((sin(max(c(0,asin(sqrt(data_cens10_r2$survival_initial[i])) - qnorm(0.975)*se_arcsine))))^2), 2)
  data_cens10_r2$arcsine_up[i]  = round(((sin(min(c(pi/2,asin(sqrt(data_cens10_r2$survival_initial[i])) + qnorm(0.975)*se_arcsine))))^2), 4)
  
  #truncation
  if(data_cens10_r2$id_low[i] <= 0) {
    data_cens10_r2$id_low[i] = 0
    data_cens10_r2$truncated_id[i] = 1} 
  
  if(data_cens10_r2$log_low[i] <= 0)  {
    data_cens10_r2$log_low[i] = 0
    data_cens10_r2$truncated_log[i] = 1} 
  
  if(data_cens10_r2$cloglog_low[i] <= 0)  {
    data_cens10_r2$cloglog_low[i] = 0
    data_cens10_r2$truncated_cloglog[i] = 1}
  
  if(data_cens10_r2$logit_low[i] <= 0)  {
    data_cens10_r2$logit_low[i] = 0
    data_cens10_r2$truncated_logit[i] = 1}
  
  if(data_cens10_r2$arcsine_low[i] <= 0)  {
    data_cens10_r2$arcsine_low[i] = 0
    data_cens10_r2$truncated_arcsine[i] = 1}
  
  
  if(data_cens10_r2$id_up[i] >= 1) {
    data_cens10_r2$id_up[i] = 1
    data_cens10_r2$truncated_id[i] = 1} 
  
  if(data_cens10_r2$log_up[i] >= 1)  {
    data_cens10_r2$log_up[i] = 1
    data_cens10_r2$truncated_log[i] = 1} 
  
  if(data_cens10_r2$cloglog_up[i] >= 1)  {
    data_cens10_r2$cloglog_up[i] = 1
    data_cens10_r2$truncated_cloglog[i] = 1} 
  
  if(data_cens10_r2$logit_up[i] >= 1)  {
    data_cens10_r2$logit_up[i] = 1
    data_cens10_r2$truncated_logit[i] = 1} 
  
  if(data_cens10_r2$arcsine_up[i] >= 1)  {
    data_cens10_r2$arcsine_up[i] = 1
    data_cens10_r2$truncated_arcsine[i] = 1} 
  
  # Absolute asymmetry 
  data_cens10_r2$AA_id[i]      = data_cens10_r2$id_up[i] + data_cens10_r2$id_low[i]- 2*data_cens10_r2$survival_rounded[i]
  data_cens10_r2$AA_log[i]     = data_cens10_r2$log_up[i] + data_cens10_r2$log_low[i]- 2*data_cens10_r2$survival_rounded[i]
  data_cens10_r2$AA_cloglog[i] = data_cens10_r2$cloglog_up[i] + data_cens10_r2$cloglog_low[i]- 2*data_cens10_r2$survival_rounded[i]
  data_cens10_r2$AA_logit[i]   = data_cens10_r2$logit_up[i] + data_cens10_r2$logit_low[i]- 2*data_cens10_r2$survival_rounded[i]
  data_cens10_r2$AA_arcsine[i] = data_cens10_r2$arcsine_up[i] + data_cens10_r2$arcsine_low[i]- 2*data_cens10_r2$survival_rounded[i]
  
  # Relative asymmetry 
  data_cens10_r2$RA_id[i]      = (data_cens10_r2$id_up[i] - data_cens10_r2$survival_rounded[i]) / (data_cens10_r2$survival_rounded[i] - data_cens10_r2$id_low[i])
  data_cens10_r2$RA_log[i]     = (data_cens10_r2$log_up[i] - data_cens10_r2$survival_rounded[i]) / (data_cens10_r2$survival_rounded[i] - data_cens10_r2$log_low[i])
  data_cens10_r2$RA_cloglog[i] = (data_cens10_r2$cloglog_up[i] - data_cens10_r2$survival_rounded[i]) / (data_cens10_r2$survival_rounded[i] - data_cens10_r2$cloglog_low[i])
  data_cens10_r2$RA_logit[i]   = (data_cens10_r2$logit_up[i] - data_cens10_r2$survival_rounded[i]) / (data_cens10_r2$survival_rounded[i] - data_cens10_r2$logit_low[i])
  data_cens10_r2$RA_arcsine[i] = (data_cens10_r2$arcsine_up[i] - data_cens10_r2$survival_rounded[i]) / (data_cens10_r2$survival_rounded[i] - data_cens10_r2$arcsine_low[i])
  
  
  # Extraction transformation is identity
  
  if(data_cens10_r2$id_low[i] == 0){
    se_true_id = (data_cens10_r2$id_up[i] - data_cens10_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens10_r2$id_up[i] == 1) {
    se_true_id = (data_cens10_r2$survival_rounded[i] - data_cens10_r2$id_low[i])/qnorm(0.975)
    
  } else {se_true_id = (data_cens10_r2$id_up[i] - data_cens10_r2$id_low[i])/(2*qnorm(0.975))}
  data_cens10_r2$true_id[i] = (se_true_id - se_id)/se_id
  
  
  if(data_cens10_r2$log_low[i] == 0){
    se_id_instead_of_log = (data_cens10_r2$log_up[i] - data_cens10_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens10_r2$log_up[i] == 1){
    se_id_instead_of_log = (data_cens10_r2$survival_rounded[i] - data_cens10_r2$log_low[i])/qnorm(0.975)
    
  } else {se_id_instead_of_log = (data_cens10_r2$log_up[i] - data_cens10_r2$log_low[i])/(2*qnorm(0.975))}
  data_cens10_r2$err_id_instead_of_log[i] = (se_id_instead_of_log - se_id)/se_id
  
  
  if(data_cens10_r2$cloglog_low[i] == 0){
    se_id_instead_of_cloglog = (data_cens10_r2$cloglog_up[i] - data_cens10_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens10_r2$cloglog_up[i] == 1){
    se_id_instead_of_cloglog = (data_cens10_r2$survival_rounded[i] - data_cens10_r2$cloglog_low[i])/qnorm(0.975)
    
  } else {se_id_instead_of_cloglog = (data_cens10_r2$cloglog_up[i] - data_cens10_r2$cloglog_low[i])/(2*qnorm(0.975))}
  data_cens10_r2$err_id_instead_of_cloglog[i] = (se_id_instead_of_cloglog - se_id)/se_id
  
  
  if(data_cens10_r2$logit_low[i] == 0){
    se_id_instead_of_logit = (data_cens10_r2$logit_up[i] - data_cens10_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens10_r2$logit_up[i] == 1){
    se_id_instead_of_logit = (data_cens10_r2$survival_rounded[i] - data_cens10_r2$logit_low[i])/qnorm(0.975)
    
  } else {se_id_instead_of_logit = (data_cens10_r2$logit_up[i] - data_cens10_r2$logit_low[i])/(2*qnorm(0.975))}
  data_cens10_r2$err_id_instead_of_logit[i] = (se_id_instead_of_logit - se_id)/se_id
  
  
  if(data_cens10_r2$arcsine_low[i] == 0){
    se_id_instead_of_arcsine = (data_cens10_r2$arcsine_up[i] - data_cens10_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens10_r2$arcsine_up[i] == 1){
    se_id_instead_of_arcsine = (data_cens10_r2$survival_rounded[i] - data_cens10_r2$arcsine_low[i])/qnorm(0.975)
    
  } else {se_id_instead_of_arcsine = (data_cens10_r2$arcsine_up[i] - data_cens10_r2$arcsine_low[i])/(2*qnorm(0.975))}
  data_cens10_r2$err_id_instead_of_arcsine[i] = (se_id_instead_of_arcsine - se_id)/se_id
  
  
  
  # Extraction transformation is logarithm
  
  if (data_cens10_r2$log_low[i] == 0) {
    se_true_log = data_cens10_r2$survival_rounded[i]*((log(data_cens10_r2$log_up[i]) - log(data_cens10_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens10_r2$log_up[i] == 1) {
    se_true_log = data_cens10_r2$survival_rounded[i]*((log(data_cens10_r2$survival_rounded[i]) - log(data_cens10_r2$log_low[i]))/(qnorm(0.975)))
    
  } else {se_true_log = data_cens10_r2$survival_rounded[i]*((log(data_cens10_r2$log_up[i]) - log(data_cens10_r2$log_low[i]))/(2*qnorm(0.975)))}
  data_cens10_r2$true_log[i] = (se_true_log - se_id)/se_id
  
  
  if (data_cens10_r2$id_low[i] == 0) {
    se_log_instead_of_id = data_cens10_r2$survival_rounded[i]*((log(data_cens10_r2$id_up[i]) - log(data_cens10_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens10_r2$id_up[i] == 1) {
    se_log_instead_of_id = data_cens10_r2$survival_rounded[i]*((log(data_cens10_r2$survival_rounded[i]) - log(data_cens10_r2$id_low[i]))/(qnorm(0.975)))
    
  } else {se_log_instead_of_id = data_cens10_r2$survival_rounded[i]*((log(data_cens10_r2$id_up[i]) - log(data_cens10_r2$id_low[i]))/(2*qnorm(0.975)))}
  data_cens10_r2$err_log_instead_of_id[i] = (se_log_instead_of_id - se_id)/se_id
  
  
  if (data_cens10_r2$cloglog_low[i] == 0) {
    se_log_instead_of_cloglog = data_cens10_r2$survival_rounded[i]*((log(data_cens10_r2$cloglog_up[i]) - log(data_cens10_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens10_r2$cloglog_up[i] == 1) {
    se_log_instead_of_cloglog = data_cens10_r2$survival_rounded[i]*((log(data_cens10_r2$survival_rounded[i]) - log(data_cens10_r2$cloglog_low[i]))/qnorm(0.975))
    
  } else {se_log_instead_of_cloglog = data_cens10_r2$survival_rounded[i]*((log(data_cens10_r2$cloglog_up[i]) - log(data_cens10_r2$cloglog_low[i]))/(2*qnorm(0.975)))}
  data_cens10_r2$err_log_instead_of_cloglog[i] = (se_log_instead_of_cloglog - se_id)/se_id
  
  
  if (data_cens10_r2$logit_low[i] == 0) {
    se_log_instead_of_logit = data_cens10_r2$survival_rounded[i]*((log(data_cens10_r2$logit_up[i]) - log(data_cens10_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens10_r2$logit_up[i] == 1) {
    se_log_instead_of_logit = data_cens10_r2$survival_rounded[i]*((log(data_cens10_r2$survival_rounded[i]) - log(data_cens10_r2$logit_low[i]))/qnorm(0.975))
    
  } else {se_log_instead_of_logit = data_cens10_r2$survival_rounded[i]*((log(data_cens10_r2$logit_up[i]) - log(data_cens10_r2$logit_low[i]))/(2*qnorm(0.975)))}
  data_cens10_r2$err_log_instead_of_logit[i] = (se_log_instead_of_logit - se_id)/se_id
  
  
  if (data_cens10_r2$arcsine_low[i] == 0) {
    se_log_instead_of_arcsine = data_cens10_r2$survival_rounded[i]*((log(data_cens10_r2$arcsine_up[i]) - log(data_cens10_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens10_r2$arcsine_up[i] == 1) {
    se_log_instead_of_arcsine = data_cens10_r2$survival_rounded[i]*((log(data_cens10_r2$survival_rounded[i]) - log(data_cens10_r2$arcsine_low[i]))/qnorm(0.975))
    
  } else {se_log_instead_of_arcsine = data_cens10_r2$survival_rounded[i]*((log(data_cens10_r2$arcsine_up[i]) - log(data_cens10_r2$arcsine_low[i]))/(2*qnorm(0.975)))}
  data_cens10_r2$err_log_instead_of_arcsine[i] = (se_log_instead_of_arcsine - se_id)/se_id
  
  
  
  # Extraction transformation is complementary log-log
  if (data_cens10_r2$cloglog_low[i] == 0) {
    se_true_cloglog = (- data_cens10_r2$survival_rounded[i]*log(data_cens10_r2$survival_rounded[i]))*((log(-log(data_cens10_r2$survival_rounded[i])) - log(-log(data_cens10_r2$cloglog_up[i])))/qnorm(0.975))
    
  } else if (data_cens10_r2$cloglog_up[i] ==1) {
    se_true_cloglog = (- data_cens10_r2$survival_rounded[i]*log(data_cens10_r2$survival_rounded[i]))*((log(-log(data_cens10_r2$cloglog_low[i])) - log(-log(data_cens10_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else {se_true_cloglog = (- data_cens10_r2$survival_rounded[i]*log(data_cens10_r2$survival_rounded[i]))*((log(-log(data_cens10_r2$cloglog_low[i])) - log(-log(data_cens10_r2$cloglog_up[i])))/(2*qnorm(0.975)))}
  data_cens10_r2$true_cloglog[i] = (se_true_cloglog - se_id)/se_id
  
  
  if (data_cens10_r2$id_low[i] == 0) {
    se_cloglog_instead_of_id = - data_cens10_r2$survival_rounded[i]*log(data_cens10_r2$survival_rounded[i])*((log(-log(data_cens10_r2$survival_rounded[i])) - log(-log(data_cens10_r2$id_up[i])))/(qnorm(0.975)))
    
  } else if (data_cens10_r2$id_up[i] == 1) {
    se_cloglog_instead_of_id = (- data_cens10_r2$survival_rounded[i]*log(data_cens10_r2$survival_rounded[i]))*((log(-log(data_cens10_r2$id_low[i])) - log(-log(data_cens10_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else{se_cloglog_instead_of_id = (- data_cens10_r2$survival_rounded[i]*log(data_cens10_r2$survival_rounded[i]))*((log(-log(data_cens10_r2$id_low[i])) - log(-log(data_cens10_r2$id_up[i])))/(2*qnorm(0.975)))}
  data_cens10_r2$err_cloglog_instead_of_id[i] = (se_cloglog_instead_of_id - se_id)/se_id
  
  
  if (data_cens10_r2$log_low[i] == 0) {
    se_cloglog_instead_of_log = (- data_cens10_r2$survival_rounded[i]*log(data_cens10_r2$survival_rounded[i]))*((log(-log(data_cens10_r2$survival_rounded[i])) - log(-log(data_cens10_r2$log_up[i])))/qnorm(0.975))
    
  } else if (data_cens10_r2$log_up[i] ==1) {
    se_cloglog_instead_of_log = (- data_cens10_r2$survival_rounded[i]*log(data_cens10_r2$survival_rounded[i]))*((log(-log(data_cens10_r2$log_low[i])) - log(-log(data_cens10_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else{se_cloglog_instead_of_log = (- data_cens10_r2$survival_rounded[i]*log(data_cens10_r2$survival_rounded[i]))*((log(-log(data_cens10_r2$log_low[i])) - log(-log(data_cens10_r2$log_up[i])))/(2*qnorm(0.975)))}
  data_cens10_r2$err_cloglog_instead_of_log[i] = (se_cloglog_instead_of_log - se_id)/se_id
  
  
  if (data_cens10_r2$logit_low[i] == 0) {
    se_cloglog_instead_of_logit = (- data_cens10_r2$survival_rounded[i]*log(data_cens10_r2$survival_rounded[i]))*((log(-log(data_cens10_r2$survival_rounded[i])) - log(-log(data_cens10_r2$logit_up[i])))/qnorm(0.975))
    
  } else if (data_cens10_r2$logit_up[i] ==1) {
    se_cloglog_instead_of_logit = (- data_cens10_r2$survival_rounded[i]*log(data_cens10_r2$survival_rounded[i]))*((log(-log(data_cens10_r2$logit_low[i])) - log(-log(data_cens10_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else{se_cloglog_instead_of_logit = (- data_cens10_r2$survival_rounded[i]*log(data_cens10_r2$survival_rounded[i]))*((log(-log(data_cens10_r2$logit_low[i])) - log(-log(data_cens10_r2$logit_up[i])))/(2*qnorm(0.975)))}
  data_cens10_r2$err_cloglog_instead_of_logit[i] = (se_cloglog_instead_of_logit - se_id)/se_id
  
  
  if (data_cens10_r2$arcsine_low[i] == 0) {
    se_cloglog_instead_of_arcsine = (- data_cens10_r2$survival_rounded[i]*log(data_cens10_r2$survival_rounded[i]))*((log(-log(data_cens10_r2$survival_rounded[i])) - log(-log(data_cens10_r2$arcsine_up[i])))/qnorm(0.975))
    
  } else if (data_cens10_r2$arcsine_up[i] ==1) {
    se_cloglog_instead_of_arcsine = (- data_cens10_r2$survival_rounded[i]*log(data_cens10_r2$survival_rounded[i]))*((log(-log(data_cens10_r2$arcsine_low[i])) - log(-log(data_cens10_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else{se_cloglog_instead_of_arcsine = (- data_cens10_r2$survival_rounded[i]*log(data_cens10_r2$survival_rounded[i]))*((log(-log(data_cens10_r2$arcsine_low[i])) - log(-log(data_cens10_r2$arcsine_up[i])))/(2*qnorm(0.975)))}
  data_cens10_r2$err_cloglog_instead_of_arcsine[i] = (se_cloglog_instead_of_arcsine - se_id)/se_id
  
  
  
  # Extraction transformation is logit
  if (data_cens10_r2$logit_low[i] == 0) {
    se_true_logit = ((logit(data_cens10_r2$logit_up[i], adjust = F) - logit(data_cens10_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i])
    
  } else if (data_cens10_r2$logit_up[i] ==1) {
    se_true_logit = ((logit(data_cens10_r2$survival_rounded[i], adjust = F) - logit(data_cens10_r2$logit_low[i], adjust = F))/qnorm(0.975))*data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i])
    
  } else {se_true_logit = ((logit(data_cens10_r2$logit_up[i], adjust = F) - logit(data_cens10_r2$logit_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i])}
  data_cens10_r2$true_logit[i] = (se_true_logit - se_id)/se_id
  
  
  if (data_cens10_r2$id_low[i] == 0) {
    se_logit_instead_of_id = ((logit(data_cens10_r2$id_up[i], adjust = F) - logit(data_cens10_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i])
    
  } else if (data_cens10_r2$id_up[i] == 1) {
    se_logit_instead_of_id = ((logit(data_cens10_r2$survival_rounded[i], adjust = F) - logit(data_cens10_r2$id_low[i], adjust = F))/qnorm(0.975))*data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i])
    
  } else{se_logit_instead_of_id = ((logit(data_cens10_r2$id_up[i], adjust = F) - logit(data_cens10_r2$id_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i])}
  data_cens10_r2$err_logit_instead_of_id[i] = (se_logit_instead_of_id - se_id)/se_id
  
  
  if (data_cens10_r2$log_low[i] == 0) {
    se_logit_instead_of_log = ((logit(data_cens10_r2$log_up[i], adjust = F) - logit(data_cens10_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i])
    
  } else if (data_cens10_r2$log_up[i] ==1) {
    se_logit_instead_of_log = ((logit(data_cens10_r2$survival_rounded[i], adjust = F) - logit(data_cens10_r2$log_low[i], adjust = F))/qnorm(0.975))*data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i])
    
  } else{se_logit_instead_of_log = ((logit(data_cens10_r2$log_up[i], adjust = F) - logit(data_cens10_r2$log_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i])}
  data_cens10_r2$err_logit_instead_of_log[i] = (se_logit_instead_of_log - se_id)/se_id
  
  
  if (data_cens10_r2$cloglog_low[i] == 0) {
    se_logit_instead_of_cloglog = ((logit(data_cens10_r2$cloglog_up[i], adjust = F) - logit(data_cens10_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i])
    
  } else if (data_cens10_r2$cloglog_up[i] ==1) {
    se_logit_instead_of_cloglog = ((logit(data_cens10_r2$survival_rounded[i], adjust = F) - logit(data_cens10_r2$cloglog_low[i], adjust = F))/qnorm(0.975))*data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i])
    
  } else{se_logit_instead_of_cloglog = ((logit(data_cens10_r2$cloglog_up[i], adjust = F) - logit(data_cens10_r2$cloglog_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i])}
  data_cens10_r2$err_logit_instead_of_cloglog[i] = (se_logit_instead_of_cloglog - se_id)/se_id
  
  
  if (data_cens10_r2$arcsine_low[i] == 0) {
    se_logit_instead_of_arcsine = ((logit(data_cens10_r2$arcsine_up[i], adjust = F) - logit(data_cens10_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i])
    
  } else if (data_cens10_r2$arcsine_up[i] ==1) {
    se_logit_instead_of_arcsine = ((logit(data_cens10_r2$survival_rounded[i], adjust = F) - logit(data_cens10_r2$arcsine_low[i], adjust = F))/qnorm(0.975))*data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i])
    
  } else{se_logit_instead_of_arcsine = ((logit(data_cens10_r2$arcsine_up[i], adjust = F) - logit(data_cens10_r2$arcsine_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i])}
  data_cens10_r2$err_logit_instead_of_arcsine[i] = (se_logit_instead_of_arcsine - se_id)/se_id
  
  
  # Extraction transformation is arcsine
  if (data_cens10_r2$arcsine_low[i] == 0) {
    se_true_arcsine = ((asin(sqrt(data_cens10_r2$arcsine_up[i])) - asin(sqrt(data_cens10_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i]))
    
  } else if (data_cens10_r2$arcsine_up[i] == 1) {
    se_true_arcsine = ((asin(sqrt(data_cens10_r2$survival_rounded[i])) - asin(sqrt(data_cens10_r2$arcsine_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i]))
    
  } else {se_true_arcsine = ((asin(sqrt(data_cens10_r2$arcsine_up[i])) - asin(sqrt(data_cens10_r2$arcsine_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i]))}
  data_cens10_r2$true_arcsine[i] = (se_true_arcsine - se_id)/se_id
  
  
  if (data_cens10_r2$id_low[i] == 0) {
    se_arcsine_instead_of_id = ((asin(sqrt(data_cens10_r2$id_up[i])) - asin(sqrt(data_cens10_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i]))
    
  } else if (data_cens10_r2$id_up[i] == 1) {
    se_arcsine_instead_of_id = ((asin(sqrt(data_cens10_r2$survival_rounded[i])) - asin(sqrt(data_cens10_r2$id_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i]))
    
  } else{se_arcsine_instead_of_id = ((asin(sqrt(data_cens10_r2$id_up[i])) - asin(sqrt(data_cens10_r2$id_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i]))}
  data_cens10_r2$err_arcsine_instead_of_id[i] = (se_arcsine_instead_of_id - se_id)/se_id
  
  
  if (data_cens10_r2$log_low[i] == 0) {
    se_arcsine_instead_of_log = ((asin(sqrt(data_cens10_r2$log_up[i])) - asin(sqrt(data_cens10_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i]))
    
  } else if (data_cens10_r2$log_up[i] ==1) {
    se_arcsine_instead_of_log = ((asin(sqrt(data_cens10_r2$survival_rounded[i])) - asin(sqrt(data_cens10_r2$log_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i]))
    
  } else{se_arcsine_instead_of_log = ((asin(sqrt(data_cens10_r2$log_up[i])) - asin(sqrt(data_cens10_r2$log_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i]))}
  data_cens10_r2$err_arcsine_instead_of_log[i] = (se_arcsine_instead_of_log - se_id)/se_id
  
  
  if (data_cens10_r2$cloglog_low[i] == 0) {
    se_arcsine_instead_of_cloglog = ((asin(sqrt(data_cens10_r2$cloglog_up[i])) - asin(sqrt(data_cens10_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i]))
    
  } else if (data_cens10_r2$arcsine_up[i] ==1) {
    se_arcsine_instead_of_cloglog = ((asin(sqrt(data_cens10_r2$survival_rounded[i])) - asin(sqrt(data_cens10_r2$cloglog_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i]))
    
  } else{se_arcsine_instead_of_cloglog = ((asin(sqrt(data_cens10_r2$cloglog_up[i])) - asin(sqrt(data_cens10_r2$cloglog_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i]))}
  data_cens10_r2$err_arcsine_instead_of_cloglog[i] = (se_arcsine_instead_of_cloglog - se_id)/se_id
  
  
  if (data_cens10_r2$logit_low[i] == 0) {
    se_arcsine_instead_of_logit = ((asin(sqrt(data_cens10_r2$logit_up[i])) - asin(sqrt(data_cens10_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i]))
    
  } else if (data_cens10_r2$logit_up[i] ==1) {
    se_arcsine_instead_of_logit = ((asin(sqrt(data_cens10_r2$survival_rounded[i])) - asin(sqrt(data_cens10_r2$logit_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i]))
    
  } else{se_arcsine_instead_of_logit = ((asin(sqrt(data_cens10_r2$logit_up[i])) - asin(sqrt(data_cens10_r2$logit_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens10_r2$survival_rounded[i]*(1-data_cens10_r2$survival_rounded[i]))}
  data_cens10_r2$err_arcsine_instead_of_logit[i] = (se_arcsine_instead_of_logit - se_id)/se_id
  
  
}
# Transform extraction errors into percentages by multiplying them by 100 
for(i in names(data_cens10_r2)) {
  if (startsWith(i, "err_") | startsWith(i, "true_")) {
    data_cens10_r2[[i]] <- data_cens10_r2[[i]] * 100
  }
}
save(data_cens10_r2, file = "your_directory/datacens10_extr_erorr_r2.Rdata")



# 25% censoring  ----

for (i in 1:nrow(data_cens25_r2)) {
  
  print(i)
  ## standard error on different scales
  data_cens25_r2$se_id_initial[i] = sqrt((data_cens25_r2$survival_initial[i]*(1-data_cens25_r2$survival_initial[i]))/data_cens25_r2$sample_size[i])
  data_cens25_r2$se_id_infl[i] = data_cens25_r2$se_id_initial[i]* data_cens25_r2$inflation_coeff[i]
  
  se_id = data_cens25_r2$se_id_infl[i]
  se_log = se_id/data_cens25_r2$survival_initial[i]
  se_cloglog = - se_id / (data_cens25_r2$survival_initial[i]*log(data_cens25_r2$survival_initial[i]))
  se_logit = se_id / (data_cens25_r2$survival_initial[i]*(1-data_cens25_r2$survival_initial[i]))
  se_arcsine =  se_id / (2*sqrt((data_cens25_r2$survival_initial[i])*(1-data_cens25_r2$survival_initial[i])))
  
  
  ## 95% confidence interval of survival computed using the different methods
  # If the lower bound is equal to 0 or the upper bound equal to 1 whether due to rounding or truncation, we consider that the interval is truncated/trimmed
  #Thus, we don't need the bounds not rounded for the following steps, we only need the rounded bounds
  
  data_cens25_r2$id_low[i] = round(data_cens25_r2$survival_initial[i] - qnorm(0.975) * se_id, 2)
  data_cens25_r2$id_up[i]  = round(data_cens25_r2$survival_initial[i] + qnorm(0.975) * se_id, 2)
  
  data_cens25_r2$log_low[i] = round(data_cens25_r2$survival_initial[i] * exp(-qnorm(0.975)*se_log), 2)
  data_cens25_r2$log_up[i]  = round(data_cens25_r2$survival_initial[i] * exp(qnorm(0.975)*se_log), 2)
  
  data_cens25_r2$cloglog_low[i] = round(data_cens25_r2$survival_initial[i] ^ (exp(qnorm(0.975)*se_cloglog)), 2)
  data_cens25_r2$cloglog_up[i]  = round(data_cens25_r2$survival_initial[i] ^ (exp(-qnorm(0.975)*se_cloglog)), 2)
  
  data_cens25_r2$logit_low[i] = round(1/(exp(-(logit(data_cens25_r2$survival_initial[i], adjust = F) - qnorm(0.975)*se_logit)) + 1), 2)
  data_cens25_r2$logit_up[i]  = round(1/(exp(-(logit(data_cens25_r2$survival_initial[i], adjust = F) + qnorm(0.975)*se_logit)) + 1), 2)
  
  data_cens25_r2$arcsine_low[i] = round(((sin(max(c(0,asin(sqrt(data_cens25_r2$survival_initial[i])) - qnorm(0.975)*se_arcsine))))^2), 2)
  data_cens25_r2$arcsine_up[i]  = round(((sin(min(c(pi/2,asin(sqrt(data_cens25_r2$survival_initial[i])) + qnorm(0.975)*se_arcsine))))^2), 2)
  
  
  #truncation
  if(data_cens25_r2$id_low[i] <= 0) {
    data_cens25_r2$id_low[i] = 0
    data_cens25_r2$truncated_id[i] = 1} 
  
  if(data_cens25_r2$log_low[i] <= 0)  {
    data_cens25_r2$log_low[i] = 0
    data_cens25_r2$truncated_log[i] = 1} 
  
  if(data_cens25_r2$cloglog_low[i] <= 0)  {
    data_cens25_r2$cloglog_low[i] = 0
    data_cens25_r2$truncated_cloglog[i] = 1}
  
  if(data_cens25_r2$logit_low[i] <= 0)  {
    data_cens25_r2$logit_low[i] = 0
    data_cens25_r2$truncated_logit[i] = 1}
  
  if(data_cens25_r2$arcsine_low[i] <= 0)  {
    data_cens25_r2$arcsine_low[i] = 0
    data_cens25_r2$truncated_arcsine[i] = 1}
  
  
  if(data_cens25_r2$id_up[i] >= 1) {
    data_cens25_r2$id_up[i] = 1
    data_cens25_r2$truncated_id[i] = 1} 
  
  if(data_cens25_r2$log_up[i] >= 1)  {
    data_cens25_r2$log_up[i] = 1
    data_cens25_r2$truncated_log[i] = 1} 
  
  if(data_cens25_r2$cloglog_up[i] >= 1)  {
    data_cens25_r2$cloglog_up[i] = 1
    data_cens25_r2$truncated_cloglog[i] = 1} 
  
  if(data_cens25_r2$logit_up[i] >= 1)  {
    data_cens25_r2$logit_up[i] = 1
    data_cens25_r2$truncated_logit[i] = 1} 
  
  if(data_cens25_r2$arcsine_up[i] >= 1)  {
    data_cens25_r2$arcsine_up[i] = 1
    data_cens25_r2$truncated_arcsine[i] = 1} 
  
  # Absolute asymmetry 
  data_cens25_r2$AA_id[i]      = data_cens25_r2$id_up[i] + data_cens25_r2$id_low[i]- 2*data_cens25_r2$survival_rounded[i]
  data_cens25_r2$AA_log[i]     = data_cens25_r2$log_up[i] + data_cens25_r2$log_low[i]- 2*data_cens25_r2$survival_rounded[i]
  data_cens25_r2$AA_cloglog[i] = data_cens25_r2$cloglog_up[i] + data_cens25_r2$cloglog_low[i]- 2*data_cens25_r2$survival_rounded[i]
  data_cens25_r2$AA_logit[i]   = data_cens25_r2$logit_up[i] + data_cens25_r2$logit_low[i]- 2*data_cens25_r2$survival_rounded[i]
  data_cens25_r2$AA_arcsine[i] = data_cens25_r2$arcsine_up[i] + data_cens25_r2$arcsine_low[i]- 2*data_cens25_r2$survival_rounded[i]
  
  # Relative asymmetry 
  data_cens25_r2$RA_id[i]      = (data_cens25_r2$id_up[i] - data_cens25_r2$survival_rounded[i]) / (data_cens25_r2$survival_rounded[i] - data_cens25_r2$id_low[i])
  data_cens25_r2$RA_log[i]     = (data_cens25_r2$log_up[i] - data_cens25_r2$survival_rounded[i]) / (data_cens25_r2$survival_rounded[i] - data_cens25_r2$log_low[i])
  data_cens25_r2$RA_cloglog[i] = (data_cens25_r2$cloglog_up[i] - data_cens25_r2$survival_rounded[i]) / (data_cens25_r2$survival_rounded[i] - data_cens25_r2$cloglog_low[i])
  data_cens25_r2$RA_logit[i]   = (data_cens25_r2$logit_up[i] - data_cens25_r2$survival_rounded[i]) / (data_cens25_r2$survival_rounded[i] - data_cens25_r2$logit_low[i])
  data_cens25_r2$RA_arcsine[i] = (data_cens25_r2$arcsine_up[i] - data_cens25_r2$survival_rounded[i]) / (data_cens25_r2$survival_rounded[i] - data_cens25_r2$arcsine_low[i])
  
  
  # Extraction transformation is identity
  
  if(data_cens25_r2$id_low[i] == 0){
    se_true_id = (data_cens25_r2$id_up[i] - data_cens25_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens25_r2$id_up[i] == 1) {
    se_true_id = (data_cens25_r2$survival_rounded[i] - data_cens25_r2$id_low[i])/qnorm(0.975)
    
  } else {se_true_id = (data_cens25_r2$id_up[i] - data_cens25_r2$id_low[i])/(2*qnorm(0.975))}
  data_cens25_r2$true_id[i] = (se_true_id - se_id)/se_id
  
  
  if(data_cens25_r2$log_low[i] == 0){
    se_id_instead_of_log = (data_cens25_r2$log_up[i] - data_cens25_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens25_r2$log_up[i] == 1){
    se_id_instead_of_log = (data_cens25_r2$survival_rounded[i] - data_cens25_r2$log_low[i])/qnorm(0.975)
    
  } else {se_id_instead_of_log = (data_cens25_r2$log_up[i] - data_cens25_r2$log_low[i])/(2*qnorm(0.975))}
  data_cens25_r2$err_id_instead_of_log[i] = (se_id_instead_of_log - se_id)/se_id
  
  
  if(data_cens25_r2$cloglog_low[i] == 0){
    se_id_instead_of_cloglog = (data_cens25_r2$cloglog_up[i] - data_cens25_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens25_r2$cloglog_up[i] == 1){
    se_id_instead_of_cloglog = (data_cens25_r2$survival_rounded[i] - data_cens25_r2$cloglog_low[i])/qnorm(0.975)
    
  } else {se_id_instead_of_cloglog = (data_cens25_r2$cloglog_up[i] - data_cens25_r2$cloglog_low[i])/(2*qnorm(0.975))}
  data_cens25_r2$err_id_instead_of_cloglog[i] = (se_id_instead_of_cloglog - se_id)/se_id
  
  
  if(data_cens25_r2$logit_low[i] == 0){
    se_id_instead_of_logit = (data_cens25_r2$logit_up[i] - data_cens25_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens25_r2$logit_up[i] == 1){
    se_id_instead_of_logit = (data_cens25_r2$survival_rounded[i] - data_cens25_r2$logit_low[i])/qnorm(0.975)
    
  } else {se_id_instead_of_logit = (data_cens25_r2$logit_up[i] - data_cens25_r2$logit_low[i])/(2*qnorm(0.975))}
  data_cens25_r2$err_id_instead_of_logit[i] = (se_id_instead_of_logit - se_id)/se_id
  
  
  if(data_cens25_r2$arcsine_low[i] == 0){
    se_id_instead_of_arcsine = (data_cens25_r2$arcsine_up[i] - data_cens25_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens25_r2$arcsine_up[i] == 1){
    se_id_instead_of_arcsine = (data_cens25_r2$survival_rounded[i] - data_cens25_r2$arcsine_low[i])/qnorm(0.975)
    
  } else {se_id_instead_of_arcsine = (data_cens25_r2$arcsine_up[i] - data_cens25_r2$arcsine_low[i])/(2*qnorm(0.975))}
  data_cens25_r2$err_id_instead_of_arcsine[i] = (se_id_instead_of_arcsine - se_id)/se_id
  
  
  
  # Extraction transformation is logarithm
  
  if (data_cens25_r2$log_low[i] == 0) {
    se_true_log = data_cens25_r2$survival_rounded[i]*((log(data_cens25_r2$log_up[i]) - log(data_cens25_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens25_r2$log_up[i] == 1) {
    se_true_log = data_cens25_r2$survival_rounded[i]*((log(data_cens25_r2$survival_rounded[i]) - log(data_cens25_r2$log_low[i]))/(qnorm(0.975)))
    
  } else {se_true_log = data_cens25_r2$survival_rounded[i]*((log(data_cens25_r2$log_up[i]) - log(data_cens25_r2$log_low[i]))/(2*qnorm(0.975)))}
  data_cens25_r2$true_log[i] = (se_true_log - se_id)/se_id
  
  
  if (data_cens25_r2$id_low[i] == 0) {
    se_log_instead_of_id = data_cens25_r2$survival_rounded[i]*((log(data_cens25_r2$id_up[i]) - log(data_cens25_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens25_r2$id_up[i] == 1) {
    se_log_instead_of_id = data_cens25_r2$survival_rounded[i]*((log(data_cens25_r2$survival_rounded[i]) - log(data_cens25_r2$id_low[i]))/(qnorm(0.975)))
    
  } else {se_log_instead_of_id = data_cens25_r2$survival_rounded[i]*((log(data_cens25_r2$id_up[i]) - log(data_cens25_r2$id_low[i]))/(2*qnorm(0.975)))}
  data_cens25_r2$err_log_instead_of_id[i] = (se_log_instead_of_id - se_id)/se_id
  
  
  if (data_cens25_r2$cloglog_low[i] == 0) {
    se_log_instead_of_cloglog = data_cens25_r2$survival_rounded[i]*((log(data_cens25_r2$cloglog_up[i]) - log(data_cens25_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens25_r2$cloglog_up[i] == 1) {
    se_log_instead_of_cloglog = data_cens25_r2$survival_rounded[i]*((log(data_cens25_r2$survival_rounded[i]) - log(data_cens25_r2$cloglog_low[i]))/qnorm(0.975))
    
  } else {se_log_instead_of_cloglog = data_cens25_r2$survival_rounded[i]*((log(data_cens25_r2$cloglog_up[i]) - log(data_cens25_r2$cloglog_low[i]))/(2*qnorm(0.975)))}
  data_cens25_r2$err_log_instead_of_cloglog[i] = (se_log_instead_of_cloglog - se_id)/se_id
  
  
  if (data_cens25_r2$logit_low[i] == 0) {
    se_log_instead_of_logit = data_cens25_r2$survival_rounded[i]*((log(data_cens25_r2$logit_up[i]) - log(data_cens25_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens25_r2$logit_up[i] == 1) {
    se_log_instead_of_logit = data_cens25_r2$survival_rounded[i]*((log(data_cens25_r2$survival_rounded[i]) - log(data_cens25_r2$logit_low[i]))/qnorm(0.975))
    
  } else {se_log_instead_of_logit = data_cens25_r2$survival_rounded[i]*((log(data_cens25_r2$logit_up[i]) - log(data_cens25_r2$logit_low[i]))/(2*qnorm(0.975)))}
  data_cens25_r2$err_log_instead_of_logit[i] = (se_log_instead_of_logit - se_id)/se_id
  
  
  if (data_cens25_r2$arcsine_low[i] == 0) {
    se_log_instead_of_arcsine = data_cens25_r2$survival_rounded[i]*((log(data_cens25_r2$arcsine_up[i]) - log(data_cens25_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens25_r2$arcsine_up[i] == 1) {
    se_log_instead_of_arcsine = data_cens25_r2$survival_rounded[i]*((log(data_cens25_r2$survival_rounded[i]) - log(data_cens25_r2$arcsine_low[i]))/qnorm(0.975))
    
  } else {se_log_instead_of_arcsine = data_cens25_r2$survival_rounded[i]*((log(data_cens25_r2$arcsine_up[i]) - log(data_cens25_r2$arcsine_low[i]))/(2*qnorm(0.975)))}
  data_cens25_r2$err_log_instead_of_arcsine[i] = (se_log_instead_of_arcsine - se_id)/se_id
  
  
  
  # Extraction transformation is complementary log-log
  if (data_cens25_r2$cloglog_low[i] == 0) {
    se_true_cloglog = (- data_cens25_r2$survival_rounded[i]*log(data_cens25_r2$survival_rounded[i]))*((log(-log(data_cens25_r2$survival_rounded[i])) - log(-log(data_cens25_r2$cloglog_up[i])))/qnorm(0.975))
    
  } else if (data_cens25_r2$cloglog_up[i] ==1) {
    se_true_cloglog = (- data_cens25_r2$survival_rounded[i]*log(data_cens25_r2$survival_rounded[i]))*((log(-log(data_cens25_r2$cloglog_low[i])) - log(-log(data_cens25_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else {se_true_cloglog = (- data_cens25_r2$survival_rounded[i]*log(data_cens25_r2$survival_rounded[i]))*((log(-log(data_cens25_r2$cloglog_low[i])) - log(-log(data_cens25_r2$cloglog_up[i])))/(2*qnorm(0.975)))}
  data_cens25_r2$true_cloglog[i] = (se_true_cloglog - se_id)/se_id
  
  
  if (data_cens25_r2$id_low[i] == 0) {
    se_cloglog_instead_of_id = - data_cens25_r2$survival_rounded[i]*log(data_cens25_r2$survival_rounded[i])*((log(-log(data_cens25_r2$survival_rounded[i])) - log(-log(data_cens25_r2$id_up[i])))/(qnorm(0.975)))
    
  } else if (data_cens25_r2$id_up[i] == 1) {
    se_cloglog_instead_of_id = (- data_cens25_r2$survival_rounded[i]*log(data_cens25_r2$survival_rounded[i]))*((log(-log(data_cens25_r2$id_low[i])) - log(-log(data_cens25_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else{se_cloglog_instead_of_id = (- data_cens25_r2$survival_rounded[i]*log(data_cens25_r2$survival_rounded[i]))*((log(-log(data_cens25_r2$id_low[i])) - log(-log(data_cens25_r2$id_up[i])))/(2*qnorm(0.975)))}
  data_cens25_r2$err_cloglog_instead_of_id[i] = (se_cloglog_instead_of_id - se_id)/se_id
  
  
  if (data_cens25_r2$log_low[i] == 0) {
    se_cloglog_instead_of_log = (- data_cens25_r2$survival_rounded[i]*log(data_cens25_r2$survival_rounded[i]))*((log(-log(data_cens25_r2$survival_rounded[i])) - log(-log(data_cens25_r2$log_up[i])))/qnorm(0.975))
    
  } else if (data_cens25_r2$log_up[i] ==1) {
    se_cloglog_instead_of_log = (- data_cens25_r2$survival_rounded[i]*log(data_cens25_r2$survival_rounded[i]))*((log(-log(data_cens25_r2$log_low[i])) - log(-log(data_cens25_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else{se_cloglog_instead_of_log = (- data_cens25_r2$survival_rounded[i]*log(data_cens25_r2$survival_rounded[i]))*((log(-log(data_cens25_r2$log_low[i])) - log(-log(data_cens25_r2$log_up[i])))/(2*qnorm(0.975)))}
  data_cens25_r2$err_cloglog_instead_of_log[i] = (se_cloglog_instead_of_log - se_id)/se_id
  
  
  if (data_cens25_r2$logit_low[i] == 0) {
    se_cloglog_instead_of_logit = (- data_cens25_r2$survival_rounded[i]*log(data_cens25_r2$survival_rounded[i]))*((log(-log(data_cens25_r2$survival_rounded[i])) - log(-log(data_cens25_r2$logit_up[i])))/qnorm(0.975))
    
  } else if (data_cens25_r2$logit_up[i] ==1) {
    se_cloglog_instead_of_logit = (- data_cens25_r2$survival_rounded[i]*log(data_cens25_r2$survival_rounded[i]))*((log(-log(data_cens25_r2$logit_low[i])) - log(-log(data_cens25_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else{se_cloglog_instead_of_logit = (- data_cens25_r2$survival_rounded[i]*log(data_cens25_r2$survival_rounded[i]))*((log(-log(data_cens25_r2$logit_low[i])) - log(-log(data_cens25_r2$logit_up[i])))/(2*qnorm(0.975)))}
  data_cens25_r2$err_cloglog_instead_of_logit[i] = (se_cloglog_instead_of_logit - se_id)/se_id
  
  
  if (data_cens25_r2$arcsine_low[i] == 0) {
    se_cloglog_instead_of_arcsine = (- data_cens25_r2$survival_rounded[i]*log(data_cens25_r2$survival_rounded[i]))*((log(-log(data_cens25_r2$survival_rounded[i])) - log(-log(data_cens25_r2$arcsine_up[i])))/qnorm(0.975))
    
  } else if (data_cens25_r2$arcsine_up[i] ==1) {
    se_cloglog_instead_of_arcsine = (- data_cens25_r2$survival_rounded[i]*log(data_cens25_r2$survival_rounded[i]))*((log(-log(data_cens25_r2$arcsine_low[i])) - log(-log(data_cens25_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else{se_cloglog_instead_of_arcsine = (- data_cens25_r2$survival_rounded[i]*log(data_cens25_r2$survival_rounded[i]))*((log(-log(data_cens25_r2$arcsine_low[i])) - log(-log(data_cens25_r2$arcsine_up[i])))/(2*qnorm(0.975)))}
  data_cens25_r2$err_cloglog_instead_of_arcsine[i] = (se_cloglog_instead_of_arcsine - se_id)/se_id
  
  
  
  # Extraction transformation is logit
  if (data_cens25_r2$logit_low[i] == 0) {
    se_true_logit = ((logit(data_cens25_r2$logit_up[i], adjust = F) - logit(data_cens25_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i])
    
  } else if (data_cens25_r2$logit_up[i] ==1) {
    se_true_logit = ((logit(data_cens25_r2$survival_rounded[i], adjust = F) - logit(data_cens25_r2$logit_low[i], adjust = F))/qnorm(0.975))*data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i])
    
  } else {se_true_logit = ((logit(data_cens25_r2$logit_up[i], adjust = F) - logit(data_cens25_r2$logit_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i])}
  data_cens25_r2$true_logit[i] = (se_true_logit - se_id)/se_id
  
  
  if (data_cens25_r2$id_low[i] == 0) {
    se_logit_instead_of_id = ((logit(data_cens25_r2$id_up[i], adjust = F) - logit(data_cens25_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i])
    
  } else if (data_cens25_r2$id_up[i] == 1) {
    se_logit_instead_of_id = ((logit(data_cens25_r2$survival_rounded[i], adjust = F) - logit(data_cens25_r2$id_low[i], adjust = F))/qnorm(0.975))*data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i])
    
  } else{se_logit_instead_of_id = ((logit(data_cens25_r2$id_up[i], adjust = F) - logit(data_cens25_r2$id_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i])}
  data_cens25_r2$err_logit_instead_of_id[i] = (se_logit_instead_of_id - se_id)/se_id
  
  
  if (data_cens25_r2$log_low[i] == 0) {
    se_logit_instead_of_log = ((logit(data_cens25_r2$log_up[i], adjust = F) - logit(data_cens25_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i])
    
  } else if (data_cens25_r2$log_up[i] ==1) {
    se_logit_instead_of_log = ((logit(data_cens25_r2$survival_rounded[i], adjust = F) - logit(data_cens25_r2$log_low[i], adjust = F))/qnorm(0.975))*data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i])
    
  } else{se_logit_instead_of_log = ((logit(data_cens25_r2$log_up[i], adjust = F) - logit(data_cens25_r2$log_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i])}
  data_cens25_r2$err_logit_instead_of_log[i] = (se_logit_instead_of_log - se_id)/se_id
  
  
  if (data_cens25_r2$cloglog_low[i] == 0) {
    se_logit_instead_of_cloglog = ((logit(data_cens25_r2$cloglog_up[i], adjust = F) - logit(data_cens25_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i])
    
  } else if (data_cens25_r2$cloglog_up[i] ==1) {
    se_logit_instead_of_cloglog = ((logit(data_cens25_r2$survival_rounded[i], adjust = F) - logit(data_cens25_r2$cloglog_low[i], adjust = F))/qnorm(0.975))*data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i])
    
  } else{se_logit_instead_of_cloglog = ((logit(data_cens25_r2$cloglog_up[i], adjust = F) - logit(data_cens25_r2$cloglog_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i])}
  data_cens25_r2$err_logit_instead_of_cloglog[i] = (se_logit_instead_of_cloglog - se_id)/se_id
  
  
  if (data_cens25_r2$arcsine_low[i] == 0) {
    se_logit_instead_of_arcsine = ((logit(data_cens25_r2$arcsine_up[i], adjust = F) - logit(data_cens25_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i])
    
  } else if (data_cens25_r2$arcsine_up[i] ==1) {
    se_logit_instead_of_arcsine = ((logit(data_cens25_r2$survival_rounded[i], adjust = F) - logit(data_cens25_r2$arcsine_low[i], adjust = F))/qnorm(0.975))*data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i])
    
  } else{se_logit_instead_of_arcsine = ((logit(data_cens25_r2$arcsine_up[i], adjust = F) - logit(data_cens25_r2$arcsine_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i])}
  data_cens25_r2$err_logit_instead_of_arcsine[i] = (se_logit_instead_of_arcsine - se_id)/se_id
  
  
  # Extraction transformation is arcsine
  if (data_cens25_r2$arcsine_low[i] == 0) {
    se_true_arcsine = ((asin(sqrt(data_cens25_r2$arcsine_up[i])) - asin(sqrt(data_cens25_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i]))
    
  } else if (data_cens25_r2$arcsine_up[i] == 1) {
    se_true_arcsine = ((asin(sqrt(data_cens25_r2$survival_rounded[i])) - asin(sqrt(data_cens25_r2$arcsine_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i]))
    
  } else {se_true_arcsine = ((asin(sqrt(data_cens25_r2$arcsine_up[i])) - asin(sqrt(data_cens25_r2$arcsine_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i]))}
  data_cens25_r2$true_arcsine[i] = (se_true_arcsine - se_id)/se_id
  
  
  if (data_cens25_r2$id_low[i] == 0) {
    se_arcsine_instead_of_id = ((asin(sqrt(data_cens25_r2$id_up[i])) - asin(sqrt(data_cens25_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i]))
    
  } else if (data_cens25_r2$id_up[i] == 1) {
    se_arcsine_instead_of_id = ((asin(sqrt(data_cens25_r2$survival_rounded[i])) - asin(sqrt(data_cens25_r2$id_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i]))
    
  } else{se_arcsine_instead_of_id = ((asin(sqrt(data_cens25_r2$id_up[i])) - asin(sqrt(data_cens25_r2$id_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i]))}
  data_cens25_r2$err_arcsine_instead_of_id[i] = (se_arcsine_instead_of_id - se_id)/se_id
  
  
  if (data_cens25_r2$log_low[i] == 0) {
    se_arcsine_instead_of_log = ((asin(sqrt(data_cens25_r2$log_up[i])) - asin(sqrt(data_cens25_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i]))
    
  } else if (data_cens25_r2$log_up[i] ==1) {
    se_arcsine_instead_of_log = ((asin(sqrt(data_cens25_r2$survival_rounded[i])) - asin(sqrt(data_cens25_r2$log_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i]))
    
  } else{se_arcsine_instead_of_log = ((asin(sqrt(data_cens25_r2$log_up[i])) - asin(sqrt(data_cens25_r2$log_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i]))}
  data_cens25_r2$err_arcsine_instead_of_log[i] = (se_arcsine_instead_of_log - se_id)/se_id
  
  
  if (data_cens25_r2$cloglog_low[i] == 0) {
    se_arcsine_instead_of_cloglog = ((asin(sqrt(data_cens25_r2$cloglog_up[i])) - asin(sqrt(data_cens25_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i]))
    
  } else if (data_cens25_r2$arcsine_up[i] ==1) {
    se_arcsine_instead_of_cloglog = ((asin(sqrt(data_cens25_r2$survival_rounded[i])) - asin(sqrt(data_cens25_r2$cloglog_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i]))
    
  } else{se_arcsine_instead_of_cloglog = ((asin(sqrt(data_cens25_r2$cloglog_up[i])) - asin(sqrt(data_cens25_r2$cloglog_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i]))}
  data_cens25_r2$err_arcsine_instead_of_cloglog[i] = (se_arcsine_instead_of_cloglog - se_id)/se_id
  
  
  if (data_cens25_r2$logit_low[i] == 0) {
    se_arcsine_instead_of_logit = ((asin(sqrt(data_cens25_r2$logit_up[i])) - asin(sqrt(data_cens25_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i]))
    
  } else if (data_cens25_r2$logit_up[i] ==1) {
    se_arcsine_instead_of_logit = ((asin(sqrt(data_cens25_r2$survival_rounded[i])) - asin(sqrt(data_cens25_r2$logit_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i]))
    
  } else{se_arcsine_instead_of_logit = ((asin(sqrt(data_cens25_r2$logit_up[i])) - asin(sqrt(data_cens25_r2$logit_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens25_r2$survival_rounded[i]*(1-data_cens25_r2$survival_rounded[i]))}
  data_cens25_r2$err_arcsine_instead_of_logit[i] = (se_arcsine_instead_of_logit - se_id)/se_id
  
  
}
# Transform extraction errors into percentages by multiplying them by 100 
for(i in names(data_cens25_r2)) {
  if (startsWith(i, "err_") | startsWith(i, "true_")) {
    data_cens25_r2[[i]] <- data_cens25_r2[[i]] * 100
  }
}
save(data_cens25_r2, file = "your_directory/datacens25_extr_erorr_r2.Rdata")


# 50% censoring  ----

for (i in 1:nrow(data_cens50_r2)) {
  
  print(i)
  ## standard error on different scales
  data_cens50_r2$se_id_initial[i] = sqrt((data_cens50_r2$survival_initial[i]*(1-data_cens50_r2$survival_initial[i]))/data_cens50_r2$sample_size[i])
  data_cens50_r2$se_id_infl[i] = data_cens50_r2$se_id_initial[i]* data_cens50_r2$inflation_coeff[i]
  
  se_id = data_cens50_r2$se_id_infl[i]
  se_log = se_id/data_cens50_r2$survival_initial[i]
  se_cloglog = - se_id / (data_cens50_r2$survival_initial[i]*log(data_cens50_r2$survival_initial[i]))
  se_logit = se_id / (data_cens50_r2$survival_initial[i]*(1-data_cens50_r2$survival_initial[i]))
  se_arcsine =  se_id / (2*sqrt((data_cens50_r2$survival_initial[i])*(1-data_cens50_r2$survival_initial[i])))
  
  
  ## 95% confidence interval of survival computed using the different methods
  # If the lower bound is equal to 0 or the upper bound equal to 1 whether due to rounding or truncation, we consider that the interval is truncated/trimmed
  #Thus, we don't need the bounds not rounded for the following steps, we only need the rounded bounds
  
  data_cens50_r2$id_low[i] = round(data_cens50_r2$survival_initial[i] - qnorm(0.975) * se_id, 2)
  data_cens50_r2$id_up[i]  = round(data_cens50_r2$survival_initial[i] + qnorm(0.975) * se_id, 2)
  
  data_cens50_r2$log_low[i] = round(data_cens50_r2$survival_initial[i] * exp(-qnorm(0.975)*se_log), 2)
  data_cens50_r2$log_up[i]  = round(data_cens50_r2$survival_initial[i] * exp(qnorm(0.975)*se_log), 2)
  
  data_cens50_r2$cloglog_low[i] = round(data_cens50_r2$survival_initial[i] ^ (exp(qnorm(0.975)*se_cloglog)), 2)
  data_cens50_r2$cloglog_up[i]  = round(data_cens50_r2$survival_initial[i] ^ (exp(-qnorm(0.975)*se_cloglog)), 2)
  
  data_cens50_r2$logit_low[i] = round(1/(exp(-(logit(data_cens50_r2$survival_initial[i], adjust = F) - qnorm(0.975)*se_logit)) + 1), 2)
  data_cens50_r2$logit_up[i]  = round(1/(exp(-(logit(data_cens50_r2$survival_initial[i], adjust = F) + qnorm(0.975)*se_logit)) + 1), 2)
  
  data_cens50_r2$arcsine_low[i] = round(((sin(max(c(0,asin(sqrt(data_cens50_r2$survival_initial[i])) - qnorm(0.975)*se_arcsine))))^2), 2)
  data_cens50_r2$arcsine_up[i]  = round(((sin(min(c(pi/2,asin(sqrt(data_cens50_r2$survival_initial[i])) + qnorm(0.975)*se_arcsine))))^2), 2)
  
  
  #truncation
  if(data_cens50_r2$id_low[i] <= 0) {
    data_cens50_r2$id_low[i] = 0
    data_cens50_r2$truncated_id[i] = 1} 
  
  if(data_cens50_r2$log_low[i] <= 0)  {
    data_cens50_r2$log_low[i] = 0
    data_cens50_r2$truncated_log[i] = 1} 
  
  if(data_cens50_r2$cloglog_low[i] <= 0)  {
    data_cens50_r2$cloglog_low[i] = 0
    data_cens50_r2$truncated_cloglog[i] = 1}
  
  if(data_cens50_r2$logit_low[i] <= 0)  {
    data_cens50_r2$logit_low[i] = 0
    data_cens50_r2$truncated_logit[i] = 1}
  
  if(data_cens50_r2$arcsine_low[i] <= 0)  {
    data_cens50_r2$arcsine_low[i] = 0
    data_cens50_r2$truncated_arcsine[i] = 1}
  
  
  if(data_cens50_r2$id_up[i] >= 1) {
    data_cens50_r2$id_up[i] = 1
    data_cens50_r2$truncated_id[i] = 1} 
  
  if(data_cens50_r2$log_up[i] >= 1)  {
    data_cens50_r2$log_up[i] = 1
    data_cens50_r2$truncated_log[i] = 1} 
  
  if(data_cens50_r2$cloglog_up[i] >= 1)  {
    data_cens50_r2$cloglog_up[i] = 1
    data_cens50_r2$truncated_cloglog[i] = 1} 
  
  if(data_cens50_r2$logit_up[i] >= 1)  {
    data_cens50_r2$logit_up[i] = 1
    data_cens50_r2$truncated_logit[i] = 1} 
  
  if(data_cens50_r2$arcsine_up[i] >= 1)  {
    data_cens50_r2$arcsine_up[i] = 1
    data_cens50_r2$truncated_arcsine[i] = 1} 
  
  # Absolute asymmetry 
  data_cens50_r2$AA_id[i]      = data_cens50_r2$id_up[i] + data_cens50_r2$id_low[i]- 2*data_cens50_r2$survival_rounded[i]
  data_cens50_r2$AA_log[i]     = data_cens50_r2$log_up[i] + data_cens50_r2$log_low[i]- 2*data_cens50_r2$survival_rounded[i]
  data_cens50_r2$AA_cloglog[i] = data_cens50_r2$cloglog_up[i] + data_cens50_r2$cloglog_low[i]- 2*data_cens50_r2$survival_rounded[i]
  data_cens50_r2$AA_logit[i]   = data_cens50_r2$logit_up[i] + data_cens50_r2$logit_low[i]- 2*data_cens50_r2$survival_rounded[i]
  data_cens50_r2$AA_arcsine[i] = data_cens50_r2$arcsine_up[i] + data_cens50_r2$arcsine_low[i]- 2*data_cens50_r2$survival_rounded[i]
  
  # Relative asymmetry 
  data_cens50_r2$RA_id[i]      = (data_cens50_r2$id_up[i] - data_cens50_r2$survival_rounded[i]) / (data_cens50_r2$survival_rounded[i] - data_cens50_r2$id_low[i])
  data_cens50_r2$RA_log[i]     = (data_cens50_r2$log_up[i] - data_cens50_r2$survival_rounded[i]) / (data_cens50_r2$survival_rounded[i] - data_cens50_r2$log_low[i])
  data_cens50_r2$RA_cloglog[i] = (data_cens50_r2$cloglog_up[i] - data_cens50_r2$survival_rounded[i]) / (data_cens50_r2$survival_rounded[i] - data_cens50_r2$cloglog_low[i])
  data_cens50_r2$RA_logit[i]   = (data_cens50_r2$logit_up[i] - data_cens50_r2$survival_rounded[i]) / (data_cens50_r2$survival_rounded[i] - data_cens50_r2$logit_low[i])
  data_cens50_r2$RA_arcsine[i] = (data_cens50_r2$arcsine_up[i] - data_cens50_r2$survival_rounded[i]) / (data_cens50_r2$survival_rounded[i] - data_cens50_r2$arcsine_low[i])
  
  
  # Extraction transformation is identity
  
  if(data_cens50_r2$id_low[i] == 0){
    se_true_id = (data_cens50_r2$id_up[i] - data_cens50_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens50_r2$id_up[i] == 1) {
    se_true_id = (data_cens50_r2$survival_rounded[i] - data_cens50_r2$id_low[i])/qnorm(0.975)
    
  } else {se_true_id = (data_cens50_r2$id_up[i] - data_cens50_r2$id_low[i])/(2*qnorm(0.975))}
  data_cens50_r2$true_id[i] = (se_true_id - se_id)/se_id
  
  
  if(data_cens50_r2$log_low[i] == 0){
    se_id_instead_of_log = (data_cens50_r2$log_up[i] - data_cens50_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens50_r2$log_up[i] == 1){
    se_id_instead_of_log = (data_cens50_r2$survival_rounded[i] - data_cens50_r2$log_low[i])/qnorm(0.975)
    
  } else {se_id_instead_of_log = (data_cens50_r2$log_up[i] - data_cens50_r2$log_low[i])/(2*qnorm(0.975))}
  data_cens50_r2$err_id_instead_of_log[i] = (se_id_instead_of_log - se_id)/se_id
  
  
  if(data_cens50_r2$cloglog_low[i] == 0){
    se_id_instead_of_cloglog = (data_cens50_r2$cloglog_up[i] - data_cens50_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens50_r2$cloglog_up[i] == 1){
    se_id_instead_of_cloglog = (data_cens50_r2$survival_rounded[i] - data_cens50_r2$cloglog_low[i])/qnorm(0.975)
    
  } else {se_id_instead_of_cloglog = (data_cens50_r2$cloglog_up[i] - data_cens50_r2$cloglog_low[i])/(2*qnorm(0.975))}
  data_cens50_r2$err_id_instead_of_cloglog[i] = (se_id_instead_of_cloglog - se_id)/se_id
  
  
  if(data_cens50_r2$logit_low[i] == 0){
    se_id_instead_of_logit = (data_cens50_r2$logit_up[i] - data_cens50_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens50_r2$logit_up[i] == 1){
    se_id_instead_of_logit = (data_cens50_r2$survival_rounded[i] - data_cens50_r2$logit_low[i])/qnorm(0.975)
    
  } else {se_id_instead_of_logit = (data_cens50_r2$logit_up[i] - data_cens50_r2$logit_low[i])/(2*qnorm(0.975))}
  data_cens50_r2$err_id_instead_of_logit[i] = (se_id_instead_of_logit - se_id)/se_id
  
  
  if(data_cens50_r2$arcsine_low[i] == 0){
    se_id_instead_of_arcsine = (data_cens50_r2$arcsine_up[i] - data_cens50_r2$survival_rounded[i])/qnorm(0.975)
    
  } else if(data_cens50_r2$arcsine_up[i] == 1){
    se_id_instead_of_arcsine = (data_cens50_r2$survival_rounded[i] - data_cens50_r2$arcsine_low[i])/qnorm(0.975)
    
  } else {se_id_instead_of_arcsine = (data_cens50_r2$arcsine_up[i] - data_cens50_r2$arcsine_low[i])/(2*qnorm(0.975))}
  data_cens50_r2$err_id_instead_of_arcsine[i] = (se_id_instead_of_arcsine - se_id)/se_id
  
  
  
  # Extraction transformation is logarithm
  
  if (data_cens50_r2$log_low[i] == 0) {
    se_true_log = data_cens50_r2$survival_rounded[i]*((log(data_cens50_r2$log_up[i]) - log(data_cens50_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens50_r2$log_up[i] == 1) {
    se_true_log = data_cens50_r2$survival_rounded[i]*((log(data_cens50_r2$survival_rounded[i]) - log(data_cens50_r2$log_low[i]))/(qnorm(0.975)))
    
  } else {se_true_log = data_cens50_r2$survival_rounded[i]*((log(data_cens50_r2$log_up[i]) - log(data_cens50_r2$log_low[i]))/(2*qnorm(0.975)))}
  data_cens50_r2$true_log[i] = (se_true_log - se_id)/se_id
  
  
  if (data_cens50_r2$id_low[i] == 0) {
    se_log_instead_of_id = data_cens50_r2$survival_rounded[i]*((log(data_cens50_r2$id_up[i]) - log(data_cens50_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens50_r2$id_up[i] == 1) {
    se_log_instead_of_id = data_cens50_r2$survival_rounded[i]*((log(data_cens50_r2$survival_rounded[i]) - log(data_cens50_r2$id_low[i]))/(qnorm(0.975)))
    
  } else {se_log_instead_of_id = data_cens50_r2$survival_rounded[i]*((log(data_cens50_r2$id_up[i]) - log(data_cens50_r2$id_low[i]))/(2*qnorm(0.975)))}
  data_cens50_r2$err_log_instead_of_id[i] = (se_log_instead_of_id - se_id)/se_id
  
  
  if (data_cens50_r2$cloglog_low[i] == 0) {
    se_log_instead_of_cloglog = data_cens50_r2$survival_rounded[i]*((log(data_cens50_r2$cloglog_up[i]) - log(data_cens50_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens50_r2$cloglog_up[i] == 1) {
    se_log_instead_of_cloglog = data_cens50_r2$survival_rounded[i]*((log(data_cens50_r2$survival_rounded[i]) - log(data_cens50_r2$cloglog_low[i]))/qnorm(0.975))
    
  } else {se_log_instead_of_cloglog = data_cens50_r2$survival_rounded[i]*((log(data_cens50_r2$cloglog_up[i]) - log(data_cens50_r2$cloglog_low[i]))/(2*qnorm(0.975)))}
  data_cens50_r2$err_log_instead_of_cloglog[i] = (se_log_instead_of_cloglog - se_id)/se_id
  
  
  if (data_cens50_r2$logit_low[i] == 0) {
    se_log_instead_of_logit = data_cens50_r2$survival_rounded[i]*((log(data_cens50_r2$logit_up[i]) - log(data_cens50_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens50_r2$logit_up[i] == 1) {
    se_log_instead_of_logit = data_cens50_r2$survival_rounded[i]*((log(data_cens50_r2$survival_rounded[i]) - log(data_cens50_r2$logit_low[i]))/qnorm(0.975))
    
  } else {se_log_instead_of_logit = data_cens50_r2$survival_rounded[i]*((log(data_cens50_r2$logit_up[i]) - log(data_cens50_r2$logit_low[i]))/(2*qnorm(0.975)))}
  data_cens50_r2$err_log_instead_of_logit[i] = (se_log_instead_of_logit - se_id)/se_id
  
  
  if (data_cens50_r2$arcsine_low[i] == 0) {
    se_log_instead_of_arcsine = data_cens50_r2$survival_rounded[i]*((log(data_cens50_r2$arcsine_up[i]) - log(data_cens50_r2$survival_rounded[i]))/(qnorm(0.975)))
    
  } else if (data_cens50_r2$arcsine_up[i] == 1) {
    se_log_instead_of_arcsine = data_cens50_r2$survival_rounded[i]*((log(data_cens50_r2$survival_rounded[i]) - log(data_cens50_r2$arcsine_low[i]))/qnorm(0.975))
    
  } else {se_log_instead_of_arcsine = data_cens50_r2$survival_rounded[i]*((log(data_cens50_r2$arcsine_up[i]) - log(data_cens50_r2$arcsine_low[i]))/(2*qnorm(0.975)))}
  data_cens50_r2$err_log_instead_of_arcsine[i] = (se_log_instead_of_arcsine - se_id)/se_id
  
  
  
  # Extraction transformation is complementary log-log
  if (data_cens50_r2$cloglog_low[i] == 0) {
    se_true_cloglog = (- data_cens50_r2$survival_rounded[i]*log(data_cens50_r2$survival_rounded[i]))*((log(-log(data_cens50_r2$survival_rounded[i])) - log(-log(data_cens50_r2$cloglog_up[i])))/qnorm(0.975))
    
  } else if (data_cens50_r2$cloglog_up[i] ==1) {
    se_true_cloglog = (- data_cens50_r2$survival_rounded[i]*log(data_cens50_r2$survival_rounded[i]))*((log(-log(data_cens50_r2$cloglog_low[i])) - log(-log(data_cens50_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else {se_true_cloglog = (- data_cens50_r2$survival_rounded[i]*log(data_cens50_r2$survival_rounded[i]))*((log(-log(data_cens50_r2$cloglog_low[i])) - log(-log(data_cens50_r2$cloglog_up[i])))/(2*qnorm(0.975)))}
  data_cens50_r2$true_cloglog[i] = (se_true_cloglog - se_id)/se_id
  
  
  if (data_cens50_r2$id_low[i] == 0) {
    se_cloglog_instead_of_id = - data_cens50_r2$survival_rounded[i]*log(data_cens50_r2$survival_rounded[i])*((log(-log(data_cens50_r2$survival_rounded[i])) - log(-log(data_cens50_r2$id_up[i])))/(qnorm(0.975)))
    
  } else if (data_cens50_r2$id_up[i] == 1) {
    se_cloglog_instead_of_id = (- data_cens50_r2$survival_rounded[i]*log(data_cens50_r2$survival_rounded[i]))*((log(-log(data_cens50_r2$id_low[i])) - log(-log(data_cens50_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else{se_cloglog_instead_of_id = (- data_cens50_r2$survival_rounded[i]*log(data_cens50_r2$survival_rounded[i]))*((log(-log(data_cens50_r2$id_low[i])) - log(-log(data_cens50_r2$id_up[i])))/(2*qnorm(0.975)))}
  data_cens50_r2$err_cloglog_instead_of_id[i] = (se_cloglog_instead_of_id - se_id)/se_id
  
  
  if (data_cens50_r2$log_low[i] == 0) {
    se_cloglog_instead_of_log = (- data_cens50_r2$survival_rounded[i]*log(data_cens50_r2$survival_rounded[i]))*((log(-log(data_cens50_r2$survival_rounded[i])) - log(-log(data_cens50_r2$log_up[i])))/qnorm(0.975))
    
  } else if (data_cens50_r2$log_up[i] ==1) {
    se_cloglog_instead_of_log = (- data_cens50_r2$survival_rounded[i]*log(data_cens50_r2$survival_rounded[i]))*((log(-log(data_cens50_r2$log_low[i])) - log(-log(data_cens50_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else{se_cloglog_instead_of_log = (- data_cens50_r2$survival_rounded[i]*log(data_cens50_r2$survival_rounded[i]))*((log(-log(data_cens50_r2$log_low[i])) - log(-log(data_cens50_r2$log_up[i])))/(2*qnorm(0.975)))}
  data_cens50_r2$err_cloglog_instead_of_log[i] = (se_cloglog_instead_of_log - se_id)/se_id
  
  
  if (data_cens50_r2$logit_low[i] == 0) {
    se_cloglog_instead_of_logit = (- data_cens50_r2$survival_rounded[i]*log(data_cens50_r2$survival_rounded[i]))*((log(-log(data_cens50_r2$survival_rounded[i])) - log(-log(data_cens50_r2$logit_up[i])))/qnorm(0.975))
    
  } else if (data_cens50_r2$logit_up[i] ==1) {
    se_cloglog_instead_of_logit = (- data_cens50_r2$survival_rounded[i]*log(data_cens50_r2$survival_rounded[i]))*((log(-log(data_cens50_r2$logit_low[i])) - log(-log(data_cens50_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else{se_cloglog_instead_of_logit = (- data_cens50_r2$survival_rounded[i]*log(data_cens50_r2$survival_rounded[i]))*((log(-log(data_cens50_r2$logit_low[i])) - log(-log(data_cens50_r2$logit_up[i])))/(2*qnorm(0.975)))}
  data_cens50_r2$err_cloglog_instead_of_logit[i] = (se_cloglog_instead_of_logit - se_id)/se_id
  
  
  if (data_cens50_r2$arcsine_low[i] == 0) {
    se_cloglog_instead_of_arcsine = (- data_cens50_r2$survival_rounded[i]*log(data_cens50_r2$survival_rounded[i]))*((log(-log(data_cens50_r2$survival_rounded[i])) - log(-log(data_cens50_r2$arcsine_up[i])))/qnorm(0.975))
    
  } else if (data_cens50_r2$arcsine_up[i] ==1) {
    se_cloglog_instead_of_arcsine = (- data_cens50_r2$survival_rounded[i]*log(data_cens50_r2$survival_rounded[i]))*((log(-log(data_cens50_r2$arcsine_low[i])) - log(-log(data_cens50_r2$survival_rounded[i])))/qnorm(0.975))
    
  } else{se_cloglog_instead_of_arcsine = (- data_cens50_r2$survival_rounded[i]*log(data_cens50_r2$survival_rounded[i]))*((log(-log(data_cens50_r2$arcsine_low[i])) - log(-log(data_cens50_r2$arcsine_up[i])))/(2*qnorm(0.975)))}
  data_cens50_r2$err_cloglog_instead_of_arcsine[i] = (se_cloglog_instead_of_arcsine - se_id)/se_id
  
  
  
  # Extraction transformation is logit
  if (data_cens50_r2$logit_low[i] == 0) {
    se_true_logit = ((logit(data_cens50_r2$logit_up[i], adjust = F) - logit(data_cens50_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i])
    
  } else if (data_cens50_r2$logit_up[i] ==1) {
    se_true_logit = ((logit(data_cens50_r2$survival_rounded[i], adjust = F) - logit(data_cens50_r2$logit_low[i], adjust = F))/qnorm(0.975))*data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i])
    
  } else {se_true_logit = ((logit(data_cens50_r2$logit_up[i], adjust = F) - logit(data_cens50_r2$logit_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i])}
  data_cens50_r2$true_logit[i] = (se_true_logit - se_id)/se_id
  
  
  if (data_cens50_r2$id_low[i] == 0) {
    se_logit_instead_of_id = ((logit(data_cens50_r2$id_up[i], adjust = F) - logit(data_cens50_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i])
    
  } else if (data_cens50_r2$id_up[i] == 1) {
    se_logit_instead_of_id = ((logit(data_cens50_r2$survival_rounded[i], adjust = F) - logit(data_cens50_r2$id_low[i], adjust = F))/qnorm(0.975))*data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i])
    
  } else{se_logit_instead_of_id = ((logit(data_cens50_r2$id_up[i], adjust = F) - logit(data_cens50_r2$id_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i])}
  data_cens50_r2$err_logit_instead_of_id[i] = (se_logit_instead_of_id - se_id)/se_id
  
  
  if (data_cens50_r2$log_low[i] == 0) {
    se_logit_instead_of_log = ((logit(data_cens50_r2$log_up[i], adjust = F) - logit(data_cens50_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i])
    
  } else if (data_cens50_r2$log_up[i] ==1) {
    se_logit_instead_of_log = ((logit(data_cens50_r2$survival_rounded[i], adjust = F) - logit(data_cens50_r2$log_low[i], adjust = F))/qnorm(0.975))*data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i])
    
  } else{se_logit_instead_of_log = ((logit(data_cens50_r2$log_up[i], adjust = F) - logit(data_cens50_r2$log_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i])}
  data_cens50_r2$err_logit_instead_of_log[i] = (se_logit_instead_of_log - se_id)/se_id
  
  
  if (data_cens50_r2$cloglog_low[i] == 0) {
    se_logit_instead_of_cloglog = ((logit(data_cens50_r2$cloglog_up[i], adjust = F) - logit(data_cens50_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i])
    
  } else if (data_cens50_r2$cloglog_up[i] ==1) {
    se_logit_instead_of_cloglog = ((logit(data_cens50_r2$survival_rounded[i], adjust = F) - logit(data_cens50_r2$cloglog_low[i], adjust = F))/qnorm(0.975))*data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i])
    
  } else{se_logit_instead_of_cloglog = ((logit(data_cens50_r2$cloglog_up[i], adjust = F) - logit(data_cens50_r2$cloglog_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i])}
  data_cens50_r2$err_logit_instead_of_cloglog[i] = (se_logit_instead_of_cloglog - se_id)/se_id
  
  
  if (data_cens50_r2$arcsine_low[i] == 0) {
    se_logit_instead_of_arcsine = ((logit(data_cens50_r2$arcsine_up[i], adjust = F) - logit(data_cens50_r2$survival_rounded[i], adjust = F))/qnorm(0.975))*data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i])
    
  } else if (data_cens50_r2$arcsine_up[i] ==1) {
    se_logit_instead_of_arcsine = ((logit(data_cens50_r2$survival_rounded[i], adjust = F) - logit(data_cens50_r2$arcsine_low[i], adjust = F))/qnorm(0.975))*data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i])
    
  } else{se_logit_instead_of_arcsine = ((logit(data_cens50_r2$arcsine_up[i], adjust = F) - logit(data_cens50_r2$arcsine_low[i], adjust = F))/(2*qnorm(0.975)))*data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i])}
  data_cens50_r2$err_logit_instead_of_arcsine[i] = (se_logit_instead_of_arcsine - se_id)/se_id
  
  
  # Extraction transformation is arcsine
  if (data_cens50_r2$arcsine_low[i] == 0) {
    se_true_arcsine = ((asin(sqrt(data_cens50_r2$arcsine_up[i])) - asin(sqrt(data_cens50_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i]))
    
  } else if (data_cens50_r2$arcsine_up[i] == 1) {
    se_true_arcsine = ((asin(sqrt(data_cens50_r2$survival_rounded[i])) - asin(sqrt(data_cens50_r2$arcsine_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i]))
    
  } else {se_true_arcsine = ((asin(sqrt(data_cens50_r2$arcsine_up[i])) - asin(sqrt(data_cens50_r2$arcsine_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i]))}
  data_cens50_r2$true_arcsine[i] = (se_true_arcsine - se_id)/se_id
  
  
  if (data_cens50_r2$id_low[i] == 0) {
    se_arcsine_instead_of_id = ((asin(sqrt(data_cens50_r2$id_up[i])) - asin(sqrt(data_cens50_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i]))
    
  } else if (data_cens50_r2$id_up[i] == 1) {
    se_arcsine_instead_of_id = ((asin(sqrt(data_cens50_r2$survival_rounded[i])) - asin(sqrt(data_cens50_r2$id_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i]))
    
  } else{se_arcsine_instead_of_id = ((asin(sqrt(data_cens50_r2$id_up[i])) - asin(sqrt(data_cens50_r2$id_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i]))}
  data_cens50_r2$err_arcsine_instead_of_id[i] = (se_arcsine_instead_of_id - se_id)/se_id
  
  
  if (data_cens50_r2$log_low[i] == 0) {
    se_arcsine_instead_of_log = ((asin(sqrt(data_cens50_r2$log_up[i])) - asin(sqrt(data_cens50_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i]))
    
  } else if (data_cens50_r2$log_up[i] ==1) {
    se_arcsine_instead_of_log = ((asin(sqrt(data_cens50_r2$survival_rounded[i])) - asin(sqrt(data_cens50_r2$log_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i]))
    
  } else{se_arcsine_instead_of_log = ((asin(sqrt(data_cens50_r2$log_up[i])) - asin(sqrt(data_cens50_r2$log_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i]))}
  data_cens50_r2$err_arcsine_instead_of_log[i] = (se_arcsine_instead_of_log - se_id)/se_id
  
  
  if (data_cens50_r2$cloglog_low[i] == 0) {
    se_arcsine_instead_of_cloglog = ((asin(sqrt(data_cens50_r2$cloglog_up[i])) - asin(sqrt(data_cens50_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i]))
    
  } else if (data_cens50_r2$arcsine_up[i] ==1) {
    se_arcsine_instead_of_cloglog = ((asin(sqrt(data_cens50_r2$survival_rounded[i])) - asin(sqrt(data_cens50_r2$cloglog_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i]))
    
  } else{se_arcsine_instead_of_cloglog = ((asin(sqrt(data_cens50_r2$cloglog_up[i])) - asin(sqrt(data_cens50_r2$cloglog_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i]))}
  data_cens50_r2$err_arcsine_instead_of_cloglog[i] = (se_arcsine_instead_of_cloglog - se_id)/se_id
  
  
  if (data_cens50_r2$logit_low[i] == 0) {
    se_arcsine_instead_of_logit = ((asin(sqrt(data_cens50_r2$logit_up[i])) - asin(sqrt(data_cens50_r2$survival_rounded[i])))/(qnorm(0.975)))*2*sqrt(data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i]))
    
  } else if (data_cens50_r2$logit_up[i] ==1) {
    se_arcsine_instead_of_logit = ((asin(sqrt(data_cens50_r2$survival_rounded[i])) - asin(sqrt(data_cens50_r2$logit_low[i])))/(qnorm(0.975)))*2*sqrt(data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i]))
    
  } else{se_arcsine_instead_of_logit = ((asin(sqrt(data_cens50_r2$logit_up[i])) - asin(sqrt(data_cens50_r2$logit_low[i])))/(2*qnorm(0.975)))*2*sqrt(data_cens50_r2$survival_rounded[i]*(1-data_cens50_r2$survival_rounded[i]))}
  data_cens50_r2$err_arcsine_instead_of_logit[i] = (se_arcsine_instead_of_logit - se_id)/se_id
  
  
}
# Transform extraction errors into percentages by multiplying them by 100 
for(i in names(data_cens50_r2)) {
  if (startsWith(i, "err_") | startsWith(i, "true_")) {
    data_cens50_r2[[i]] <- data_cens50_r2[[i]] * 100
  }
}

save(data_cens50_r2, file = "your_directory/datacens50_extr_erorr_r2.Rdata")
