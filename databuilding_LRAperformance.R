# Packages ----
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


# Load datasets ----
load("your_directory/datacens0_extr_erorr_r4.Rdata")
load("your_directory/datacens10_extr_erorr_r4.Rdata")
load("your_directory/datacens25_extr_erorr_r4.Rdata")
load("your_directory/datacens50_extr_erorr_r4.Rdata")

load("your_directory/datacens0_extr_erorr_r3.Rdata")
load("your_directory/datacens10_extr_erorr_r3.Rdata")
load("your_directory/datacens25_extr_erorr_r3.Rdata")
load("your_directory/datacens50_extr_erorr_r3.Rdata")

load("your_directory/datacens0_extr_erorr_r2.Rdata")
load("your_directory/datacens10_extr_erorr_r2.Rdata")
load("your_directory/datacens25_extr_erorr_r2.Rdata")
load("your_directory/datacens50_extr_erorr_r2.Rdata")


# Create one dataframe for each f and censoring rate (transformation used by the authors a.k.a fcomputation)
# One dataframe contains the 95%CI computed with the corresponding fcomputation then the LRA is computed for each fE (a.k.a fextraction)

var_names_df_LRA_fA_id <- c("sample_size", "survival_initial", "survival_rounded_r2", "survival_rounded_r3", "survival_rounded_r4", 
                            "id_low_r2", "id_low_r3", "id_low_r4", "id_up_r2", "id_up_r3", "id_up_r4","f_author", 
                            "LRA_fE_id_r2", "LRA_fE_id_r3", "LRA_fE_id_r4",
                            "LRA_fE_log_r2", "LRA_fE_log_r3", "LRA_fE_log_r4",
                            "LRA_fE_cloglog_r2", "LRA_fE_cloglog_r3", "LRA_fE_cloglog_r4",
                            "LRA_fE_logit_r2", "LRA_fE_logit_r3", "LRA_fE_logit_r4",
                            "LRA_fE_arcsine_r2", "LRA_fE_arcsine_r3", "LRA_fE_arcsine_r4",
                            "fE_givenby_LRA_r2", "fE_givenby_LRA_r3", "fE_givenby_LRA_r4",
                            "success_LRA_r2", "success_LRA_r3", "success_LRA_r4",
                            "true_id_r2", "true_id_r3", "true_id_r4",
                            "err_log_instead_of_id_r2", "err_log_instead_of_id_r3", "err_log_instead_of_id_r4",
                            "err_cloglog_instead_of_id_r2", "err_cloglog_instead_of_id_r3", "err_cloglog_instead_of_id_r4", 
                            "err_logit_instead_of_id_r2", "err_logit_instead_of_id_r3", "err_logit_instead_of_id_r4", 
                            "err_arcsine_instead_of_id_r2", "err_arcsine_instead_of_id_r3", "err_arcsine_instead_of_id_r4")

var_names_df_LRA_fA_log <- c("sample_size", "survival_initial", "survival_rounded_r2", "survival_rounded_r3", "survival_rounded_r4", 
                             "log_low_r2", "log_low_r3", "log_low_r4", "log_up_r2", "log_up_r3", "log_up_r4","f_author", 
                             "LRA_fE_id_r2", "LRA_fE_id_r3", "LRA_fE_id_r4",
                             "LRA_fE_log_r2", "LRA_fE_log_r3", "LRA_fE_log_r4",
                             "LRA_fE_cloglog_r2", "LRA_fE_cloglog_r3", "LRA_fE_cloglog_r4",
                             "LRA_fE_logit_r2", "LRA_fE_logit_r3", "LRA_fE_logit_r4",
                             "LRA_fE_arcsine_r2", "LRA_fE_arcsine_r3", "LRA_fE_arcsine_r4",
                             "fE_givenby_LRA_r2", "fE_givenby_LRA_r3", "fE_givenby_LRA_r4",
                             "success_LRA_r2", "success_LRA_r3", "success_LRA_r4",
                             "true_log_r2", "true_log_r3", "true_log_r4",
                             "err_id_instead_of_log_r2", "err_id_instead_of_log_r3", "err_id_instead_of_log_r4",
                             "err_cloglog_instead_of_log_r2", "err_cloglog_instead_of_log_r3", "err_cloglog_instead_of_log_r4", 
                             "err_logit_instead_of_log_r2", "err_logit_instead_of_log_r3", "err_logit_instead_of_log_r4", 
                             "err_arcsine_instead_of_log_r2", "err_arcsine_instead_of_log_r3", "err_arcsine_instead_of_log_r4")

var_names_df_LRA_fA_cloglog <- c("sample_size", "survival_initial", "survival_rounded_r2", "survival_rounded_r3", "survival_rounded_r4", 
                                 "cloglog_low_r2", "cloglog_low_r3", "cloglog_low_r4", "cloglog_up_r2", "cloglog_up_r3", "cloglog_up_r4","f_author", 
                                 "LRA_fE_id_r2", "LRA_fE_id_r3", "LRA_fE_id_r4",
                                 "LRA_fE_log_r2", "LRA_fE_log_r3", "LRA_fE_log_r4",
                                 "LRA_fE_cloglog_r2", "LRA_fE_cloglog_r3", "LRA_fE_cloglog_r4",
                                 "LRA_fE_logit_r2", "LRA_fE_logit_r3", "LRA_fE_logit_r4",
                                 "LRA_fE_arcsine_r2", "LRA_fE_arcsine_r3", "LRA_fE_arcsine_r4",
                                 "fE_givenby_LRA_r2", "fE_givenby_LRA_r3", "fE_givenby_LRA_r4",
                                 "success_LRA_r2", "success_LRA_r3", "success_LRA_r4",
                                 "true_cloglog_r2", "true_cloglog_r3", "true_cloglog_r4",
                                 "err_id_instead_of_cloglog_r2", "err_id_instead_of_cloglog_r3", "err_id_instead_of_cloglog_r4",
                                 "err_log_instead_of_cloglog_r2", "err_log_instead_of_cloglog_r3", "err_log_instead_of_cloglog_r4", 
                                 "err_logit_instead_of_cloglog_r2", "err_logit_instead_of_cloglog_r3", "err_logit_instead_of_cloglog_r4", 
                                 "err_arcsine_instead_of_cloglog_r2", "err_arcsine_instead_of_cloglog_r3", "err_arcsine_instead_of_cloglog_r4")

var_names_df_LRA_fA_logit <- c("sample_size", "survival_initial", "survival_rounded_r2", "survival_rounded_r3", "survival_rounded_r4", 
                               "logit_low_r2", "logit_low_r3", "logit_low_r4", "logit_up_r2", "logit_up_r3", "logit_up_r4","f_author", 
                               "LRA_fE_id_r2", "LRA_fE_id_r3", "LRA_fE_id_r4",
                               "LRA_fE_log_r2", "LRA_fE_log_r3", "LRA_fE_log_r4",
                               "LRA_fE_cloglog_r2", "LRA_fE_cloglog_r3", "LRA_fE_cloglog_r4",
                               "LRA_fE_logit_r2", "LRA_fE_logit_r3", "LRA_fE_logit_r4",
                               "LRA_fE_arcsine_r2", "LRA_fE_arcsine_r3", "LRA_fE_arcsine_r4",
                               "fE_givenby_LRA_r2", "fE_givenby_LRA_r3", "fE_givenby_LRA_r4",
                               "success_LRA_r2", "success_LRA_r3", "success_LRA_r4",
                               "true_logit_r2", "true_logit_r3", "true_logit_r4",
                               "err_id_instead_of_logit_r2", "err_id_instead_of_logit_r3", "err_id_instead_of_logit_r4",
                               "err_log_instead_of_logit_r2", "err_log_instead_of_logit_r3", "err_log_instead_of_logit_r4", 
                               "err_cloglog_instead_of_logit_r2", "err_cloglog_instead_of_logit_r3", "err_cloglog_instead_of_logit_r4", 
                               "err_arcsine_instead_of_logit_r2", "err_arcsine_instead_of_logit_r3", "err_arcsine_instead_of_logit_r4")

var_names_df_LRA_fA_arcsine <- c("sample_size", "survival_initial", "survival_rounded_r2", "survival_rounded_r3", "survival_rounded_r4", 
                                 "arcsine_low_r2", "arcsine_low_r3", "arcsine_low_r4", "arcsine_up_r2", "arcsine_up_r3", "arcsine_up_r4","f_author", 
                                 "LRA_fE_id_r2", "LRA_fE_id_r3", "LRA_fE_id_r4",
                                 "LRA_fE_log_r2", "LRA_fE_log_r3", "LRA_fE_log_r4",
                                 "LRA_fE_cloglog_r2", "LRA_fE_cloglog_r3", "LRA_fE_cloglog_r4",
                                 "LRA_fE_logit_r2", "LRA_fE_logit_r3", "LRA_fE_logit_r4",
                                 "LRA_fE_arcsine_r2", "LRA_fE_arcsine_r3", "LRA_fE_arcsine_r4",
                                 "fE_givenby_LRA_r2", "fE_givenby_LRA_r3", "fE_givenby_LRA_r4",
                                 "success_LRA_r2", "success_LRA_r3", "success_LRA_r4",
                                 "true_arcsine_r2", "true_arcsine_r3", "true_arcsine_r4",
                                 "err_id_instead_of_arcsine_r2", "err_id_instead_of_arcsine_r3", "err_id_instead_of_arcsine_r4",
                                 "err_log_instead_of_arcsine_r2", "err_log_instead_of_arcsine_r3", "err_log_instead_of_arcsine_r4", 
                                 "err_cloglog_instead_of_arcsine_r2", "err_cloglog_instead_of_arcsine_r3", "err_cloglog_instead_of_arcsine_r4", 
                                 "err_logit_instead_of_arcsine_r2", "err_logit_instead_of_arcsine_r3", "err_logit_instead_of_arcsine_r4")


df_LRA_fA_id_cens0      = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_id), dimnames = list(NULL, var_names_df_LRA_fA_id)))
df_LRA_fA_log_cens0     = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_log), dimnames = list(NULL, var_names_df_LRA_fA_log)))
df_LRA_fA_cloglog_cens0 = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_cloglog), dimnames = list(NULL, var_names_df_LRA_fA_cloglog)))
df_LRA_fA_logit_cens0   = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_logit), dimnames = list(NULL, var_names_df_LRA_fA_logit)))
df_LRA_fA_arcsine_cens0 = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_arcsine), dimnames = list(NULL, var_names_df_LRA_fA_arcsine)))

df_LRA_fA_id_cens10      = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_id), dimnames = list(NULL, var_names_df_LRA_fA_id)))
df_LRA_fA_log_cens10     = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_log), dimnames = list(NULL, var_names_df_LRA_fA_log)))
df_LRA_fA_cloglog_cens10 = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_cloglog), dimnames = list(NULL, var_names_df_LRA_fA_cloglog)))
df_LRA_fA_logit_cens10   = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_logit), dimnames = list(NULL, var_names_df_LRA_fA_logit)))
df_LRA_fA_arcsine_cens10 = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_arcsine), dimnames = list(NULL, var_names_df_LRA_fA_arcsine)))

df_LRA_fA_id_cens25      = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_id), dimnames = list(NULL, var_names_df_LRA_fA_id)))
df_LRA_fA_log_cens25     = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_log), dimnames = list(NULL, var_names_df_LRA_fA_log)))
df_LRA_fA_cloglog_cens25 = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_cloglog), dimnames = list(NULL, var_names_df_LRA_fA_cloglog)))
df_LRA_fA_logit_cens25   = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_logit), dimnames = list(NULL, var_names_df_LRA_fA_logit)))
df_LRA_fA_arcsine_cens25 = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_arcsine), dimnames = list(NULL, var_names_df_LRA_fA_arcsine)))

df_LRA_fA_id_cens50      = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_id), dimnames = list(NULL, var_names_df_LRA_fA_id)))
df_LRA_fA_log_cens50     = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_log), dimnames = list(NULL, var_names_df_LRA_fA_log)))
df_LRA_fA_cloglog_cens50 = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_cloglog), dimnames = list(NULL, var_names_df_LRA_fA_cloglog)))
df_LRA_fA_logit_cens50   = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_logit), dimnames = list(NULL, var_names_df_LRA_fA_logit)))
df_LRA_fA_arcsine_cens50 = as.data.frame(matrix(0, nrow = nrow(data_cens0_r4), ncol = length(var_names_df_LRA_fA_arcsine), dimnames = list(NULL, var_names_df_LRA_fA_arcsine)))


# Fill the new datasets ----
## Censoring = 0% ----
df_LRA_fA_id_cens0[, c("sample_size", "survival_initial", "survival_rounded_r2","id_low_r2", "id_up_r2",
                       "true_id_r2", "err_log_instead_of_id_r2", "err_cloglog_instead_of_id_r2", "err_logit_instead_of_id_r2", "err_arcsine_instead_of_id_r2",  
                       
                       "survival_rounded_r3", "id_low_r3", "id_up_r3", 
                       "true_id_r3", "err_log_instead_of_id_r3", "err_cloglog_instead_of_id_r3", "err_logit_instead_of_id_r3", "err_arcsine_instead_of_id_r3", 
                       
                       "survival_rounded_r4", "id_low_r4",  "id_up_r4", 
                       "true_id_r4", "err_log_instead_of_id_r4", "err_cloglog_instead_of_id_r4", "err_logit_instead_of_id_r4","err_arcsine_instead_of_id_r4")] <-
  
  c(data_cens0_r2[,c("sample_size","survival_initial", "survival_rounded", "id_low", "id_up", 
                   "true_id", "err_log_instead_of_id", "err_cloglog_instead_of_id", "err_logit_instead_of_id", "err_arcsine_instead_of_id")],
    
    data_cens0_r3[,c("survival_rounded", "id_low", "id_up", 
                     "true_id", "err_log_instead_of_id", "err_cloglog_instead_of_id", "err_logit_instead_of_id", "err_arcsine_instead_of_id")],
    
    data_cens0_r4[,c("survival_rounded", "id_low", "id_up", 
                     "true_id", "err_log_instead_of_id", "err_cloglog_instead_of_id", "err_logit_instead_of_id", "err_arcsine_instead_of_id")])
 


df_LRA_fA_log_cens0[, c("sample_size", "survival_initial", "survival_rounded_r2","log_low_r2", "log_up_r2",
                       "true_log_r2", "err_id_instead_of_log_r2", "err_cloglog_instead_of_log_r2", "err_logit_instead_of_log_r2", "err_arcsine_instead_of_log_r2",  
                       
                       "survival_rounded_r3", "log_low_r3", "log_up_r3", 
                       "true_log_r3", "err_id_instead_of_log_r3", "err_cloglog_instead_of_log_r3", "err_logit_instead_of_log_r3", "err_arcsine_instead_of_log_r3", 
                       
                       "survival_rounded_r4", "log_low_r4",  "log_up_r4", 
                       "true_log_r4", "err_id_instead_of_log_r4", "err_cloglog_instead_of_log_r4", "err_logit_instead_of_log_r4","err_arcsine_instead_of_log_r4")] <-
  
  c(data_cens0_r2[,c("sample_size","survival_initial", "survival_rounded", "log_low", "log_up", 
                     "true_log", "err_id_instead_of_log", "err_cloglog_instead_of_log", "err_logit_instead_of_log", "err_arcsine_instead_of_log")],
    
    data_cens0_r3[,c("survival_rounded", "log_low", "log_up", 
                     "true_log", "err_id_instead_of_log", "err_cloglog_instead_of_log", "err_logit_instead_of_log", "err_arcsine_instead_of_log")],
    
    data_cens0_r4[,c("survival_rounded", "log_low", "log_up", 
                     "true_log", "err_id_instead_of_log", "err_cloglog_instead_of_log", "err_logit_instead_of_log", "err_arcsine_instead_of_log")])



df_LRA_fA_cloglog_cens0[, c("sample_size", "survival_initial", "survival_rounded_r2","cloglog_low_r2", "cloglog_up_r2",
                       "true_cloglog_r2", "err_id_instead_of_cloglog_r2", "err_log_instead_of_cloglog_r2", "err_logit_instead_of_cloglog_r2", "err_arcsine_instead_of_cloglog_r2",  
                       
                       "survival_rounded_r3", "cloglog_low_r3", "cloglog_up_r3", 
                       "true_cloglog_r3", "err_id_instead_of_cloglog_r3", "err_log_instead_of_cloglog_r3", "err_logit_instead_of_cloglog_r3", "err_arcsine_instead_of_cloglog_r3", 
                       
                       "survival_rounded_r4", "cloglog_low_r4",  "cloglog_up_r4", 
                       "true_cloglog_r4", "err_id_instead_of_cloglog_r4", "err_log_instead_of_cloglog_r4", "err_logit_instead_of_cloglog_r4","err_arcsine_instead_of_cloglog_r4")] <-
  
  c(data_cens0_r2[,c("sample_size","survival_initial", "survival_rounded", "cloglog_low", "cloglog_up", 
                     "true_cloglog", "err_id_instead_of_cloglog", "err_log_instead_of_cloglog", "err_logit_instead_of_cloglog", "err_arcsine_instead_of_cloglog")],
    
    data_cens0_r3[,c("survival_rounded", "cloglog_low", "cloglog_up", 
                     "true_cloglog", "err_id_instead_of_cloglog", "err_log_instead_of_cloglog", "err_logit_instead_of_cloglog", "err_arcsine_instead_of_cloglog")],
    
    data_cens0_r4[,c("survival_rounded", "cloglog_low", "cloglog_up", 
                     "true_cloglog", "err_id_instead_of_cloglog", "err_log_instead_of_cloglog", "err_logit_instead_of_cloglog", "err_arcsine_instead_of_cloglog")])



df_LRA_fA_logit_cens0[, c("sample_size", "survival_initial", "survival_rounded_r2","logit_low_r2", "logit_up_r2",
                            "true_logit_r2", "err_id_instead_of_logit_r2", "err_log_instead_of_logit_r2", "err_cloglog_instead_of_logit_r2", "err_arcsine_instead_of_logit_r2",  
                            
                            "survival_rounded_r3", "logit_low_r3", "logit_up_r3", 
                            "true_logit_r3", "err_id_instead_of_logit_r3", "err_log_instead_of_logit_r3", "err_cloglog_instead_of_logit_r3", "err_arcsine_instead_of_logit_r3", 
                            
                            "survival_rounded_r4", "logit_low_r4",  "logit_up_r4", 
                            "true_logit_r4", "err_id_instead_of_logit_r4", "err_log_instead_of_logit_r4", "err_cloglog_instead_of_logit_r4","err_arcsine_instead_of_logit_r4")] <-
  
  c(data_cens0_r2[,c("sample_size","survival_initial", "survival_rounded", "logit_low", "logit_up", 
                     "true_logit", "err_id_instead_of_logit", "err_log_instead_of_logit", "err_cloglog_instead_of_logit", "err_arcsine_instead_of_logit")],
    
    data_cens0_r3[,c("survival_rounded", "logit_low", "logit_up", 
                     "true_logit", "err_id_instead_of_logit", "err_log_instead_of_logit", "err_cloglog_instead_of_logit", "err_arcsine_instead_of_logit")],
    
    data_cens0_r4[,c("survival_rounded", "logit_low", "logit_up", 
                     "true_logit", "err_id_instead_of_logit", "err_log_instead_of_logit", "err_cloglog_instead_of_logit", "err_arcsine_instead_of_logit")])



df_LRA_fA_arcsine_cens0[, c("sample_size", "survival_initial", "survival_rounded_r2","arcsine_low_r2", "arcsine_up_r2",
                          "true_arcsine_r2", "err_id_instead_of_arcsine_r2", "err_log_instead_of_arcsine_r2", "err_cloglog_instead_of_arcsine_r2", "err_logit_instead_of_arcsine_r2",  
                          
                          "survival_rounded_r3", "arcsine_low_r3", "arcsine_up_r3", 
                          "true_arcsine_r3", "err_id_instead_of_arcsine_r3", "err_log_instead_of_arcsine_r3", "err_cloglog_instead_of_arcsine_r3", "err_logit_instead_of_arcsine_r3", 
                          
                          "survival_rounded_r4", "arcsine_low_r4",  "arcsine_up_r4", 
                          "true_arcsine_r4", "err_id_instead_of_arcsine_r4", "err_log_instead_of_arcsine_r4", "err_cloglog_instead_of_arcsine_r4","err_logit_instead_of_arcsine_r4")] <-
  
  c(data_cens0_r2[,c("sample_size","survival_initial", "survival_rounded", "arcsine_low", "arcsine_up", 
                     "true_arcsine", "err_id_instead_of_arcsine", "err_log_instead_of_arcsine", "err_cloglog_instead_of_arcsine", "err_logit_instead_of_arcsine")],
    
    data_cens0_r3[,c("survival_rounded", "arcsine_low", "arcsine_up", 
                     "true_arcsine", "err_id_instead_of_arcsine", "err_log_instead_of_arcsine", "err_cloglog_instead_of_arcsine", "err_logit_instead_of_arcsine")],
    
    data_cens0_r4[,c("survival_rounded", "arcsine_low", "arcsine_up", 
                     "true_arcsine", "err_id_instead_of_arcsine", "err_log_instead_of_arcsine", "err_cloglog_instead_of_arcsine", "err_logit_instead_of_arcsine")])


## Censoring = 10% ----

df_LRA_fA_id_cens10[, c("sample_size", "survival_initial", "survival_rounded_r2","id_low_r2", "id_up_r2",
                        "true_id_r2", "err_log_instead_of_id_r2", "err_cloglog_instead_of_id_r2", "err_logit_instead_of_id_r2", "err_arcsine_instead_of_id_r2",  
                        
                        "survival_rounded_r3", "id_low_r3", "id_up_r3", 
                        "true_id_r3", "err_log_instead_of_id_r3", "err_cloglog_instead_of_id_r3", "err_logit_instead_of_id_r3", "err_arcsine_instead_of_id_r3", 
                        
                        "survival_rounded_r4", "id_low_r4",  "id_up_r4", 
                        "true_id_r4", "err_log_instead_of_id_r4", "err_cloglog_instead_of_id_r4", "err_logit_instead_of_id_r4","err_arcsine_instead_of_id_r4")] <-
  
  c(data_cens10_r2[,c("sample_size","survival_initial", "survival_rounded", "id_low", "id_up", 
                      "true_id", "err_log_instead_of_id", "err_cloglog_instead_of_id", "err_logit_instead_of_id", "err_arcsine_instead_of_id")],
    
    data_cens10_r3[,c("survival_rounded", "id_low", "id_up", 
                      "true_id", "err_log_instead_of_id", "err_cloglog_instead_of_id", "err_logit_instead_of_id", "err_arcsine_instead_of_id")],
    
    data_cens10_r4[,c("survival_rounded", "id_low", "id_up", 
                      "true_id", "err_log_instead_of_id", "err_cloglog_instead_of_id", "err_logit_instead_of_id", "err_arcsine_instead_of_id")])



df_LRA_fA_log_cens10[, c("sample_size", "survival_initial", "survival_rounded_r2","log_low_r2", "log_up_r2",
                         "true_log_r2", "err_id_instead_of_log_r2", "err_cloglog_instead_of_log_r2", "err_logit_instead_of_log_r2", "err_arcsine_instead_of_log_r2",  
                         
                         "survival_rounded_r3", "log_low_r3", "log_up_r3", 
                         "true_log_r3", "err_id_instead_of_log_r3", "err_cloglog_instead_of_log_r3", "err_logit_instead_of_log_r3", "err_arcsine_instead_of_log_r3", 
                         
                         "survival_rounded_r4", "log_low_r4",  "log_up_r4", 
                         "true_log_r4", "err_id_instead_of_log_r4", "err_cloglog_instead_of_log_r4", "err_logit_instead_of_log_r4","err_arcsine_instead_of_log_r4")] <-
  
  c(data_cens10_r2[,c("sample_size","survival_initial", "survival_rounded", "log_low", "log_up", 
                      "true_log", "err_id_instead_of_log", "err_cloglog_instead_of_log", "err_logit_instead_of_log", "err_arcsine_instead_of_log")],
    
    data_cens10_r3[,c("survival_rounded", "log_low", "log_up", 
                      "true_log", "err_id_instead_of_log", "err_cloglog_instead_of_log", "err_logit_instead_of_log", "err_arcsine_instead_of_log")],
    
    data_cens10_r4[,c("survival_rounded", "log_low", "log_up", 
                      "true_log", "err_id_instead_of_log", "err_cloglog_instead_of_log", "err_logit_instead_of_log", "err_arcsine_instead_of_log")])



df_LRA_fA_cloglog_cens10[, c("sample_size", "survival_initial", "survival_rounded_r2","cloglog_low_r2", "cloglog_up_r2",
                             "true_cloglog_r2", "err_id_instead_of_cloglog_r2", "err_log_instead_of_cloglog_r2", "err_logit_instead_of_cloglog_r2", "err_arcsine_instead_of_cloglog_r2",  
                             
                             "survival_rounded_r3", "cloglog_low_r3", "cloglog_up_r3", 
                             "true_cloglog_r3", "err_id_instead_of_cloglog_r3", "err_log_instead_of_cloglog_r3", "err_logit_instead_of_cloglog_r3", "err_arcsine_instead_of_cloglog_r3", 
                             
                             "survival_rounded_r4", "cloglog_low_r4",  "cloglog_up_r4", 
                             "true_cloglog_r4", "err_id_instead_of_cloglog_r4", "err_log_instead_of_cloglog_r4", "err_logit_instead_of_cloglog_r4","err_arcsine_instead_of_cloglog_r4")] <-
  
  c(data_cens10_r2[,c("sample_size","survival_initial", "survival_rounded", "cloglog_low", "cloglog_up", 
                      "true_cloglog", "err_id_instead_of_cloglog", "err_log_instead_of_cloglog", "err_logit_instead_of_cloglog", "err_arcsine_instead_of_cloglog")],
    
    data_cens10_r3[,c("survival_rounded", "cloglog_low", "cloglog_up", 
                      "true_cloglog", "err_id_instead_of_cloglog", "err_log_instead_of_cloglog", "err_logit_instead_of_cloglog", "err_arcsine_instead_of_cloglog")],
    
    data_cens10_r4[,c("survival_rounded", "cloglog_low", "cloglog_up", 
                      "true_cloglog", "err_id_instead_of_cloglog", "err_log_instead_of_cloglog", "err_logit_instead_of_cloglog", "err_arcsine_instead_of_cloglog")])



df_LRA_fA_logit_cens10[, c("sample_size", "survival_initial", "survival_rounded_r2","logit_low_r2", "logit_up_r2",
                           "true_logit_r2", "err_id_instead_of_logit_r2", "err_log_instead_of_logit_r2", "err_cloglog_instead_of_logit_r2", "err_arcsine_instead_of_logit_r2",  
                           
                           "survival_rounded_r3", "logit_low_r3", "logit_up_r3", 
                           "true_logit_r3", "err_id_instead_of_logit_r3", "err_log_instead_of_logit_r3", "err_cloglog_instead_of_logit_r3", "err_arcsine_instead_of_logit_r3", 
                           
                           "survival_rounded_r4", "logit_low_r4",  "logit_up_r4", 
                           "true_logit_r4", "err_id_instead_of_logit_r4", "err_log_instead_of_logit_r4", "err_cloglog_instead_of_logit_r4","err_arcsine_instead_of_logit_r4")] <-
  
  c(data_cens10_r2[,c("sample_size","survival_initial", "survival_rounded", "logit_low", "logit_up", 
                      "true_logit", "err_id_instead_of_logit", "err_log_instead_of_logit", "err_cloglog_instead_of_logit", "err_arcsine_instead_of_logit")],
    
    data_cens10_r3[,c("survival_rounded", "logit_low", "logit_up", 
                      "true_logit", "err_id_instead_of_logit", "err_log_instead_of_logit", "err_cloglog_instead_of_logit", "err_arcsine_instead_of_logit")],
    
    data_cens10_r4[,c("survival_rounded", "logit_low", "logit_up", 
                      "true_logit", "err_id_instead_of_logit", "err_log_instead_of_logit", "err_cloglog_instead_of_logit", "err_arcsine_instead_of_logit")])



df_LRA_fA_arcsine_cens10[, c("sample_size", "survival_initial", "survival_rounded_r2","arcsine_low_r2", "arcsine_up_r2",
                             "true_arcsine_r2", "err_id_instead_of_arcsine_r2", "err_log_instead_of_arcsine_r2", "err_cloglog_instead_of_arcsine_r2", "err_logit_instead_of_arcsine_r2",  
                             
                             "survival_rounded_r3", "arcsine_low_r3", "arcsine_up_r3", 
                             "true_arcsine_r3", "err_id_instead_of_arcsine_r3", "err_log_instead_of_arcsine_r3", "err_cloglog_instead_of_arcsine_r3", "err_logit_instead_of_arcsine_r3", 
                             
                             "survival_rounded_r4", "arcsine_low_r4",  "arcsine_up_r4", 
                             "true_arcsine_r4", "err_id_instead_of_arcsine_r4", "err_log_instead_of_arcsine_r4", "err_cloglog_instead_of_arcsine_r4","err_logit_instead_of_arcsine_r4")] <-
  
  c(data_cens10_r2[,c("sample_size","survival_initial", "survival_rounded", "arcsine_low", "arcsine_up", 
                      "true_arcsine", "err_id_instead_of_arcsine", "err_log_instead_of_arcsine", "err_cloglog_instead_of_arcsine", "err_logit_instead_of_arcsine")],
    
    data_cens10_r3[,c("survival_rounded", "arcsine_low", "arcsine_up", 
                      "true_arcsine", "err_id_instead_of_arcsine", "err_log_instead_of_arcsine", "err_cloglog_instead_of_arcsine", "err_logit_instead_of_arcsine")],
    
    data_cens10_r4[,c("survival_rounded", "arcsine_low", "arcsine_up", 
                      "true_arcsine", "err_id_instead_of_arcsine", "err_log_instead_of_arcsine", "err_cloglog_instead_of_arcsine", "err_logit_instead_of_arcsine")])




## Censoring = 25% ----

df_LRA_fA_id_cens25[, c("sample_size", "survival_initial", "survival_rounded_r2","id_low_r2", "id_up_r2",
                       "true_id_r2", "err_log_instead_of_id_r2", "err_cloglog_instead_of_id_r2", "err_logit_instead_of_id_r2", "err_arcsine_instead_of_id_r2",  
                       
                       "survival_rounded_r3", "id_low_r3", "id_up_r3", 
                       "true_id_r3", "err_log_instead_of_id_r3", "err_cloglog_instead_of_id_r3", "err_logit_instead_of_id_r3", "err_arcsine_instead_of_id_r3", 
                       
                       "survival_rounded_r4", "id_low_r4",  "id_up_r4", 
                       "true_id_r4", "err_log_instead_of_id_r4", "err_cloglog_instead_of_id_r4", "err_logit_instead_of_id_r4","err_arcsine_instead_of_id_r4")] <-
  
  c(data_cens25_r2[,c("sample_size","survival_initial", "survival_rounded", "id_low", "id_up", 
                     "true_id", "err_log_instead_of_id", "err_cloglog_instead_of_id", "err_logit_instead_of_id", "err_arcsine_instead_of_id")],
    
    data_cens25_r3[,c("survival_rounded", "id_low", "id_up", 
                     "true_id", "err_log_instead_of_id", "err_cloglog_instead_of_id", "err_logit_instead_of_id", "err_arcsine_instead_of_id")],
    
    data_cens25_r4[,c("survival_rounded", "id_low", "id_up", 
                     "true_id", "err_log_instead_of_id", "err_cloglog_instead_of_id", "err_logit_instead_of_id", "err_arcsine_instead_of_id")])



df_LRA_fA_log_cens25[, c("sample_size", "survival_initial", "survival_rounded_r2","log_low_r2", "log_up_r2",
                        "true_log_r2", "err_id_instead_of_log_r2", "err_cloglog_instead_of_log_r2", "err_logit_instead_of_log_r2", "err_arcsine_instead_of_log_r2",  
                        
                        "survival_rounded_r3", "log_low_r3", "log_up_r3", 
                        "true_log_r3", "err_id_instead_of_log_r3", "err_cloglog_instead_of_log_r3", "err_logit_instead_of_log_r3", "err_arcsine_instead_of_log_r3", 
                        
                        "survival_rounded_r4", "log_low_r4",  "log_up_r4", 
                        "true_log_r4", "err_id_instead_of_log_r4", "err_cloglog_instead_of_log_r4", "err_logit_instead_of_log_r4","err_arcsine_instead_of_log_r4")] <-
  
  c(data_cens25_r2[,c("sample_size","survival_initial", "survival_rounded", "log_low", "log_up", 
                     "true_log", "err_id_instead_of_log", "err_cloglog_instead_of_log", "err_logit_instead_of_log", "err_arcsine_instead_of_log")],
    
    data_cens25_r3[,c("survival_rounded", "log_low", "log_up", 
                     "true_log", "err_id_instead_of_log", "err_cloglog_instead_of_log", "err_logit_instead_of_log", "err_arcsine_instead_of_log")],
    
    data_cens25_r4[,c("survival_rounded", "log_low", "log_up", 
                     "true_log", "err_id_instead_of_log", "err_cloglog_instead_of_log", "err_logit_instead_of_log", "err_arcsine_instead_of_log")])



df_LRA_fA_cloglog_cens25[, c("sample_size", "survival_initial", "survival_rounded_r2","cloglog_low_r2", "cloglog_up_r2",
                            "true_cloglog_r2", "err_id_instead_of_cloglog_r2", "err_log_instead_of_cloglog_r2", "err_logit_instead_of_cloglog_r2", "err_arcsine_instead_of_cloglog_r2",  
                            
                            "survival_rounded_r3", "cloglog_low_r3", "cloglog_up_r3", 
                            "true_cloglog_r3", "err_id_instead_of_cloglog_r3", "err_log_instead_of_cloglog_r3", "err_logit_instead_of_cloglog_r3", "err_arcsine_instead_of_cloglog_r3", 
                            
                            "survival_rounded_r4", "cloglog_low_r4",  "cloglog_up_r4", 
                            "true_cloglog_r4", "err_id_instead_of_cloglog_r4", "err_log_instead_of_cloglog_r4", "err_logit_instead_of_cloglog_r4","err_arcsine_instead_of_cloglog_r4")] <-
  
  c(data_cens25_r2[,c("sample_size","survival_initial", "survival_rounded", "cloglog_low", "cloglog_up", 
                     "true_cloglog", "err_id_instead_of_cloglog", "err_log_instead_of_cloglog", "err_logit_instead_of_cloglog", "err_arcsine_instead_of_cloglog")],
    
    data_cens25_r3[,c("survival_rounded", "cloglog_low", "cloglog_up", 
                     "true_cloglog", "err_id_instead_of_cloglog", "err_log_instead_of_cloglog", "err_logit_instead_of_cloglog", "err_arcsine_instead_of_cloglog")],
    
    data_cens25_r4[,c("survival_rounded", "cloglog_low", "cloglog_up", 
                     "true_cloglog", "err_id_instead_of_cloglog", "err_log_instead_of_cloglog", "err_logit_instead_of_cloglog", "err_arcsine_instead_of_cloglog")])



df_LRA_fA_logit_cens25[, c("sample_size", "survival_initial", "survival_rounded_r2","logit_low_r2", "logit_up_r2",
                          "true_logit_r2", "err_id_instead_of_logit_r2", "err_log_instead_of_logit_r2", "err_cloglog_instead_of_logit_r2", "err_arcsine_instead_of_logit_r2",  
                          
                          "survival_rounded_r3", "logit_low_r3", "logit_up_r3", 
                          "true_logit_r3", "err_id_instead_of_logit_r3", "err_log_instead_of_logit_r3", "err_cloglog_instead_of_logit_r3", "err_arcsine_instead_of_logit_r3", 
                          
                          "survival_rounded_r4", "logit_low_r4",  "logit_up_r4", 
                          "true_logit_r4", "err_id_instead_of_logit_r4", "err_log_instead_of_logit_r4", "err_cloglog_instead_of_logit_r4","err_arcsine_instead_of_logit_r4")] <-
  
  c(data_cens25_r2[,c("sample_size","survival_initial", "survival_rounded", "logit_low", "logit_up", 
                     "true_logit", "err_id_instead_of_logit", "err_log_instead_of_logit", "err_cloglog_instead_of_logit", "err_arcsine_instead_of_logit")],
    
    data_cens25_r3[,c("survival_rounded", "logit_low", "logit_up", 
                     "true_logit", "err_id_instead_of_logit", "err_log_instead_of_logit", "err_cloglog_instead_of_logit", "err_arcsine_instead_of_logit")],
    
    data_cens25_r4[,c("survival_rounded", "logit_low", "logit_up", 
                     "true_logit", "err_id_instead_of_logit", "err_log_instead_of_logit", "err_cloglog_instead_of_logit", "err_arcsine_instead_of_logit")])



df_LRA_fA_arcsine_cens25[, c("sample_size", "survival_initial", "survival_rounded_r2","arcsine_low_r2", "arcsine_up_r2",
                            "true_arcsine_r2", "err_id_instead_of_arcsine_r2", "err_log_instead_of_arcsine_r2", "err_cloglog_instead_of_arcsine_r2", "err_logit_instead_of_arcsine_r2",  
                            
                            "survival_rounded_r3", "arcsine_low_r3", "arcsine_up_r3", 
                            "true_arcsine_r3", "err_id_instead_of_arcsine_r3", "err_log_instead_of_arcsine_r3", "err_cloglog_instead_of_arcsine_r3", "err_logit_instead_of_arcsine_r3", 
                            
                            "survival_rounded_r4", "arcsine_low_r4",  "arcsine_up_r4", 
                            "true_arcsine_r4", "err_id_instead_of_arcsine_r4", "err_log_instead_of_arcsine_r4", "err_cloglog_instead_of_arcsine_r4","err_logit_instead_of_arcsine_r4")] <-
  
  c(data_cens25_r2[,c("sample_size","survival_initial", "survival_rounded", "arcsine_low", "arcsine_up", 
                     "true_arcsine", "err_id_instead_of_arcsine", "err_log_instead_of_arcsine", "err_cloglog_instead_of_arcsine", "err_logit_instead_of_arcsine")],
    
    data_cens25_r3[,c("survival_rounded", "arcsine_low", "arcsine_up", 
                     "true_arcsine", "err_id_instead_of_arcsine", "err_log_instead_of_arcsine", "err_cloglog_instead_of_arcsine", "err_logit_instead_of_arcsine")],
    
    data_cens25_r4[,c("survival_rounded", "arcsine_low", "arcsine_up", 
                     "true_arcsine", "err_id_instead_of_arcsine", "err_log_instead_of_arcsine", "err_cloglog_instead_of_arcsine", "err_logit_instead_of_arcsine")])



## Censoring = 50% ----

df_LRA_fA_id_cens50[, c("sample_size", "survival_initial", "survival_rounded_r2","id_low_r2", "id_up_r2",
                        "true_id_r2", "err_log_instead_of_id_r2", "err_cloglog_instead_of_id_r2", "err_logit_instead_of_id_r2", "err_arcsine_instead_of_id_r2",  
                        
                        "survival_rounded_r3", "id_low_r3", "id_up_r3", 
                        "true_id_r3", "err_log_instead_of_id_r3", "err_cloglog_instead_of_id_r3", "err_logit_instead_of_id_r3", "err_arcsine_instead_of_id_r3", 
                        
                        "survival_rounded_r4", "id_low_r4",  "id_up_r4", 
                        "true_id_r4", "err_log_instead_of_id_r4", "err_cloglog_instead_of_id_r4", "err_logit_instead_of_id_r4","err_arcsine_instead_of_id_r4")] <-
  
  c(data_cens50_r2[,c("sample_size","survival_initial", "survival_rounded", "id_low", "id_up", 
                      "true_id", "err_log_instead_of_id", "err_cloglog_instead_of_id", "err_logit_instead_of_id", "err_arcsine_instead_of_id")],
    
    data_cens50_r3[,c("survival_rounded", "id_low", "id_up", 
                      "true_id", "err_log_instead_of_id", "err_cloglog_instead_of_id", "err_logit_instead_of_id", "err_arcsine_instead_of_id")],
    
    data_cens50_r4[,c("survival_rounded", "id_low", "id_up", 
                      "true_id", "err_log_instead_of_id", "err_cloglog_instead_of_id", "err_logit_instead_of_id", "err_arcsine_instead_of_id")])



df_LRA_fA_log_cens50[, c("sample_size", "survival_initial", "survival_rounded_r2","log_low_r2", "log_up_r2",
                         "true_log_r2", "err_id_instead_of_log_r2", "err_cloglog_instead_of_log_r2", "err_logit_instead_of_log_r2", "err_arcsine_instead_of_log_r2",  
                         
                         "survival_rounded_r3", "log_low_r3", "log_up_r3", 
                         "true_log_r3", "err_id_instead_of_log_r3", "err_cloglog_instead_of_log_r3", "err_logit_instead_of_log_r3", "err_arcsine_instead_of_log_r3", 
                         
                         "survival_rounded_r4", "log_low_r4",  "log_up_r4", 
                         "true_log_r4", "err_id_instead_of_log_r4", "err_cloglog_instead_of_log_r4", "err_logit_instead_of_log_r4","err_arcsine_instead_of_log_r4")] <-
  
  c(data_cens50_r2[,c("sample_size","survival_initial", "survival_rounded", "log_low", "log_up", 
                      "true_log", "err_id_instead_of_log", "err_cloglog_instead_of_log", "err_logit_instead_of_log", "err_arcsine_instead_of_log")],
    
    data_cens50_r3[,c("survival_rounded", "log_low", "log_up", 
                      "true_log", "err_id_instead_of_log", "err_cloglog_instead_of_log", "err_logit_instead_of_log", "err_arcsine_instead_of_log")],
    
    data_cens50_r4[,c("survival_rounded", "log_low", "log_up", 
                      "true_log", "err_id_instead_of_log", "err_cloglog_instead_of_log", "err_logit_instead_of_log", "err_arcsine_instead_of_log")])



df_LRA_fA_cloglog_cens50[, c("sample_size", "survival_initial", "survival_rounded_r2","cloglog_low_r2", "cloglog_up_r2",
                             "true_cloglog_r2", "err_id_instead_of_cloglog_r2", "err_log_instead_of_cloglog_r2", "err_logit_instead_of_cloglog_r2", "err_arcsine_instead_of_cloglog_r2",  
                             
                             "survival_rounded_r3", "cloglog_low_r3", "cloglog_up_r3", 
                             "true_cloglog_r3", "err_id_instead_of_cloglog_r3", "err_log_instead_of_cloglog_r3", "err_logit_instead_of_cloglog_r3", "err_arcsine_instead_of_cloglog_r3", 
                             
                             "survival_rounded_r4", "cloglog_low_r4",  "cloglog_up_r4", 
                             "true_cloglog_r4", "err_id_instead_of_cloglog_r4", "err_log_instead_of_cloglog_r4", "err_logit_instead_of_cloglog_r4","err_arcsine_instead_of_cloglog_r4")] <-
  
  c(data_cens50_r2[,c("sample_size","survival_initial", "survival_rounded", "cloglog_low", "cloglog_up", 
                      "true_cloglog", "err_id_instead_of_cloglog", "err_log_instead_of_cloglog", "err_logit_instead_of_cloglog", "err_arcsine_instead_of_cloglog")],
    
    data_cens50_r3[,c("survival_rounded", "cloglog_low", "cloglog_up", 
                      "true_cloglog", "err_id_instead_of_cloglog", "err_log_instead_of_cloglog", "err_logit_instead_of_cloglog", "err_arcsine_instead_of_cloglog")],
    
    data_cens50_r4[,c("survival_rounded", "cloglog_low", "cloglog_up", 
                      "true_cloglog", "err_id_instead_of_cloglog", "err_log_instead_of_cloglog", "err_logit_instead_of_cloglog", "err_arcsine_instead_of_cloglog")])



df_LRA_fA_logit_cens50[, c("sample_size", "survival_initial", "survival_rounded_r2","logit_low_r2", "logit_up_r2",
                           "true_logit_r2", "err_id_instead_of_logit_r2", "err_log_instead_of_logit_r2", "err_cloglog_instead_of_logit_r2", "err_arcsine_instead_of_logit_r2",  
                           
                           "survival_rounded_r3", "logit_low_r3", "logit_up_r3", 
                           "true_logit_r3", "err_id_instead_of_logit_r3", "err_log_instead_of_logit_r3", "err_cloglog_instead_of_logit_r3", "err_arcsine_instead_of_logit_r3", 
                           
                           "survival_rounded_r4", "logit_low_r4",  "logit_up_r4", 
                           "true_logit_r4", "err_id_instead_of_logit_r4", "err_log_instead_of_logit_r4", "err_cloglog_instead_of_logit_r4","err_arcsine_instead_of_logit_r4")] <-
  
  c(data_cens50_r2[,c("sample_size","survival_initial", "survival_rounded", "logit_low", "logit_up", 
                      "true_logit", "err_id_instead_of_logit", "err_log_instead_of_logit", "err_cloglog_instead_of_logit", "err_arcsine_instead_of_logit")],
    
    data_cens50_r3[,c("survival_rounded", "logit_low", "logit_up", 
                      "true_logit", "err_id_instead_of_logit", "err_log_instead_of_logit", "err_cloglog_instead_of_logit", "err_arcsine_instead_of_logit")],
    
    data_cens50_r4[,c("survival_rounded", "logit_low", "logit_up", 
                      "true_logit", "err_id_instead_of_logit", "err_log_instead_of_logit", "err_cloglog_instead_of_logit", "err_arcsine_instead_of_logit")])



df_LRA_fA_arcsine_cens50[, c("sample_size", "survival_initial", "survival_rounded_r2","arcsine_low_r2", "arcsine_up_r2",
                             "true_arcsine_r2", "err_id_instead_of_arcsine_r2", "err_log_instead_of_arcsine_r2", "err_cloglog_instead_of_arcsine_r2", "err_logit_instead_of_arcsine_r2",  
                             
                             "survival_rounded_r3", "arcsine_low_r3", "arcsine_up_r3", 
                             "true_arcsine_r3", "err_id_instead_of_arcsine_r3", "err_log_instead_of_arcsine_r3", "err_cloglog_instead_of_arcsine_r3", "err_logit_instead_of_arcsine_r3", 
                             
                             "survival_rounded_r4", "arcsine_low_r4",  "arcsine_up_r4", 
                             "true_arcsine_r4", "err_id_instead_of_arcsine_r4", "err_log_instead_of_arcsine_r4", "err_cloglog_instead_of_arcsine_r4","err_logit_instead_of_arcsine_r4")] <-
  
  c(data_cens50_r2[,c("sample_size","survival_initial", "survival_rounded", "arcsine_low", "arcsine_up", 
                      "true_arcsine", "err_id_instead_of_arcsine", "err_log_instead_of_arcsine", "err_cloglog_instead_of_arcsine", "err_logit_instead_of_arcsine")],
    
    data_cens50_r3[,c("survival_rounded", "arcsine_low", "arcsine_up", 
                      "true_arcsine", "err_id_instead_of_arcsine", "err_log_instead_of_arcsine", "err_cloglog_instead_of_arcsine", "err_logit_instead_of_arcsine")],
    
    data_cens50_r4[,c("survival_rounded", "arcsine_low", "arcsine_up", 
                      "true_arcsine", "err_id_instead_of_arcsine", "err_log_instead_of_arcsine", "err_cloglog_instead_of_arcsine", "err_logit_instead_of_arcsine")])

df_LRA_fA_id_cens0$f_author <- df_LRA_fA_id_cens10$f_author <- df_LRA_fA_id_cens25$f_author <- df_LRA_fA_id_cens50$f_author <- rep("id", nrow(df_LRA_fA_id_cens50))
df_LRA_fA_log_cens0$f_author <- df_LRA_fA_log_cens10$f_author <- df_LRA_fA_log_cens25$f_author <- df_LRA_fA_log_cens50$f_author <- rep("log", nrow(df_LRA_fA_log_cens50))
df_LRA_fA_cloglog_cens0$f_author <- df_LRA_fA_cloglog_cens10$f_author <- df_LRA_fA_cloglog_cens25$f_author <- df_LRA_fA_cloglog_cens50$f_author <- rep("cloglog", nrow(df_LRA_fA_cloglog_cens50))
df_LRA_fA_logit_cens0$f_author <- df_LRA_fA_logit_cens10$f_author <- df_LRA_fA_logit_cens25$f_author <- df_LRA_fA_logit_cens50$f_author <- rep("logit", nrow(df_LRA_fA_logit_cens50))
df_LRA_fA_arcsine_cens0$f_author <- df_LRA_fA_arcsine_cens10$f_author <- df_LRA_fA_arcsine_cens25$f_author <- df_LRA_fA_arcsine_cens50$f_author <- rep("arcsine", nrow(df_LRA_fA_arcsine_cens50))


# LRA computation for each fE ----

list_df_LRA <- as.list(ls(pattern = "^df_LRA_fA_"))
fE = c("id", "log", "cloglog", "logit", "arcsine")
cens = c("cens0", "cens10", "cens25", "cens50")
rounding = c("r2", "r3", "r4")


for (i in fE) {
  for (j in cens) {
    datafr <- get(paste0("df_LRA_fA_", i, "_", j)) 
    for (k in rounding) {
      
      # fE = id
      diff_up_fE_id  = datafr[paste0(i,"_up_",k)] - datafr[paste0("survival_rounded_",k)]  # distance between the upper bound and the survival
      diff_low_fE_id = datafr[paste0("survival_rounded_",k)] - datafr[paste0(i,"_low_",k)] # distance between the lower bound and the survival
      datafr[paste0("LRA_fE_id_",k)] = log(diff_up_fE_id/diff_low_fE_id)
      
      # fE = log
      diff_up_fE_log  = log(datafr[paste0(i,"_up_",k)]) - log(datafr[paste0("survival_rounded_",k)])
      diff_low_fE_log = log(datafr[paste0("survival_rounded_",k)]) - log(datafr[paste0(i,"_low_",k)])
      datafr[paste0("LRA_fE_log_",k)] = log(diff_up_fE_log/diff_low_fE_log)
      
      # fE = cloglog                                           
      diff_up_fE_cloglog  = log(-log(datafr[paste0("survival_rounded_",k)])) - log(-log(datafr[paste0(i,"_up_",k)]))
      diff_low_fE_cloglog = log(-log(datafr[paste0(i,"_low_",k)])) - log(-log(datafr[paste0("survival_rounded_",k)]))
      datafr[paste0("LRA_fE_cloglog_",k)] = log(diff_up_fE_cloglog/diff_low_fE_cloglog)
      
      # fE = logit
      diff_up_fE_logit  = logit(datafr[paste0("survival_rounded_",k)], adjust = F) - logit(datafr[paste0(i,"_up_",k)], adjust = F)
      diff_low_fE_logit = logit(datafr[paste0(i,"_low_",k)], adjust = F) - logit(datafr[paste0("survival_rounded_",k)], adjust = F)
      datafr[paste0("LRA_fE_logit_",k)] = log(diff_up_fE_logit/diff_low_fE_logit)
      
      # fE = arcsine
      diff_up_fE_arcsine  = asin(sqrt(datafr[paste0("survival_rounded_",k)])) - asin(sqrt(datafr[paste0(i,"_up_",k)]))
      diff_low_fE_arcsine = asin(sqrt(datafr[paste0(i,"_low_",k)])) - asin(sqrt(datafr[paste0("survival_rounded_",k)]))
      datafr[paste0("LRA_fE_arcsine_",k)] = log(diff_up_fE_arcsine/diff_low_fE_arcsine)
      
      
    }
    assign(paste0("df_LRA_fA_", i,"_", j), as.data.frame(datafr))
  }
}

# LRA performance ----
#Find the transformation leading to the smallest LRA

for (i in fE) {
  for (j in cens) {
    
    datafr <- get(paste0("df_LRA_fA_", i, "_", j)) 
    
    for (k in rounding) {
        datafr[, paste0("fE_givenby_LRA_",k)] <-  ifelse(apply(abs(datafr[,paste0(c("LRA_fE_id_", "LRA_fE_log_","LRA_fE_cloglog_", "LRA_fE_logit_", "LRA_fE_arcsine_"), k)]), 1, min) == abs(datafr[,paste0("LRA_fE_id_",k)]), "id",
                                                  ifelse(apply(abs(datafr[,paste0(c("LRA_fE_id_", "LRA_fE_log_","LRA_fE_cloglog_", "LRA_fE_logit_", "LRA_fE_arcsine_"), k)]), 1, min) == abs(datafr[,paste0("LRA_fE_log_",k)]), "log",
                                                  ifelse(apply(abs(datafr[,paste0(c("LRA_fE_id_", "LRA_fE_log_","LRA_fE_cloglog_", "LRA_fE_logit_", "LRA_fE_arcsine_"), k)]), 1, min) == abs(datafr[,paste0("LRA_fE_cloglog_",k)]),"cloglog",
                                                  ifelse(apply(abs(datafr[,paste0(c("LRA_fE_id_", "LRA_fE_log_","LRA_fE_cloglog_", "LRA_fE_logit_", "LRA_fE_arcsine_"), k)]), 1, min) == abs(datafr[,paste0("LRA_fE_logit_",k)]), "logit",
                                                  ifelse(apply(abs(datafr[,paste0(c("LRA_fE_id_", "LRA_fE_log_","LRA_fE_cloglog_", "LRA_fE_logit_", "LRA_fE_arcsine_"), k)]), 1, min) == abs(datafr[,paste0("LRA_fE_arcsine_",k)]), "arcsine", NA)))))
        
        datafr[,paste0("success_LRA_",k)] <- ifelse(datafr[,paste0("fE_givenby_LRA_",k)] == datafr[,"f_author"], 1,0)
        
      }
      assign(paste0("df_LRA_fA_", i,"_", j), data.frame(datafr))
    }
    
   
}


# Remove truncated cases

for (i in fE) {
  for (j in cens) {
    
    datafr <- get(paste0("df_LRA_fA_", i, "_", j)) 
    
    for (k in rounding) {
      for (l in 1:nrow(datafr)) {
        
        if (datafr[l, paste0(fE, "_up_", rounding)] == 1 | datafr[l, paste0(fE, "_low_", rounding)] == 0 | 
            datafr[l, paste0("survival_rounded_", rounding)] == datafr[l, paste0(fE, "_up_", rounding)] |
            datafr[l, paste0("survival_rounded_", rounding)] == datafr[l, paste0(fE, "_low_", rounding)]|
            datafr[l, paste0(fE, "_up_", rounding)] == datafr[l, paste0(fE, "_low_", rounding)]){
          datafr$fE_givenby_LRA[l] = "not applicable"
          datafr$success_LRA[i] = "not applicable"}
        
      }
      
    }assign(paste0("df_LRA_fA_", i, "_", j), as.data.frame(df_log), envir = .GlobalEnv)
  }
}
  



#Extraction error after LRA implementation ----
#Create a new extraction_error_afterLRA for all dataframes
list_df_LRA_fA <- ls(pattern = "^df_LRA_fA_")

for (df_name in list_df_LRA_fA) {
  df <- get(df_name)
  
  df$extraction_error_afterLRA_r2 <- numeric(nrow(df))
  df$extraction_error_afterLRA_r3 <- numeric(nrow(df))
  df$extraction_error_afterLRA_r4 <- numeric(nrow(df))
  
  assign(df_name, as.data.frame(df))
}

df_LRA_id_names <- ls(pattern = "^df_LRA_fA_id_") # Dataframes names starting with df_LRA_fA_id in my environment
for (k in df_LRA_id_names) {
  df_id <- get(k)
  for (i in rounding) {
    df_id[, paste0("extraction_error_afterLRA_",i)]= ifelse(df_id[, paste0("fE_givenby_LRA_",i)] =="id", df_id[, paste0("true_id_",i)],
                                                      ifelse(df_id[, paste0("fE_givenby_LRA_",i)] =="log", df_id[, paste0("err_log_instead_of_id_",i)],
                                                      ifelse(df_id[, paste0("fE_givenby_LRA_",i)] =="cloglog", df_id[, paste0("err_cloglog_instead_of_id_",i)],
                                                      ifelse(df_id[, paste0("fE_givenby_LRA_",i)] =="logit", df_id[, paste0("err_logit_instead_of_id_",i)],
                                                      ifelse(df_id[, paste0("fE_givenby_LRA_",i)] =="arcsine", df_id[, paste0("err_arcsine_instead_of_id_",i)], NA)))))
  } 
  assign(k, as.data.frame(df_id), envir = .GlobalEnv)
}
  


df_LRA_log_names <- ls(pattern = "^df_LRA_fA_log_") # Dataframes names starting with df_LRA_fA_log in my environment
for (k in df_LRA_log_names) {
  df_log <- get(k)
  for (i in rounding) {
      df_log[, paste0("extraction_error_afterLRA_",i)]=  ifelse(df_log[, paste0("fE_givenby_LRA_",i)] =="log", df_log[,paste0("true_log_",i)],
                                                         ifelse(df_log[, paste0("fE_givenby_LRA_",i)] =="id", df_log[,paste0("err_id_instead_of_log_",i)],
                                                         ifelse(df_log[, paste0("fE_givenby_LRA_",i)] =="cloglog", df_log[,paste0("err_cloglog_instead_of_log_",i)],
                                                         ifelse(df_log[, paste0("fE_givenby_LRA_",i)] =="logit", df_log[,paste0("err_logit_instead_of_log_",i)],
                                                         ifelse(df_log[, paste0("fE_givenby_LRA_",i)] =="arcsine", df_log[,paste0("err_arcsine_instead_of_log_",i)], NA)))))
    }
  assign(k, as.data.frame(df_log), envir = .GlobalEnv)
  }



df_LRA_cloglog_names <- ls(pattern = "^df_LRA_fA_cloglog_") # Dataframes names starting with df_LRA_fA_cloglog in my environment
for (k in df_LRA_cloglog_names) {
  df_cloglog <- get(k)
  for (i in rounding) {
      df_cloglog[, paste0("extraction_error_afterLRA_",i)]= ifelse(df_cloglog[, paste0("fE_givenby_LRA_",i)] =="cloglog", df_cloglog[,paste0("true_cloglog_",i)],
                                                             ifelse(df_cloglog[, paste0("fE_givenby_LRA_",i)] =="log", df_cloglog[,paste0("err_log_instead_of_cloglog_",i)],
                                                             ifelse(df_cloglog[, paste0("fE_givenby_LRA_",i)] =="id", df_cloglog[,paste0("err_id_instead_of_cloglog_",i)],
                                                             ifelse(df_cloglog[, paste0("fE_givenby_LRA_",i)] =="logit", df_cloglog[,paste0("err_logit_instead_of_cloglog_",i)],
                                                             ifelse(df_cloglog[, paste0("fE_givenby_LRA_",i)] =="arcsine", df_cloglog[,paste0("err_arcsine_instead_of_cloglog_",i)], NA)))))
  
  }
  assign(k, as.data.frame(df_cloglog), envir = .GlobalEnv)
}


df_LRA_logit_names <- ls(pattern = "^df_LRA_fA_logit_") # Dataframes names starting with df_LRA_fA_logit in my environment
for (k in df_LRA_logit_names) {
  df_logit <- get(k)
  for (i in rounding) {
      df_logit[, paste0("extraction_error_afterLRA_",i)]= ifelse(df_logit[, paste0("fE_givenby_LRA_",i)] =="logit", df_logit[,paste0("true_logit_",i)],
                                                           ifelse(df_logit[, paste0("fE_givenby_LRA_",i)] =="log", df_logit[,paste0("err_log_instead_of_logit_",i)],
                                                           ifelse(df_logit[, paste0("fE_givenby_LRA_",i)] =="id", df_logit[,paste0("err_id_instead_of_logit_",i)],
                                                           ifelse(df_logit[, paste0("fE_givenby_LRA_",i)] =="cloglog", df_logit[,paste0("err_cloglog_instead_of_logit_",i)],
                                                           ifelse(df_logit[, paste0("fE_givenby_LRA_",i)] =="arcsine", df_logit[,paste0("err_arcsine_instead_of_logit_",i)], NA)))))
    }
  
  
  assign(k, as.data.frame(df_logit), envir = .GlobalEnv)
}


df_LRA_arcsine_names <- ls(pattern = "^df_LRA_fA_arcsine_") # Dataframes names starting with df_LRA_fA_arcsine in my environment
for (k in df_LRA_arcsine_names) {
  df_arcsine <- get(k)
  for (i in rounding) {
      df_arcsine[, paste0("extraction_error_afterLRA_",i)]= ifelse(df_arcsine[, paste0("fE_givenby_LRA_",i)] =="arcsine", df_arcsine[,paste0("true_arcsine_",i)],
                                                             ifelse(df_arcsine[, paste0("fE_givenby_LRA_",i)] =="log", df_arcsine[,paste0("err_log_instead_of_arcsine_",i)],
                                                             ifelse(df_arcsine[, paste0("fE_givenby_LRA_",i)] =="id", df_arcsine[,paste0("err_id_instead_of_arcsine_",i)],
                                                             ifelse(df_arcsine[, paste0("fE_givenby_LRA_",i)] =="cloglog", df_arcsine[,paste0("err_cloglog_instead_of_arcsine_",i)],
                                                             ifelse(df_arcsine[, paste0("fE_givenby_LRA_",i)] =="logit", df_arcsine[,paste0("err_logit_instead_of_arcsine_",i)], NA)))))
    }
  
  
  assign(k, as.data.frame(df_arcsine), envir = .GlobalEnv)
}


## Create categories of extraction error after using LRA ----
list_df_LRA_fA <- ls(pattern = "^df_LRA_fA_")
for (df_name in list_df_LRA_fA) {
  df <- get(df_name)
  
  df$extraction_error_afterLRA_r2_cat <- numeric(nrow(df))
  df$extraction_error_afterLRA_r3_cat <- numeric(nrow(df))
  df$extraction_error_afterLRA_r4_cat <- numeric(nrow(df))
  
  assign(df_name, as.data.frame(df))
}


for (df_name in ls(pattern = "^df_LRA_fA")) {
  df <- get(df_name)
  for (i in rounding) {
    
  df[, paste0("extraction_error_afterLRA_",i,"_cat")]  = ifelse(is.na(df[,paste0("extraction_error_afterLRA_",i)]), "not applicable",
                                                      ifelse(df[,paste0("extraction_error_afterLRA_",i)] < -50, 1,
                                                      ifelse(df[,paste0("extraction_error_afterLRA_",i)] < -25, 2,
                                                      ifelse(df[,paste0("extraction_error_afterLRA_",i)] < -10, 3,
                                                      ifelse(df[,paste0("extraction_error_afterLRA_",i)] <  -5, 4,
                                                      ifelse(df[,paste0("extraction_error_afterLRA_",i)] < - 1, 5,
                                                      ifelse(df[,paste0("extraction_error_afterLRA_",i)] <   1, 6,
                                                      ifelse(df[,paste0("extraction_error_afterLRA_",i)] <   5, 7,
                                                      ifelse(df[,paste0("extraction_error_afterLRA_",i)] <  10, 8,
                                                      ifelse(df[,paste0("extraction_error_afterLRA_",i)] <  25, 9,
                                                      ifelse(df[,paste0("extraction_error_afterLRA_",i)] <= 50, 10, 11)))))))))))
  
    
  df[,paste0("extraction_error_afterLRA_",i,"_cat")] <- as.character(df[,paste0("extraction_error_afterLRA_",i,"_cat")])
  assign(df_name, as.data.frame(df), envir = .GlobalEnv)
}
}

