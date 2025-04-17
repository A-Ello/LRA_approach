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
library(meta)
library(metafor)
library(officer)
library(flextable)


# Dataset ----

load("your_directory/ExampleMA_ReRTSyst.Rdata")

max(ExampleMA_ReRTSyst$TimeToPFS)  # 43.25671
min(ExampleMA_ReRTSyst$TimeToPFS)  # 0.05590062
#The time is in months

# Survival for each study ----

  ## Bergman 2020 ----
time_Bergman <- ExampleMA_ReRTSyst$TimeToPFS[ExampleMA_ReRTSyst$Studies == "Bergman 2020"]
event_Bergman <- ExampleMA_ReRTSyst$EventPFS[ExampleMA_ReRTSyst$Studies == "Bergman 2020"]

KMfit_Bergman_plain   <- survfit(Surv(time_Bergman, event_Bergman) ~ 1, conf.type="plain")
sum_KMfit_Bergman_plain <- summary(KMfit_Bergman_plain, times = 0:40)

KMfit_Bergman_log     <- survfit(Surv(time_Bergman, event_Bergman) ~ 1, conf.type="log")
sum_KMfit_Bergman_log <- summary(KMfit_Bergman_log, times = 0:40)

KMfit_Bergman_cloglog <- survfit(Surv(time_Bergman, event_Bergman) ~ 1, conf.type="log-log")
sum_KMfit_Bergman_cloglog <- summary(KMfit_Bergman_cloglog, times = 0:40)

KMfit_Bergman_logit   <- survfit(Surv(time_Bergman, event_Bergman) ~ 1, conf.type="logit")
sum_KMfit_Bergman_logit <- summary(KMfit_Bergman_logit, times = 0:40)

KMfit_Bergman_arcsine <- survfit(Surv(time_Bergman, event_Bergman) ~ 1, conf.type="arcsin")
sum_KMfit_Bergman_arcsine <- summary(KMfit_Bergman_arcsine, times = 0:40)


  ## Bovi 2020 ----
time_Bovi <- ExampleMA_ReRTSyst$TimeToPFS[ExampleMA_ReRTSyst$Studies == "Bovi 2020"]
event_Bovi <- ExampleMA_ReRTSyst$EventPFS[ExampleMA_ReRTSyst$Studies == "Bovi 2020"]

KMfit_Bovi_plain   <- survfit(Surv(time_Bovi, event_Bovi) ~ 1, conf.type="plain")
sum_KMfit_Bovi_plain <- summary(KMfit_Bovi_plain, times = 0:40)

KMfit_Bovi_log     <- survfit(Surv(time_Bovi, event_Bovi) ~ 1, conf.type="log")
sum_KMfit_Bovi_log <- summary(KMfit_Bovi_log, times = 0:40)

KMfit_Bovi_cloglog <- survfit(Surv(time_Bovi, event_Bovi) ~ 1, conf.type="log-log")
sum_KMfit_Bovi_cloglog <- summary(KMfit_Bovi_cloglog, times = 0:40)

KMfit_Bovi_logit   <- survfit(Surv(time_Bovi, event_Bovi) ~ 1, conf.type="logit")
sum_KMfit_Bovi_logit <- summary(KMfit_Bovi_logit, times = 0:40)

KMfit_Bovi_arcsine <- survfit(Surv(time_Bovi, event_Bovi) ~ 1, conf.type="arcsin")
sum_KMfit_Bovi_arcsine <- summary(KMfit_Bovi_arcsine, times = 0:40)


  ## Yasuda 2018 ---- 
time_Yasuda <- ExampleMA_ReRTSyst$TimeToPFS[ExampleMA_ReRTSyst$Studies == "Yasuda 2018"]
event_Yasuda <- ExampleMA_ReRTSyst$EventPFS[ExampleMA_ReRTSyst$Studies == "Yasuda 2018"]

KMfit_Yasuda_plain   <- survfit(Surv(time_Yasuda, event_Yasuda) ~ 1, conf.type="plain")
sum_KMfit_Yasuda_plain <- summary(KMfit_Yasuda_plain, times = 0:40)

KMfit_Yasuda_log     <- survfit(Surv(time_Yasuda, event_Yasuda) ~ 1, conf.type="log")
sum_KMfit_Yasuda_log <- summary(KMfit_Yasuda_log, times = 0:40)

KMfit_Yasuda_cloglog <- survfit(Surv(time_Yasuda, event_Yasuda) ~ 1, conf.type="log-log")
sum_KMfit_Yasuda_cloglog <- summary(KMfit_Yasuda_cloglog, times = 0:40)

KMfit_Yasuda_logit   <- survfit(Surv(time_Yasuda, event_Yasuda) ~ 1, conf.type="logit")
sum_KMfit_Yasuda_logit <- summary(KMfit_Yasuda_logit, times = 0:40)

KMfit_Yasuda_arcsine <- survfit(Surv(time_Yasuda, event_Yasuda) ~ 1, conf.type="arcsin")
sum_KMfit_Yasuda_arcsine <- summary(KMfit_Yasuda_arcsine, times = 0:40)


  ## Kim 2015 ----
time_Kim <- ExampleMA_ReRTSyst$TimeToPFS[ExampleMA_ReRTSyst$Studies == "Kim 2015"]
event_Kim <- ExampleMA_ReRTSyst$EventPFS[ExampleMA_ReRTSyst$Studies == "Kim 2015"]

KMfit_Kim_plain   <- survfit(Surv(time_Kim, event_Kim) ~ 1, conf.type="plain")
sum_KMfit_Kim_plain <- summary(KMfit_Kim_plain, times = 0:40)

KMfit_Kim_log     <- survfit(Surv(time_Kim, event_Kim) ~ 1, conf.type="log")
sum_KMfit_Kim_log <- summary(KMfit_Kim_log, times = 0:40)

KMfit_Kim_cloglog <- survfit(Surv(time_Kim, event_Kim) ~ 1, conf.type="log-log")
sum_KMfit_Kim_cloglog <- summary(KMfit_Kim_cloglog, times = 0:40)

KMfit_Kim_logit   <- survfit(Surv(time_Kim, event_Kim) ~ 1, conf.type="logit")
sum_KMfit_Kim_logit <- summary(KMfit_Kim_logit, times = 0:40)

KMfit_Kim_arcsine <- survfit(Surv(time_Kim, event_Kim) ~ 1, conf.type="arcsin")
sum_KMfit_Kim_arcsine <- summary(KMfit_Kim_arcsine, times = 0:40)



  ## Tsien 2022 ----
time_Tsien <- ExampleMA_ReRTSyst$TimeToPFS[ExampleMA_ReRTSyst$Studies == "Tsien 2022"]
event_Tsien <- ExampleMA_ReRTSyst$EventPFS[ExampleMA_ReRTSyst$Studies == "Tsien 2022"]

KMfit_Tsien_plain   <- survfit(Surv(time_Tsien, event_Tsien) ~ 1, conf.type="plain")
sum_KMfit_Tsien_plain <- summary(KMfit_Tsien_plain, times = 0:40)

KMfit_Tsien_log     <- survfit(Surv(time_Tsien, event_Tsien) ~ 1, conf.type="log")
sum_KMfit_Tsien_log <- summary(KMfit_Tsien_log, times = 0:40)

KMfit_Tsien_cloglog <- survfit(Surv(time_Tsien, event_Tsien) ~ 1, conf.type="log-log")
sum_KMfit_Tsien_cloglog <- summary(KMfit_Tsien_cloglog, times = 0:40)

KMfit_Tsien_logit   <- survfit(Surv(time_Tsien, event_Tsien) ~ 1, conf.type="logit")
sum_KMfit_Tsien_logit <- summary(KMfit_Tsien_logit, times = 0:40)

KMfit_Tsien_arcsine <- survfit(Surv(time_Tsien, event_Tsien) ~ 1, conf.type="arcsin")
sum_KMfit_Tsien_arcsine <- summary(KMfit_Tsien_arcsine, times = 0:40)


# Plot all survival curves
plot(KMfit_Bovi_log, col = "blue", lty = 1, lwd = 1.5, conf.int = FALSE)
lines(KMfit_Tsien_log, col = "red", lty = 1, lwd = 1.5, conf.int = FALSE)
lines(KMfit_Bergman_log, col = "green", lty = 1, lwd = 1.5, conf.int = FALSE)
lines(KMfit_Kim_log, col = "purple", lty = 1, lwd = 1.5, conf.int = FALSE)
lines(KMfit_Yasuda_log, col = "orange", lty = 1, lwd = 1.5, conf.int = FALSE)

legend("topright", legend = c("Bovi", "Tsien", "Bergman", "Kim", "Yasuda"), # Labels for each curve
       col = c("blue", "red", "green", "purple", "orange"), lty = 1, lwd = 1.5, cex = 0.8)



# Selected survival times ----
selected_times <- c(4, 6,10)
fE = c("id", "log", "cloglog", "logit", "arcsine")

#Put them in a dataframe and round to 3 decimal places
names <- c("study", "survival", "true_se_id", "id_up", "id_low", 
           "log_up", "log_low", "cloglog_up", "cloglog_low", 
           "logit_up", "logit_low", "arcsine_up", "arcsine_low", 
           
           "se_correct_extr_id", "se_log_instead_of_id", "se_cloglog_instead_of_id", "se_logit_instead_of_id", 
           "se_arcsine_instead_of_id", "se_LRA_instead_of_id",
           
           "se_correct_extr_log", "se_id_instead_of_log", "se_cloglog_instead_of_log", "se_logit_instead_of_log", 
           "se_arcsine_instead_of_log", "se_LRA_instead_of_log",
           
           "se_correct_extr_cloglog", "se_id_instead_of_cloglog", "se_log_instead_of_cloglog", "se_logit_instead_of_cloglog", 
           "se_arcsine_instead_of_cloglog", "se_LRA_instead_of_cloglog",
           
           "se_correct_extr_logit", "se_id_instead_of_logit", "se_log_instead_of_logit", "se_cloglog_instead_of_logit", 
           "se_arcsine_instead_of_logit", "se_LRA_instead_of_logit",
           
           "se_correct_extr_arcsine", "se_id_instead_of_arcsine", "se_log_instead_of_arcsine", "se_cloglog_instead_of_arcsine", 
           "se_logit_instead_of_arcsine", "se_LRA_instead_of_arcsine")


for (i in selected_times) {
  assign(paste0("data_t",i), matrix(0, nrow = length(unique(ExampleMA_ReRTSyst$Studies)), 
                                    ncol = length(names), dimnames = list(NULL, names)))
}

for (i in selected_times) {
  
  mat <- get(paste0("data_t",i))
  

  mat[,"study"] <- unique(ExampleMA_ReRTSyst$Studies)
  
  mat[,"survival"] <- round(c(sum_KMfit_Bergman_log$surv[sum_KMfit_Bergman_log$time==i], 
                          sum_KMfit_Bovi_log$surv[sum_KMfit_Bovi_log$time==i],
                          sum_KMfit_Yasuda_log$surv[sum_KMfit_Yasuda_log$time==i],
                          sum_KMfit_Kim_log$surv[sum_KMfit_Kim_log$time==i],
                          sum_KMfit_Tsien_log$surv[sum_KMfit_Tsien_log$time==i]), 3)

  mat[,"true_se_id"] <- c(sum_KMfit_Bergman_plain$std.err[sum_KMfit_Bergman_plain$time==i], 
                          sum_KMfit_Bovi_plain$std.err[sum_KMfit_Bovi_plain$time==i],
                          sum_KMfit_Yasuda_plain$std.err[sum_KMfit_Yasuda_plain$time==i],
                          sum_KMfit_Kim_plain$std.err[sum_KMfit_Kim_plain$time==i],
                          sum_KMfit_Tsien_plain$std.err[sum_KMfit_Tsien_plain$time==i])

  mat[,"id_up"] <- round(c(sum_KMfit_Bergman_plain$upper[sum_KMfit_Bergman_plain$time==i], 
                       sum_KMfit_Bovi_plain$upper[sum_KMfit_Bovi_plain$time==i],
                       sum_KMfit_Yasuda_plain$upper[sum_KMfit_Yasuda_plain$time==i],
                       sum_KMfit_Kim_plain$upper[sum_KMfit_Kim_plain$time==i],
                       sum_KMfit_Tsien_plain$upper[sum_KMfit_Tsien_plain$time==i]), 3)

  mat[,"id_low"] <- round(c(sum_KMfit_Bergman_plain$lower[sum_KMfit_Bergman_plain$time==i], 
                       sum_KMfit_Bovi_plain$lower[sum_KMfit_Bovi_plain$time==i],
                       sum_KMfit_Yasuda_plain$lower[sum_KMfit_Yasuda_plain$time==i],
                       sum_KMfit_Kim_plain$lower[sum_KMfit_Kim_plain$time==i],
                       sum_KMfit_Tsien_plain$lower[sum_KMfit_Tsien_plain$time==i]), 3)

  mat[,"log_up"] <- round(c(sum_KMfit_Bergman_log$upper[sum_KMfit_Bergman_log$time==i], 
                       sum_KMfit_Bovi_log$upper[sum_KMfit_Bovi_log$time==i],
                       sum_KMfit_Yasuda_log$upper[sum_KMfit_Yasuda_log$time==i],
                       sum_KMfit_Kim_log$upper[sum_KMfit_Kim_log$time==i],
                       sum_KMfit_Tsien_log$upper[sum_KMfit_Tsien_log$time==i]), 3)

  mat[,"log_low"] <- round(c(sum_KMfit_Bergman_log$lower[sum_KMfit_Bergman_log$time==i], 
                        sum_KMfit_Bovi_log$lower[sum_KMfit_Bovi_log$time==i],
                        sum_KMfit_Yasuda_log$lower[sum_KMfit_Yasuda_log$time==i],
                        sum_KMfit_Kim_log$lower[sum_KMfit_Kim_log$time==i],
                        sum_KMfit_Tsien_log$lower[sum_KMfit_Tsien_log$time==i]), 3)

  mat[,"cloglog_up"] <- round(c(sum_KMfit_Bergman_cloglog$upper[sum_KMfit_Bergman_cloglog$time==i], 
                        sum_KMfit_Bovi_cloglog$upper[sum_KMfit_Bovi_cloglog$time==i],
                        sum_KMfit_Yasuda_cloglog$upper[sum_KMfit_Yasuda_cloglog$time==i],
                        sum_KMfit_Kim_cloglog$upper[sum_KMfit_Kim_cloglog$time==i],
                        sum_KMfit_Tsien_cloglog$upper[sum_KMfit_Tsien_cloglog$time==i]), 3)

  mat[,"cloglog_low"] <- round(c(sum_KMfit_Bergman_cloglog$lower[sum_KMfit_Bergman_cloglog$time==i], 
                         sum_KMfit_Bovi_cloglog$lower[sum_KMfit_Bovi_cloglog$time==i],
                         sum_KMfit_Yasuda_cloglog$lower[sum_KMfit_Yasuda_cloglog$time==i],
                         sum_KMfit_Kim_cloglog$lower[sum_KMfit_Kim_cloglog$time==i],
                         sum_KMfit_Tsien_cloglog$lower[sum_KMfit_Tsien_cloglog$time==i]), 3)

  mat[,"logit_up"] <- round(c(sum_KMfit_Bergman_logit$upper[sum_KMfit_Bergman_logit$time==i], 
                            sum_KMfit_Bovi_logit$upper[sum_KMfit_Bovi_logit$time==i],
                            sum_KMfit_Yasuda_logit$upper[sum_KMfit_Yasuda_logit$time==i],
                            sum_KMfit_Kim_logit$upper[sum_KMfit_Kim_logit$time==i],
                            sum_KMfit_Tsien_logit$upper[sum_KMfit_Tsien_logit$time==i]), 3)

  mat[,"logit_low"] <- round(c(sum_KMfit_Bergman_logit$lower[sum_KMfit_Bergman_logit$time==i], 
                             sum_KMfit_Bovi_logit$lower[sum_KMfit_Bovi_logit$time==i],
                             sum_KMfit_Yasuda_logit$lower[sum_KMfit_Yasuda_logit$time==i],
                             sum_KMfit_Kim_logit$lower[sum_KMfit_Kim_logit$time==i],
                             sum_KMfit_Tsien_logit$lower[sum_KMfit_Tsien_logit$time==i]), 3)

  mat[,"arcsine_up"] <- round(c(sum_KMfit_Bergman_arcsine$upper[sum_KMfit_Bergman_arcsine$time==i], 
                          sum_KMfit_Bovi_arcsine$upper[sum_KMfit_Bovi_arcsine$time==i],
                          sum_KMfit_Yasuda_arcsine$upper[sum_KMfit_Yasuda_arcsine$time==i],
                          sum_KMfit_Kim_arcsine$upper[sum_KMfit_Kim_arcsine$time==i],
                          sum_KMfit_Tsien_arcsine$upper[sum_KMfit_Tsien_arcsine$time==i]), 3)

  mat[,"arcsine_low"] <- round(c(sum_KMfit_Bergman_arcsine$lower[sum_KMfit_Bergman_arcsine$time==i], 
                           sum_KMfit_Bovi_arcsine$lower[sum_KMfit_Bovi_arcsine$time==i],
                           sum_KMfit_Yasuda_arcsine$lower[sum_KMfit_Yasuda_arcsine$time==i],
                           sum_KMfit_Kim_arcsine$lower[sum_KMfit_Kim_arcsine$time==i],
                           sum_KMfit_Tsien_arcsine$lower[sum_KMfit_Tsien_arcsine$time==i]), 3)

 assign(paste0("data_t",i), as.data.frame(mat))
}

# Convertir toutes les colonnes en numérique sauf study
data_t4[, -which(colnames(mat) == "study")] <- lapply(data_t4[, -which(colnames(mat) == "study")], as.numeric)
data_t6[, -which(colnames(mat) == "study")] <- lapply(data_t6[, -which(colnames(mat) == "study")], as.numeric)
data_t10[, -which(colnames(mat) == "study")] <- lapply(data_t10[, -which(colnames(mat) == "study")], as.numeric)


# Extraction of standard error ----
# All the SE are on the scale of the survival

    ## Known functions ----
for (i in selected_times) {
  
  mat <- get(paste0("data_t",i))
  
  colnames(mat)[colnames(mat) == "se_correct_extr_id"] <- "se_id_instead_of_id"
  colnames(mat)[colnames(mat) == "se_correct_extr_log"] <- "se_log_instead_of_log"
  colnames(mat)[colnames(mat) == "se_correct_extr_cloglog"] <- "se_cloglog_instead_of_cloglog"
  colnames(mat)[colnames(mat) == "se_correct_extr_logit"] <- "se_logit_instead_of_logit"
  colnames(mat)[colnames(mat) == "se_correct_extr_arcsine"] <- "se_arcsine_instead_of_arcsine"
    
  for (compute in fE) {
 
  # Extraction with id
  mat[,paste0("se_id_instead_of_", compute)]  <- (mat[, paste0(compute,"_up")] - mat[,paste0(compute,"_low")])/(2*qnorm(0.975))
  
  # Extraction with log
  mat[,paste0("se_log_instead_of_", compute)]  <- ((log(mat[, paste0(compute,"_up")]) - log(mat[, paste0(compute,"_low")]))/(2*qnorm(0.975)))*mat[, "survival"]
  
  # Extraction with cloglog
  mat[,paste0("se_cloglog_instead_of_", compute)] <- -((log(-log(mat[, paste0(compute,"_low")])) - log(-log(mat[, paste0(compute,"_up")])))/(2*qnorm(0.975)))*mat[, "survival"]*log(mat[, "survival"])
  
  # Extraction with logit
  mat[,paste0("se_logit_instead_of_",compute)]    <- ((logit(mat[, paste0(compute,"_up")], adjust = F) - logit(mat[, paste0(compute,"_low")], adjust = F))/(2*qnorm(0.975)))*mat[, "survival"]*(1-mat[, "survival"])
  
  # Extraction with arcsine
  mat[,paste0("se_arcsine_instead_of_",compute)]  <- ((asin(sqrt(mat[, paste0(compute,"_up")])) - asin(sqrt(mat[, paste0(compute,"_low")])))/(2*qnorm(0.975)))*2*sqrt(mat[, "survival"]*(1-mat[, "survival"]))
   
  
  # When truncated
   #upper = 1
  up1 = which(mat[,paste0(compute,"_up")]==1)
  mat[up1, paste0("se_id_instead_of_",compute)] = (mat[up1, "survival"] - mat[up1, paste0(compute,"_low")])/(qnorm(0.975)) 
  
  mat[up1, paste0("se_log_instead_of_",compute)] = ((log(mat[up1, "survival"]) - log(mat[up1 ,paste0(compute,"_low")]))/(qnorm(0.975)))*mat[up1, "survival"]
  
  mat[up1, paste0("se_cloglog_instead_of_",compute)] = -((log(-log(mat[up1, paste0(compute,"_low")])) - log(-log(mat[up1, "survival"])))/(qnorm(0.975)))*mat[up1, "survival"]*log(mat[up1, "survival"])
  
  mat[up1, paste0("se_logit_instead_of_",compute)] = ((logit(mat[up1, "survival"], adjust = F) - logit(mat[up1, paste0(compute,"_low")], adjust = F))/(qnorm(0.975)))*mat[up1, "survival"]*(1-mat[up1, "survival"])
  
  mat[up1, paste0("se_arcsine_instead_of_",compute)] =  ((asin(sqrt(mat[up1, "survival"])) - asin(sqrt(mat[up1, paste0(compute,"_low")])))/(qnorm(0.975)))*2*sqrt(mat[up1, "survival"]*(1-mat[up1, "survival"]))
  
  
  #low=0
  low0 = which(mat[,paste0(compute,"_low")]==0)
  mat[low0, paste0("se_id_instead_of_",compute)] = (mat[low0, paste0(compute,"_up")] - mat[low0, "survival"])/(qnorm(0.975)) 
  
  mat[low0, paste0("se_log_instead_of_",compute)] = ((log(mat[low0,  paste0(compute,"_up")]) - log(mat[low0 ,"survival"]))/(qnorm(0.975)))*mat[low0, "survival"]
  
  mat[low0, paste0("se_cloglog_instead_of_",compute)] = -((log(-log(mat[low0,  "survival"])) - log(-log(mat[low0,  paste0(compute,"_up")])))/(qnorm(0.975)))*mat[low0, "survival"]*log(mat[low0, "survival"])
  
  mat[low0, paste0("se_logit_instead_of_",compute)] = ((logit(mat[low0,  paste0(compute,"_up")], adjust = F) - logit(mat[low0, "survival"], adjust = F))/(qnorm(0.975)))*mat[low0, "survival"]*(1-mat[low0, "survival"])
  
  mat[low0, paste0("se_arcsine_instead_of_",compute)] =  ((asin(sqrt(mat[low0,  paste0(compute,"_up")])) - asin(sqrt(mat[low0, "survival"])))/(qnorm(0.975)))*2*sqrt(mat[low0, "survival"]*(1-mat[low0, "survival"]))
  
  
  }
  
  colnames(mat)[colnames(mat) == "se_id_instead_of_id"]           <- "se_correct_extr_id"
  colnames(mat)[colnames(mat) == "se_log_instead_of_log"]         <- "se_correct_extr_log"
  colnames(mat)[colnames(mat) == "se_cloglog_instead_of_cloglog"] <- "se_correct_extr_cloglog"
  colnames(mat)[colnames(mat) == "se_logit_instead_of_logit"]     <- "se_correct_extr_logit"
  colnames(mat)[colnames(mat) == "se_arcsine_instead_of_arcsine"] <- "se_correct_extr_arcsine"
  
    assign(paste0("data_t",i), as.data.frame(mat))
    }

 
    ## LRA ----

for (j in selected_times) {
  
  datafr <- get(paste0("data_t",j))
  
  for (i in fE) {
    
      # fE = id
      diff_up_fE_id  = datafr[, paste0(i,"_up")] - datafr[,"survival"]  # distance between the upper bound and the survival
      diff_low_fE_id = datafr[, "survival"] - datafr[,paste0(i,"_low")] # distance between the lower bound and the survival
      LRA_fe_id = log(diff_up_fE_id/diff_low_fE_id)
      
      # fE = log
      diff_up_fE_log  = log(datafr[,paste0(i,"_up")]) - log(datafr[,paste0("survival")])
      diff_low_fE_log = log(datafr[,paste0("survival")]) - log(datafr[,paste0(i,"_low")])
      LRA_fe_log = log(diff_up_fE_log/diff_low_fE_log)
      
      # fE = cloglog                                           
      diff_up_fE_cloglog  = log(-log(datafr[,paste0("survival")])) - log(-log(datafr[,paste0(i,"_up")]))
      diff_low_fE_cloglog = log(-log(datafr[,paste0(i,"_low")])) - log(-log(datafr[,paste0("survival")]))
      LRA_fe_cloglog = log(diff_up_fE_cloglog/diff_low_fE_cloglog)
      
      # fE = logit
      diff_up_fE_logit  = logit(datafr[, paste0("survival")], adjust = F) - logit(datafr[,paste0(i,"_up")], adjust = F)
      diff_low_fE_logit = logit(datafr[, paste0(i,"_low")], adjust = F) - logit(datafr[,paste0("survival")], adjust = F)
      LRA_fe_logit = log(diff_up_fE_logit/diff_low_fE_logit)
      
      # fE = arcsine
      diff_up_fE_arcsine  = asin(sqrt(datafr[,paste0("survival")])) - asin(sqrt(datafr[,paste0(i,"_up")]))
      diff_low_fE_arcsine = asin(sqrt(datafr[,paste0(i,"_low")])) - asin(sqrt(datafr[,paste0("survival")]))
      LRA_fe_arcsine = log(diff_up_fE_arcsine/diff_low_fE_arcsine)
      
      se_extr_with_id      = (datafr[,paste0(i,"_up")] - datafr[,paste0(i,"_low")])/(2*qnorm(0.975))
      
      se_extr_with_log     = datafr[,paste0("survival")]*((log(datafr[,paste0(i,"_up")]) - log(datafr[,paste0(i,"_low")]))/(2*qnorm(0.975)))
      
      se_extr_with_cloglog = (- datafr[,paste0("survival")]*log(datafr[,paste0("survival")]))*((log(-log(datafr[,paste0(i,"_low")])) - log(-log(datafr[,paste0(i,"_up")])))/(2*qnorm(0.975)))
     
      se_extr_with_logit   = ((logit(datafr[,paste0(i,"_up")], adjust = F) - logit(datafr[,paste0(i,"_low")], adjust = F))/(2*qnorm(0.975)))*datafr[,paste0("survival")]*(1-datafr[,paste0("survival")])
     
      se_extr_with_arcsine = ((asin(sqrt(datafr[,paste0(i,"_up")])) - asin(sqrt(datafr[,paste0(i,"_low")])))/(2*qnorm(0.975)))*2*sqrt(datafr[,paste0("survival")]*(1-datafr[,paste0("survival")]))
      
      datafr[, paste0("se_LRA_instead_of_",i)] <-ifelse(apply(abs(cbind(LRA_fe_id, LRA_fe_log, LRA_fe_cloglog, LRA_fe_logit, LRA_fe_arcsine)), 1, min) == abs(LRA_fe_id), se_extr_with_id,
                                                 ifelse(apply(abs(cbind(LRA_fe_id, LRA_fe_log, LRA_fe_cloglog, LRA_fe_logit, LRA_fe_arcsine)), 1, min) == abs(LRA_fe_log), se_extr_with_log,
                                                 ifelse(apply(abs(cbind(LRA_fe_id, LRA_fe_log, LRA_fe_cloglog, LRA_fe_logit, LRA_fe_arcsine)), 1, min) == abs(LRA_fe_cloglog), se_extr_with_cloglog,
                                                 ifelse(apply(abs(cbind(LRA_fe_id, LRA_fe_log, LRA_fe_cloglog, LRA_fe_logit, LRA_fe_arcsine)), 1, min) == abs(LRA_fe_logit), se_extr_with_logit,
                                                 ifelse(apply(abs(cbind(LRA_fe_id, LRA_fe_log, LRA_fe_cloglog, LRA_fe_logit, LRA_fe_arcsine)), 1, min) == abs(LRA_fe_arcsine), se_extr_with_arcsine, NA)))))
  
      pb = which(datafr[,paste0(i,"_up")]==1 | datafr[,paste0(i,"_low")]==0)
    datafr[pb, paste0("se_LRA_instead_of_",i)] <- NA
    
    }
    assign(paste0("data_t",j), as.data.frame(datafr))
  }
    

# Checkpoint ----

for (i in fE) {
  
 print(data_t10$survival > data_t10[,paste0(fE, "_low")])
 print(data_t10$survival < data_t10[,paste0(fE, "_up")])
 print(data_t10[,paste0(fE, "_up")] > data_t10[,paste0(fE, "_low")])
  
} 

# Checkpoint
data_t10$se_correct_extr_arcsine - data_t10$se_LRA_instead_of_arcsine
data_t10$se_correct_extr_id - data_t10$se_LRA_instead_of_id
data_t10$se_correct_extr_log - data_t10$se_LRA_instead_of_log
data_t10$se_correct_extr_logit - data_t10$se_LRA_instead_of_logit
data_t10$se_correct_extr_cloglog - data_t10$se_LRA_instead_of_cloglog



# Meta analyses on cloglog scale ----
# Take timepoint 10 (t10)
comp_function = c("id", "log", "cloglog", "logit", "arcsine")
extr_function = c("id", "log", "cloglog", "logit", "arcsine", "LRA")

row_names = c(paste0("comp_id_extr_",extr_function), paste0("comp_log_extr_",extr_function),
              paste0("comp_cloglog_extr_",extr_function), paste0("comp_logit_extr_",extr_function),
              paste0("comp_arcsine_extr_",extr_function))
col_names <- c("pooled_est", "lower_bound", "upper_bound")

results_MA_fixed_t10 <- matrix(0, nrow= 30, ncol = 3, dimnames = list(row_names,col_names ))
results_MA_random_t10 <- matrix(0, nrow= 30, ncol = 3, dimnames = list(row_names,col_names ))

#Change columns names for it to match the loop
colnames(data_t10)[colnames(data_t10) == "se_correct_extr_id"] <- "se_id_instead_of_id"
colnames(data_t10)[colnames(data_t10) == "se_correct_extr_log"] <- "se_log_instead_of_log"
colnames(data_t10)[colnames(data_t10) == "se_correct_extr_cloglog"] <- "se_cloglog_instead_of_cloglog"
colnames(data_t10)[colnames(data_t10) == "se_correct_extr_logit"] <- "se_logit_instead_of_logit"
colnames(data_t10)[colnames(data_t10) == "se_correct_extr_arcsine"] <- "se_arcsine_instead_of_arcsine"

for (comp in comp_function) {
  
  for (extr in extr_function) {
    surv = data_t10[,"survival"]
    se   = data_t10[,paste0("se_", extr, "_instead_of_",comp)]
    
    results_MA_fixed_t10[paste0("comp_", comp, "_extr_", extr), "pooled_est"] =
      exp(-exp(summary(rma.uni(log(-log(surv)), vi=(-se/(surv*log(surv)))^2, method="FE"))$b))
    
    #Inverser les bound parce que cloglog
    results_MA_fixed_t10[paste0("comp_", comp, "_extr_", extr), "lower_bound"] =
      exp(-exp(summary(rma.uni(log(-log(surv)), vi=(-se/(surv*log(surv)))^2, method="FE"))$ci.ub))
    
    results_MA_fixed_t10[paste0("comp_", comp, "_extr_", extr), "upper_bound"] =
      exp(-exp(summary(rma.uni(log(-log(surv)), vi=(-se/(surv*log(surv)))^2, method="FE"))$ci.lb))
  
    
    
    results_MA_random_t10[paste0("comp_", comp, "_extr_", extr), "pooled_est"] =
      exp(-exp(summary(rma.uni(log(-log(surv)), vi=(-se/(surv*log(surv)))^2, method="REML"))$b))
    
    #Inverser les bound parce que cloglog
    results_MA_random_t10[paste0("comp_", comp, "_extr_", extr), "lower_bound"] =
      exp(-exp(summary(rma.uni(log(-log(surv)), vi=(-se/(surv*log(surv)))^2, method="REML"))$ci.ub))
    
    results_MA_random_t10[paste0("comp_", comp, "_extr_", extr), "upper_bound"] =
      exp(-exp(summary(rma.uni(log(-log(surv)), vi=(-se/(surv*log(surv)))^2, method="REML"))$ci.lb))
    
    if(sum(is.na(se))>0) {
      
      results_MA_fixed_t10[paste0("comp_", comp, "_extr_", extr), "pooled_est"]   =
      results_MA_fixed_t10[paste0("comp_", comp, "_extr_", extr), "lower_bound"]  =
      results_MA_fixed_t10[paste0("comp_", comp, "_extr_", extr), "upper_bound"]  =
      results_MA_random_t10[paste0("comp_", comp, "_extr_", extr), "pooled_est"]  =
      results_MA_random_t10[paste0("comp_", comp, "_extr_", extr), "lower_bound"] =
      results_MA_random_t10[paste0("comp_", comp, "_extr_", extr), "upper_bound"] = NA
    }
  } 
}


#Change columns names back to original
colnames(data_t10)[colnames(data_t10) == "se_id_instead_of_id"]           <- "se_correct_extr_id"
colnames(data_t10)[colnames(data_t10) == "se_log_instead_of_log"]         <- "se_correct_extr_log"
colnames(data_t10)[colnames(data_t10) == "se_cloglog_instead_of_cloglog"] <- "se_correct_extr_cloglog"
colnames(data_t10)[colnames(data_t10) == "se_logit_instead_of_logit"]     <- "se_correct_extr_logit"
colnames(data_t10)[colnames(data_t10) == "se_arcsine_instead_of_arcsine"] <- "se_correct_extr_arcsine"

View(results_MA_random_t10)
View(results_MA_fixed_t10)



# Create the tables ----
results_MA_fixed_t10_output = as.matrix(paste0(sprintf("%.4f",results_MA_fixed_t10[,"pooled_est"]), " (", sprintf("%.4f",results_MA_fixed_t10[,"lower_bound"]),
                                     " to ", sprintf("%.4f",results_MA_fixed_t10[,"upper_bound"]),")"))

results_MA_random_t10_output = as.matrix(paste0(sprintf("%.4f",results_MA_random_t10[,"pooled_est"]), " (", sprintf("%.4f",results_MA_random_t10[,"lower_bound"]),
                                     " to ", sprintf("%.4f",results_MA_random_t10[,"upper_bound"]),")"))

rownames(results_MA_random_t10_output) <- row_names
rownames(results_MA_fixed_t10_output) <- row_names

outputable_fixed <- matrix(0, nrow = 5, ncol=6, dimnames = list(comp_function, extr_function))

for (i in rownames(results_MA_fixed_t10_output)) {
  pattern_comp <- str_extract(i,"(?<=comp_)[^_]+(?=_extr)")
  pattern_extr <- str_extract(i,"(?<=extr_).*")
  
  outputable_fixed[pattern_comp,pattern_extr] = results_MA_fixed_t10_output[i,]
}


outputable_random <- matrix(0, nrow = 5, ncol=6, dimnames = list(comp_function, extr_function))

for (i in rownames(results_MA_random_t10_output)) {
  pattern_comp <- str_extract(i,"(?<=comp_)[^_]+(?=_extr)")
  pattern_extr <- str_extract(i,"(?<=extr_).*")
  
  outputable_random[pattern_comp,pattern_extr] = results_MA_random_t10_output[i,]
}




outputable_random_flex <- flextable(as.data.frame(cbind(rownames(outputable_random),outputable_random)))
outputable_random_flex <- read_docx() %>% 
  body_add_flextable(outputable_random_flex)
print(outputable_random_flex, 
      target = "your_directory/tableau_random_effects.docx")



outputable_fixed_flex <- flextable(as.data.frame(cbind(rownames(outputable_fixed),outputable_fixed)))
outputable_fixed_flex <- read_docx() %>% 
  body_add_flextable(outputable_fixed_flex)

print(outputable_fixed_flex, 
      target = "your_directory/tableau_fixed_effects.docx")

