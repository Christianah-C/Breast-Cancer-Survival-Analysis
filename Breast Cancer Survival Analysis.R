#Import the data
data <- read.csv("C:\\Users\\Christianah.O_BROOKS\\Downloads\\Breast Cancer METABRIC.csv", header = TRUE, stringsAsFactors = FALSE)
head(data)

#Data Cleaning
Clean_data <- na.omit(data)
Clean_data
View(Clean_data)
cat("Original number of rows:", nrow(data), "\n")
cat("Number of rows after removing missing values:", nrow(Clean_data), "\n")

#Load libraries
library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(survminer)
library(tidyr)

#Kaplan-Mier Analysis
Clean_data$PatientStatus <- ifelse(Clean_data$status == "Living", 1, 0)
km <- with(Clean_data, Surv(time, PatientStatus))
head(km, 80)

km_fit <- survfit(Surv(time, PatientStatus) ~ 1, data=Clean_data)
summary(km_fit, times = c(1,30,60,90*(1:10)))

km_trt_fit <- survfit(Surv(time, PatientStatus) ~ Chemotherapy, data=Clean_data)
autoplot(km_trt_fit)

km_trt_fit <- survfit(Surv(time, PatientStatus) ~ Hormone.Therapy, data=Clean_data)
autoplot(km_trt_fit)

km_trt_fit <- survfit(Surv(time, PatientStatus) ~ Radio.Therapy, data=Clean_data)
autoplot(km_trt_fit)

km_trt_fit <- survfit(Surv(time, PatientStatus) ~ Type.of.Breast.Surgery, data=Clean_data)
autoplot(km_trt_fit)

#Cox Proportional Hazard Model
cox_model <- coxph(Surv(time, PatientStatus) ~ Age + Tumor.Size + ER.Status + 
                     PR.Status + HER2.Status + Chemotherapy + Tumor.Stage, data = Clean_data)
summary(cox_model)

cox_fit <- survfit(cox_model)
autoplot(cox_fit)

#Aalen’s additive regression model
aa_fit <-aareg(Surv(time, PatientStatus) ~ Age + Tumor.Size + ER.Status + 
                     PR.Status + HER2.Status + Chemotherapy + Tumor.Stage, data = Clean_data)
summary(aa_fit)
autoplot(aa_fit)


#Random forest model
r_fit <- ranger(Surv(time, PatientStatus) ~ Age + Tumor.Size + ER.Status + 
                  PR.Status + HER2.Status + Chemotherapy + Tumor.Stage,
                data = Clean_data,
                mtry = 4,
                importance = "permutation",
                splitrule = "extratrees",
                verbose = TRUE)
r_fit

# Average the survival models
death_times <- r_fit$unique.death.times 
surv_prob <- data.frame(r_fit$survival)
avg_prob <- sapply(surv_prob,mean)

# Plot the survival models for each patient
plot(r_fit$unique.death.times,r_fit$survival[1,], 
     type = "l", 
     ylim = c(0,1),
     col = "red",
     xlab = "Days",
     ylab = "survival",
     main = "Patient Survival Curves")

cols <- colors()
for (n in sample(c(2:dim(Clean_data)[1]), 20)){
  lines(r_fit$unique.death.times, r_fit$survival[n,], type = "l", col = cols[n])
}

lines(death_times, avg_prob, lwd = 2)
legend(500, 0.7, legend = c('Average = black'))

vi <- data.frame(sort(round(r_fit$variable.importance, 4), decreasing = TRUE))
names(vi) <- "importance"
head(vi)

cat("Prediction Error = 1 - Harrell's c-index = ", r_fit$prediction.error)


#ggplot eyeball comparism of models
kmi <- rep("KM",length(km_fit$time))
km_df <- data.frame(km_fit$time,km_fit$surv,kmi)
names(km_df) <- c("Time","Surv","Model")

coxi <- rep("Cox",length(cox_fit$time))
cox_df <- data.frame(cox_fit$time,cox_fit$surv,coxi)
names(cox_df) <- c("Time","Surv","Model")

rfi <- rep("RF",length(r_fit$unique.death.times))
rf_df <- data.frame(r_fit$unique.death.times,avg_prob,rfi)
names(rf_df) <- c("Time","Surv","Model")

plot_df <- rbind(km_df,cox_df,rf_df)

p <- ggplot(plot_df, aes(x = Time, y = Surv, color = Model))
p + geom_line()


# Exponential survival model
s <- with(Clean_data,Surv(time,PatientStatus))
fKM <- survfit(s ~ Type.of.Breast.Surgery,data=Clean_data)
sExp <- survreg(s ~ as.factor(Type.of.Breast.Surgery),dist='exp',data=Clean_data)
summary(sExp)

pred.Type.of.Breast.Surgery1 = predict(sExp, newdata=list(Type.of.Breast.Surgery="Breast Conserving"),type="quantile",p=seq(.01,.99,by=.01))
pred.Type.of.Breast.Surgery2 = predict(sExp, newdata=list(Type.of.Breast.Surgery="Mastectomy"),type="quantile",p=seq(.01,.99,by=.01))

df = data.frame(y=seq(.99,.01,by=-.01), Type.of.Breast.Surgery1=pred.Type.of.Breast.Surgery1, Type.of.Breast.Surgery2=pred.Type.of.Breast.Surgery2)
df_long = gather(df, key= "Type.of.Breast.Surgery", value="time", -y)

p = ggsurvplot(fKM, data = Clean_data, risk.table = T)
p$plot = p$plot + geom_line(data=df_long, aes(x=time, y=y, group=Type.of.Breast.Surgery))
p$plot

p$table

#Weibull survival model
s <- with(Clean_data,Surv(time,PatientStatus))
fKM <- survfit(s ~  Type.of.Breast.Surgery,data=Clean_data)
sWei <- survreg(s ~ as.factor(Type.of.Breast.Surgery), dist = 'weibull', data = Clean_data)
summary(sWei)

pred.Type.of.Breast.Surgery1 = predict(sWei, newdata=list(Type.of.Breast.Surgery="Breast Conserving"),type="quantile",p=seq(.01,.99,by=.01))
pred.Type.of.Breast.Surgery2 = predict(sWei, newdata=list(Type.of.Breast.Surgery="Mastectomy"),type="quantile",p=seq(.01,.99,by=.01))

df = data.frame(y=seq(.99,.01,by=-.01), Type.of.Breast.Surgery1=pred.Type.of.Breast.Surgery1, Type.of.Breast.Surgery2=pred.Type.of.Breast.Surgery2)
df_long = gather(df, key= "Type.of.Breast.Surgery", value="time", -y)

p = ggsurvplot(fKM, data = Clean_data, risk.table = T)
p$plot = p$plot + geom_line(data=df_long, aes(x=time, y=y, group=Type.of.Breast.Surgery))
p$plot
  
p$table      
