# Download necessary pacakges and load libraries
install.packages("GGally")
install.packages("tidyverse")
install.packages("car")
install.packages("mice")         # includes functions for multiple imputation
install.packages("ggmice")       # visualizations for missing data
install.packages("reshape2")     # for converting wide to long format
install.packages("leaps")        # for model selection
install.packages("olsrr")        # Package to assess influential observations
install.packages("multcomp")     # for contrasts and multiple comparisons
install.packages("lmtest")       # for breusch-pagan test
install.packages("glmnet")      # for lasso, ridge, and elastic nets
install.packages("readxl")

library(GGally)
library(tidyverse)  # for ggplot2, dplyr, and others
library(car)
library(mice)
library(ggmice)
library(reshape2)
library(leaps)
library(olsrr)
library(multcomp)
library(lmtest)
library(glmnet)
library(readxl)

# Read and explore dataset
copd <- read_xlsx("/Users/annaywj/Desktop/copd_data.xlsx")
head(copd)

# EDA I 
## summary table
summary(copd) 

## Check missing values
options(repr.plot.width=20, repr.plot.height=10)
plot_pattern(copd)

## Column hay_fever has values 0,1,3, convert into a single categorical variable
head(copd$hay_fever)
copd$fever_level <- as.factor(copd$hay_fever)
levels(copd$fever_level) <- c("No Hay Fever", "Mild Hay Fever", "Several Hay Fever")
head(copd$fever_level)

## Convert all other categorical variables into factors
copd$race <- as.factor(copd$race)
copd$gender <- as.factor(copd$gender)
copd$asthma <- as.factor(copd$asthma)
copd$bronchitis_attack <- as.factor(copd$bronchitis_attack)
copd$pneumonia <- as.factor(copd$pneumonia)
copd$chronic_bronchitis <- as.factor(copd$chronic_bronchitis)
copd$emphysema <- as.factor(copd$emphysema)
copd$copd <- as.factor(copd$copd)
copd$sleep_apnea <- as.factor(copd$sleep_apnea)
copd$smoking_status <- as.factor(copd$smoking_status)

# EDA II
## Plotting before fitting regression (correlations)
num <- select_if(copd, is.numeric)
cor <- cor(num, use = "complete.obs")
melted <- melt(cor)
melted$value <- round(melted$value, 2)

options(repr.plot.width=20, repr.plot.height=12)
ggplot(melted, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  scale_fill_gradient2(low = "#619CFF", high = "#F8766D", mid = "white",
                       limit = c(-1, 1)) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size = 15),
        axis.text.y = element_text(hjust=1, size = 15))

# to narrow down more
options(repr.plot.width=20, repr.plot.height=10)

# removing non-numeric and ID/date columns
copd2 <- na.omit(dplyr::select(copd, -sid, -visit_year, -visit_date))
copd2
# Check for multicollinearity
cor <- cor(select_if(copd2, is.numeric))
high_cor <- which(abs(cor) > 0.70, arr.ind = TRUE)

# Create predictors and response for the best subset selection by BIC
response <- copd2$pct_emphysema
predictors <- model.matrix(pct_emphysema ~ ., data = copd2)[,-1]

# Plot best subsets selection
best_subsets <- regsubsets(predictors, y = response, nbest = 1, nvmax = NULL, really.big = TRUE)
plot(best_subsets, scale = "bic")

# Pairwise scatter plot and correlations
for (i in c(2:5, 7:ncol(copd2))) {
  if (i == 22) {
    next
  }
  print(ggpairs(copd2[,c(i, 22)]))
}

# narrow down copd2 to copd3 based on plots above
copd3 <- na.omit(dplyr::select(copd2, -gender, -height_cm, -sysBP, -hr, -hay_fever, -pneumonia, -emphysema,
                               -race, -weight_kg, -diasBP, -asthma, -chronic_bronchitis, -copd, -SmokStartAge, 
                               -fever_level, -CigPerDaySmokAvg, -Duration_Smoking))
print(names(copd3))

# Dig into correlation between pct_emphysema and other variables \
## visit_age
options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = copd3, aes(visit_age, pct_emphysema)) +
  geom_point() + geom_smooth()

## O2_hours_day
options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = copd3, aes(O2_hours_day, pct_emphysema)) +
  geom_point() + geom_smooth()

fit_linear <- lm(pct_emphysema ~ poly(O2_hours_day, 1), data = copd3)

options(repr.plot.width=8, repr.plot.height=8)
newdata <- data.frame(O2_hours_day = seq(min(copd3$O2_hours_day),
                                         max(copd3$O2_hours_day), 1))
plot(copd3$O2_hours_day, copd3$pct_emphysema)
lines(newdata$O2_hours_day, predict(fit_linear, newdata = newdata),
      col = "blue", lwd = 2)

## BMI
options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = copd3, aes(bmi, pct_emphysema)) +
  geom_point() + geom_smooth()

copd3$bmi25 <- as.factor(as.numeric(copd3$bmi > 25))
levels(copd3$bmi25) <- c("<25", ">25")

options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = copd3, aes(bmi25, pct_emphysema)) +
  geom_boxplot()

## Bronchitis Attack
options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = copd3, aes(bronchitis_attack, pct_emphysema)) +
  geom_boxplot()

## Sleep Apnea
options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = copd3, aes(sleep_apnea, pct_emphysema)) +
  geom_boxplot()

## Smoking Status
options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = copd3, aes(smoking_status, pct_emphysema)) +
  geom_boxplot()

## Total Lung Capacity
options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = copd3, aes(total_lung_capacity, pct_emphysema)) +
  geom_point() + geom_smooth()

## Functional Residual Capacity
options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = copd3, aes(functional_residual_capacity, pct_emphysema)) +
  geom_point() + geom_smooth()

fit_quadratic <- lm(pct_emphysema ~ poly(functional_residual_capacity, 2), data = copd3)

options(repr.plot.width=8, repr.plot.height=8)
newdata <- data.frame(functional_residual_capacity = seq(min(copd3$functional_residual_capacity),
                                                         max(copd3$functional_residual_capacity), 1))
plot(copd3$functional_residual_capacity, copd3$pct_emphysema)
lines(newdata$functional_residual_capacity, predict(fit_quadratic, newdata = newdata),
      col = "red", lwd = 2)

## Percentage of Gastrapping
options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = copd3, aes(pct_gastrapping, pct_emphysema)) +
  geom_point() + geom_smooth()

fit_quadratic <- lm(pct_emphysema ~ poly(pct_gastrapping, 2), data = copd3)

options(repr.plot.width=8, repr.plot.height=8)
newdata <- data.frame(pct_gastrapping = seq(min(copd3$pct_gastrapping),
                                            max(copd3$pct_gastrapping), 1))
plot(copd3$pct_gastrapping, copd3$pct_emphysema)
lines(newdata$pct_gastrapping, predict(fit_quadratic, newdata = newdata),
      col = "red", lwd = 2)

## Insp_Meanatt
options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = copd3, aes(insp_meanatt, pct_emphysema)) +
  geom_point() + geom_smooth()

fit_linear <- lm(pct_emphysema ~ poly(insp_meanatt, 1), data = copd3)

options(repr.plot.width=8, repr.plot.height=8)
newdata <- data.frame(insp_meanatt = seq(min(copd3$insp_meanatt),
                                         max(copd3$insp_meanatt), 1))
plot(copd3$insp_meanatt, copd3$pct_emphysema)
lines(newdata$insp_meanatt, predict(fit_linear, newdata = newdata),
      col = "red", lwd = 2)

## Exp_Meanatt
options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = copd3, aes(exp_meanatt, pct_emphysema)) +
  geom_point() + geom_smooth()

fit_linear <- lm(pct_emphysema ~ poly(exp_meanatt, 1), data = copd3)

options(repr.plot.width=8, repr.plot.height=8)
newdata <- data.frame(exp_meanatt = seq(min(copd3$exp_meanatt),
                                        max(copd3$exp_meanatt), 1))
plot(copd3$exp_meanatt, copd3$pct_emphysema)
lines(newdata$exp_meanatt, predict(fit_linear, newdata = newdata),
      col = "red", lwd = 2)

## FEV1 vs FVC Ratio
options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = copd3, aes(FEV1_FVC_ratio, pct_emphysema)) +
  geom_point() + geom_smooth()

copd3$FEV1_FVC_ratio <- as.factor(as.numeric(copd3$FEV1_FVC_ratio > 0.5))
levels(copd3$FEV1_FVC_ratio) <- c("<0.5", ">0.5")

options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = copd3, aes(FEV1_FVC_ratio, pct_emphysema)) +
  geom_boxplot()

## FEV1
options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = copd3, aes(FEV1, pct_emphysema)) +
  geom_point() + geom_smooth()

copd3$FEV1 <- as.factor(as.numeric(copd3$FEV1 > 1.5))
levels(copd3$FEV1) <- c("<1.5", ">1.5")

options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = copd3, aes(FEV1, pct_emphysema)) +
  geom_boxplot()

## FVC
options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = copd3, aes(FVC, pct_emphysema)) +
  geom_point() + geom_smooth()

fit_linear <- lm(pct_emphysema ~ poly(FVC, 1), data = copd3)
fit_quadratic <- lm(pct_emphysema ~ poly(FVC, 2), data = copd3)

options(repr.plot.width=8, repr.plot.height=8)
newdata <- data.frame(FVC = seq(min(copd3$FVC),
                                         max(copd3$FVC), 1))
plot(copd3$FVC, copd3$pct_emphysema)
lines(newdata$FVC, predict(fit_quadratic, newdata = newdata),
      col = "red", lwd = 2)
lines(newdata$FVC, predict(fit_linear, newdata = newdata),
      col = "blue", lwd = 2)

## FEV1 Phase2
options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = copd3, aes(FEV1_phase2, pct_emphysema)) +
  geom_point() + geom_smooth()

copd3$FEV1_phase2 <- as.factor(as.numeric(copd3$FEV1_phase2 > 1.5))
levels(copd3$FEV1_phase2) <- c("<1.5", ">1.5")

options(repr.plot.width=8, repr.plot.height=8)
ggplot(data = copd3, aes(FEV1_phase2, pct_emphysema)) +
  geom_boxplot()



#Regression Analysis - Multivariable Regression and Model Building
## Fit regression
fit <- lm(pct_emphysema ~
            visit_age +
            O2_hours_day +
            bmi +
            smoking_status +
            poly(total_lung_capacity, 2) +
            poly(functional_residual_capacity, 2) +
            poly(pct_gastrapping, 2) +
            insp_meanatt +
            exp_meanatt +
            FEV1_FVC_ratio +
            FEV1 +
            FVC +
            FEV1_phase2,
          data = copd3)
summary(fit)

## Best Subset Selection to choose better predictors
options(repr.plot.width=20, repr.plot.height=10)
response     <- copd3$pct_emphysema
preds        <- model.matrix(fit)
best_subsets <-  regsubsets(preds, y = response,
                            nbest = 100,    # save the best "nbest" for each number of variables
                            nvmax = 20,    # maximum number of variables allowed in the model
                            really.big=T, # for larger datasets
                            intercept = FALSE)
plot(best_subsets, scale = "bic")

## Step-wise Selection to narrow down more based on AIC
fit_final <- step(fit)
## Choose the final model
fit_final <- lm(pct_emphysema ~
                  visit_age +
                  O2_hours_day +
                  poly(total_lung_capacity, 2) +
                  poly(functional_residual_capacity, 2) +
                  poly(pct_gastrapping, 2) +
                  FEV1_FVC_ratio,
                data = copd3)
summary(fit_final)

# Regression Analysis - Residual Diagnostics 
## Check assumption violations
options(repr.plot.width=15, repr.plot.height=5)
par(mfrow = c(1,3))
plot(fit_final, which = 1:3)

## Check nomality
ks.test(fit_final$residuals, 'pnorm')

## Check for serial correlation
durbinWatsonTest(fit_final)

## Check for heteroskedasticity
bptest(fit_final)
## Weighted Least Squares (WLS)
abs_resid  <- abs(residuals(fit_final))
fitted     <- fitted(fit_final)
fit_resid  <- lm(abs_resid ~ fitted)

shat <- fitted(fit_resid)
w <- 1 / shat^2
fit_wls <- lm(summary(fit_final)$call,
              weights = w, data = copd3)

## Breusch Pagan Test
bptest(fit_wls)

## Plot to check violations again
options(repr.plot.width=15, repr.plot.height=5)
par(mfrow = c(1,3))
plot(fit_wls, which = 1:3)

## Compare the differences
summary(fit_final)
summary(fit_wls)

# Check for presence of multicollinearity
vif(fit_final)
# Some variabels need further investigation
fit <- lm(pct_emphysema ~ total_lung_capacity + I(total_lung_capacity^2), data = copd3)
poly_tlc <- poly(copd3$total_lung_capacity, 2)
fit <- lm(pct_emphysema ~ poly_tlc[,1] + poly_tlc[,2], data = copd3)
summary(fit)
print(vif(fit))

poly_gas <- poly(copd3$pct_gastrapping, 2)
fit <- lm(pct_emphysema ~ poly_gas[,1] + poly_gas[,2], data = copd3)
summary(fit)
fit <- lm(pct_emphysema ~ pct_gastrapping + I(pct_gastrapping^2), data = copd3)
summary(fit)
print(vif(fit))

# Regression Analysis & Diagnosis without missing values 
## check and show any missing values
copd[copd == -1] <- NA
na_count <- colSums(is.na(copd))
print(na_count)

## Remove missing value
copd_clean <- na.omit(copd)

## Convert factors
head(copd_clean$hay_fever)
copd_clean$fever_level <- as.factor(copd_clean$hay_fever)
levels(copd_clean$fever_level) <- c("No Hay Fever", "Mild Hay Fever", "Several Hay Fever")
head(copd_clean$fever_level)

copd_clean$race <- as.factor(copd_clean$race)
copd_clean$gender <- as.factor(copd_clean$gender)
copd_clean$asthma <- as.factor(copd_clean$asthma)
copd_clean$bronchitis_attack <- as.factor(copd_clean$bronchitis_attack)
copd_clean$pneumonia <- as.factor(copd_clean$pneumonia)
copd_clean$chronic_bronchitis <- as.factor(copd_clean$chronic_bronchitis)
copd_clean$emphysema <- as.factor(copd_clean$emphysema)
copd_clean$copd <- as.factor(copd_clean$copd)
copd_clean$sleep_apnea <- as.factor(copd_clean$sleep_apnea)
copd_clean$smoking_status <- as.factor(copd_clean$smoking_status)

## Narrow down variables like before
copd2_clean <- na.omit(dplyr::select(copd_clean, -sid, -visit_year, -visit_date))
copd3_clean <- na.omit(dplyr::select(copd2, -gender, -height_cm, -sysBP, -hr, -hay_fever, -pneumonia, -emphysema,
                                     -race, -weight_kg, -diasBP, -asthma, -chronic_bronchitis, -copd, -SmokStartAge,
                                     -fever_level, -CigPerDaySmokAvg, -Duration_Smoking))
summary(copd3_clean)

## Fit the model without missing values
fit_final_clean <- lm(pct_emphysema ~
                        visit_age +
                        O2_hours_day +
                        poly(total_lung_capacity, 2) +
                        poly(functional_residual_capacity, 2) +
                        poly(pct_gastrapping, 2) +
                        FEV1_FVC_ratio,
                      data = copd3_clean)
summary(fit_final_clean)

## Check Assumption Violations
options(repr.plot.width=15, repr.plot.height=5)
par(mfrow = c(1,3))
plot(fit_final_clean, which = 1:3)

## Check normality
ks.test(fit_final_clean$residuals, 'pnorm')

## Check Serial Correlation
durbinWatsonTest(fit_final_clean)

## Breusch Pagan Test
bptest(fit_final_clean)

## Transform Reponse Variable
fit_final_clean2 <- lm(log(pct_emphysema) ~
                         visit_age +
                         O2_hours_day +
                         poly(total_lung_capacity, 2) +
                         poly(functional_residual_capacity, 2) +
                         poly(pct_gastrapping, 2) +
                         FEV1_FVC_ratio,
                       data = copd3_clean)
summary(fit_final_clean2)

options(repr.plot.width=15, repr.plot.height=5)
par(mfrow = c(1,3))
plot(fit_final_clean2, which = 1:3)

## Check for presence of outliers
options(repr.plot.width=15, repr.plot.height=5)
par(mfrow = c(1,3))
plot(fit_final, which = 4:6)

options(repr.plot.width=15, repr.plot.height=8)
ols_plot_resid_lev(fit_final)

options(repr.plot.width=15, repr.plot.height=8)
ols_plot_cooksd_bar(fit_final)

options(repr.plot.width=20, repr.plot.height=8)
ols_plot_dfbetas(fit_final)

### Try to remove one outlier to see the difference
fit1 <- lm(summary(fit_final)$call, data = copd3[-7,])
summary(fit1)
### compare to fit_final
summary(fit_final)

## Determine the necessity of transformaiton using AIC and ANOVA
fit_notran <- lm(pct_emphysema ~  
                   O2_hours_day +
                   poly(total_lung_capacity, 2) +
                   poly(functional_residual_capacity, 2) +
                   poly(pct_gastrapping, 2) +
                   FEV1_FVC_ratio,
                 data = copd3)

options(repr.plot.width=15, repr.plot.height=5)
par(mfrow = c(1,3))

AIC(fit_final)
AIC(fit_notran)

anova(fit_notran, fit_final)

summary(fit_notran)

## Added Variable Plots
options(repr.plot.width=20, repr.plot.height=20)
avPlots(fit_notran)

## Contrast to test the different variables
poly_mort <- poly(copd3$pct_gastrapping, 2)
gastrapping0 <- gastrapping80 <- model.matrix(fit_final)[1,]

gastrapping0[2:3] <- poly(0, 2, coefs = attr(poly_mort, "coefs"))
gastrapping80[2:3] <- poly(80, 2, coefs = attr(poly_mort, "coefs"))

K <- matrix(gastrapping0 - gastrapping80, 1)
K
cont <- glht(fit_final, linfct = K)
summary(cont)

poly_mort <- poly(copd3$functional_residual_capacity, 2)
residual0 <- residual8 <- model.matrix(fit_final)[1,]

residual0[4:5] <- poly(0, 2, coefs = attr(poly_mort, "coefs"))
residual8[4:5] <- poly(8, 2, coefs = attr(poly_mort, "coefs"))

K <- matrix(residual0 - residual8, 1)
K
cont <- glht(fit_final, linfct = K)
summary(cont)

poly_mort <- poly(copd3$total_lung_capacity, 2)
capacity0 <- capacity12 <- model.matrix(fit_final)[1,]

capacity0[6:7] <- poly(0, 2, coefs = attr(poly_mort, "coefs"))
capacity12[6:7] <- poly(12, 2, coefs = attr(poly_mort, "coefs"))

K <- matrix(capacity0 - capacity12, 1)
K
cont <- glht(fit_final, linfct = K)
summary(cont)


## Confidence Intervals of the Model using Bonferroni
n <- length(coef(fit_final)) - 1
level <- 1 - (0.05/n)
confint(fit_notran, level = level)
