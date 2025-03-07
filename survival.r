# Load and install libraries
#install.packages("survminer")
#install.packages("caret")
#install.packages("glmnet")
library(glmnet)
library(caret)
library(survival)
library(survminer)
library(dplyr)
library(car)


# Read the CSV files

df <- read.csv("dataset/tumor_dataset.csv")

set.seed(221)  # For reproducibility
trainIndex <- createDataPartition(df$event, p = 0.8, list = FALSE)

# Split the data
train <- df[trainIndex, ]
test <- df[-trainIndex, ]


#recode
train <- train %>% select(-num_stage)
train$fac_in_subcohort <- ifelse(train$fac_in_subcohort == 'True' , 0, 1)
train$fac_instit <- ifelse(train$fac_instit == 'Favourable' , 0, 1)
train$fac_histol <- ifelse(train$fac_histol == 'Favourable' , 0, 1)
train$fac_study <- ifelse(train$fac_study == 3 , 0, 1)

test <- test %>% select(-num_stage)
test$fac_in_subcohort <- ifelse(test$fac_in_subcohort == 'True' , 0, 1)
test$fac_instit <- ifelse(test$fac_instit == 'Favourable' , 0, 1)
test$fac_histol <- ifelse(test$fac_histol == 'Favourable' , 0, 1)
test$fac_study <- ifelse(test$fac_study == 3 , 0, 1)


# Fit the Kaplan-Meier survival curve on the training set


surv_obj <- Surv(train$time, train$event) # creates a survival object using time and event indicator, default use Kaplan-Meier 
km_fit <- survfit(surv_obj ~ 1, data = train)
# we would have 'group' instead of '1', if there were groups, with group being a variable 
ggsurvplot(km_fit, 
           data = train,
           conf.int = TRUE, 
           risk.table = TRUE, # number of subjects at risk at different time points
           #pval = TRUE,    # when multiple groups automatically performs a log-rank test and displays the p-value, which tests the null hypothesis that there is no difference between the survival curves 
           
           ggtheme = theme_survminer(),
           title = "Kaplan-Meier Survival Curve (Data Frame)",
           risk.table.height = 0.35,  # Adjust the height of the risk table
           #xlim = c(2000, 2014),
           xlab = "Time (months)",
           break.time.by = 1000, # Set the x-axis limits
           break.y.by = 0.25)  # Set the y-axis steps


# For in_subcohort
km_fit <- survfit(surv_obj ~ fac_in_subcohort, data = train) # creates survival curves 
ggsurvplot(km_fit, 
           data = train,
           conf.int = TRUE, 
           risk.table = TRUE, # number of subjects at risk at different time points
           pval = TRUE,    # when multiple groups automatically performs a log-rank test and displays the p-value, which tests the null hypothesis that there is no difference between the survival curves 
           
           ggtheme = theme_survminer(),
           title = "Kaplan-Meier Survival Curve (Data Frame)",
           risk.table.height = 0.35,  # Adjust the height of the risk table
           #xlim = c(2000, 2014),
           xlab = "Time (months)",
           break.time.by = 1000, # Set the x-axis limits
           break.y.by = 0.25)  # Set the y-axis steps


# For instit
km_fit <- survfit(surv_obj ~ fac_instit, data = train)
ggsurvplot(km_fit, 
           data = train,
           conf.int = TRUE, 
           risk.table = TRUE, # number of subjects at risk at different time points
           pval = TRUE,    # when multiple groups automatically performs a log-rank test and displays the p-value, which tests the null hypothesis that there is no difference between the survival curves 
           
           ggtheme = theme_survminer(),
           title = "Kaplan-Meier Survival Curve (Data Frame)",
           risk.table.height = 0.35,  # Adjust the height of the risk table
           #xlim = c(2000, 2014),
           xlab = "Time (months)",
           break.time.by = 1000, # Set the x-axis limits
           break.y.by = 0.25)  # Set the y-axis steps


# For histol
km_fit <- survfit(surv_obj ~ fac_histol, data = train)
ggsurvplot(km_fit, 
           data = train,
           conf.int = TRUE, 
           risk.table = TRUE, # number of subjects at risk at different time points
           pval = TRUE,    # when multiple groups automatically performs a log-rank test and displays the p-value, which tests the null hypothesis that there is no difference between the survival curves 
           
           ggtheme = theme_survminer(),
           title = "Kaplan-Meier Survival Curve (Data Frame)",
           risk.table.height = 0.35,  # Adjust the height of the risk table
           #xlim = c(2000, 2014),
           xlab = "Time (months)",
           break.time.by = 1000, # Set the x-axis limits
           break.y.by = 0.25)  # Set the y-axis steps


# For stage
km_fit <- survfit(surv_obj ~ fac_stage, data = train)
ggsurvplot(km_fit, 
           data = train,
           conf.int = TRUE, 
           risk.table = TRUE, # number of subjects at risk at different time points
           pval = TRUE,    # when multiple groups automatically performs a log-rank test and displays the p-value, which tests the null hypothesis that there is no difference between the survival curves 
           
           ggtheme = theme_survminer(),
           title = "Kaplan-Meier Survival Curve (Data Frame)",
           risk.table.height = 0.35,  # Adjust the height of the risk table
           #xlim = c(2000, 2014),
           xlab = "Time (months)",
           break.time.by = 1000, # Set the x-axis limits
           break.y.by = 0.25)  # Set the y-axis steps

# For study
km_fit <- survfit(surv_obj ~ fac_study, data = train)
ggsurvplot(km_fit, 
           data = train,
           conf.int = TRUE, 
           risk.table = TRUE, # number of subjects at risk at different time points
           pval = TRUE,    # when multiple groups automatically performs a log-rank test and displays the p-value, which tests the null hypothesis that there is no difference between the survival curves 
           
           ggtheme = theme_survminer(),
           title = "Kaplan-Meier Survival Curve (Data Frame)",
           risk.table.height = 0.35,  # Adjust the height of the risk table
           #xlim = c(2000, 2014),
           xlab = "Time (months)",
           break.time.by = 1000, # Set the x-axis limits
           break.y.by = 0.25)  # Set the y-axis steps


# Fit a Cox proportional hazards model

# First we're going to choose our variables

#Identify strong predictors 

predictor_vars <- colnames(train)[4:ncol(train)]

univariate_results <- lapply(predictor_vars, function(x) { # univariate Cox regression for each predictor to identify strong predictors     #lapply always give list
  summary(coxph(as.formula(paste("Surv(time, event) ~", x)), data = train)) # as.formula converts a string into a formula object. coxph need a a formula object
})

print(univariate_results)

significant_predictors <- c("num_age", "fac_stage", "fac_study", "fac_instit","fac_histol") #  p-value <0.05 from the Wald test  # sapply Simplifies the result to a vector or matrix if possible.

print(significant_predictors)
#All the categorial variable that have significant differences into them are significant predictors

# We're going to test now their multicolinearity 

cor_matrix <- cor(train[, significant_predictors])
View(cor_matrix)

# fac_histol and fac_instit are correlated. We are going to keep the better predictor.
summary(coxph(Surv(time, event) ~ fac_instit, data = train))
summary(coxph(Surv(time, event) ~ fac_histol, data = train))           

# Elles ont la même significativité, on va choisir fac_histol because in the context of the study it should be a more accurate variable 


low_cor_significant_predictors <- significant_predictors[significant_predictors != "fac_instit"] # I remove one predictor from every pair with high correlation
print(low_cor_significant_predictors)

low_cor_significant_predictors <- c("num_age" ,"fac_stage"  , "fac_histol")

# Cox model
cox_model <- coxph(Surv(time, event) ~  num_age + strata(fac_stage) + strata(fac_histol), data = train)
# Summarize the model
summary(cox_model)


# Evaluate the model on the test set
test$predicted_risk <- predict(cox_model, newdata = test, type = "risk")
surv_obj_test <- Surv(test$time, test$event)
concordance_test <- concordance(surv_obj_test ~ test$predicted_risk)
print(concordance_test)

# Check the proportional hazards assumption
ph_test <- cox.zph(cox_model)
print(ph_test)
plot(ph_test)

# Cox model
cox_model <- coxph(Surv(train$time, train$event) ~ tt(num_age) + strata(fac_stage) + strata(fac_histol) , data = train, tt = function(x, t, ...) (-x) * log(t) )
# Summarize the model
summary(cox_model)

# Evaluate the model on the test set
test$predicted_risk <- predict(cox_model, newdata = test, type = "risk")
surv_obj_test <- Surv(test$time, test$event)
concordance_test <- concordance(surv_obj_test ~ test$predicted_risk)
print(concordance_test)








x <- model.matrix(Surv(time, event) ~ ., data = train)[, -1]
y <- Surv(train$time, train$event)
lasso_model <- cv.glmnet(x, y, family = "cox", alpha = 1)
ridge_model <- cv.glmnet(x, y, family = "cox", alpha = 0)
elastic_net_model <- cv.glmnet(x, y, family = "cox", alpha = 0.5)

# Evaluate the model on the test set
test <- test %>% select(-predicted_risk)
test_matrix <- model.matrix(~ ., data = test %>% select(-event, -time) )[, -1]
test$predicted_risk <- predict(lasso_model, newx = test_matrix, type = "link")
surv_obj_test <- Surv(test$time, test$event)
concordance_test <- concordance(surv_obj_test ~ test$predicted_risk)
print(concordance_test)

# Evaluate the model on the test set
test <- test %>% select(-predicted_risk)
test_matrix <- model.matrix(~ ., data = test %>% select(-event, -time) )[, -1]
test$predicted_risk <- predict(ridge_model, newx = test_matrix, type = "link")
surv_obj_test <- Surv(test$time, test$event)
concordance_test <- concordance(surv_obj_test ~ test$predicted_risk)
print(concordance_test)

# Evaluate the model on the test set
test <- test %>% select(-predicted_risk)
test_matrix <- model.matrix(~ ., data = test %>% select(-event, -time) )[, -1]
test$predicted_risk <- predict(elastic_net_model, newx = test_matrix, type = "link")
surv_obj_test <- Surv(test$time, test$event)
concordance_test <- concordance(surv_obj_test ~ test$predicted_risk)
print(concordance_test)

