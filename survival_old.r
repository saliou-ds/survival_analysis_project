# Load and install libraries
#install.packages("survminer")
#install.packages("FactoMineR")
#install.packages("explor")
library(FactoMineR)
library(explor)
library(survival)
library(survminer)
library(dplyr)
library(car)

# Read the CSV files

train <- read.csv("dataset_paper/financial_train.csv")
validation <- read.csv("dataset_paper/financial_validation.csv")
test <- read.csv("dataset_paper/financial_test.csv")

#recode
train <- train %>% mutate(status_label = ifelse(status_label == "alive", 0, 1)) # recode(status_label, "alive" = 0, "failed" = 1) doesnt work anymore, dont know why (because conflict with packages / same fonction name)
validation <- validation %>% mutate(status_label = ifelse(status_label == "alive", 0, 1))
test <- test %>% mutate(status_label = ifelse(status_label == "alive", 0, 1))

train <- train %>% mutate(fyear = fyear - min(fyear) + 1)
validation <- validation %>% mutate(fyear = fyear - min(fyear) + 1)
test <- test %>% mutate(fyear = fyear - min(fyear) + 1)

train <- rbind(train, validation)

# Fit the Kaplan-Meier survival curve on the training set
surv_obj <- Surv(train$fyear, train$status_label) # creates a survival object using time and event indicator, default use Kaplan-Meier 
km_fit <- survfit(surv_obj ~ 1, data = train) # creates survival curves 
                                            # we would have 'group' instead of '1', if there were groups, with group being a variable 

ggsurvplot(km_fit, 
           data = train,
           conf.int = TRUE, 
           risk.table = TRUE, # number of subjects at risk at different time points
           #pval = TRUE,    # when multiple groups automatically performs a log-rank test and displays the p-value, which tests the null hypothesis that there is no difference between the survival curves 
           
           ggtheme = theme_survminer(),
           title = "Kaplan-Meier Survival Curve (Training Data)",
           risk.table.height = 0.35,  # Adjust the height of the risk table
           #xlim = c(2000, 2014),
           break.time.by = 2, # Set the x-axis limits
           break.y.by = 0.25)  # Set the y-axis steps


# Fit a Cox proportional hazards model

# First we're going to choose our variables

#Identify strong predictors 

predictor_vars <- colnames(train)[4:ncol(train)]

univariate_results <- lapply(predictor_vars, function(x) { # univariate Cox regression for each predictor to identify strong predictors     #lapply always give list
  summary(coxph(as.formula(paste("Surv(fyear, status_label) ~", x)), data = train)) # as.formula converts a string into a formula object. coxph need a a formula object
})

print(univariate_results[1:5])

significant_predictors <- predictor_vars[sapply(univariate_results, function(x) x$coefficients[,"Pr(>|z|)"] < 0.05)] #  p-value <0.05 from the Wald test  # sapply Simplifies the result to a vector or matrix if possible.

# We're going to test now their multicolinearity 

cor_matrix <- cor(train[, significant_predictors])
View(cor_matrix)

high_cor_pairs <- abs(cor_matrix) > 0.7
diag(high_cor_pairs) <- FALSE
high_cor_pairs <- which(high_cor_pairs, arr.ind = TRUE) # each row represents the emplacement (row and column) of a pair where cor_matrice > 0.7

low_cor_significant_predictors <- significant_predictors[- unique(high_cor_pairs[, "col"]) ] # I remove one predictor from every pair with high correlation

# Cox model
cox_model <- coxph(Surv(train$fyear, train$status_label) ~ ., data = train[, low_cor_significant_predictors])
# Summarize the model
summary(cox_model)

#X1_net_income alone was significant but with X3_net_income  it's not so 

cox_model <- coxph(Surv(fyear, status_label) ~ X3_net_income, data = train)
summary(cox_model)


# Check the proportional hazards assumption
ph_test <- cox.zph(cox_model)
print(ph_test)
plot(ph_test)

# Evaluate the model on the test set
test$predicted_risk <- predict(cox_model, newdata = test, type = "risk")
surv_obj_test <- Surv(test$fyear, test$status_label)
concordance_test <- concordance(surv_obj_test ~ test$predicted_risk)
print(concordance_test)


# keep all significant_predictors
res.pca <- PCA(train[, significant_predictors], graph = FALSE)
#explor(res.pca)

# Extract the principal components
pca_components <- as.data.frame(res.pca$ind$coord[,c("Dim.1","Dim.4","Dim.5")])

# Add the principal components to the combined training set
train_pca <- cbind(train, pca_components)

# Fit the Cox model using the principal components
cox_model_pca <- coxph(Surv(fyear, status_label) ~ ., data = train_pca[, c("fyear", "status_label", colnames(pca_components))])
summary(cox_model_pca)

# Check the proportional hazards assumption
ph_test_pca <- cox.zph(cox_model_pca)
print(ph_test_pca)
plot(ph_test_pca)

# Evaluate the model on the test set
test_pca <- predict(res.pca, newdata = test[, significant_predictors])
test_pca <- as.data.frame(test_pca$coord)
test_pca <- cbind(test, test_pca)


test_pca$predicted_risk <- predict(cox_model_pca, newdata = test_pca, type = "risk")
surv_obj_test <- Surv(test_pca$fyear, test_pca$status_label)
concordance_test_pca <- concordance(surv_obj_test ~ test_pca$predicted_risk)
print(concordance_test_pca)

# By using PCA, it's a little better but not enough


######################################################

# Use of log-abs for my low_cor_significant_predictors

# Read the CSV files
train <- read.csv("dataset_paper/financial_train.csv")
validation <- read.csv("dataset_paper/financial_validation.csv")
test <- read.csv("dataset_paper/financial_test.csv")

# Recode status_label using ifelse
train <- train %>% mutate(status_label = ifelse(status_label == "alive", 0, 1))
validation <- validation %>% mutate(status_label = ifelse(status_label == "alive", 0, 1))
test <- test %>% mutate(status_label = ifelse(status_label == "alive", 0, 1))

# Apply log(abs(x) + 1) transformation to all variables from columns 4 to the last column
train[, 4:ncol(train)] <- log(abs(train[, 4:ncol(train)]) + 1)
validation[, 4:ncol(validation)] <- log(abs(validation[, 4:ncol(validation)]) + 1)
test[, 4:ncol(test)] <- log(abs(test[, 4:ncol(test)]) + 1)

train <- rbind(train, validation)

# Fit the Cox model
cox_model <- coxph(Surv(train$fyear, train$status_label) ~ ., data = train[, low_cor_significant_predictors])
summary(cox_model)

# Check the proportional hazards assumption
ph_test <- cox.zph(cox_model)
print(ph_test)
plot(ph_test)

# Evaluate the model on the test set
test$predicted_risk <- predict(cox_model, newdata = test, type = "risk")
surv_obj_test <- Surv(test$fyear, test$status_label)
concordance_test <- concordance(surv_obj_test ~ test$predicted_risk)
print(concordance_test)

