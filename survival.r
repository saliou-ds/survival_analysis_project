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

df <- read.csv("dataset/Financial Distress.csv")

#recode
df$Financial.Distress <- ifelse(df$Financial.Distress < -0.5, 1, df$Financial.Distress)


# Fit the Kaplan-Meier survival curve on the training set
surv_obj <- Surv(df$Time, df$Financial.Distress) # creates a survival object using time and event indicator, default use Kaplan-Meier 
km_fit <- survfit(surv_obj ~ 1, data = df) # creates survival curves 
# we would have 'group' instead of '1', if there were groups, with group being a variable 

ggsurvplot(km_fit, 
           data = df,
           conf.int = TRUE, 
           risk.table = TRUE, # number of subjects at risk at different time points
           #pval = TRUE,    # when multiple groups automatically performs a log-rank test and displays the p-value, which tests the null hypothesis that there is no difference between the survival curves 
           
           ggtheme = theme_survminer(),
           title = "Kaplan-Meier Survival Curve (Data Frame)",
           risk.table.height = 0.35,  # Adjust the height of the risk table
           #xlim = c(2000, 2014),
           break.time.by = 2, # Set the x-axis limits
           break.y.by = 0.25)  # Set the y-axis steps


# Fit a Cox proportional hazards model

# First we're going to choose our variables

#Identify strong predictors 

predictor_vars <- colnames(df)[4:ncol(df)]

univariate_results <- lapply(predictor_vars, function(x) { # univariate Cox regression for each predictor to identify strong predictors     #lapply always give list
  summary(coxph(as.formula(paste("Surv(Time, Financial.Distress) ~", x)), data = df)) # as.formula converts a string into a formula object. coxph need a a formula object
})

print(univariate_results[1:5])

significant_predictors <- predictor_vars[sapply(univariate_results, function(x) x$coefficients[,"Pr(>|z|)"] < 0.05)] #  p-value <0.05 from the Wald test  # sapply Simplifies the result to a vector or matrix if possible.

# We're going to test now their multicolinearity 

cor_matrix <- cor(df[, significant_predictors])
View(cor_matrix)

high_cor_pairs <- abs(cor_matrix) > 0.7
diag(high_cor_pairs) <- FALSE
high_cor_pairs <- which(high_cor_pairs, arr.ind = TRUE) # each row represents the emplacement (row and column) of a pair where cor_matrice > 0.7

low_cor_significant_predictors <- significant_predictors[- unique(high_cor_pairs[, "col"]) ] # I remove one predictor from every pair with high correlation

# Standardize Predictors
df[, low_cor_significant_predictors] <- scale(df[, low_cor_significant_predictors])

# Cox model
cox_model <- coxph(Surv(df$Time, df$Financial.Distress) ~ ., data = df[, low_cor_significant_predictors])
# Summarize the model
summary(cox_model)


# Check the proportional hazards assumption
ph_test <- cox.zph(cox_model)
print(ph_test)
plot(ph_test)


#X1_net_income alone was significant but with X3_net_income  it's not so 
cox_summary <- summary(cox_model)
final_predictors <- rownames(cox_summary$coefficients)[cox_summary$coefficients[, "Pr(>|z|)"] < 0.05]








ph_test_results <- data.frame(
  predictor = rownames(ph_test$table),
  p_value = ph_test$table[, "p"]
)

# Filter predictors with p-value > 0.05
final_predictors <- ph_test_results$predictor[ph_test_results$p_value > 0.05]

# Remove the 'GLOBAL' row if present
final_predictors <- final_predictors[final_predictors != "GLOBAL"]

# Print the final predictors
print(final_predictors)











# Cox model
cox_model <- coxph(Surv(df$Time, df$Financial.Distress) ~ ., data = df[, final_predictors])
# Summarize the model
summary(cox_model)

# Check the proportional hazards assumption
ph_test <- cox.zph(cox_model)
print(ph_test)
plot(ph_test)



