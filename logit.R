## author: yao

# Load necessary libraries
library(stats)

# Load data from CSV file
phenodata_50660_42861_all_clean <- read.csv('./phenodata_50660_42861_all_clean.csv')

# Subset the data to include only Smoking values of 0 and 3
phenodata_50660_42861_all_clean_only0and3 <- subset(phenodata_50660_42861_all_clean, Smoking %in% c(0, 3))

# Convert 'Smoking' column to factor for better interpretation
phenodata_50660_42861_all_clean_only0and3$Smoking <- as.factor(phenodata_50660_42861_all_clean_only0and3$Smoking)

# predictor variable
X <- phenodata_50660_42861_all_clean_only0and3$Predicted_minus_Chronological.lr.top500

# Dependent variable
y <- phenodata_50660_42861_all_clean_only0and3$Binary_Smoking

# Fit the logistic regression model
model <- glm(y ~ X, family = "binomial")

# Print the regression summary
summary(model)

McFadden_Pseudo_R2 = 1 - model$deviance / model$null.deviance # works for glm
