##### author: 姚博瀚

#### Part 1
### Download gse40279 and run analysis by using "MEAL" package (from original code 4.R)
library(MEAL)
library(minfi)
library(limma)
library(ggplot2)

### Installation of GEOquery
## Approach 1: From my github
# Install remotes from CRAN
install.packages("remotes")
# Download modified GEOquery package from my github 
# by using function(install_github()) from 'remotes' package. 
library(remotes)
install_github("curryhank08/GEOquery_with_modifiable_timeout_seconds", force = TRUE)

## Approach 2: From github of GEOquery's author
# Install GEOquery with fix of timeout bug
library(remotes)
install_github("seandavi/GEOquery")

# Load modified GEOquery or GEOquery 2.70 version
library(GEOquery)
# Setting the max timeout_seconds
options(timeout=100000)
# Check the input timeout_seconds
getOption("timeout")

# Download GSE40279 by a fuction getGEO() from GEOquery package.
gse40279 <- getGEO("GSE40279", GSEMatrix = TRUE, AnnotGPL = TRUE)
gse40279_matrix <- gse40279[[1]]

data <- exprs(gse40279_matrix)

# Create age categories
age <- pData(gse40279_matrix)$characteristics_ch1

# Remove "age (y):" and convert to numeric
age <- sub("^\\s*age \\(y\\): ", "", age)
age <- as.numeric(age)

# The ^ character denotes the start of the string,
# \\s* matches any number of leading whitespace characters,
# and "age \\(y\\): " matches the exact string "age (y): ". 

# Assign age values to a new column in pData of gse40279_matrix
pData(gse40279_matrix)$age <- age

# Define age categories based on specific age ranges
age_categories <- cut(age,
                      breaks = c(0, 30, 65, Inf),
                      labels = c("Young", "Middle", "Old"),
                      include.lowest = TRUE)

# Assign age categories to the pData of gse40279_matrix
pData(gse40279_matrix)$age_category <- age_categories

# Run MEAL pipeline on the categorized data
res <- runPipeline(set = gse40279_matrix,
                   variable_names = "age",
                   betas = TRUE,
                   analyses = c("DiffMean", "DiffVar"))

# Extract the result of the DiffMean analysis
result_Meal <- getProbeResults(res, rid = 1, 
                               fNames = c("UCSC_RefGene_Name", "RANGE_START", "CHR", "ID"))
result_Meal$CHR <- as.numeric(result_Meal$CHR)

# Save result_Meal data frame in local for the use in python
# Specify the file path to save the CSV file
file_path <- "./result_Meal.csv"

# Use write.csv() to export the data frame to a CSV file
write.csv(result_Meal, file = file_path, row.names = TRUE)

# Print a message indicating successful export
cat("Data frame successfully exported to:", file_path, "\n")

# Remove rows with missing values
# result_Meal_clean <- na.omit(result_Meal)

#### Part 1-2
## manhattan plot for all
library(qqman)
# function from qqman to plot manhattan 
manhattan(result_Meal, 
          main = "Manhattan Plot for gse40279 (Analysis of DiffMean from regressing age on beta-value for each CpG)",
          cex = 0.6,
          ylim = c(0, 200),
          chr="CHR", 
          bp="RANGE_START", 
          snp= "ID", 
          p="P.Value" )

### Part 2
## Load gse40279 and run analysis by using "limma" packages (from original code limma.R)
library(limma)

data <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex)}

x1 <- pData(gse40279_matrix)$age
design <- model.matrix(~x1)
fit <- lmFit(gse40279_matrix, design)
fit <- eBayes(fit)

result_limma <- topTable(fit, coef="x1", number = Inf, adjust.method = "BH", sort.by = "P")
result_limma$CHR <- as.numeric(result_limma$CHR)

# Remove rows with missing values
# result_limma_clean <- na.omit(result_limma)

## Save result_Meal data frame in local for the use in python
# Specify the file path to save the CSV file
file_path <- "./result_limma.csv"

# Use write.csv() to export the data frame to a CSV file
write.csv(result_limma, file = file_path, row.names = TRUE)

# Print a message indicating successful export
cat("Data frame successfully exported to:", file_path, "\n")


## manhattan plot for result_limma
library(qqman)

# function from qqman to plot manhattan 
manhattan(result_limma, 
          main = "Manhattan Plot of all CpGs after regressing age on beta-value",
          cex = 0.3,
          ylim = c(0, 130),
          chr="CHR", 
          bp="RANGE_START", 
          snp= "ID", 
          p="P.Value",
          genomewideline = FALSE,
          suggestiveline = -log10(1e-05))


#### Below code is not for this code but I remain it for latter usage
### Part 3
## Merge the results of analysis from MEAL and limma, and then create a scatter plot (from original code 6.R)
library(ggplot2)

# Merge results
res_limma_YO_c_p <- data.frame(name = row.names(result_limma_YO), limma_p = result_limma_YO$P.Value)
res_MEAL_YO_c_p <- data.frame(name = row.names(result_Meal_sub_YO_clean), MEAL_p = result_Meal_sub_YO_clean$P.Value)

limma_Meal_YO_p <- merge(res_limma_YO_c_p,
                         res_MEAL_YO_c_p,
                         by.x = "name",
                         by.y = "name")

# Create a scatter plot with x-axis: p-value from limma and y-axis: p-value from MEAL.
plot(limma_Meal_YO_p$limma_p, limma_Meal_YO_p$MEAL_p, 
     xlab = "limma_P.value", ylab = "MEAL_P.value", 
     main = "P.value from Analysis of DiffMean (Young - Old) on MEAL and limma", 
     pch = 20, col = "#8bc34a", cex = 1)

# Mean-Difference Plot from limma
as.numeric(limma_Meal_YO_p$limma_p)
as.numeric(limma_Meal_YO_p$MEAL_p)
x1 <- log(limma_Meal_YO_p$limma_p)
x2 <- log(limma_Meal_YO_p$MEAL_p)
oldpar <- par(mfrow=c(1,2))
plot(x1,x2,xlab="loge(limma_p)",ylab="loge(MEAL_p)",main="Scatter Plot")
mdplot(cbind(x1,x2),bg.pch=1,bg.cex=1,main = "Mean-Difference Plot")
par(oldpar)
