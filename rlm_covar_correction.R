#Set the working directory as needed
getwd()
setwd("./")

# # Load necessary libraries
library(MASS)
library(MSnSet.utils)
library(ggplot2)

# Read CSV file (protein quant)
data <- read.csv("./input/Quant_file.csv", row.names = 1)
data <-data[,-c(1,2)]
data <- log2(data)

#Seperate the data that has no missing value
df <- data[complete.cases(data),]

#Keep the data that had missing value
exclude<- data[!complete.cases(data),]
exclude <- t(exclude)
exclude <- exclude[order(rownames(exclude)), ]

# Working on the data that has no missing value
data <- t(df)
data <- data[order(rownames(data)), ]

# Read CSV file for metadata
metadata <- read.csv("./input/meta_file.csv", row.names = 1)  # Assuming TMT.channel is in the first column
metadata <- metadata[order(rownames(metadata)), ]
hist(metadata$Cell_line)

# Make data frame for corrected data to save
corrected_data <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = ncol(data)))
colnames(corrected_data) <- colnames(data)
rownames(corrected_data) <- rownames(data)

# Fit multivariate linear model for each protein and obtain residuals
for (protein in colnames(data)) {
  # Extract variables for the model
  y <- data[, protein]
  cell_line <- metadata$Cell_line
  
  # Create a data frame for the model
  model_data <- data.frame(response_variable = y, cell_line=cell_line, maxit=30)
  
  # Fit multivariate linear model
  model <- MASS::rlm(response_variable ~ cell_line, data = model_data)

    # Obtain the residuals and store in corrected_data
  corrected_data[, protein] <- residuals(model) + coef(model)[1]
}

# Merge corrected_data with excluded_data
complete_correcteddata <- cbind(corrected_data, exclude)
complete_correcteddata <- complete_correcteddata[order(rownames(complete_correcteddata)), ]
write.csv(complete_correcteddata,"./output/Covar_corrected_quantfile.csv")

#########################################################################################
#Visualize in PCA
##---before correction-----###
data <- t(data)
data <- as.matrix(data)

msnset <- MSnSet(exprs = data)
annotated_data <- AnnotatedDataFrame(metadata)
phenoData(msnset) <- annotated_data

plot_pca(msnset, phenotype = "Cell_type", legend_title = "Cell_type", show_ellipse = T, components = c(1,2))
plot_pca(msnset, phenotype = "Cell_type", legend_title = "Cell_type", show_ellipse = T, components = c(1,3) )
ggsave("./output/Example_before_correction.png",width = 4, height = 4, dpi = 300 )

##---after correction-----###
corrected_data <- t(corrected_data)
corrected_data <- as.matrix(corrected_data)

msnset <- MSnSet(exprs = corrected_data)
annotated_data <- AnnotatedDataFrame(metadata)
phenoData(msnset) <- annotated_data

plot_pca(msnset, phenotype = "Cell_type", legend_title = "Cell_type", show_ellipse = T, components = c(1,2))
plot_pca(msnset, phenotype = "Cell_type", legend_title = "Cell_type", show_ellipse = T, components = c(1,3) )
ggsave("./output/Example_covar_correction.png",width = 4, height = 4, dpi = 300 )
#########################################################################################
#########################################################################################