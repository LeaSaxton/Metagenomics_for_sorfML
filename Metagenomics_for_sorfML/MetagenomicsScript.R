#Author: Lea Saxton
#Date: 30/09/2023
#Purpose: Assessing the correlation between the bacterial species from
# a metagenomics dataset and what was found on the surface preparation using FTIR for different food products

#loading the libraries
library(randomForest)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(factoextra)
library(FactoMineR)
library(caret)
library(Boruta)
library(ggplot2)
library(ggfortify)
library(corrplot)

rm(list=ls())

#loading the datasets
path_NGS <- "Mussels_NGS_all.xlsx"
sheets_NGS <- openxlsx::getSheetNames(path_NGS) #accessing the first sheet
data_frame_NGS <- lapply(sheets_NGS[1], openxlsx::read.xlsx, xlsxFile=path_NGS)
dataSet_NGS <- data_frame_NGS[[1]]
path_FTIR <- "FTIR mussels_AUA.xlsx"
sheets_FTIR <- openxlsx::getSheetNames(path_FTIR)
data_frame_FTIR_chilensis <- lapply(sheets_FTIR[1], openxlsx::read.xlsx, xlsxFile=path_FTIR)
data_frame_FTIR_galSansShell <- lapply(sheets_FTIR[2], openxlsx::read.xlsx, xlsxFile=path_FTIR)
data_frame_FTIR_galAvecShell <- lapply(sheets_FTIR[3], openxlsx::read.xlsx, xlsxFile=path_FTIR)
dataSet_chilensis <- data_frame_FTIR_chilensis[[1]]
dataSet_galSansShell <- data_frame_FTIR_galSansShell[[1]]
dataSet_galAvecShell <- data_frame_FTIR_galAvecShell[[1]]

dataPretreatment <- function(dataSet_features,dataSet_bacteria){
        #remove the rows containing NA
        dataSet_features <- na.omit(dataSet_features)
        #removing irrelevant columns from the dataset
        columns_to_delete <- c("codename", "Batch", "batch", "day.of.storage","Day.of.storage", "temperature","Temperature", "TVC.(LOG10).cfu/g","TVC(LOG10).cfu/g", "origin")
        dataSet_features <- dataSet_features[,!(names(dataSet_features)%in%columns_to_delete)]
        # Count occurrences of each unique value in the "NGS.coding" column
        ngs_coding_counts <- table(dataSet_features$NGS.coding)
        # Filter dataSet_bacteria based on matching genus names in dataSet_features
        filtered_dataSet_NGS <- dataSet_bacteria[dataSet_bacteria$genus %in% names(ngs_coding_counts), ]
        # Multiply rows by corresponding counts from ngs_coding_counts
        matched_genus <- filtered_dataSet_NGS$genus
        counts <- unname(ngs_coding_counts[matched_genus])
        # Duplicate rows based on counts
        duplicated_rows <- lapply(1:length(counts), function(i) dataSet_bacteria[dataSet_bacteria$genus == matched_genus[i], ])
        dataSet_NGS_final <- do.call(rbind, lapply(1:length(counts), function(i) duplicated_rows[[i]][rep(1, counts[i]), ]))
        # Create a list to store the modified datasets
        modified_datasets <- list(
                dataSet_features = dataSet_features,
                dataSet_NGS_final = dataSet_NGS_final
        )
        return(modified_datasets)
}

Pretreatment_chilensis <- dataPretreatment(dataSet_chilensis, dataSet_NGS)
dataSet_chilensis_final <- Pretreatment_chilensis$dataSet_features
dataSet_NGS_final_chilensis <- Pretreatment_chilensis$dataSet_NGS_final
Pretreatment_galSansShell <- dataPretreatment(dataSet_galSansShell, dataSet_NGS)
dataSet_galSansShell_final <- Pretreatment_galSansShell$dataSet_features
dataSet_NGS_final_Sans <- Pretreatment_galSansShell$dataSet_NGS_final
Pretreatment_galAvecShell <- dataPretreatment(dataSet_galAvecShell, dataSet_NGS)
dataSet_galAvecShell_final <- Pretreatment_galAvecShell$dataSet_features
dataSet_NGS_final_Avec <- Pretreatment_galAvecShell$dataSet_NGS_final


results_df<- data.frame()

GettingCorrelations <- function(dataSet_features_final, dataSet_NGS_final, name_product){

        # Loop through each target variable
        for (target_col in colnames(dataSet_NGS_final)[colnames(dataSet_NGS_final) != "genus"]) {
                cat("loop for: \n")
                print(target_col)
                # Extract the target variable
                target <- dataSet_NGS_final[[target_col]]
                genus <- dataSet_NGS_final[["genus"]]
                subset <- target
                subset <- as.data.frame(subset)
                colnames(subset) <- target_col
                # Merge features and target for the current target variable
                #keep only in the dataSet the rows with the NGS.coding found in the metadata file
                merged_data <- cbind(dataSet_features_final, subset)
                column <- "NGS.coding"
                merged_data <- merged_data[,!(names(merged_data)%in%column)]
                #Reducing the number of features for each dataset
                NGS.coding <- merged_data[, 1]
                merged_data <- merged_data[,-1]
                # Perform Boruta search
                boruta_output <- Boruta(as.formula(paste(target_col, "~ .")), data = na.omit(merged_data), doTrace = 0)
                boruta_signif <- getSelectedAttributes(boruta_output, withTentative = TRUE)
                cat("Number of selected attributes:\n")
                print(length(boruta_signif))
                directory_path <- paste0(name_product)
                if (length(boruta_signif)>0){
                        dataSet <- merged_data[, (names(merged_data)%in%boruta_signif)]

                        #Data Visualisation:
                        if(all(sapply(dataSet, is.numeric))) {
                                headers <- as.numeric(colnames(dataSet))
                                num_columns <- length(headers)
                                group_size <- 5
                                # Calculate the number of full groups and the remaining columns
                                full_groups <- floor(num_columns / group_size)
                                remaining_columns <- num_columns %% group_size
                                # Create a list to store the ranges
                                ranges <- list()
                                # Populate full groups
                                for (i in 1:full_groups) {
                                        start <- headers[(i - 1) * group_size + 1]
                                        end <- headers[i * group_size]
                                        ranges[[paste(start, "-", end, sep = "")]] <- colnames(dataSet)[headers >= start & headers <= end]
                                }
                                # Handle remaining columns
                                if (remaining_columns > 0) {
                                        start <- headers[num_columns - remaining_columns + 1]
                                        end <- headers[num_columns]
                                        ranges[[paste(start, "-", end, sep = "")]] <- colnames(dataSet)[headers >= start & headers <= end]
                                }
                                # Create a list to store average values for each range
                                average_measures <- matrix(NA, nrow = nrow(dataSet), ncol = length(ranges))
                                # Calculate average for each range and assign column names
                                for (i in 1:length(ranges)) {
                                        range <- ranges[[i]]
                                        range_data <- dataSet[, range, drop = FALSE]
                                        range_average <- rowMeans(range_data, na.rm = TRUE)
                                        average_measures[, i] <- range_average
                                }

                                # Extract the original range names from the 'ranges' list
                                #range_names <- unlist(sapply(ranges, function(range) paste(range, collapse = "-")))
                                range_names <- names(ranges)
                                # Assign original range names directly to the columns of average_measures matrix
                                colnames(average_measures) <- range_names
                                print(average_measures)
                                # Create range labels for the x-axis
                                range_labels <- range_names
                                # Create the directory if it doesn't exist
                                if (!dir.exists(directory_path)) {
                                        dir.create(directory_path)
                                }
                                # Specify the file path for saving the plot
                                file_path <- file.path(directory_path, paste0(target_col, "_barplot.png"))
                                # Calculate the width based on the number of range_labels
                                width_coefficient <- length(range_labels) * 100
                                png_width <- ifelse(width_coefficient < 1000, 1000, width_coefficient)  # Ensure a minimum width of 1000
                                # Create the PNG file with the calculated width
                                png(file_path, width = png_width, height = 1000)
                                par(mar = c(20, 6, 4, 2) + 0.1)  # Set the margins (bottom, left, top, right)
                                barplot(average_measures,
                                        #names.arg = range_labels,
                                        col = "skyblue",
                                        main = "Average Measures at Different Wavelengths",
                                        xlab = "Wavelengths", ylab = "Average Measures", las =2)
                                dev.off()
                        } else {
                                print("Error: Not all elements in the list are numeric.")
                        }

                        dataSet <- cbind(dataSet, subset)
                        # Calculate correlations between features and the target variable
                        correlations <- cor(dataSet, target)
                        correlations <- correlations[correlations[,1]>0.5]
                        # Add the results to the data frame
                        results_df <- as.data.frame(correlations)
                }
                stop()
        }
        write.csv(results_df, file = paste0(name_product, "_correlations.csv"), row.names = FALSE)
}

Correlation_Chilensis <- GettingCorrelations(dataSet_chilensis_final, dataSet_NGS_final_chilensis, "Chilensis")
corr_final_chilensis <- mean(results_df_chilensis$correlation)#Average for the correlations between each bacterial species and features measured with FTIR


Correlation_GalSansShell <- GettingCorrelations(dataSet_galSansShell_final, dataSet_NGS_final_Sans, "GalSansShell")
results_df_galSans <- read.csv("GalSansShell_correlations.csv", header = TRUE)
corr_final_GalSansShell <- mean(results_df_galSans$correlation)



Correlation_GalAvecShell <- GettingCorrelations(dataSet_galAvecShell_final, dataSet_NGS_final_Avec, "GalAvecShell")
results_df_galAvec <- read.csv("GalAvecShell_correlations.csv", header=TRUE)
corr_final_GalAvecShell <- mean(results_df_galAvec$correlation)

cat("Correlation Chilensis Mussels: \n")
print(corr_final_chilensis)
cat("Correlation Galloprovincialis without shell mussels: \n ")
print(corr_final_GalSansShell)
cat("Correlation Galloprovincialis with shell mussels: \n")
print(corr_final_GalAvecShell)

