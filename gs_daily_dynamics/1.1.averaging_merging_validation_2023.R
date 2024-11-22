#########################################################################################
# Data preparation before analysis - slac1 2023
#########################################################################################

# Base path where all measurement folders are located
base_path <- "/local/workdir/hb435/Research/CROPPS/experimental_designs_2023/slac1_HB/measurments"

# List of folders to process
folders <- c("07192023", "07282023", "08062023", "08112023", "08262023", "09062023", "09152023")

# Load necessary libraries
library(dplyr)
library(tidyr)

#########################################################################################
# Logging Function
#########################################################################################

# Define a function to log messages to 'validation_log.txt'
log_message <- function(message, is_warning = FALSE, log_file = "validation_log.txt") {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  if (is_warning) {
    formatted_message <- paste0("[", timestamp, "] WARNING: ", message, "\n")
  } else {
    formatted_message <- paste0("[", timestamp, "] ", message, "\n")
  }
  cat(formatted_message, file = log_file, append = TRUE)
}

#########################################################################################
# Process each folder
#########################################################################################

for (folder in folders) {
  # Set the working directory
  folder_path <- file.path(base_path, folder)
  setwd(folder_path)
  
  # Initialize or clear the validation log for this folder
  log_file <- file.path(folder_path, "validation_log.txt")
  if (file.exists(log_file)) {
    file.remove(log_file)
  }
  
  # Find the file ending with 'Z_1.csv'
  files <- list.files(pattern = "Z_1\\.csv$")
  
  if (length(files) == 0) {
    message("No files ending with 'Z_1.csv' found in folder ", folder)
    next
  }
  
  # Assume the first matching file is the one to process
  file_to_read <- files[1]
  
  # Read the phenotype data file
  licor <- read.csv(file_to_read, header = FALSE, skip = 1)
  
  # Assign row 2 as headers and remove it from the data frame
  names(licor) <- as.character(unlist(licor[1,]))
  licor <- licor[-1, ]
  licor <- licor[-1, ]
  
  # Convert numeric columns from the 10th column onwards to numeric type
  num_cols <- sapply(licor[, 10:ncol(licor)], is.character)
  licor[, 10:ncol(licor)][, num_cols] <- lapply(licor[, 10:ncol(licor)][, num_cols], as.numeric)
  
  # Validation Step 1: Check column names and structure
  log_message("Validation 1: Column names and structure", log_file = log_file)
  expected_cols <- c("Obs#", "Time", "Date", "configName", "configAuthor", "remark", "Observation", 
                     "gsw", "gbw", "gtw", "E_apparent", "VPcham", "VPref", "VPleaf", "VPDleaf", 
                     "H2O_r", "H2O_s", "Qamb", "Tleaf")
  missing_cols <- setdiff(expected_cols, names(licor))
  if (length(missing_cols) > 0) {
    log_message(paste("Missing columns:", paste(missing_cols, collapse = ", ")), is_warning = TRUE, log_file = log_file)
  } else {
    log_message("All expected columns are present.", log_file = log_file)
  }
  
  # Subset the dataframe to keep only the specified columns
  cols_to_keep <- expected_cols
  licor <- licor[, cols_to_keep]
  
  # Validation Step 2: Check 'Observation' column integrity
  log_message("Validation 2: 'Observation' column integrity", log_file = log_file)
  licor$Observation <- as.numeric(licor$Observation)
  missing_obs <- sum(is.na(licor$Observation))
  if (missing_obs > 0) {
    log_message(paste("Missing Observation values:", missing_obs), is_warning = TRUE, log_file = log_file)
  } else {
    log_message("No missing values in 'Observation' column.", log_file = log_file)
  }
  
  duplicate_obs <- licor$Observation[duplicated(licor$Observation)]
  if (length(duplicate_obs) > 0) {
    log_message(paste("Duplicate Observation values:", paste(unique(duplicate_obs), collapse = ", ")), is_warning = TRUE, log_file = log_file)
  } else {
    log_message("No duplicate Observation values.", log_file = log_file)
  }
  
  #########################################################################################
  # 1. Create average for every 4 reads (technical replicates)
  #########################################################################################
  
  # Create a new dataframe for averages
  average_df <- licor[seq(1, nrow(licor), by = 4), 1:9]
  
  # Compute averages for each 4 rows from column 8 onwards
  for (i in 8:ncol(licor)) {
    average_values <- tapply(licor[, i], (seq(nrow(licor)) - 1) %/% 4, mean, na.rm = TRUE)
    average_df[, names(licor)[i]] <- average_values
  }
  
  # Validation Step 3: Check row counts after averaging
  log_message("Validation 3: Row count after averaging", log_file = log_file)
  expected_rows <- ceiling(nrow(licor) / 4)
  actual_rows <- nrow(average_df)
  if (abs(actual_rows - expected_rows) > 1) {
    log_message("Row count after averaging does not match expected value.", is_warning = TRUE, log_file = log_file)
  }
  log_message(paste("Expected rows:", expected_rows, "Actual rows:", actual_rows), log_file = log_file)
  
  # Write to new CSV file
  # Generate output filenames based on the input file
  file_base <- sub("\\.csv$", "", file_to_read)
  means_file <- paste0(file_base, '_means.csv')
  write.csv(average_df, means_file, row.names = FALSE)
  
  #########################################################################################
  # 2. Add the mean table to the headmap data
  #########################################################################################
  
  # Read the headmap file
  headmap_file <- "/local/workdir/hb435/Research/CROPPS/experimental_designs_2023/slac1_HB/headmap.csv"
  headmap <- read.csv(headmap_file)
  
  # Use only the uneven rows from headmap since each plot is a two-row plot with different plot numbers
  headmap <- headmap[seq(1, nrow(headmap), by = 2), ]
  
  # Define the columns to add from headmap
  cols_to_add <- c("rep", "plot_num", "range", "row", "geno")
  
  # Repeat the headmap data frame as many times as necessary
  repeated_headmap <- headmap[rep(seq_len(nrow(headmap)), length.out = nrow(average_df)), cols_to_add]
  
  # Validation Step 4: Check headmap repetition
  log_message("Validation 4: Headmap repetition", log_file = log_file)
  expected_repeats <- nrow(average_df)
  actual_repeats <- nrow(repeated_headmap)
  if (actual_repeats != expected_repeats) {
    log_message("Mismatch in headmap repetition.", is_warning = TRUE, log_file = log_file)
  }
  log_message(paste("Expected headmap rows:", expected_repeats, "Actual rows:", actual_repeats), log_file = log_file)
  
  # Bind the columns from headmap to average_df
  average_df <- cbind(repeated_headmap, average_df)
  
  # Write to new CSV file
  w_headmap_file <- paste0(file_base, '_w_headmap.csv')
  write.csv(average_df, w_headmap_file, row.names = FALSE)
  
  #########################################################################################
  # 3. Add window data
  #########################################################################################
  
  # Initialize the window column
  average_df$window <- 1
  
  # Update the "window" column based on range and row
  for (i in 2:nrow(average_df)) {
    if (average_df$range[i] == 23 && average_df$row[i] == 25) {
      average_df$window[i] <- average_df$window[i - 1] + 1
    } else {
      average_df$window[i] <- average_df$window[i - 1]
    }
  }
  
  # Validation Step 5: Verify window assignment
  log_message("Validation 5: Window assignment", log_file = log_file)
  expected_windows <- length(unique(average_df$window))
  actual_windows <- max(average_df$window)
  if (actual_windows != expected_windows) {
    log_message("Mismatch in window assignment.", is_warning = TRUE, log_file = log_file)
  }
  log_message(paste("Expected windows:", expected_windows, "Actual windows:", actual_windows), log_file = log_file)
  
  # Separate allele and genetic background from geno column
  average_df$allele <- substr(average_df$geno, nchar(average_df$geno) - 1, nchar(average_df$geno))
  average_df$geno_ID <- substr(average_df$geno, 1, 5)
  average_df <- subset(average_df, select = -geno)
  
  # Write the final data to CSV
  final_file <- paste0(file_base, '_w_headmap_final.csv')
  write.csv(average_df, final_file, row.names = FALSE)
  
}
