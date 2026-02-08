# Combine output tables from 2-compartment simulation batches across multiple
# parameter sets and dates into single merged tables per file type.

library(data.table)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Initialize variables
cohort <- NULL
params_list <- NULL
dates_list <- NULL
additional_pattern <- ""
output_suffix <- "combined"

base_dir <- "OUT_files/Tables/"
path_to_tables <- base_dir

# Parse arguments
i <- 1
while (i <= length(args)) {
  if (args[i] == "--cohort") {
    cohort <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--params") {
    params_list <- strsplit(args[i + 1], ",")[[1]]
    i <- i + 2
  } else if (args[i] == "--dates") {
    dates_list <- strsplit(args[i + 1], ",")[[1]]
    i <- i + 2
  } else if (args[i] == "--pattern") {
    additional_pattern <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--output-suffix") {
    output_suffix <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--path") {
    path_to_tables <- args[i + 1]
    i <- i + 2
  } else {
    cat("Unknown argument:", args[i], "\n")
    quit(status = 1)
  }
}

if (length(params_list) != length(dates_list)) {
  cat("Error: Number of params and dates must match\n")
  quit(status = 1)
}

# Print configuration
cat("  Cohort:", cohort, "\n")
cat("  Parameters:", paste(params_list, collapse = ", "), "\n")
cat("  Dates:", paste(dates_list, collapse = ", "), "\n")
cat("  Additional pattern:", additional_pattern, "\n")
cat("  Output suffix:", output_suffix, "\n")

# Build search patterns for each parameter set
search_patterns <- c()
for (i in seq_along(params_list)) {
  base_pattern <- paste0(cohort, "_", params_list[i])
  date_pattern <- paste0("(sim)?", dates_list[i])
  pattern <- paste0(base_pattern, "_", date_pattern)

  if (additional_pattern != "") {
    pattern <- paste0(pattern, "_", additional_pattern)
  }
  search_patterns <- c(search_patterns, pattern)
}

# Get all files matching the patterns
all_files <- c()
for (i in seq_along(params_list)) {
  base_pattern <- paste0(cohort, "_", params_list[i])
  date_variations <- c(dates_list[i], paste0("sim", dates_list[i]))

  for (date_var in date_variations) {
    pattern <- paste0(base_pattern, "_", date_var)
    files <- list.files(path_to_tables, pattern = pattern, full.names = FALSE)
    all_files <- c(all_files, files)
  }

  cat("Files found for parameter set", params_list[i], ":",
      length(list.files(path_to_tables, pattern = base_pattern, full.names = FALSE)), "\n")
}

all_files <- unique(all_files)

# Extract the file type prefix from a filename by removing the cohort/param/date portion
extract_file_type <- function(filename, cohort, params_list, dates_list) {
  base_name <- sub("\\.txt$", "", filename)

  for (i in seq_along(params_list)) {
    param_pattern <- paste0(cohort, "_", params_list[i])
    date_variations <- c(dates_list[i], paste0("sim", dates_list[i]))

    for (date_var in date_variations) {
      full_pattern <- paste0(param_pattern, "_", date_var)
      if (grepl(full_pattern, base_name)) {
        pattern_pos <- regexpr(full_pattern, base_name)
        if (pattern_pos > 1) {
          file_type <- substr(base_name, 1, pattern_pos - 2)
          return(file_type)
        }
      }
    }
  }
  return(NA)
}

file_types <- unique(sapply(all_files, extract_file_type,
                            cohort = cohort, params_list = params_list, dates_list = dates_list))
file_types <- file_types[!is.na(file_types)]

cat("File types identified:\n")
for (ft in file_types) {
  type_count <- sum(sapply(all_files, function(f) {
    extracted_type <- extract_file_type(f, cohort, params_list, dates_list)
    !is.na(extracted_type) && extracted_type == ft
  }))
  cat("  -", ft, "(", type_count, "files )\n")
}

# Combine tables of the same type across batches
combine_tables_by_type <- function(file_type, all_files, patterns, output_suffix, path_to_tables) {
  cat("Processing file type:", file_type, "\n")

  type_files <- c()
  for (file in all_files) {
    if (startsWith(file, paste0(file_type, "_"))) {
      type_files <- c(type_files, file)
    }
  }

  if (length(type_files) == 0) {
    cat("  No files found for type:", file_type, "\n")
    return()
  }

  cat("  Files to combine:", length(type_files), "\n")
  for (f in type_files) {
    cat("    -", f, "\n")
  }

  combined_table <- NULL

  for (file in type_files) {
    cat("  Reading:", file, "\n")
    tryCatch({
      table_data <- read.table(file.path(path_to_tables, file),
                               sep = '\t',
                               header = TRUE,
                               stringsAsFactors = FALSE)

      if (is.null(combined_table)) {
        combined_table <- table_data
      } else {
        combined_table <- rbind(combined_table, table_data)
      }
    }, error = function(e) {
      cat("    Error reading file:", file, "-", e$message, "\n")
    })
  }

  if (!is.null(combined_table)) {
    param_part <- paste(params_list, collapse = "_")
    date_part <- paste(dates_list, collapse = "_")

    output_filename <- paste0(file_type, "_", cohort, "_", param_part, "_",
                              date_part, "_", output_suffix, ".txt")

    output_path <- file.path(path_to_tables, output_filename)
    write.table(combined_table, output_path, sep = '\t',
                row.names = FALSE, col.names = TRUE, quote = FALSE)

    cat("  Output written to:", output_filename, "\n")
  }
}

# Process each file type
for (file_type in file_types) {
  combine_tables_by_type(file_type, all_files, search_patterns, output_suffix, path_to_tables)
}

