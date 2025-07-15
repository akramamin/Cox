# =============================================================================
# ENHANCED WARD CLUSTERING WITH TRAJECTORY ANALYSIS
# =============================================================================
# Author : Enhanced Analysis Pipeline
# Created: 2024
# Purpose: Comprehensive patient clustering workflow combining optimal baseline 
#          clustering with advanced trajectory analysis and enhanced features
# 
# Key Enhancements:
# - Ward's method with Jaccard distance for robust baseline clustering
# - Enhanced trajectory analysis with new mutation detection
# - Improved feature engineering (VAF slopes, volatility, mutation dynamics)
# - Advanced BMT filtering with Excel serial number support
# - Comprehensive statistical analysis and visualization
# - Cox regression analysis for survival outcomes
# - Optimized package management with graceful fallbacks

# =============================================================================
# SECTION 1: ENHANCED SETUP AND CONFIGURATION
# =============================================================================

# Enhanced Package Management with Graceful Fallbacks -------------------------
required_packages <- c("data.table", "tidyverse", "pheatmap", "RColorBrewer", 
                       "cluster", "NbClust", "readxl", "lubridate", "ggplot2", 
                       "reshape2", "survival", "survminer", "dendextend", "proxy")

# Function to install packages with fallback options
install_package_safe <- function(pkg) {
  tryCatch({
    install.packages(pkg, repos = "https://cran.r-project.org/")
    return(TRUE)
  }, error = function(e1) {
    message(paste("Standard install failed for", pkg, "- trying user library..."))
    tryCatch({
      if (!dir.exists(Sys.getenv("R_LIBS_USER"))) {
        dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)
      }
      install.packages(pkg, lib = Sys.getenv("R_LIBS_USER"), repos = "https://cran.r-project.org/")
      .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
      return(TRUE)
    }, error = function(e2) {
      message(paste("Failed to install", pkg, "- some features may be limited"))
      return(FALSE)
    })
  })
}

# Function to load packages with graceful fallback
load_package_safe <- function(pkg, required = TRUE) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    library(pkg, character.only = TRUE)
    return(TRUE)
  } else {
    if (required) {
      message(paste("Installing required package:", pkg))
      if (install_package_safe(pkg)) {
        library(pkg, character.only = TRUE)
        return(TRUE)
      } else {
        stop(paste("Cannot proceed without", pkg, "package"))
      }
    } else {
      message(paste("Optional package", pkg, "not available - some features may be limited"))
      return(FALSE)
    }
  }
}

# Enhanced package loading with status tracking
message("=== ENHANCED WARD CLUSTERING WITH TRAJECTORY ANALYSIS ===")
message("Checking and loading required packages...")

package_status <- list()

# Essential packages - must be loaded
essential_packages <- c("data.table", "cluster")
for (pkg in essential_packages) {
  package_status[[pkg]] <- load_package_safe(pkg, required = TRUE)
}

# High-priority packages - try to install if missing
high_priority_packages <- c("survival", "ggplot2")
for (pkg in high_priority_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing high-priority package:", pkg))
    install_package_safe(pkg)
  }
  package_status[[pkg]] <- load_package_safe(pkg, required = FALSE)
}

# Optional packages - load if available, don't install if missing
optional_packages <- c("readxl", "tidyverse", "pheatmap", "RColorBrewer", "NbClust", 
                      "lubridate", "ggplot2", "reshape2", "survival", 
                      "survminer", "dendextend", "proxy")

for (pkg in optional_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    tryCatch({
      library(pkg, character.only = TRUE)
      package_status[[pkg]] <- TRUE
      message(paste("✓ Loaded", pkg))
    }, error = function(e) {
      package_status[[pkg]] <- FALSE
      message(paste("⚠ Could not load", pkg))
    })
  } else {
    package_status[[pkg]] <- FALSE
    message(paste("⚠ Optional package", pkg, "not available"))
  }
}

# Report package status
message("\nPackage Status Summary:")
for (pkg in names(package_status)) {
  status <- if (package_status[[pkg]]) "✓ Available" else "✗ Missing"
  message(paste("  ", pkg, ":", status))
}

# =============================================================================
# SECTION 2: ENHANCED UTILITY FUNCTIONS
# =============================================================================

# Function 1: Enhanced MRN and NGS Date Extraction
extract_mrn_ngs_enhanced <- function(excel_file, sheet_name = "Sheet2", output_file = "MRN_NGS_dates.csv") {
  message("Extracting MRN and NGS date columns with enhanced processing...")
  
  if (!package_status[["readxl"]]) {
    stop("readxl package required for Excel processing")
  }
  
  df <- as.data.table(read_excel(excel_file, sheet = sheet_name, .name_repair = 'minimal'))
  
  # Enhanced column detection
  mrn_idx <- which(tolower(trimws(names(df))) %in% c('mrn', 'patient id', 'patient_id', 'id'))
  ngs_idx <- which(grepl('ngs.*date|date.*ngs', tolower(names(df))))
  
  if (length(ngs_idx) == 0) {
    stop('No NGS date columns found')
  }
  
  sel_idx <- c(mrn_idx, ngs_idx)
  sub <- df[, ..sel_idx]
  
  # Enhanced column naming
  if (length(mrn_idx) > 0) {
    setnames(sub, 1, 'MRN')
  } else {
    sub[, MRN := NA]
    setcolorder(sub, c('MRN', setdiff(names(sub), 'MRN')))
  }
  
  # Create meaningful NGS column names
  ngs_names <- paste0('NGS_', seq_len(ncol(sub) - 1), '_date')
  setnames(sub, 2:ncol(sub), ngs_names)
  
  fwrite(sub, output_file)
  message(paste('Enhanced extraction saved', nrow(sub), 'rows and', ncol(sub), 'cols to', output_file))
  
  return(sub)
}

# Function 2: Enhanced Wide-to-Long Conversion with Better Error Handling
convert_wide_to_long_enhanced <- function(excel_file, sheet_name = "Sheet2", output_file = "CMML_Serial_long.csv") {
  message(paste('Enhanced wide-to-long conversion for', sheet_name, 'from', excel_file))
  
  if (!package_status[["readxl"]]) {
    stop("readxl package required for Excel processing")
  }
  
  dt <- as.data.table(read_excel(excel_file, sheet = sheet_name, .name_repair = 'minimal'))
  
  # Enhanced column identification
  mrn_idx <- which(tolower(trimws(names(dt))) %in% c('mrn', 'patient id', 'patient_id', 'id'))
  if (length(mrn_idx) == 0) {
    # Try broader search
    mrn_idx <- which(grepl('mrn|patient|id', tolower(names(dt))))
    if (length(mrn_idx) == 0) stop('No patient identifier column found')
  }
  
  ngs_idx <- which(grepl('ngs.*date|date.*ngs', tolower(names(dt))))
  if (length(ngs_idx) == 0) stop('No "NGS date" columns found')
  
  ngs_idx <- sort(ngs_idx)
  seg_end <- c(ngs_idx[-1] - 1, ncol(dt))
  
  segments <- vector('list', length(ngs_idx))
  
  for (i in seq_along(ngs_idx)) {
    col_range <- ngs_idx[i]:seg_end[i]
    seg <- dt[, c(mrn_idx[1], col_range), with = FALSE]
    setnames(seg, 1:2, c('MRN', 'NGS_date'))
    
    # Enhanced gene column naming with conflict resolution
    if (ncol(seg) > 2) {
      gnam <- names(seg)
      gnam[3:length(gnam)] <- make.unique(gnam[3:length(gnam)], sep = '_v')
      setnames(seg, gnam)
    }
    
    # Enhanced filtering - remove empty and invalid dates
    seg <- seg[!is.na(NGS_date) & trimws(as.character(NGS_date)) != '']
    segments[[i]] <- seg
    message(paste('Segment', i, ':', nrow(seg), 'rows,', ncol(seg), 'cols'))
  }
  
  # Enhanced merging with better column handling
  long_dt <- rbindlist(segments, use.names = TRUE, fill = TRUE)
  setcolorder(long_dt, c('MRN', 'NGS_date', setdiff(names(long_dt), c('MRN', 'NGS_date'))))
  
  fwrite(long_dt, output_file)
  message(paste('Enhanced long-form data saved to', output_file, '(', nrow(long_dt), 'rows)'))
  
  return(long_dt)
}

# Function 3: Enhanced Normalization Functions
normalize_l2 <- function(x) {
  nrm <- sqrt(sum(x^2))
  if (nrm == 0) return(x)
  x / nrm
}

# Function 4: Enhanced Date Parsing with Excel Serial Support
parse_date_enhanced <- function(date_vec) {
  vec <- as.character(date_vec)
  parsed_dates <- rep(NA_real_, length(vec))
  
  # Handle Excel serial numbers
  is_serial <- grepl("^[0-9]+(\\.[0-9]+)?$", vec) & !is.na(vec) & vec != ""
  if (any(is_serial)) {
    serial_nums <- as.numeric(vec[is_serial])
    # Filter out obviously invalid serials (too small or too large)
    valid_serials <- serial_nums > 1 & serial_nums < 100000
    if (any(valid_serials)) {
      parsed_dates[is_serial][valid_serials] <- as.numeric(as.Date(serial_nums[valid_serials], origin = "1899-12-30"))
    }
  }
  
  # Handle regular date formats
  if (package_status[["lubridate"]] && any(!is_serial & !is.na(vec) & vec != "")) {
    non_serial_idx <- which(!is_serial & !is.na(vec) & vec != "")
    date_attempts <- suppressWarnings(lubridate::parse_date_time(vec[non_serial_idx], 
                                                               orders = c("ymd", "mdy", "dmy", "ymd HMS", "mdy HMS", "dmy HMS")))
    parsed_dates[non_serial_idx] <- as.numeric(as.Date(date_attempts))
  }
  
  final_dates <- as.Date(parsed_dates, origin = "1970-01-01")
  
  # Handle Excel serial 0 (1899-12-30) as missing
  final_dates[final_dates == as.Date("1899-12-30")] <- NA
  
  return(final_dates)
}

# =============================================================================
# SECTION 3: ENHANCED CONFIGURATION AND DATA LOADING
# =============================================================================

# Enhanced Configuration Parameters
clinical_excel <- "CMML_Project_2.xlsx"   # Primary workbook
clinical_sheet <- "Sheet1"                # Clinical data sheet
serial_sheet   <- "Sheet2"                # Serial NGS data sheet
mutation_file  <- "CMML_Sheet2.csv"       # Baseline mutations (CSV backup)
output_prefix  <- "enhanced_ward"         # Output file prefix
set.seed(123)                             # Reproducibility

message("=== ENHANCED DATA LOADING AND PREPARATION ===")

# Enhanced Clinical Data Loading -----------------------------------------------
clinical_csv <- "Sheet1.csv"
if (file.exists(clinical_csv)) {
  message("Loading clinical data from CSV...")
  df_clinical <- fread(clinical_csv, stringsAsFactors = FALSE)
  
  # Enhanced date parsing for clinical data
  date_columns <- c("last follow up", "Date", "BMT dates")
  for (col in date_columns) {
    if (col %in% names(df_clinical)) {
      df_clinical[[col]] <- parse_date_enhanced(df_clinical[[col]])
    }
  }
  
} else if (package_status[["readxl"]] && file.exists(clinical_excel)) {
  message("Loading clinical data from Excel...")
  df_clinical <- read_excel(clinical_excel, sheet = clinical_sheet, .name_repair = "minimal")
  df_clinical <- as.data.table(df_clinical)
  
  # Enhanced date parsing
  date_columns <- c("last follow up", "Date", "BMT dates")
  for (col in date_columns) {
    if (col %in% names(df_clinical)) {
      df_clinical[[col]] <- parse_date_enhanced(df_clinical[[col]])
    }
  }
} else {
  stop("No clinical data source available. Need either Sheet1.csv or Excel file with readxl package.")
}

message(paste("Clinical data loaded:", nrow(df_clinical), "patients x", ncol(df_clinical), "variables"))

# Enhanced Diagnosis Column Detection ------------------------------------------
diag_candidates <- grep("diagnosis|subtype|disease", names(df_clinical), ignore.case = TRUE, value = TRUE)
diag_col <- if (length(diag_candidates) > 0) diag_candidates[1] else NULL
if (!is.null(diag_col)) {
  message(paste("Diagnosis column identified:", diag_col))
}

# Enhanced Survival Data Preparation -------------------------------------------
if (all(c("OS months", "Survival status") %in% names(df_clinical))) {
  message("Processing overall survival data...")
  df_clinical$`OS months` <- suppressWarnings(as.numeric(df_clinical$`OS months`))
  
  # Enhanced survival status parsing
  stat_col <- df_clinical$`Survival status`
  event_parsed <- rep(NA_integer_, length(stat_col))
  
  if (is.numeric(stat_col)) {
    event_parsed <- ifelse(is.na(stat_col), NA_integer_, ifelse(stat_col > 0, 1L, 0L))
  } else {
    stat_chr <- tolower(trimws(as.character(stat_col)))
    # Death indicators
    event_parsed[stat_chr %in% c("dead", "deceased", "expired", "death", "died", "1", "yes", "y", "true")] <- 1L
    # Alive indicators
    event_parsed[stat_chr %in% c("alive", "living", "0", "no", "n", "false", "censored", "unknown", "na", "")] <- 0L
  }
  df_clinical$event <- event_parsed
  
  n_events <- sum(df_clinical$event == 1, na.rm = TRUE)
  n_censored <- sum(df_clinical$event == 0, na.rm = TRUE)
  message(paste("OS events:", n_events, "deaths,", n_censored, "censored"))
}

# Enhanced Progression-Free Survival Processing --------------------------------
if ("Progression Free Survival Months" %in% names(df_clinical)) {
  message("Processing progression-free survival data...")
  df_clinical$`Progression Free Survival Months` <- suppressWarnings(as.numeric(df_clinical$`Progression Free Survival Months`))
  
  # Enhanced PFS event detection
  pfs_event_cols <- c("Leukemia progression", "Progression of MDS", "Disease progression", "Progression")
  pfs_col_found <- NULL
  
  for (col in pfs_event_cols) {
    if (col %in% names(df_clinical)) {
      pfs_col_found <- col
      break
    }
  }
  
  if (!is.null(pfs_col_found)) {
    pfs_status <- tolower(trimws(as.character(df_clinical[[pfs_col_found]])))
    df_clinical$pfs_event <- ifelse(pfs_status %in% c("yes", "y", "1", "true", "progression"), 1, 0)
    n_pfs_events <- sum(df_clinical$pfs_event == 1, na.rm = TRUE)
    message(paste("PFS events:", n_pfs_events, "progressions detected"))
  }
}

# Enhanced BMT Date Processing with Advanced Filtering -------------------------
if ("BMT dates" %in% names(df_clinical)) {
  message("Processing BMT dates with enhanced filtering...")
  df_clinical$bmt_date_parsed <- parse_date_enhanced(df_clinical$`BMT dates`)
  
  n_bmt_patients <- sum(!is.na(df_clinical$bmt_date_parsed))
  message(paste("BMT dates parsed:", n_bmt_patients, "patients with valid BMT dates"))
  
  if (n_bmt_patients > 0) {
    # Quality check on BMT dates
    bmt_range <- range(df_clinical$bmt_date_parsed, na.rm = TRUE)
    message(paste("BMT date range:", as.character(bmt_range[1]), "to", as.character(bmt_range[2])))
  }
} else {
  message("No BMT dates column found - all NGS data will be included in trajectory analysis")
  df_clinical$bmt_date_parsed <- NA
}

# Enhanced Mutation Data Loading -----------------------------------------------
message("Loading mutation data with enhanced processing...")

if (file.exists(mutation_file)) {
  message("Loading from existing mutation CSV...")
  df_genes_raw <- fread(mutation_file, header = TRUE, stringsAsFactors = FALSE)
  df_genes <- df_genes_raw
  
} else {
  message("Extracting mutation data from source files...")
  
  # Try Sheet2.csv first
  sheet2_csv <- "Sheet2.csv"
  if (file.exists(sheet2_csv)) {
    message("Extracting baseline mutations from Sheet2.csv...")
    sheet2_data <- fread(sheet2_csv, stringsAsFactors = FALSE)
    
    # Enhanced NGS_1 column detection
    ngs1_cols <- grep("^NGS.*1|^NGS_1", names(sheet2_data), value = TRUE)
    if (length(ngs1_cols) > 0) {
      # Include patient identifier if available
      id_cols <- grep("^MRN|^Patient|^ID", names(sheet2_data), value = TRUE)
      if (length(id_cols) > 0) {
        df_genes <- sheet2_data[, c(id_cols[1], ngs1_cols), with = FALSE]
      } else {
        df_genes <- sheet2_data[, ngs1_cols, with = FALSE]
      }
      
      fwrite(df_genes, mutation_file)
      message(paste("Extracted and saved", ncol(df_genes), "baseline mutation columns"))
    } else {
      stop("No NGS_1 columns found in Sheet2 data")
    }
    
  } else if (package_status[["readxl"]] && file.exists(clinical_excel)) {
    message("Extracting from Excel workbook...")
    df_genes_from_excel <- read_excel(clinical_excel, sheet = serial_sheet, .name_repair = "minimal")
    ngs1_cols <- grep("NGS.*1|NGS_1", names(df_genes_from_excel), value = TRUE)
    
    if (length(ngs1_cols) > 0) {
      df_genes <- as.data.table(df_genes_from_excel[, ngs1_cols])
      fwrite(df_genes, mutation_file)
      message(paste("Extracted", ncol(df_genes), "baseline mutation columns from Excel"))
    } else {
      stop("Cannot identify NGS_1 columns in Excel data")
    }
  } else {
    stop("No mutation data source available")
  }
}

# Enhanced data alignment validation
if (nrow(df_genes) != nrow(df_clinical)) {
  warning(paste("Row count mismatch: clinical (", nrow(df_clinical), ") vs mutation (", nrow(df_genes), ")"))
  min_rows <- min(nrow(df_genes), nrow(df_clinical))
  df_genes <- df_genes[1:min_rows, ]
  df_clinical <- df_clinical[1:min_rows, ]
  message(paste("Aligned datasets to", min_rows, "rows"))
}

message(paste("Final dataset:", nrow(df_clinical), "patients ready for analysis"))

# =============================================================================
# SECTION 4: ENHANCED BASELINE CLUSTERING WITH WARD'S METHOD
# =============================================================================

message("=== ENHANCED BASELINE CLUSTERING ANALYSIS ===")

# Enhanced Binary Mutation Matrix Construction ---------------------------------
message("Constructing enhanced binary mutation matrix...")

# Remove duplicate columns and clean gene data
cols_to_keep <- colnames(df_genes)[!duplicated(colnames(df_genes))]
df_genes <- df_genes[, ..cols_to_keep]

# Enhanced gene extraction for NGS_1 data
cols_ngs1 <- grep("^NGS.*1|^NGS_1", colnames(df_genes), value = TRUE)
df_genes_base <- df_genes[, ..cols_ngs1]

# Enhanced gene name extraction with better parsing
all_cols <- colnames(df_genes_base)
unique_genes <- sort(unique(sapply(all_cols, function(col) {
  # Remove NGS_1 prefix and extract gene name
  clean_col <- gsub("^NGS.*1[._\\s]*", "", col)
  gene <- strsplit(clean_col, "[._\\s]")[[1]][1]
  
  # Handle NA and empty gene names
  if (is.na(gene) || nchar(trimws(gene)) == 0) {
    return(NA)
  }
  
  gene <- trimws(gene)
  
  # Filter out technical annotations
  if (!startsWith(gene, "p.") && 
      !startsWith(gene, "c.") && 
      !grepl("^(date|pos|vaf|nm|type)$", gene, ignore.case = TRUE)) {
    return(gene)
  } else {
    return(NA)
  }
})))
unique_genes <- unique_genes[!is.na(unique_genes) & unique_genes != ""]

message(paste("Identified", length(unique_genes), "unique genes for analysis"))

# Enhanced mutation matrix construction with better aggregation
mutation_matrix <- sapply(unique_genes, function(gene) {
  # Find all columns for this gene
  gene_pattern <- paste0("^NGS.*1[._\\s]*", gene, "([._\\s]|$)")
  gene_cols <- grep(gene_pattern, colnames(df_genes_base), value = TRUE)
  
  # Exclude technical columns (VAF, position, etc.)
  mut_cols <- gene_cols[!grepl("(pos|vaf|nm|type|date)", gene_cols, ignore.case = TRUE)]
  
  if (length(mut_cols) == 0) return(rep(0, nrow(df_genes_base)))
  
  # Enhanced numeric conversion and aggregation
  gene_dat <- as.data.frame(df_genes_base[, ..mut_cols])
  gene_dat[] <- lapply(gene_dat, function(x) {
    # Convert to numeric, handling various formats
    num_x <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(x))))
    ifelse(is.na(num_x), 0, num_x)
  })
  
  # Return maximum value across all columns for this gene (mutation present if any column > 0)
  apply(gene_dat, 1, max, na.rm = TRUE)
})

mutation_matrix <- as.data.frame(mutation_matrix)
rownames(mutation_matrix) <- 1:nrow(mutation_matrix)

# Enhanced validation and filtering
if (nrow(mutation_matrix) == 0 || ncol(mutation_matrix) == 0) {
  stop("ERROR: Invalid mutation matrix dimensions")
}

# Enhanced variance filtering with informative messages
variances <- apply(mutation_matrix, 2, var, na.rm = TRUE)
keep_cols <- !(variances == 0 | is.na(variances))
n_removed <- sum(!keep_cols)

if (sum(keep_cols) == 0) {
  stop("ERROR: No genes with variance detected")
}

mutation_matrix <- mutation_matrix[, keep_cols, drop = FALSE]
message(paste("Mutation matrix:", nrow(mutation_matrix), "patients x", ncol(mutation_matrix), "genes"))
if (n_removed > 0) {
  message(paste("Removed", n_removed, "genes with zero variance"))
}

# Enhanced Binary Matrix for Jaccard Distance ----------------------------------
# Convert to strict binary matrix (0/1) for Jaccard distance calculation
binary_matrix <- as.matrix(mutation_matrix)
storage.mode(binary_matrix) <- "numeric"
binary_matrix <- ifelse(is.na(binary_matrix) | binary_matrix == 0, 0, 1)

# Track patients with mutations
row_sums <- rowSums(binary_matrix)
kept_idx <- which(row_sums > 0)
excluded_patients <- setdiff(1:nrow(binary_matrix), kept_idx)

if (length(excluded_patients) > 0) {
  write.csv(excluded_patients, "enhanced_excluded_patients_no_mutations.csv")
  message(paste("Excluded", length(excluded_patients), "patients with no mutations"))
}

binary_matrix_filtered <- binary_matrix[kept_idx, , drop = FALSE]
message(paste("Clustering matrix:", nrow(binary_matrix_filtered), "patients with mutations"))

# Enhanced Jaccard Distance Calculation ----------------------------------------
message("Computing enhanced Jaccard distance matrix...")

if (package_status[["proxy"]]) {
  jaccard_dist <- dist(binary_matrix_filtered, method = "Jaccard")
  message("Using proxy package for Jaccard distance")
} else {
  message("Computing Jaccard distance manually...")
  n_samples <- nrow(binary_matrix_filtered)
  jaccard_matrix <- matrix(0, n_samples, n_samples)
  
  # Enhanced manual Jaccard calculation with progress indication
  for (i in 1:(n_samples-1)) {
    if (i %% 10 == 0) message(paste("Processing patient", i, "of", n_samples))
    for (j in (i+1):n_samples) {
      x <- binary_matrix_filtered[i, ]
      y <- binary_matrix_filtered[j, ]
      intersection <- sum(x & y)
      union <- sum(x | y)
      jaccard_dist_val <- if (union == 0) 0 else 1 - intersection/union
      jaccard_matrix[i, j] <- jaccard_dist_val
      jaccard_matrix[j, i] <- jaccard_dist_val
    }
  }
  jaccard_dist <- as.dist(jaccard_matrix)
}

# Enhanced validation of distance matrix
if (any(is.na(jaccard_dist))) {
  warning("NA values detected in distance matrix")
  jaccard_dist[is.na(jaccard_dist)] <- 1  # Set NA distances to maximum
}

# Enhanced Ward's Hierarchical Clustering --------------------------------------
message("Performing enhanced Ward's clustering...")

if (nrow(binary_matrix_filtered) < 2) {
  stop("ERROR: Insufficient samples for clustering")
}

# Use Ward.D2 method (appropriate for distance matrices)
hc_ward <- hclust(jaccard_dist, method = "ward.D2")

# Enhanced Dendrogram Visualization
png("enhanced_ward_dendrogram.png", width = 1200, height = 800)
plot(hc_ward, main = "Enhanced Ward's Clustering (Jaccard Distance)", 
     xlab = "Patients", ylab = "Distance", cex = 0.7)
dev.off()

# Enhanced Optimal k Determination with Multiple Methods -----------------------
message("Determining optimal number of clusters...")

# Method 1: Silhouette analysis
max_k <- min(10, nrow(binary_matrix_filtered) - 1)
sil_scores <- sapply(2:max_k, function(k) {
  cluster_assignments <- cutree(hc_ward, k = k)
  sil <- silhouette(cluster_assignments, jaccard_dist)
  mean(sil[, 3])
})

opt_k_silhouette <- which.max(sil_scores) + 1
message(paste("Optimal k by silhouette:", opt_k_silhouette))

# Method 2: Elbow method using within-cluster sum of squares
wss_scores <- sapply(1:max_k, function(k) {
  cluster_assignments <- cutree(hc_ward, k = k)
  total_wss <- 0
  for (i in 1:k) {
    points <- which(cluster_assignments == i)
    if (length(points) > 1) {
      # All pairwise distances within the cluster
      dist_mat <- as.matrix(jaccard_dist)[points, points]
      total_wss <- total_wss + sum(dist_mat[upper.tri(dist_mat)])
    }
  }
  return(total_wss)
})

# Enhanced visualization of cluster selection criteria
png("enhanced_cluster_selection.png", width = 1200, height = 800)
par(mfrow = c(1, 2))

# Silhouette plot
plot(2:max_k, sil_scores, type = "b", pch = 19, col = "blue",
     main = "Silhouette Analysis", xlab = "Number of Clusters (k)", 
     ylab = "Average Silhouette Width")
abline(v = opt_k_silhouette, col = "red", lty = 2)
text(opt_k_silhouette, max(sil_scores), paste("k =", opt_k_silhouette), 
     pos = 4, col = "red", font = 2)

# Elbow plot
plot(1:max_k, wss_scores, type = "b", pch = 19, col = "darkgreen",
     main = "Elbow Method", xlab = "Number of Clusters (k)", 
     ylab = "Within-cluster Sum of Squares")

dev.off()

# Use silhouette-based optimal k (or override here if needed)
opt_k_final <- opt_k_silhouette
final_clusters <- cutree(hc_ward, k = opt_k_final)

message(paste("Final clustering: k =", opt_k_final))
cluster_table <- table(final_clusters)
message("Cluster sizes:")
print(cluster_table)

# Enhanced Silhouette Analysis for Final Clustering ----------------------------
sil_final <- silhouette(final_clusters, jaccard_dist)
avg_sil <- mean(sil_final[, 3])

png("enhanced_silhouette_final.png", width = 1200, height = 800)
plot(sil_final, col = rainbow(opt_k_final), border = NA,
     main = paste("Enhanced Silhouette Analysis (k =", opt_k_final, ")"))
mtext(paste("Average silhouette width:", round(avg_sil, 3)), 
      side = 3, line = 0.5, cex = 1.1)
dev.off()

message(paste("Average silhouette width:", round(avg_sil, 3)))

# Map clusters back to all patients (including those with no mutations)
baseline_cluster <- rep(NA_integer_, nrow(mutation_matrix))
baseline_cluster[kept_idx] <- final_clusters

# Add cluster information to clinical data
df_clinical$baseline_cluster <- baseline_cluster

# Enhanced output files
fwrite(df_clinical, paste0("CMML_clinical_", output_prefix, "_clustered.csv"))
message("Enhanced baseline clustering results saved")

# =============================================================================
# SECTION 5: ENHANCED TRAJECTORY ANALYSIS
# =============================================================================

message("=== ENHANCED TRAJECTORY ANALYSIS ===")

# Enhanced Serial NGS Data Loading ---------------------------------------------
serial_long_csv <- "CMML_Serial_long.csv"

if (file.exists(serial_long_csv)) {
  message("Loading serial NGS data from long CSV...")
  df_serial_raw <- fread(serial_long_csv, stringsAsFactors = FALSE)
  
} else if (file.exists("Sheet2.csv")) {
  message("Converting Sheet2.csv to long format for trajectory analysis...")
  
  # Load Sheet2.csv directly and convert to long format
  sheet2_data <- fread("Sheet2.csv", stringsAsFactors = FALSE)
  
  # Enhanced column identification
  mrn_idx <- which(grepl("mrn|patient|id", tolower(names(sheet2_data))))
  if (length(mrn_idx) == 0) mrn_idx <- 1  # Use first column as fallback
  
  # Find NGS date columns
  ngs_date_idx <- which(grepl('ngs.*date|date.*ngs', tolower(names(sheet2_data))))
  if (length(ngs_date_idx) == 0) {
    message("No NGS date columns found - skipping trajectory analysis")
    df_serial_raw <- NULL
  } else {
    message(paste("Found", length(ngs_date_idx), "NGS date columns for trajectory analysis"))
    
    # Convert to long format manually
    segments <- vector('list', length(ngs_date_idx))
    
    for (i in seq_along(ngs_date_idx)) {
      # Define column range for this NGS timepoint
      start_col <- ngs_date_idx[i]
      end_col <- if (i < length(ngs_date_idx)) ngs_date_idx[i+1] - 1 else ncol(sheet2_data)
      
      col_range <- start_col:end_col
      seg <- sheet2_data[, c(mrn_idx[1], col_range), with = FALSE]
      
      # Standardize column names
      old_names <- names(seg)
      new_names <- c("MRN", "NGS_date", paste0("Gene_", seq_len(length(old_names)-2)))
      if (length(old_names) > 2) {
        # Try to extract actual gene names from column headers
        gene_cols <- old_names[3:length(old_names)]
        clean_gene_names <- sapply(gene_cols, function(x) {
          # Remove NGS prefix and extract gene name
          clean_name <- gsub("^NGS.*[0-9][._\\s]*", "", x)
          clean_name <- strsplit(clean_name, "[._\\s]")[[1]][1]
          if (is.na(clean_name) || nchar(trimws(clean_name)) == 0) {
            return(x)  # Return original if can't clean
          }
          return(trimws(clean_name))
        })
        new_names[3:length(new_names)] <- make.unique(clean_gene_names, sep = "_v")
      }
      setnames(seg, old_names, new_names)
      
      # Remove empty dates
      seg <- seg[!is.na(NGS_date) & trimws(as.character(NGS_date)) != '']
      
      segments[[i]] <- seg
      message(paste('NGS timepoint', i, ':', nrow(seg), 'patients with data'))
    }
    
    # Combine all timepoints
    df_serial_raw <- rbindlist(segments, use.names = TRUE, fill = TRUE)
    setcolorder(df_serial_raw, c('MRN', 'NGS_date', setdiff(names(df_serial_raw), c('MRN', 'NGS_date'))))
    
    # Save the long format data
    fwrite(df_serial_raw, serial_long_csv)
    message(paste('Long-form trajectory data saved to', serial_long_csv, '(', nrow(df_serial_raw), 'rows)'))
  }
  
} else if (package_status[["readxl"]] && file.exists(clinical_excel)) {
  message("Converting Excel Sheet2 to long format...")
  df_serial_raw <- convert_wide_to_long_enhanced(clinical_excel, serial_sheet, serial_long_csv)
  
} else {
  message("No serial NGS data found - skipping trajectory analysis")
  df_serial_raw <- NULL
}

if (!is.null(df_serial_raw) && nrow(df_serial_raw) > 0) {
  
  # Enhanced column identification and cleaning
  orig_names <- names(df_serial_raw)
  clean_names <- gsub("[^A-Za-z0-9_]", "_", orig_names)
  clean_names <- make.unique(clean_names, sep = "_")
  colnames(df_serial_raw) <- clean_names
  
  # Enhanced column detection
  if ("MRN" %in% names(df_serial_raw)) {
    id_col <- "MRN"
  } else {
    id_candidates <- names(df_serial_raw)[grepl("Patient|Subject|ID|MRN", names(df_serial_raw), ignore.case = TRUE)]
    id_col <- if (length(id_candidates) > 0) id_candidates[1] else NULL
  }
  
  date_candidates <- names(df_serial_raw)[grepl("date|Date", names(df_serial_raw), ignore.case = TRUE)]
  vaf_candidates <- names(df_serial_raw)[grepl("VAF|vaf", names(df_serial_raw), ignore.case = TRUE)]
  
  if (is.null(id_col)) {
    df_serial_raw$PatientID <- paste0("Patient_", 1:nrow(df_serial_raw))
    id_col <- "PatientID"
  }
  
  if (length(date_candidates) == 0) {
    message("No date columns found - skipping trajectory analysis")
    df_serial_raw <- NULL
  } else {
    date_col <- date_candidates[1]
  }
  
  if (length(vaf_candidates) == 0) {
    message("No VAF columns found - skipping trajectory analysis") 
    df_serial_raw <- NULL
  }
}

if (!is.null(df_serial_raw) && nrow(df_serial_raw) > 0) {
  
  message(paste("Processing trajectory data with columns:", id_col, "|", date_col, "|", length(vaf_candidates), "VAF columns"))
  
  # Enhanced data melting to long format
  dt_serial <- as.data.table(df_serial_raw[, c(id_col, date_col, vaf_candidates), with = FALSE])
  setnames(dt_serial, c(id_col, date_col), c("Patient", "Date"))
  
  # Convert all VAF columns to numeric before melting
  vaf_cols <- vaf_candidates
  for (col in vaf_cols) {
    if (col %in% names(dt_serial)) {
      dt_serial[[col]] <- suppressWarnings(as.numeric(as.character(dt_serial[[col]])))
    }
  }
  
  long_serial <- melt(dt_serial,
                      id.vars = c("Patient", "Date"),
                      variable.name = "Gene",
                      value.name = "VAF",
                      na.rm = TRUE)
  
  long_serial <- as.data.table(long_serial)
  
  # Enhanced date and VAF parsing
  long_serial[, Date := parse_date_enhanced(Date)]
  long_serial[, VAF := as.numeric(VAF)]
  
  # Remove invalid data
  long_serial <- long_serial[!is.na(VAF) & !is.na(Date) & VAF >= 0]
  
  # Enhanced BMT filtering
  if (any(!is.na(df_clinical$bmt_date_parsed))) {
    message("Applying enhanced BMT filtering...")
    
    # Create patient-BMT mapping
    patient_id_candidates <- names(df_clinical)[grepl("MRN|Patient|Subject|ID", names(df_clinical), ignore.case = TRUE)]
    
    if (length(patient_id_candidates) > 0) {
      patient_id_col <- patient_id_candidates[1]
      bmt_subset <- df_clinical[!is.na(bmt_date_parsed)]
      bmt_mapping <- data.table(Patient = as.character(bmt_subset[[patient_id_col]]), 
                                BMT_Date = bmt_subset$bmt_date_parsed)
    } else {
      bmt_mapping <- data.table(Patient = as.character(1:nrow(df_clinical)), 
                                BMT_Date = df_clinical$bmt_date_parsed)
      bmt_mapping <- bmt_mapping[!is.na(BMT_Date)]
    }
    
    if (nrow(bmt_mapping) > 0) {
      # Convert patient IDs to character for consistent matching
      long_serial[, Patient := as.character(Patient)]
      
      # Merge and filter
      long_serial_with_bmt <- merge(long_serial, bmt_mapping, by = "Patient", all.x = TRUE)
      n_before_bmt <- nrow(long_serial)
      long_serial <- long_serial_with_bmt[is.na(BMT_Date) | Date < BMT_Date]
      long_serial[, BMT_Date := NULL]
      
      n_after_bmt <- nrow(long_serial)
      message(paste("BMT filtering: removed", n_before_bmt - n_after_bmt, "post-BMT data points"))
    }
  }
  
  # Enhanced patient trajectory identification
  pat_time_counts <- long_serial[!is.na(Date), uniqueN(Date), by = Patient]
  multi_tp_pats <- pat_time_counts[V1 >= 2]$Patient
  
  message(paste("Patients with ≥2 timepoints:", length(multi_tp_pats)))
  
  if (length(multi_tp_pats) >= 2) {
    
    # Filter to multi-timepoint patients
    long_serial <- long_serial[Patient %in% multi_tp_pats]
    
    # Enhanced trajectory feature engineering
    message("Computing enhanced trajectory features...")
    
    # Basic trajectory features (mean VAF per patient-date)
    pat_date_vaf <- long_serial[, .(mean_vaf = mean(VAF, na.rm = TRUE)), by = .(Patient, Date)]
    
    # Enhanced new mutation detection and slope computation
    compute_enhanced_new_mutation_features <- function(patient_data) {
      setorder(patient_data, Date)
      
      if (nrow(patient_data) < 2) return(list(
        n_new_mutations = 0.0, mean_new_mut_slope = 0.0, max_new_mut_slope = 0.0,
        min_new_mut_slope = 0.0, sd_new_mut_slope = 0.0, new_mut_volatility = 0.0
      ))
      
      dates <- unique(patient_data$Date)
      if (length(dates) < 2) return(list(
        n_new_mutations = 0.0, mean_new_mut_slope = 0.0, max_new_mut_slope = 0.0,
        min_new_mut_slope = 0.0, sd_new_mut_slope = 0.0, new_mut_volatility = 0.0
      ))
      
      # Identify baseline mutations (first timepoint)
      baseline_genes <- patient_data[Date == dates[1] & VAF > 0]$Gene
      
      new_mutation_slopes <- c()
      new_mutation_vafs <- c()
      
      for (i in 2:length(dates)) {
        current_date <- dates[i]
        current_genes <- patient_data[Date == current_date & VAF > 0]$Gene
        new_genes <- setdiff(current_genes, baseline_genes)
        
        if (length(new_genes) > 0) {
          for (gene in new_genes) {
            gene_data <- patient_data[Gene == gene & Date >= current_date]
            if (nrow(gene_data) >= 2) {
              time_numeric <- as.numeric(gene_data$Date - min(gene_data$Date))
              if (length(unique(time_numeric)) > 1 && var(time_numeric) > 0) {
                slope <- lm(VAF ~ time_numeric, data = gene_data)$coefficients[2]
                if (!is.na(slope)) {
                  new_mutation_slopes <- c(new_mutation_slopes, slope)
                  new_mutation_vafs <- c(new_mutation_vafs, gene_data$VAF)
                }
              }
            }
          }
        }
      }
      
      if (length(new_mutation_slopes) > 0) {
        return(list(
          n_new_mutations = as.numeric(length(new_mutation_slopes)),
          mean_new_mut_slope = as.numeric(mean(new_mutation_slopes, na.rm = TRUE)),
          max_new_mut_slope = as.numeric(max(new_mutation_slopes, na.rm = TRUE)),
          min_new_mut_slope = as.numeric(min(new_mutation_slopes, na.rm = TRUE)),
          sd_new_mut_slope = as.numeric(sd(new_mutation_slopes, na.rm = TRUE)),
          new_mut_volatility = as.numeric(mean(abs(diff(new_mutation_vafs)), na.rm = TRUE))
        ))
      } else {
        return(list(
          n_new_mutations = 0.0, mean_new_mut_slope = 0.0, max_new_mut_slope = 0.0,
          min_new_mut_slope = 0.0, sd_new_mut_slope = 0.0, new_mut_volatility = 0.0
        ))
      }
    }
    
    # Compute enhanced new mutation features
    new_mutation_features <- long_serial[, compute_enhanced_new_mutation_features(.SD), by = Patient]
    
    # Enhanced basic trajectory features
    traj_features <- pat_date_vaf[, {
      setorder(.SD, Date)
      vafs <- mean_vaf
      if (length(vafs) < 2) {
        NULL
      } else {
        time_idx <- seq_along(vafs)
        time_days <- as.numeric(Date - min(Date))
        
        # Enhanced feature computation
        trend_temporal <- if (var(time_days) > 0) lm(vafs ~ time_days)$coefficients[2] else 0
        trend_ordinal <- lm(vafs ~ time_idx)$coefficients[2]
        
        list(
          n_timepoints = length(vafs),
          mean_vaf = mean(vafs, na.rm = TRUE),
          median_vaf = median(vafs, na.rm = TRUE),
          sd_vaf = sd(vafs, na.rm = TRUE),
          trend_temporal = trend_temporal,
          trend_ordinal = trend_ordinal,
          volatility = mean(abs(diff(vafs)), na.rm = TRUE),
          max_vaf = max(vafs, na.rm = TRUE),
          min_vaf = min(vafs, na.rm = TRUE),
          range_vaf = max(vafs, na.rm = TRUE) - min(vafs, na.rm = TRUE),
          cv_vaf = sd(vafs, na.rm = TRUE) / mean(vafs, na.rm = TRUE),
          time_span_days = max(time_days) - min(time_days)
        )
      }
    }, by = Patient]
    
    # Merge enhanced features
    traj_features_enhanced <- merge(traj_features, new_mutation_features, by = "Patient", all.x = TRUE)
    
    # Fill missing values
    new_mut_cols <- c("n_new_mutations", "mean_new_mut_slope", "max_new_mut_slope", 
                      "min_new_mut_slope", "sd_new_mut_slope", "new_mut_volatility")
    for (col in new_mut_cols) {
      traj_features_enhanced[is.na(get(col)), (col) := 0]
    }
    
    # Enhanced trajectory clustering
    traj_features_complete <- traj_features_enhanced[complete.cases(traj_features_enhanced)]
    
    if (nrow(traj_features_complete) >= 3) {
      message(paste("Clustering", nrow(traj_features_complete), "patient trajectories"))
      
      # Enhanced feature scaling
      feat_cols <- setdiff(names(traj_features_complete), "Patient")
      feat_mat <- scale(traj_features_complete[, ..feat_cols])
      
      # Enhanced optimal k determination for trajectories
      max_k_traj <- min(6, nrow(traj_features_complete) - 1)
      if (max_k_traj >= 2) {
        sil_traj <- sapply(2:max_k_traj, function(k) {
          pam_fit <- pam(feat_mat, k)
          mean(silhouette(pam_fit$clustering, dist(feat_mat))[, 3])
        })
        opt_k_traj <- which.max(sil_traj) + 1
        message(paste("Optimal trajectory clusters:", opt_k_traj))
        
        # Final trajectory clustering
        pam_final <- pam(feat_mat, opt_k_traj)
        traj_features_complete$trajectory_cluster <- pam_final$clustering
      } else {
        opt_k_traj <- 1
        traj_features_complete$trajectory_cluster <- 1
      }
      
      # Save enhanced trajectory features
      fwrite(traj_features_complete, paste0("enhanced_trajectory_features_", output_prefix, ".csv"))
      
      # Enhanced trajectory visualization
      if (package_status[["ggplot2"]]) {
        
        # Enhanced trend vs volatility plot
        p1 <- ggplot(traj_features_complete, aes(x = trend_temporal, y = volatility, 
                                                color = factor(trajectory_cluster))) +
          geom_point(size = 3, alpha = 0.8) +
          theme_minimal() +
          labs(title = "Enhanced Trajectory Analysis: Temporal Trend vs Volatility",
               x = "Temporal Trend (VAF change per day)",
               y = "VAF Volatility",
               color = "Trajectory Cluster")
        
        ggsave(paste0("enhanced_trajectory_trend_volatility_", output_prefix, ".png"), 
               p1, width = 10, height = 8, dpi = 300)
        
        # Enhanced new mutations analysis
        if (any(traj_features_complete$n_new_mutations > 0)) {
          p2 <- ggplot(traj_features_complete, 
                       aes(x = n_new_mutations, y = mean_new_mut_slope,
                           color = factor(trajectory_cluster),
                           size = time_span_days)) +
            geom_point(alpha = 0.8) +
            theme_minimal() +
            labs(title = "Enhanced New Mutations Analysis",
                 x = "Number of New Mutations",
                 y = "Mean VAF Slope of New Mutations",
                 color = "Trajectory Cluster",
                 size = "Follow-up Duration (days)")
          
          ggsave(paste0("enhanced_new_mutations_analysis_", output_prefix, ".png"), 
                 p2, width = 10, height = 8, dpi = 300)
        }
        
        # Enhanced PCA of trajectory features
        if (nrow(traj_features_complete) >= 3) {
          pca_feat_cols <- setdiff(names(traj_features_complete), c("Patient", "trajectory_cluster"))
          pca_traj <- prcomp(traj_features_complete[, ..pca_feat_cols], center = TRUE, scale. = TRUE)
          
          pca_df <- data.frame(
            PC1 = pca_traj$x[, 1],
            PC2 = pca_traj$x[, 2],
            Patient = traj_features_complete$Patient,
            Cluster = factor(traj_features_complete$trajectory_cluster),
            N_NewMut = traj_features_complete$n_new_mutations
          )
          
          p3 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, size = N_NewMut)) +
            geom_point(alpha = 0.8) +
            theme_minimal() +
            labs(title = "Enhanced PCA: Trajectory Feature Space",
                 color = "Trajectory Cluster",
                 size = "New Mutations") +
            theme(panel.background = element_rect(fill = "white"),
                  plot.background = element_rect(fill = "white"))
          
          ggsave(paste0("enhanced_pca_trajectory_", output_prefix, ".png"), 
                 p3, width = 10, height = 8, dpi = 300, bg = "white")
        }
      }
      
      message("Enhanced trajectory analysis completed successfully")
      
    } else {
      message("Insufficient patients for trajectory clustering")
      traj_features_complete <- NULL
    }
    
  } else {
    message("Insufficient patients with multiple timepoints for trajectory analysis")
    traj_features_complete <- NULL
  }
  
} else {
  message("No serial NGS data available for trajectory analysis")
  traj_features_complete <- NULL
}

# =============================================================================
# SECTION 6: ENHANCED VISUALIZATION AND ANALYSIS
# =============================================================================

message("=== ENHANCED VISUALIZATION AND ANALYSIS ===")

# Enhanced PCA visualization for baseline clusters
if (nrow(binary_matrix_filtered) >= 3) {
  message("Creating enhanced PCA visualization...")
  
  # Use normalized data for PCA
  mutation_mat_norm <- t(apply(binary_matrix_filtered, 1, normalize_l2))
  pca_base <- prcomp(mutation_mat_norm, center = TRUE, scale. = FALSE)
  
  pca_df <- data.frame(
    PC1 = pca_base$x[, 1],
    PC2 = pca_base$x[, 2],
    Patient_Index = kept_idx,
    Cluster = factor(final_clusters)
  )
  
  # Add diagnosis information if available
  if (!is.null(diag_col) && diag_col %in% names(df_clinical)) {
    pca_df$Diagnosis <- df_clinical[[diag_col]][pca_df$Patient_Index]
  }
  
  # Enhanced PCA plot
  if (package_status[["ggplot2"]]) {
    # Calculate variance explained
    var_explained <- round(100 * pca_base$sdev^2 / sum(pca_base$sdev^2), 1)
    
    p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster)) +
      geom_point(size = 3, alpha = 0.8) +
      theme_minimal() +
      labs(title = "Enhanced PCA: Baseline Mutation Profiles",
           x = paste0("PC1 (", var_explained[1], "% variance)"),
           y = paste0("PC2 (", var_explained[2], "% variance)"),
           color = "Baseline Cluster") +
      theme(panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white"))
    
    if (!is.null(diag_col) && diag_col %in% names(df_clinical)) {
      p_pca <- p_pca + aes(shape = factor(Diagnosis)) +
        labs(shape = "Diagnosis")
    }
    
    ggsave(paste0("enhanced_pca_baseline_", output_prefix, ".png"), 
           p_pca, width = 10, height = 8, dpi = 300, bg = "white")
  }
}

# Enhanced heatmap visualization with comprehensive annotations
if (package_status[["pheatmap"]]) {
  message("Creating enhanced comprehensive heatmaps...")
  
  # Prepare enhanced binary display matrix
  binary_display <- binary_matrix[, colSums(binary_matrix) > 0, drop = FALSE]
  
  # Create comprehensive row annotations with clinical data
  row_annot <- data.frame(
    Baseline_Cluster = factor(ifelse(is.na(baseline_cluster), "No_Mutations", 
                                   paste0("Cluster_", baseline_cluster)))
  )
  
  # Add trajectory cluster if available
  if (!is.null(traj_features_complete)) {
    # Map trajectory clusters back to all patients
    traj_cluster_all <- rep(NA, nrow(df_clinical))
    for (i in 1:nrow(traj_features_complete)) {
      patient_match <- which(as.character(1:nrow(df_clinical)) == as.character(traj_features_complete$Patient[i]))
      if (length(patient_match) > 0) {
        traj_cluster_all[patient_match] <- traj_features_complete$trajectory_cluster[i]
      }
    }
    row_annot$Trajectory_Cluster <- factor(ifelse(is.na(traj_cluster_all), "No_Trajectory", 
                                                paste0("Traj_", traj_cluster_all)))
  }
  
  rownames(row_annot) <- rownames(binary_display)
  
  # Enhanced color schemes for annotations
  cluster_colors <- rainbow(length(unique(baseline_cluster[!is.na(baseline_cluster)])))
  names(cluster_colors) <- paste0("Cluster_", sort(unique(baseline_cluster[!is.na(baseline_cluster)])))
  cluster_colors["No_Mutations"] <- "gray90"
  
  annotation_colors <- list(Baseline_Cluster = cluster_colors)
  
  # Add diagnosis colors if available
  if ("Diagnosis" %in% names(row_annot)) {
    diag_levels <- levels(row_annot$Diagnosis)
    n_diag <- length(diag_levels)
    if (n_diag > 1) {
      # RColorBrewer requires at least 3 colors for Set3
      if (n_diag < 3) {
        # Use a manual palette or repeat colors
        base_colors <- RColorBrewer::brewer.pal(3, "Set3")[1:n_diag]
        annotation_colors$Diagnosis <- base_colors
      } else {
        annotation_colors$Diagnosis <- RColorBrewer::brewer.pal(min(n_diag, 11), "Set3")[1:n_diag]
      }
      names(annotation_colors$Diagnosis) <- diag_levels
    }
  }
  
  # Add survival colors if available
  if ("Survival_Status" %in% names(row_annot)) {
    annotation_colors$Survival_Status <- c("Death" = "red", "Alive" = "green")
  }
  
  # Add trajectory colors if available
  if ("Trajectory_Cluster" %in% names(row_annot)) {
    traj_levels <- levels(row_annot$Trajectory_Cluster)
    traj_colors <- RColorBrewer::brewer.pal(min(length(traj_levels), 11), "Paired")[1:length(traj_levels)]
    names(traj_colors) <- traj_levels
    traj_colors["No_Trajectory"] <- "gray80"
    annotation_colors$Trajectory_Cluster <- traj_colors
  }
  
  # Create enhanced global mutation heatmap
  pheatmap(binary_display,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = FALSE,
           show_colnames = TRUE,
           fontsize_col = 6,
           fontsize_row = 6,
           fontsize = 8,
           annotation_row = row_annot,
           annotation_colors = annotation_colors,
           color = c("white", "navy"),
           legend_breaks = c(0, 1),
           legend_labels = c("No Mutation", "Mutation Present"),
           main = "Comprehensive Mutation Profile by Clinical Clusters",
           filename = paste0("enhanced_baseline_mutation_heatmap_", output_prefix, ".png"),
           width = max(12, min(25, ncol(binary_display) * 0.3)),
           height = max(10, min(20, nrow(binary_display) * 0.15)))
  
  # Create detailed per-cluster mutation heatmaps with enhanced features
  for (cl in sort(unique(baseline_cluster[!is.na(baseline_cluster)]))) {
    cluster_patients <- which(baseline_cluster == cl)
    if (length(cluster_patients) > 1) {
      cluster_matrix <- binary_matrix[cluster_patients, , drop = FALSE]
      
      # Include genes with any variation in this cluster
      variant_genes <- colSums(cluster_matrix) > 0
      
      if (sum(variant_genes) > 1) {
        cluster_matrix <- cluster_matrix[, variant_genes, drop = FALSE]
        
        # Create cluster-specific annotations
        cluster_row_annot <- row_annot[rownames(cluster_matrix), , drop = FALSE]
        
        # Add mutation burden for this cluster
        mutation_burden <- rowSums(cluster_matrix)
        # DEBUG: Print mutation_burden summary
        print("[DEBUG] mutation_burden summary:")
        print(summary(mutation_burden))
        print("[DEBUG] unique mutation_burden values:")
        print(unique(mutation_burden))
        
        breaks <- unique(quantile(mutation_burden, probs = c(0, 0.33, 0.67, 1)))
        print("[DEBUG] breaks for cut:")
        print(breaks)
        if (length(breaks) < 2) {
          # fallback: all in one bin
          cluster_row_annot$Mutation_Burden <- factor(rep("All", length(mutation_burden)))
        } else {
          cluster_row_annot$Mutation_Burden <- cut(mutation_burden, 
                                                  breaks = breaks, 
                                                  labels = c("Low", "Medium", "High")[1:(length(breaks)-1)],
                                                  include.lowest = TRUE)
        }
        
        # Enhanced color scheme for cluster-specific analysis
        cluster_annotation_colors <- annotation_colors
        cluster_annotation_colors$Mutation_Burden <- c("Low" = "lightblue", "Medium" = "orange", "High" = "darkred")
        
        pheatmap(cluster_matrix,
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 show_rownames = FALSE,
                 show_colnames = TRUE,
                 fontsize_col = 8,
                 fontsize_row = 6,
                 fontsize = 8,
                 annotation_row = cluster_row_annot,
                 annotation_colors = cluster_annotation_colors,
                 color = c("white", "darkblue"),
                 legend_breaks = c(0, 1),
                 legend_labels = c("No Mutation", "Mutation"),
                 main = paste("Cluster", cl, "Detailed Mutation Profile (n=", length(cluster_patients), ")"),
                 filename = paste0("enhanced_cluster_", cl, "_detailed_heatmap_", output_prefix, ".png"),
                 width = max(10, min(25, ncol(cluster_matrix) * 0.4)),
                 height = max(8, min(20, nrow(cluster_matrix) * 0.25)))
      }
    }
  }
  
  message("Enhanced heatmaps created successfully")
}

# Enhanced statistical summary
# =============================================================================
# SECTION 7: COMPREHENSIVE COX REGRESSION SURVIVAL ANALYSIS
# =============================================================================

message("=== COMPREHENSIVE COX REGRESSION SURVIVAL ANALYSIS ===")

# Cox Regression Analysis for Baseline Clusters
if (package_status[["survival"]] && package_status[["survminer"]] && 
    all(c("OS months", "event") %in% names(df_clinical))) {
  
  message("Performing Cox regression analysis for baseline clusters...")
  
  # Prepare survival data
  surv_data <- df_clinical[!is.na(df_clinical$`OS months`) & !is.na(df_clinical$event), ]
  surv_data$baseline_cluster_factor <- factor(surv_data$baseline_cluster)
  
  if (nrow(surv_data) > 10 && length(unique(surv_data$baseline_cluster_factor)) > 1) {
    
    # Univariate Cox regression for baseline clusters
    cox_baseline <- coxph(Surv(`OS months`, event) ~ baseline_cluster_factor, data = surv_data)
    
    message("Baseline Cluster Cox Regression Results:")
    print(summary(cox_baseline))
    
    # Create Kaplan-Meier survival curves
    if (package_status[["survminer"]]) {
      surv_fit_baseline <- survfit(Surv(`OS months`, event) ~ baseline_cluster_factor, data = surv_data)
      
      # Enhanced Kaplan-Meier plot
      km_plot <- ggsurvplot(
        surv_fit_baseline,
        data = surv_data,
        pval = TRUE,
        conf.int = TRUE,
        risk.table = TRUE,
        risk.table.col = "strata",
        linetype = "strata",
        surv.median.line = "hv",
        ggtheme = theme_minimal(),
        palette = rainbow(length(unique(surv_data$baseline_cluster_factor))),
        title = "Survival Analysis by Baseline Mutation Clusters",
        xlab = "Time (months)",
        ylab = "Overall Survival Probability",
        legend.title = "Baseline Cluster",
        legend.labs = paste("Cluster", levels(surv_data$baseline_cluster_factor))
      )
      
      # Save Kaplan-Meier plot
      ggsave(paste0("enhanced_kaplan_meier_baseline_", output_prefix, ".png"), 
             km_plot$plot, width = 12, height = 10, dpi = 300, bg = "white")
      
      # Save risk table
      ggsave(paste0("enhanced_risk_table_baseline_", output_prefix, ".png"), 
             km_plot$table, width = 12, height = 4, dpi = 300, bg = "white")
    }
    
    # Multivariable Cox regression if additional covariates available
    covariate_cols <- c()
    if (!is.null(diag_col) && diag_col %in% names(surv_data)) {
      covariate_cols <- c(covariate_cols, diag_col)
    }
    
    # Add age if available
    age_cols <- grep("age|Age", names(surv_data), value = TRUE)
    if (length(age_cols) > 0) {
      age_col <- age_cols[1]
      surv_data[[age_col]] <- suppressWarnings(as.numeric(surv_data[[age_col]]))
      if (sum(!is.na(surv_data[[age_col]])) > nrow(surv_data) * 0.5) {
        covariate_cols <- c(covariate_cols, age_col)
      }
    }
    
    # Add gender if available
    gender_cols <- grep("gender|sex|Sex|Gender", names(surv_data), value = TRUE)
    if (length(gender_cols) > 0) {
      covariate_cols <- c(covariate_cols, gender_cols[1])
    }
    
    if (length(covariate_cols) > 0) {
      message("Performing multivariable Cox regression...")
      
      # Create formula for multivariable analysis
      formula_str <- paste("Surv(`OS months`, event) ~ baseline_cluster_factor +", 
                          paste(paste0("`", covariate_cols, "`"), collapse = " + "))
      cox_formula <- as.formula(formula_str)
      
      # Before running Cox regression, check for collinearity and zero-variance predictors
      if (exists("cox_formula")) {
        # Extract variable names from formula
        vars_in_model <- all.vars(cox_formula)[-1] # remove response
        # Check for zero or near-zero variance
        nzv_vars <- vars_in_model[sapply(vars_in_model, function(v) {
          x <- surv_data[[v]]
          if (is.factor(x) || is.character(x)) {
            length(unique(x[!is.na(x)])) <= 1
          } else {
            var(x, na.rm = TRUE) == 0 || all(is.na(x))
          }
        })]
        if (length(nzv_vars) > 0) {
          message("[DEBUG] Removing zero-variance or constant variables from Cox model:", paste(nzv_vars, collapse=", "))
          # Remove from formula
          vars_to_keep <- setdiff(vars_in_model, nzv_vars)
          if (length(vars_to_keep) > 0) {
            formula_str <- paste("Surv(`OS months`, event) ~ baseline_cluster_factor +", 
                                paste(paste0("`", vars_to_keep, "`"), collapse = " + "))
            cox_formula <- as.formula(formula_str)
          } else {
            stop("No variables left for Cox regression after removing zero-variance predictors.")
          }
        }
      }
      
      tryCatch({
        cox_multivariable <- coxph(cox_formula, data = surv_data)
        message("Multivariable Cox Regression Results:")
        print(summary(cox_multivariable))
        
        # Create forest plot for multivariable results
        if (package_status[["survminer"]]) {
          forest_plot <- ggforest(cox_multivariable, data = surv_data,
                                main = "Multivariable Cox Regression - Baseline Clusters",
                                fontsize = 0.8)
          
          ggsave(paste0("enhanced_forest_plot_baseline_", output_prefix, ".png"), 
                 forest_plot, width = 12, height = 8, dpi = 300, bg = "white")
        }
        
      }, error = function(e) {
        message(paste("Multivariable Cox regression failed:", e$message))
      })
    }
    
  } else {
    message("Insufficient data for Cox regression analysis")
  }
  
  # Cox Regression Analysis for Trajectory Clusters (if available)
  if (!is.null(traj_features_complete) && nrow(traj_features_complete) > 5) {
    message("Performing Cox regression analysis for trajectory clusters...")
    # Identify patient ID column in df_clinical
    id_candidates <- names(df_clinical)[grepl("MRN|Patient|Subject|ID", names(df_clinical), ignore.case = TRUE)]
    if (length(id_candidates) > 0) {
      id_col <- id_candidates[1]
    } else {
      stop("No patient identifier column found in df_clinical")
    }
    # Ensure surv_data has Patient column
    surv_data$Patient <- as.character(surv_data[[id_col]])
    # Map traj_features_complete$Patient to the same ID values as in df_clinical
    # If traj_features_complete$Patient is numeric and matches row numbers, map to df_clinical[[id_col]]
    if (is.numeric(traj_features_complete$Patient) && all(traj_features_complete$Patient %in% seq_len(nrow(df_clinical)))) {
      traj_features_complete$Patient <- as.character(df_clinical[[id_col]][as.integer(traj_features_complete$Patient)])
    } else {
      traj_features_complete$Patient <- as.character(traj_features_complete$Patient)
    }
    # DEBUG: Print column names and sample IDs before merge
    print("[DEBUG] surv_data columns:")
    print(colnames(surv_data))
    print("[DEBUG] traj_features_complete columns:")
    print(colnames(traj_features_complete))
    print("[DEBUG] Sample surv_data$Patient:")
    print(head(surv_data$Patient))
    print("[DEBUG] Sample traj_features_complete$Patient:")
    print(head(traj_features_complete$Patient))
    # Merge trajectory data with survival data by 'Patient'
    if ("Patient" %in% colnames(surv_data) && "Patient" %in% colnames(traj_features_complete)) {
      traj_surv_data <- merge(surv_data, traj_features_complete, by = "Patient", all.x = FALSE)
    } else {
      stop("[ERROR] Cannot merge: 'Patient' column not found in both data.tables.")
    }
    if (nrow(traj_surv_data) > 5 && length(unique(traj_surv_data$trajectory_cluster)) > 1) {
      traj_surv_data$trajectory_cluster_factor <- factor(traj_surv_data$trajectory_cluster)
      
      # Univariate Cox regression for trajectory clusters
      cox_trajectory <- coxph(Surv(`OS months`, event) ~ trajectory_cluster_factor, data = traj_surv_data)
      
      message("Trajectory Cluster Cox Regression Results:")
      print(summary(cox_trajectory))
      
      # Create Kaplan-Meier survival curves for trajectory clusters
      if (package_status[["survminer"]]) {
        surv_fit_trajectory <- survfit(Surv(`OS months`, event) ~ trajectory_cluster_factor, data = traj_surv_data)
        
        km_plot_traj <- ggsurvplot(
          surv_fit_trajectory,
          data = traj_surv_data,
          pval = TRUE,
          conf.int = TRUE,
          risk.table = TRUE,
          risk.table.col = "strata",
          linetype = "strata",
          surv.median.line = "hv",
          ggtheme = theme_minimal(),
          palette = RColorBrewer::brewer.pal(length(unique(traj_surv_data$trajectory_cluster_factor)), "Set2"),
          title = "Survival Analysis by Trajectory Clusters",
          xlab = "Time (months)",
          ylab = "Overall Survival Probability",
          legend.title = "Trajectory Cluster",
          legend.labs = paste("Trajectory", levels(traj_surv_data$trajectory_cluster_factor))
        )
        
        ggsave(paste0("enhanced_kaplan_meier_trajectory_", output_prefix, ".png"), 
               km_plot_traj$plot, width = 12, height = 10, dpi = 300, bg = "white")
      }
      
      # Combined Cox model with both baseline and trajectory clusters
      tryCatch({
        cox_combined <- coxph(Surv(`OS months`, event) ~ baseline_cluster_factor + trajectory_cluster_factor, 
                            data = traj_surv_data)
        message("Combined Baseline + Trajectory Cox Regression Results:")
        print(summary(cox_combined))
        
        # Forest plot for combined model
        if (package_status[["survminer"]]) {
          forest_plot_combined <- ggforest(cox_combined, data = traj_surv_data,
                                         main = "Combined Cox Regression - Baseline + Trajectory Clusters",
                                         fontsize = 0.8)
          
          ggsave(paste0("enhanced_forest_plot_combined_", output_prefix, ".png"), 
                 forest_plot_combined, width = 12, height = 8, dpi = 300, bg = "white")
        }
        
      }, error = function(e) {
        message(paste("Combined Cox regression failed:", e$message))
      })
    }
  }
  
  # Progression-Free Survival Analysis (if available)
  if ("Progression Free Survival Months" %in% names(df_clinical) && "pfs_event" %in% names(df_clinical)) {
    message("Performing progression-free survival analysis...")
    
    pfs_data <- df_clinical[!is.na(df_clinical$`Progression Free Survival Months`) & 
                           !is.na(df_clinical$pfs_event), ]
    pfs_data$baseline_cluster_factor <- factor(pfs_data$baseline_cluster)
    
    if (nrow(pfs_data) > 10 && length(unique(pfs_data$baseline_cluster_factor)) > 1) {
      cox_pfs <- coxph(Surv(`Progression Free Survival Months`, pfs_event) ~ baseline_cluster_factor, 
                      data = pfs_data)
      
      message("Progression-Free Survival Cox Regression Results:")
      print(summary(cox_pfs))
      
      # PFS Kaplan-Meier curves
      if (package_status[["survminer"]]) {
        surv_fit_pfs <- survfit(Surv(`Progression Free Survival Months`, pfs_event) ~ baseline_cluster_factor, 
                               data = pfs_data)
        
        km_plot_pfs <- ggsurvplot(
          surv_fit_pfs,
          data = pfs_data,
          pval = TRUE,
          conf.int = TRUE,
          risk.table = TRUE,
          ggtheme = theme_minimal(),
          title = "Progression-Free Survival by Baseline Clusters",
          xlab = "Time (months)",
          ylab = "Progression-Free Survival Probability"
        )
        
        ggsave(paste0("enhanced_kaplan_meier_pfs_", output_prefix, ".png"), 
               km_plot_pfs$plot, width = 12, height = 10, dpi = 300, bg = "white")
      }
    }
  }
  
  message("Cox regression survival analysis completed")
  
} else {
  missing_packages <- c()
  if (!package_status[["survival"]]) missing_packages <- c(missing_packages, "survival")
  if (!package_status[["survminer"]]) missing_packages <- c(missing_packages, "survminer")
  
  missing_data <- c()
  if (!"OS months" %in% names(df_clinical)) missing_data <- c(missing_data, "OS months")
  if (!"event" %in% names(df_clinical)) missing_data <- c(missing_data, "event")
  
  if (length(missing_packages) > 0) {
    message(paste("Skipping Cox regression - missing packages:", paste(missing_packages, collapse = ", ")))
  }
  if (length(missing_data) > 0) {
    message(paste("Skipping Cox regression - missing data columns:", paste(missing_data, collapse = ", ")))
  }
}

message("=== ENHANCED ANALYSIS SUMMARY ===")

# Create comprehensive summary report
summary_stats <- list(
  "Analysis Type" = "Enhanced Ward Clustering with Trajectory Analysis + Cox Regression",
  "Total Patients" = nrow(df_clinical),
  "Patients with Mutations" = length(kept_idx),
  "Genes Analyzed" = ncol(mutation_matrix),
  "Optimal Baseline Clusters" = opt_k_final,
  "Average Silhouette Width" = round(avg_sil, 3),
  "Clustering Method" = "Ward.D2 with Jaccard Distance",
  "Survival Analysis" = ifelse(all(c("OS months", "event") %in% names(df_clinical)), "Cox Regression Completed", "No Survival Data"),
  "Cox Regression Features" = "Baseline Clusters, Trajectory Clusters, Multivariable Models"
)

if (!is.null(traj_features_complete)) {
  summary_stats <- c(summary_stats, list(
    "Patients with Trajectories" = nrow(traj_features_complete),
    "Trajectory Clusters" = max(traj_features_complete$trajectory_cluster),
    "Patients with New Mutations" = sum(traj_features_complete$n_new_mutations > 0),
    "Total New Mutations Detected" = sum(traj_features_complete$n_new_mutations)
  ))
}

# Save summary report
capture.output({
  cat("=== ENHANCED WARD CLUSTERING ANALYSIS REPORT ===\n\n")
  cat("Analysis completed at:", format(Sys.time()), "\n\n")
  
  for (i in 1:length(summary_stats)) {
    cat(names(summary_stats)[i], ":", summary_stats[[i]], "\n")
  }
  
  cat("\nBaseline Cluster Distribution:\n")
  print(table(baseline_cluster, useNA = "ifany"))
  
  if (!is.null(traj_features_complete)) {
    cat("\nTrajectory Cluster Distribution:\n")
    print(table(traj_features_complete$trajectory_cluster))
  }
  
  cat("\nPackage Availability:\n")
  for (pkg in names(package_status)) {
    status <- if (package_status[[pkg]]) "Available" else "Missing"
    cat("  ", pkg, ":", status, "\n")
  }
  
}, file = paste0("enhanced_analysis_summary_", output_prefix, ".txt"))

message("=== ENHANCED WARD CLUSTERING ANALYSIS COMPLETED ===")
message(paste("Total runtime:", round(proc.time()[3], 2), "seconds"))
message("All enhanced results saved with prefix:", output_prefix)
message(paste("Summary report saved to: enhanced_analysis_summary_", output_prefix, ".txt"))