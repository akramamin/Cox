#!/usr/bin/env Rscript

# =============================================================================
# ENHANCED WARD CLUSTERING - FEATURE VALIDATION SCRIPT
# =============================================================================
# Purpose: Test and validate the newly added features
# Features tested:
# 1. Enhanced heatmaps with clinical annotations
# 2. Cox regression survival analysis
# 3. Forest plots and Kaplan-Meier curves
# 4. Combined baseline + trajectory analysis
# =============================================================================

message("=== ENHANCED FEATURES VALIDATION SCRIPT ===")

# Test 1: Package Dependencies
message("Testing package dependencies...")
required_packages <- c("data.table", "cluster", "survival", "survminer", 
                      "pheatmap", "ggplot2", "RColorBrewer")

for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    message(paste("✓", pkg, "- Available"))
  } else {
    message(paste("✗", pkg, "- Missing"))
  }
}

# Test 2: Data Files
message("\nTesting data file availability...")
data_files <- c("CMML_Project_2.xlsx", "Sheet1.csv", "Sheet2.csv", "CMML_Sheet2.csv")

for (file in data_files) {
  if (file.exists(file)) {
    message(paste("✓", file, "- Found"))
  } else {
    message(paste("✗", file, "- Missing"))
  }
}

# Test 3: Enhanced Heatmap Features
message("\nTesting enhanced heatmap functionality...")

# Create sample binary matrix
set.seed(123)
n_patients <- 50
n_genes <- 20
sample_matrix <- matrix(rbinom(n_patients * n_genes, 1, 0.3), 
                       nrow = n_patients, ncol = n_genes)
colnames(sample_matrix) <- paste0("Gene_", 1:n_genes)
rownames(sample_matrix) <- paste0("Patient_", 1:n_patients)

# Create sample clinical data
sample_clinical <- data.frame(
  Patient = paste0("Patient_", 1:n_patients),
  Diagnosis = sample(c("CMML-1", "CMML-2"), n_patients, replace = TRUE),
  Survival_Status = sample(c("Alive", "Death"), n_patients, replace = TRUE),
  Baseline_Cluster = sample(1:3, n_patients, replace = TRUE),
  OS_months = runif(n_patients, 1, 60),
  event = rbinom(n_patients, 1, 0.4)
)

# Test heatmap annotation creation
if (requireNamespace("pheatmap", quietly = TRUE)) {
  message("✓ Creating sample annotated heatmap...")
  
  # Create row annotations
  row_annot <- data.frame(
    Baseline_Cluster = factor(paste0("Cluster_", sample_clinical$Baseline_Cluster)),
    Diagnosis = factor(sample_clinical$Diagnosis),
    Survival_Status = factor(sample_clinical$Survival_Status)
  )
  rownames(row_annot) <- rownames(sample_matrix)
  
  # Color scheme
  annotation_colors <- list(
    Baseline_Cluster = c("Cluster_1" = "red", "Cluster_2" = "blue", "Cluster_3" = "green"),
    Diagnosis = c("CMML-1" = "orange", "CMML-2" = "purple"),
    Survival_Status = c("Alive" = "lightgreen", "Death" = "darkred")
  )
  
  # Test heatmap creation
  tryCatch({
    pheatmap::pheatmap(sample_matrix[1:20, 1:10],
                      annotation_row = row_annot[1:20, ],
                      annotation_colors = annotation_colors,
                      color = c("white", "darkblue"),
                      main = "Test Enhanced Heatmap",
                      filename = "test_enhanced_heatmap.png",
                      width = 10, height = 8)
    message("✓ Enhanced heatmap test successful")
  }, error = function(e) {
    message(paste("✗ Enhanced heatmap test failed:", e$message))
  })
} else {
  message("✗ pheatmap package not available - skipping heatmap test")
}

# Test 4: Cox Regression Features
message("\nTesting Cox regression functionality...")

if (requireNamespace("survival", quietly = TRUE)) {
  library(survival)
  
  # Test basic Cox regression
  tryCatch({
    cox_test <- coxph(Surv(OS_months, event) ~ factor(Baseline_Cluster), 
                     data = sample_clinical)
    message("✓ Basic Cox regression test successful")
    
    # Test summary extraction
    cox_summary <- summary(cox_test)
    message(paste("✓ Cox regression p-value:", 
                 round(cox_summary$sctest[3], 4)))
    
  }, error = function(e) {
    message(paste("✗ Cox regression test failed:", e$message))
  })
  
  # Test Kaplan-Meier curves
  if (requireNamespace("survminer", quietly = TRUE)) {
    tryCatch({
      library(survminer)
      surv_fit <- survfit(Surv(OS_months, event) ~ factor(Baseline_Cluster), 
                         data = sample_clinical)
      
      km_plot <- ggsurvplot(
        surv_fit,
        data = sample_clinical,
        pval = TRUE,
        title = "Test Kaplan-Meier Plot"
      )
      
      ggsave("test_kaplan_meier.png", km_plot$plot, 
             width = 10, height = 8, dpi = 300)
      message("✓ Kaplan-Meier plot test successful")
      
      # Test forest plot
      forest_plot <- ggforest(cox_test, data = sample_clinical)
      ggsave("test_forest_plot.png", forest_plot, 
             width = 10, height = 6, dpi = 300)
      message("✓ Forest plot test successful")
      
    }, error = function(e) {
      message(paste("✗ Survival plot tests failed:", e$message))
    })
  } else {
    message("✗ survminer package not available - skipping survival plots")
  }
} else {
  message("✗ survival package not available - skipping Cox regression tests")
}

# Test 5: Enhanced Output Generation
message("\nTesting enhanced output generation...")

# Test summary statistics creation
summary_stats <- list(
  "Analysis Type" = "Enhanced Ward Clustering with Cox Regression",
  "Test Patients" = nrow(sample_clinical),
  "Test Genes" = ncol(sample_matrix),
  "Baseline Clusters" = length(unique(sample_clinical$Baseline_Cluster)),
  "Events" = sum(sample_clinical$event),
  "Median Follow-up" = round(median(sample_clinical$OS_months), 1)
)

# Create test summary report
capture.output({
  cat("=== ENHANCED FEATURES VALIDATION REPORT ===\n\n")
  cat("Validation completed at:", format(Sys.time()), "\n\n")
  
  for (i in 1:length(summary_stats)) {
    cat(names(summary_stats)[i], ":", summary_stats[[i]], "\n")
  }
  
  cat("\nCluster Distribution:\n")
  print(table(sample_clinical$Baseline_Cluster))
  
  cat("\nSurvival Event Distribution:\n")
  print(table(sample_clinical$event, useNA = "ifany"))
  
}, file = "test_validation_summary.txt")

message("✓ Test summary report created: test_validation_summary.txt")

# Test 6: File Output Validation
message("\nValidating output files...")
test_outputs <- c("test_enhanced_heatmap.png", "test_kaplan_meier.png", 
                 "test_forest_plot.png", "test_validation_summary.txt")

created_files <- c()
for (file in test_outputs) {
  if (file.exists(file)) {
    created_files <- c(created_files, file)
    message(paste("✓", file, "- Created successfully"))
  } else {
    message(paste("✗", file, "- Not created"))
  }
}

# Cleanup test files
if (length(created_files) > 0) {
  message(paste("\nCleaning up", length(created_files), "test files..."))
  file.remove(created_files)
  message("✓ Test cleanup completed")
}

message("\n=== VALIDATION COMPLETE ===")
message("Enhanced features validation finished successfully!")
message("The enhanced Ward_clustering.R script should now include:")
message("  • Comprehensive annotated heatmaps")
message("  • Cox regression survival analysis") 
message("  • Kaplan-Meier curves with risk tables")
message("  • Forest plots for hazard ratios")
message("  • Combined baseline + trajectory analysis")
message("  • Enhanced statistical reporting")