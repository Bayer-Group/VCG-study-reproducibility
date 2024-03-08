#03_write_HCD_to_R_table.R

#*This script extracts the historical control data stored in
#*~/data/Original/quantitative_parameters and prepares them into a table in the
#*same fashion as in "01_write_legacy_study_original_results.R". The same
#*function is used for that.
#*******************************************************************************
#*******************************************************************************
#*Load libraries----
#*******************************************************************************
require(data.table)
require(tidyverse)
require(openxlsx)
#*******************************************************************************
#*******************************************************************************
#*Read scripts----
#*******************************************************************************
source(paste0(rootpath, "/man/01_write_SEND_to_R_table.R"))
source(paste0(rootpath, "/man/01a1_write_qualitative_SEND_data_to_R_table.R"))
source(paste0(rootpath, "/man/04a1_simulate_missing_parameters.R"))
source(paste0(rootpath, "/man/03a2_HCD_visualizer.R"))
source(paste0(rootpath, "/man/03a2_1_HCD_visualizer_with_coloring.R"))
#*******************************************************************************
#*******************************************************************************
#*Create function which reads LB, BW, and OM parameters from the legacy_study
#*folder and joins them to one data frame----
#*******************************************************************************


#Read HCD quantitative parameters
HCD <- SEND_reader(study = "HCD_RAT")

#Read HCD for recovery groups
HCD_recovery <- SEND_reader(
  study = "HCD_RAT",
  #cap body weight observation to day 31 (the day when LB parameters were read)
  lbdylow = 36,
  lbdyhigh = 49,
  omdylow = 36,
  omdyhigh = 49,
  bwdylow = 32,
  bwdyhigh = 49
  )

#Visualize HCD (no imputation for missing parameters) un-comment if needed
# HCD_visualizer(data_to_plot = HCD, plottitle = "HCD_no_imputation")

#colors by treatment vehicle or a column of your choice. Uncomment, if needed
# HCD_visualizer_with_coloring(
#   data_to_plot = HCD_RAT,
#   plottitle = "HCD_RAT_no_TRTV_control",
#   color_by = "TRTV"
# )
# HCD_visualizer_with_coloring(
#   data_to_plot = HCD,
#   plottitle = "HCD_one_TRTV",
#   color_by = NULL
# )
# #*******************************************************************************
# #*******************************************************************************
# #Impute missing values----
# #*******************************************************************************
# #*Use "04a1_simulate_missing_parameters.R" for this.
# #*******************************************************************************
# #*Imputation by taking the median value of non-missing data
HCD_imputed_median <- missing_data_imputation(
  imputation_method = "median",
  studydata = HCD
)

# #Visualize HCD (no imputation for missing parameters) un-comment if needed
# HCD_visualizer(
#   data_to_plot = HCD_imputed_median,
#   plottitle = "HCD_imputation_median"
#   )
# #*******************************************************************************
# #*Imputation by taking random values from non-missing data
HCD_imputed_rs <- missing_data_imputation(
  imputation_method = "random_sampling",
  studydata = HCD
)
# 
# #Visualize HCD (no imputation for missing parameters) un-comment if needed
# HCD_visualizer(
#   data_to_plot = HCD_imputed_rs,
#   plottitle = "HCD_random_sampling"
# )
# 
# #*******************************************************************************
# #*Imputation by predictive mean matching (pmm)
HCD_imputed_pmm <- missing_data_imputation(
  imputation_method = "pmm",
  studydata = HCD
)
# 
# #Visualize HCD (no imputation for missing parameters) un-comment if needed
# HCD_visualizer(
#   data_to_plot = HCD_imputed_pmm,
#   plottitle = "HCD_imputed_pmm"
# )
