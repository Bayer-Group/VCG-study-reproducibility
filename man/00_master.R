#*********************************************************************************;
#* Program Identification ***********************    RCSS / Bayer AG             **;
#*;
# clear the environment first
# remove (almost) everything in the working environment
gc()
rm(list = ls(all.names = TRUE))
# name of the program
program <- "00_master.R"
#* Version:           <1 [2023_04_20]>;
#* Project:           <VCG>;
#* Author:            <Alexander Gurjanov>;
#*;
#*********************************************************************************;
#* Program Description  **********************************************************;
#*;
#*This program is the execution file of the calculation and visualization steps
#*of the remaining R-scripts.
#*The functions of the R-scripts are executed here with respect
#*to the entered values.
#*;
#*********************************************************************************;
#*********************************************************************************;
#*********************************************************************************;
#* Program Specification  ********************************************************;
#*;
#* This program reads and executes the function from all other R-scripts:
#* "01_write_legacy_study_original_results.R" extracts the original results from
#* the CCG and the respective dose groups and writes them as readable tables.
#* <input:  Electrolyte values of the VCG data set and legacy studies> ;
#* <steps:  Calculating original results of significance tests> ;
#* <output: XLSX file showing original results.> ;
#*;
#*********************************************************************************;
#*********************************************************************************;
# retrieve the current study path
rootpath <- sub("*\\/man", "",dirname(rstudioapi::getSourceEditorContext()$path))

# path to folder for results / graphs / derived data sets
path_res <- paste(rootpath,'/Results', sep='')
#* 
#* DO NOT EDIT BELOW THESE LINES
#*

# path to folder derived data sets
der <- paste(rootpath,'/data/Derived', sep='')

cat("---------------------",program,"---------------------","\n\n")
cat(date(),"\n\n")
cat("\n",version$version.string,"\n\n")

setwd(path_res)

cat("\n Results will be stored in folder:",getwd()," \n")

#* 
#* DO NOT EDIT ABOVE THESE LINES
#*

##### program starts here
#******************************************************************************
#******************************************************************************
#******************************************************************************
#Load libraries
#for fast loading of large tables
library(data.table)
#for effective data processing
library(tidyverse)
#for reading Excel sheets
library(openxlsx)
#for calculating significance tests (Dunnett's and heterogenicity Dunnett's)
library(PMCMRplus)
#for calculating significance test (Wilcox test)
library(rstatix)
#for calculating effect size
library(effsize)
#for creating nice heat-map like tables
library(gt)
#for extracting gt and plotly HTML objects as static images
library(webshot2)
#for missing data imputation
library(mice)
#cleaning corrupt parameter names
library(janitor)
#writing text in ggplot plots
library(grid)
#combine several plots into one
library(cowplot)
#combining plots with panel titles
library(gridExtra)
#calculate the power of simple tests with two groups (e.g., t-test)
library(pwr)
#split work on many cores
library(parallel)
#*******************************************************************************
#set random seed for reproducibility of results
set.seed(1)
#*******************************************************************************
#Source all scripts----
#Read all files and impute missing values
source(paste0(rootpath, "/man/03_write_HCD_to_R_table.R"))
#print table of missing data
source(paste0(rootpath, "/man/03a1_missing_data_calculator.R"))
#Resampling experiment
Sys.time()
source(paste0(rootpath, "/man/04_resampling.R"))
Sys.time()
#Visualize as heatmap-tables
source(paste0(rootpath, "/man/05_visualize_as_table.R"))
#print effect sizes
source(paste0(rootpath, "/man/06_visualize_as_plot.R"))
#print MI results
source(paste0(rootpath, "/man/06a2_MI_plot.R"))

#*******************************************************************************
# ##### program ends here
# save.image(file = paste(rootpath,"\\Programs\\",program,".RData",sep=""))
# #* 
# #* DO NOT EDIT BELOW THESE LINES
# #*
# savehistory(file = paste(rootpath,"\\Programs\\",program,".Rhistory",sep=""))
# sessionInfo()