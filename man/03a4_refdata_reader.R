#*03a4_refdata_reader.R
#*
#*This script reads the reference data which is entered as DOCX documents
#*and transforms them in to tables easily readable by R. These reference values
#*are then used to compare whether animal values are within or outside of the
#*HCD 2s reference range.
#*******************************************************************************
#*******************************************************************************
#*read libraries----
#*******************************************************************************
require(tidyverse)
require(openxlsx)
#*******************************************************************************
#*******************************************************************************
#Read mapping table for harmonization of parameter names----
#*******************************************************************************
data_mapping <- read.xlsx(
  paste0(
    rootpath, "/Data/LB_unit_mapping.xlsx"
  )
) %>%
  rename(
    "LBTEST" = "LBTEST_old",
    "LBSPEC" = "LBSPEC_old",
    "LBORRESU" = "LBORRESU_old"
  )
#*******************************************************************************
#*******************************************************************************
#Read the DOCX document----
#*******************************************************************************
#*Create a nested function which reads the respective reference tables and
#*combines the quantitative parameters into one.
refdata_reader <- function(study, param){
  #read parameters
  params <- read.xlsx(
    paste0(
      rootpath, "/Data/Original/legacy_studies/",
      study,
      "/",
      param,
      ".xlsx"
    )
  )
  
  #remove leading and trailing whitespace
  params_stripped <- params %>%
    mutate(
      across(
        .cols = intersect(
          c(
            "n", "N", "Mean", "S.D.", "Range-2s", "Range+2s", "Min",
            "Percentile2.5th", "Percentile5th", "Median", "Percentile95th",
            "Percentile97.5th", "Max", "Grading/Percent1", "Grading/Percent2",
            "Grading/Percent3", "Grading/Percent4"
          ),
          names(.data)
        ),
        .fns = ~if_else(.col %in% names(.data), trimws(.x, "l"), .x)
      ),
      #*if a "min" is present in the minimum range, the value of the "Min" column
      #*is taken instead (to avoid having negative biological numbers)
      `Range-2s` = if_else(
        `Range-2s` == "Min", as.character(Min), as.character(`Range-2s`)
      ),
      across(
        #turn all numbers into numerics
        .cols = intersect(
          c(
            "n", "N", "Mean", "S.D.", "Range-2s", "Range+2s", "Min",
            "Percentile2.5th", "Percentile5th", "Median", "Percentile95th",
            "Percentile97.5th", "Max"
          ),
          names(.data)
        ),
        .fns = ~if_else(.col %in% names(.data), as.numeric(.x), .x)
      ),
      #make LBTEST values large
      LBTEST = str_trim(toupper(LBTEST)),
      #replace protected white spaces with normal ones
      LBTEST = gsub("\u00A0", " ", LBTEST, perl = TRUE),
      #*turn units into character (bugfix for LBORRESU being treated as double
      #*when completely empty)
      LBORRESU = as.character(LBORRESU)
    )
  
  cleaned_parameters <- params_stripped %>%
    left_join(data_mapping, by = c("LBTEST", "LBSPEC", "LBORRESU")) %>%
    select(-LBTEST, -LBSPEC, -LBORRESU, -LBTESTCD_old, -meanval) %>%
    rename(
      "LBTEST" = "LBTEST_new",
      "LBORRESU" = "LBORRESU_new",
      "LBTESTCD" = "LBTESTCD_new",
      "LBSPEC" = "LBSPEC_new"
    )
  
  return(cleaned_parameters)
}