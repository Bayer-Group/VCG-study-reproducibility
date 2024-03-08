#04a1_simulate_missing_parameters.R

#*This function aims to simulate parameters if no sufficient data is present
#*in the VCG data set. This is done in the following fashion:
#*
#*Input:  historical control data (HCD) and if present, sentinel animals.
#*Steps:  - imute missing data by mean value of non-missing data
#*        - impute missing data by randomly sampling non-missing data
#*        - impute missing data by predictive mean matching (pmm)
#*        
#*Output: Imputed data


#*******************************************************************************
#Check if packages are installed------------------------------------------------
#*******************************************************************************
#*******************************************************************************
require(tidyverse)
require(mice)
require(janitor)
require(ranger)
#*******************************************************************************

#*Beginning of function---------------------------------------------------------
#*******************************************************************************
#*******************************************************************************
missing_data_imputation <- function(
    imputation_method = "median",
    studydata = HCD_selected_LBTESTCDs
)
{
  #*Split studydata by sex to avoid simulating testes in females
  studydata_by_sex <- split(studydata, studydata$SEX)
  #*pivot the LBTESTCD of the studydata to have each measured parameter as a
  #*separate column
  studydata_pivoted <- lapply(
    studydata_by_sex,
    function(x){
      x %>% 
        select(USUBJID, LBTESTCD, LBORRES, LBSPEC, LBORRESU) %>%
        pivot_wider(
        id_cols = USUBJID,
        names_from = c(LBTESTCD, LBSPEC, LBORRESU),
        values_from = LBORRES
      )
    }
  )
  
  #clean column names
  janitorized_colnames <- lapply(
    studydata_pivoted,
    function(x){
      names(x) <- clean_names(x)
    }
  )
  #*****************************************************************************
  #*****************************************************************************
  #*Simple imputation methods----
  #*****************************************************************************
  #*Here, simple imputation methods are performed. I.e., the median of all
  #*non-missing variables is taken or random samples from non-missing data is
  #*used to fill the gaps
  #*****************************************************************************
  #*create a function to calculate the median and fill the gaps
  median_filling <- function(){
    median_values <- lapply(janitorized_colnames, function(x) {
      sapply(x[-1], function(numcol) {
        if (is.numeric(numcol)) {
          median(numcol, na.rm = TRUE)
        } else {
          NA
        }
      })
    })
    
    #fill in missing values with median of non-missing values
    studydata_filled <- Map(
      function(studydata, medians){
      studydata[-1] <- Map(
        function(col, median){
        if(is.numeric(col)){
          ifelse(is.na(col), median, col)
          }else{
          col
            }
          },
        studydata[-1],
        medians
        )
      studydata
    },
    janitorized_colnames,
    median_values
    )
    
    return(studydata_filled)
  }
  #*****************************************************************************
  #*create a function to fill gaps with random values from non-missing fields
  random_filling <- function(){
    #collect all non-missing values
    random_values <- lapply(janitorized_colnames, function(x) {
      sapply(x[-1], function(numcol) {
        if (is.numeric(numcol)) {
          numcol[!is.na(numcol)]
        } else {
          NA
        }
      })
    })
    
    #fill in missing values with median of non-missing values
    studydata_filled <- Map(
      function(studydata, randoms){
        studydata[-1] <- Map(
          function(col, random){
            if(is.numeric(col)){
              ifelse(is.na(col), sample(random, 1), col)
            }else{
              col
            }
          },
          studydata[-1],
          randoms
        )
        studydata
      },
      janitorized_colnames,
      random_values
    )
    return(studydata_filled)
  }
  
  
  #*****************************************************************************
  #*****************************************************************************
  #*Advanced imputation methods----
  #*****************************************************************************
  #*Predictive mean matching
  pmm_filling <- function(){
    #create a model to impute missing data
    studydata_imputation_model <- lapply(
      janitorized_colnames,
      function(x){
        mice(
          x,
          m = 5,
          method = "pmm",
          seed = 1
        )
      }
    )
    
    #fill in missing values
    studydata_filled <- lapply(
      studydata_imputation_model,
      function(x){
        complete(x, 1)
      }
    )
  }
  
  #*****************************************************************************
  #*Bayesian linear regression
  lr_filling <- function(){
    #create a model to impute missing data
    studydata_imputation_model <- lapply(
      janitorized_colnames,
      function(x){
        mice(
          x,
          m = 5,
          method = "norm",
          seed = 1
        )
      }
    )
    
    #fill in missing values
    studydata_filled <- lapply(
      studydata_imputation_model,
      function(x){
        complete(x, 1)
      }
    )
  }
  #*****************************************************************************
  #*****************************************************************************
  #*Apply function based on selection----
  #*****************************************************************************
  studydata_filled <- if(imputation_method == "median"){
    median_filling()
  }else if(imputation_method == "random_sampling"){
    random_filling()
  }else if(imputation_method == "pmm"){
    pmm_filling()
  }else if(imputation_method == "lr"){
    lr_filling()
  }
  
  #*****************************************************************************
  #*****************************************************************************
  #*get the old column names back so that you can translate the data frames
  #*back to the old structure----
  #*****************************************************************************
  studydata_filled_with_old_names <- Map(
    setNames,
    studydata_filled,
    lapply(studydata_pivoted, names)
    )
  
  #transform the new study data back to the old structure
  studydata_filled_old_structure <- lapply(
    studydata_filled_with_old_names,
    function(x){
      x %>%
      pivot_longer(
        cols = c(everything(), -USUBJID),
        names_to = "LBTESTCD",
        values_to = "LBORRES"
      ) %>%
      separate(
        col = "LBTESTCD",
        into = c("LBTESTCD", "LBSPEC", "LBORRESU"),
        sep = "_"
      )
    }
  ) %>%
    #*gives a warning message about body weights not being able to separate
    #*into three columns. This message is suppressed. The bug is fixed in the
    #*next snippet.
    suppressWarnings()
  
  #*****************************************************************************
  #unlist results of males and females
  studydata_filled_flatten <- studydata_filled_old_structure %>%
    imap_dfr(~ .x %>% as_tibble(), .id = "SEX") %>%
    #clean body weight, organ weight, and FW values
    mutate(
      LBTESTCD = case_when(
        LBTESTCD == "BW" ~ paste(LBTESTCD, LBSPEC, sep = "_"),
        grepl("WEIGHT", LBTESTCD) ~ paste(LBTESTCD, LBSPEC, sep = "_"),
        grepl("OWBW", LBTESTCD) ~ paste(LBTESTCD, LBSPEC, sep = "_"),
        grepl("FC", LBTESTCD) ~ paste(LBTESTCD, LBSPEC, sep = "_"),
        grepl("WC", LBTESTCD) ~ paste(LBTESTCD, LBSPEC, sep = "_"),
        TRUE ~ LBTESTCD
      ),
      LBSPEC = case_when(
        grepl("BW_D", LBTESTCD) ~ "BW",
        grepl("FC_W", LBTESTCD) ~ "FW",
        grepl("WC_W", LBTESTCD) ~ "FW",
        TRUE ~ LBSPEC
      ),
      #*manual correction of "glands" and "vesicles" into "gland" and "vesicle"
      LBSPEC = gsub("GLANDS", "GLAND", LBSPEC),
      LBSPEC = gsub("VESICLES", "VESICLE", LBSPEC),
      LBORRESU = case_when(
        grepl("BW_D", LBTESTCD) ~ "g",
        grepl("WEIGHT", LBTESTCD) ~ "g",
        grepl("OWBW", LBTESTCD) ~ "%",
        grepl("FC_W", LBTESTCD) ~ "g/day",
        grepl("WC_W", LBTESTCD) ~ "g/day",
        TRUE ~ LBORRESU
      )
    )
  
  #*****************************************************************************
  #join the new values to those of the old studydata
  new_studydata <- studydata_filled_flatten %>%
    left_join(
      studydata %>% select(-LBORRES),
      by = c("USUBJID", "SEX", "LBSPEC", "LBORRESU", "LBTESTCD")
    ) %>%
    mutate(trial_set_description = "CG")
  
  return(new_studydata)
  }
