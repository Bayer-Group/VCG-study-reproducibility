#*13_power_interpretations
#*
#*This is a standalone script which means it is not executed by 00_master.R
#*
#*This script extracts the power calculations of the data and checks which of
#*the parameter values were actually well powered.
#*
#*Input:  Statistical results of CCG and resampling VCGs and their power.
#*Steps:  - Get the power.
#*        - Group by parameters and dose group and sex.
#*        - Highlight either those parameters which were well powered or under
#*          powered.
#*Output: - List of parameters which were powered or under powered.
#*******************************************************************************
#*******************************************************************************
#*beginning of function----
#*******************************************************************************
power_interpreter <- function(legacy_results, resampling_results){
  #*****************************************************************************
  #*drop BG and BW results for now. And remove control groups and values where
  #*effsize is NA
  #*****************************************************************************
  legacy_results[[3]] <- legacy_results[[3]] %>%
    filter(
      !grepl("^BG_|^BW_", LBTESTCD),
      !is.na(power)
      )
  
  legacy <- legacy_results[[3]]
  
  resampling_results[[3]] <- resampling_results[[3]] %>%
    filter(
      !grepl("^BG_|^BW_", LBTESTCD),
      !is.na(power)
    )
  
  resampling <- resampling_results[[3]]
  #*****************************************************************************
  #*drop parameters which are not within resampling list----
  #*****************************************************************************
  #make parameter ID (PID)
  legacy <- legacy %>%
    mutate(PID = paste(LBTESTCD, LBSPEC, LBORRESU, sep = "."))
  #get vector with PID from resampling
  PID_vector <- resampling %>%
    mutate(PID = paste(LBTESTCD, LBSPEC, LBORRESU, sep = ".")) %>%
    pull(PID)
  #filter by this vector
  legacy_filtered <- legacy %>%
    filter(PID %in% PID_vector) %>% select(-PID)
  #*****************************************************************************
  #*get information as described in the script description----
  #*****************************************************************************
  #are the results above 80 percent?
  above80 <- resampling %>%
    #add CCG power to the resampling table
    bind_rows(
      legacy_filtered %>% mutate(iteration = "0")
      ) %>%
    mutate(
      above80 = if_else(power >= .8, TRUE, FALSE)
    )
  
  #get the percentage of everything above 80
  absolute_percentage_above80 <- round(
    nrow(above80 %>% filter(above80)) / nrow(above80) * 100
    )
  
  #per parameter:are they enough powered?
  per_param <- above80 %>%
    group_by(LBTESTCD, LBSPEC, LBORRESU, SEX) %>%
    summarize(
      per_param_percentage = sum(above80),
      .groups = "drop"
    ) %>%
    mutate(
      params_enough_powered = if_else(per_param_percentage >= 95, TRUE, FALSE)
      )
  
  #per parameter: are there parameters enough powered in both sexes?
  per_param_both_sexes <- per_param %>%
    #if at least one sex is under powered, both aren't sufficient
    group_by(LBTESTCD, LBSPEC, LBORRESU) %>%
    summarize(
      params_enough_powered_both_sexes = sum(params_enough_powered),
      .groups = "drop"
      ) %>%
    mutate(
      params_enough_powered_both_sexes = if_else(
        params_enough_powered_both_sexes == 2, TRUE, FALSE
      )
    )
  
  #how many parameters are enough powered in both sexes?
  per_param_both_sexes_percentage <- round(
    nrow(per_param_both_sexes %>% filter(params_enough_powered_both_sexes)) /
      nrow(per_param_both_sexes) * 100
  )
  #*****************************************************************************
  #*summarize all results in a list and export----
  #*****************************************************************************
  results <- list(
    above80,
    absolute_percentage_above80,
    per_param,
    per_param_both_sexes,
    per_param_both_sexes_percentage
  )
  
  names(results) <- c(
    "above80",
    "absolute_percentage_above80",
    "per_param",
    "per_param_both_sexes",
    "per_param_both_sexes_percentage"
  )
  
  return(results)
} #end of function
#*******************************************************************************
#*******************************************************************************
#*apply functions to results----
#*******************************************************************************
#legacy study A
study_A_power <- power_interpreter(CCG_legacy_study_A, VCG_legacy_study_A)
A_params <- study_A_power[["per_param_both_sexes"]] %>%
  filter(params_enough_powered_both_sexes) %>% pull(LBTESTCD)
study_A_power[["per_param_both_sexes_percentage"]]
#*******************************************************************************
#legacy study B
study_B_power <- power_interpreter(CCG_legacy_study_B, VCG_legacy_study_B)
B_params <- study_B_power[["per_param_both_sexes"]] %>%
  filter(params_enough_powered_both_sexes) %>% pull(LBTESTCD)
study_B_power[["per_param_both_sexes_percentage"]]
#*******************************************************************************
#legacy study B_recovery
study_B_recovery_power <- power_interpreter(
  CCG_legacy_study_B_recovery, VCG_legacy_study_B_recovery
  )
B_rec_params <- study_B_recovery_power[["per_param_both_sexes"]] %>%
  filter(params_enough_powered_both_sexes) %>% pull(LBTESTCD)
study_B_recovery_power[["per_param_both_sexes_percentage"]]
#*******************************************************************************
#legacy study C
study_C_power <- power_interpreter(CCG_legacy_study_C, VCG_legacy_study_C)
C_params <- study_C_power[["per_param_both_sexes"]] %>%
  filter(params_enough_powered_both_sexes) %>% pull(LBTESTCD)
study_C_power[["per_param_both_sexes_percentage"]]
#*******************************************************************************
#Check overall parameters which are well powered across all studies
Reduce(intersect, list(A_params, B_params, C_params))
