#*12_effect_size_interpretation.R
#*
#*This is a standalone script which means it is not executed by 00_master.R
#*
#*This script interprets the results of Cliff's Delta effect sizes. The
#*following questions I wanted to answer:
#*1) Is the CCG between 5 % and 95 % percentile of VCG range? Y/N,
#*2) If so, does the majority vote for the same thing? Y/N
#*3) How big is the smallest effect leading to a significant difference?
#*4) How big is the biggest effect leading to a non significant difference?
#*
#*For this script, the effect size results and the voting results of the VCG
#*script are used.
#*
#*Input:  Effect sizes and majority votes of 04_resampling.R
#*Steps:  - Extract effect sizes per parameter.
#*        - Check if the CCG is within the 95 % percentile of the VCG effects.
#*        - Check whether the majority of votes voted for the same significance.
#*        - Calculate closest distance leading to a significant result.
#*        - Calculated largest distance leading to a non-significant result.
#*        - The last two points are done only if applicable.
#*Output: Table with interpreted results.
#*******************************************************************************
#*******************************************************************************
#*beginning of function----
#*******************************************************************************
effect_interpreter <- function(legacy_results, resampling_results){
  #*****************************************************************************
  #*drop BG results for now. And remove control groups and values where
  #*effsize is NA
  #*****************************************************************************
  legacy_results[[2]] <- legacy_results[[2]] %>%
    filter(
      !grepl("BG_", LBTESTCD),
      trial_set_description != "CG",
      !is.na(effsize)
      )
  
  resampling_results[[2]] <- resampling_results[[2]] %>%
    filter(
      !grepl("BG_", LBTESTCD),
      trial_set_description != "CG",
      !is.na(effsize)
    )
  #*****************************************************************************
  #*turn everything into a list split by doses and sex----
  #*****************************************************************************
  legacy_effect_list <- split(
    legacy_results[[2]],
    paste(
      legacy_results[[2]]$trial_set_description,
      legacy_results[[2]]$SEX,
      sep = "_"
      )
  )
  
  resampling_effect_list <- split(
    resampling_results[[2]],
    paste(
      resampling_results[[2]] %>%
        pull(trial_set_description),
      resampling_results[[2]] %>%
        pull(SEX),
      sep = "_"
    )
  )
  #*****************************************************************************
  #*drop parameters which are not within resampling list----
  #*****************************************************************************
  legacy_effect_list_filtered <- mapply(function(x,y){
    #make parameter ID (PID)
    x <- x %>%
      mutate(PID = paste(LBTESTCD, LBSPEC, LBORRESU, sep = "."))
    #get vector with PID from resampling
    PID_vector <- y %>%
      mutate(PID = paste(LBTESTCD, LBSPEC, LBORRESU, sep = ".")) %>%
      pull(PID)
    #filter by this vector
    x <- x %>% filter(PID %in% PID_vector) %>% select(-PID)
    return(x)
  },
  x = legacy_effect_list, y = resampling_effect_list, SIMPLIFY = FALSE
  )
  #*****************************************************************************
  #*get information as described in the script description----
  #*****************************************************************************
  tests <- mapply(function(x,y){
    #get information from VCG iterations
    y_summaries <- y %>%
      group_by(LBTESTCD, LBSPEC, LBORRESU) %>%
      summarize(
        #get 5 % percentile
        minval = quantile(effsize, .05),
        #get 95 % percentile
        maxval = quantile(effsize, .95),
        #get mode (i.e. most frequent significance)
        voted_sig = {
          sig <- sum(significance == "1", na.rm = TRUE)
          nonsig <- sum(significance == "0", na.rm = TRUE)
          ifelse(sig >= nonsig, "1", "0")
        },
        #get lowest effect leading to significance (if applicable)
        lowest_sig = ifelse(
          any(significance == 1),
          effsize[which.min(abs(effsize[significance == 1]))],
          NA
        ),
        largest_nonsig = ifelse(
          any(significance == 0),
          effsize[which.max(abs(effsize[significance == 0]))],
          NA
        ),
        .groups = "drop"
      )
    
    #combine the results to original CCG res
    interpretations <- x %>%
      left_join(
        y_summaries, by = c("LBTESTCD", "LBSPEC", "LBORRESU")
      ) %>%
      mutate(
        #is the CCG effect within 95 % percentile of VCG effects?
        CCG_eff_between_VCG = between(effsize, minval, maxval),
        #is the significance consistent to majority vote of VCGs?
        consistent = if_else(significance == voted_sig, TRUE, FALSE)
      )
    return(interpretations)
  },
  x = legacy_effect_list_filtered, y = resampling_effect_list, SIMPLIFY = FALSE
  )
  #*****************************************************************************
  #*flatten the results and return from function----
  #*****************************************************************************
  tests_flatten <- tests %>%
    imap_dfr(~.x %>% as_tibble())
  
  return(tests_flatten)
} #end of function
#*******************************************************************************
#*******************************************************************************
#*apply functions to results----
#*******************************************************************************
#legacy study A
study_A_effects <- effect_interpreter(CCG_legacy_study_A, VCG_legacy_study_A)
#summarize results
study_A_effects_sum <- study_A_effects %>%
  group_by(LBTESTCD, LBSPEC, LBORRESU) %>%
  summarize(
    #per parameter: are there any CCG effect sizes outside of 95 percentiles?
    any_eff_outside = sum(!CCG_eff_between_VCG),
    #per parameter: are there any discrepancies?
    any_discrepancies = sum(!consistent),
    .groups = "drop"
  )

#calculate percentages of parameters with at least one CCG outside of VCG range
study_A_effects_percentage_CCG_within_VCG <- round(
  study_A_effects %>% filter(CCG_eff_between_VCG) %>% nrow() /
  nrow(study_A_effects) * 100
  )
#how many are there that are still consistent with original result?
study_A_effects_percentage_outside_but_where_still_con <- round(
  study_A_effects %>% filter(!CCG_eff_between_VCG & consistent) %>% nrow() /
    study_A_effects %>% filter(!CCG_eff_between_VCG) %>% nrow() * 100
)
#write results as XLSX
write.xlsx(study_A_effects, "legacy_study_A_effect_table.xlsx")
#*******************************************************************************
#legacy study B
study_B_effects <- effect_interpreter(CCG_legacy_study_B, VCG_legacy_study_B)
#summarize results
study_B_effects_sum <- study_B_effects %>%
  group_by(LBTESTCD, LBSPEC, LBORRESU) %>%
  summarize(
    #per parameter: are there any CCG effect sizes outside of 95 percentiles?
    any_eff_outside = sum(!CCG_eff_between_VCG),
    #per parameter: are there any discrepancies?
    any_discrepancies = sum(!consistent),
    .groups = "drop"
  )

#calculate percentages of parameters with at least one CCG outside of VCG range
study_B_effects_percentage_CCG_within_VCG <- round(
  study_B_effects %>% filter(CCG_eff_between_VCG) %>% nrow() /
    nrow(study_B_effects) * 100
)
#how many are there that are still consistent with original result?
study_B_effects_percentage_outside_but_where_still_con <- round(
  study_B_effects %>% filter(!CCG_eff_between_VCG & consistent) %>% nrow() /
    study_B_effects %>% filter(!CCG_eff_between_VCG) %>% nrow() * 100
)
#write results as XLSX
write.xlsx(study_B_effects, "legacy_study_B_effect_table.xlsx")
#*******************************************************************************
#legacy study B_recovery
study_B_recovery_effects <- effect_interpreter(CCG_legacy_study_B_recovery, VCG_legacy_study_B_recovery)
#summarize results
study_B_recovery_effects_sum <- study_B_recovery_effects %>%
  group_by(LBTESTCD, LBSPEC, LBORRESU) %>%
  summarize(
    #per parameter: are there any CCG effect sizes outside of 95 percentiles?
    any_eff_outside = sum(!CCG_eff_between_VCG),
    #per parameter: are there any discrepancies?
    any_discrepancies = sum(!consistent),
    .groups = "drop"
  )

#calculate percentages of parameters with at least one CCG outside of VCG range
study_B_recovery_effects_percentage_CCG_within_VCG <- round(
  study_B_recovery_effects %>% filter(CCG_eff_between_VCG) %>% nrow() /
    nrow(study_B_recovery_effects) * 100
)
#how many are there that are still consistent with original result?
study_B_recovery_effects_percentage_outside_but_where_still_con <- round(
  study_B_recovery_effects %>% filter(!CCG_eff_between_VCG & consistent) %>% nrow() /
    study_B_recovery_effects %>% filter(!CCG_eff_between_VCG) %>% nrow() * 100
)
#write results as XLSX
write.xlsx(study_B_recovery_effects, "legacy_study_B_recovery_effect_table.xlsx")
#*******************************************************************************
#legacy study C
study_C_effects <- effect_interpreter(CCG_legacy_study_C, VCG_legacy_study_C)
#summarize results
study_C_effects_sum <- study_C_effects %>%
  group_by(LBTESTCD, LBSPEC, LBORRESU) %>%
  summarize(
    #per parameter: are there any CCG effect sizes outside of 95 percentiles?
    any_eff_outside = sum(!CCG_eff_between_VCG),
    #per parameter: are there any discrepancies?
    any_discrepancies = sum(!consistent),
    .groups = "drop"
  )

#calculate percentages of parameters with at least one CCG outside of VCG range
study_C_effects_percentage_CCG_within_VCG <- round(
  study_C_effects %>% filter(CCG_eff_between_VCG) %>% nrow() /
    nrow(study_C_effects) * 100
)
#how many are there that are still consistent with original result?
study_C_effects_percentage_outside_but_where_still_con <- round(
  study_C_effects %>% filter(!CCG_eff_between_VCG & consistent) %>% nrow() /
    study_C_effects %>% filter(!CCG_eff_between_VCG) %>% nrow() * 100
)
#write results as XLSX
write.xlsx(study_C_effects, "legacy_study_C_effect_table.xlsx")
