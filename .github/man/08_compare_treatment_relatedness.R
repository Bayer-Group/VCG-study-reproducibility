#*08_compare_test_substance_relatedness.R
#*
#*This is a standalone script which means it is not executed by 00_master.R
#*
#*This script collects the XLSX documents created in 05_visualize_as_table.R
#*and summarizes the results of them.
#*The XLSX files serve to support a study director in decision making process.
#*Basically, the study directors can use the XLSX file to classify whether a
#*significant finding is treatment related or not and - if not - provide reasons
#*for their decisions.
#*In this script, we summarize the decisions and the discrepancies between CCG
#*results and VCG results.
#*
#*Input:  XLSX files from legacy studies with VCGs, filled out by subject matter
#*        experts (SMEs).
#*Steps:  - call XLSX files.
#*        - compare the results of treatment-relatedness classification
#*        - summarize the results and the reasons for classification
#*Output: Summary of treatment related results which is compared between CCG
#*        and VCGs.
#*******************************************************************************
#*******************************************************************************
#*load libraries----
#*******************************************************************************
require(tidyverse)
require(openxlsx)
#*******************************************************************************
#*define a function which summarizes:
#*1) the total number of parameters
#*2) the number of parameters where at least one significant difference is found
#*3) the number of parameters where a significant finding was classified as not
#*   test-substance related
#*4) how often each argument was used for this decision

decision_summary <- function(studydata){
  #total number of parameters
  num_parameter <- tibble(value = "num parameter", n = studydata %>% nrow())
  
  #number of parameters with at least one statistically significant result
  sig_parameter <- studydata %>%
    filter(apply(., 1, function(row) any(grepl("\\*", row))))
  
  num_sig_parameter <- tibble(
    value = "num sig parameter", n = sig_parameter %>% nrow()
    )
  
  #number of parameters which were classified as not-test-substance related
  not_test_substance_related <- sig_parameter %>%
    filter(!test_substance_related)
  
  num_not_test_substance_related <- tibble(
    value = "not test-substance related", n = not_test_substance_related %>%
      nrow()
  )
  
  #how often was which decision made?
  number_reasons <- not_test_substance_related %>%
    mutate(across(everything(), as.character)) %>%
    pivot_longer(
      cols = starts_with("reason_"),
      names_to = "reason",
      values_to = "value"
    ) %>% 
    filter(!is.na(value)) %>%
    count(value) %>%
    arrange(desc(n))
  
  summarized_result <- rbind(
    num_parameter,
    num_sig_parameter,
    num_not_test_substance_related,
    number_reasons
    ) %>%
    mutate(
      percentage = round(n / num_not_test_substance_related$n * 100),
      percentage = if_else(percentage > 100, NaN, percentage)
      ) %>%
    rbind(
      tibble(
        value = "percentage for non-relevace decision",
        n = 0,
        percentage = round(num_not_test_substance_related$n / num_sig_parameter$n  * 100))
    )
}

#*******************************************************************************
#*read study test-substance related results (trr)----
#*******************************************************************************
#*Legacy study A----
#*******************************************************************************
study_A_CCG_trr <- read.xlsx(
  paste0(rootpath, "/data/derived/legacy_study_A_CCG_results_table.xlsx")
)


study_A_trr <- read.xlsx(
  paste0(rootpath, "/data/derived/legacy_study_A_VCG_results_table.xlsx")
  )
#*******************************************************************************
#*extract the ones where there were discrepancies in the statistics and check
#*whether the conclusion of the study director would lead to the same
#*outcome
study_A_discrepancies_extracted <- study_A_trr %>%
  #remove the ones where there wasn't any statistical discrepancy
  filter(!consistent)
#*******************************************************************************
#check the inconsistently significant onec
study_A_incon_sig <- study_A_discrepancies_extracted %>%
  filter(
    grepl("CCG: (?!\\*)", where_was_the_discrepancy, perl = TRUE)
  )

#*check the ones where the findings were not test_substance_related. These are
#*the trouble makers
study_A_incon_sig_TR <- study_A_incon_sig %>%
  filter(test_substance_related)

#*get all treatment related noteworthy findings, irrelevant whether they were
#*consistent with the original study or not
study_A_sig_TR <- study_A_trr %>% filter(test_substance_related)
#*******************************************************************************
#check the inconsistently non-significant ones
study_A_incon_non_sig <- study_A_discrepancies_extracted %>%
  filter(
    grepl("VCG: (?!\\*)", where_was_the_discrepancy, perl = TRUE)
  )

#*check the ones which are now treatment related in comparison to the original
#*study. These may alter the toxicological outcome
study_A_incon_non_sig_not_test_substance_related <- study_A_incon_non_sig %>%
  filter(test_substance_related_in_original_study)
#*******************************************************************************
#*******************************************************************************
#*Legacy study B----
#*******************************************************************************
study_B_CCG_trr <- read.xlsx(
  paste0(rootpath, "/data/derived/legacy_study_B_CCG_results_table.xlsx")
)

study_B_trr <- read.xlsx(
  paste0(rootpath, "/data/derived/legacy_study_B_VCG_results_table.xlsx")
)
#*******************************************************************************
#*extract the ones where there were discrepancies in the statistics and check
#*whether the conclusion of the study director would lead to the same
#*outcome
study_B_discrepancies_extracted <- study_B_trr %>%
  #remove the ones where there wasn't any statistical discrepancy
  filter(!consistent)
#*******************************************************************************
#check the inconsistently significant onec
study_B_incon_sig <- study_B_discrepancies_extracted %>%
  filter(
    grepl("CCG: (?!\\*)", where_was_the_discrepancy, perl = TRUE)
  )

#*check the ones where the findings were not test_substance_related. These are
#*the trouble makers
study_B_incon_sig_TR <- study_B_incon_sig %>%
  filter(test_substance_related)

#*get all treatment related noteworthy findings, irrelevant whether they were
#*consistent with the original study or not
study_B_sig_TR <- study_B_trr %>% filter(test_substance_related)
#*******************************************************************************
#check the inconsistently non-significant ones
study_B_incon_non_sig <- study_B_discrepancies_extracted %>%
  filter(
    grepl("VCG: (?!\\*)", where_was_the_discrepancy, perl = TRUE)
  )

#*check the ones which are now treatment related in comparison to the original
#*study. These may alter the toxicological outcome
study_B_incon_non_sig_not_test_substance_related <- study_B_incon_non_sig %>%
  filter(test_substance_related_in_original_study)
#*******************************************************************************
#*******************************************************************************
#*Legacy study C----
#*******************************************************************************
study_C_CCG_trr <- read.xlsx(
  paste0(rootpath, "/data/derived/legacy_study_C_CCG_results_table.xlsx")
)

study_C_trr <- read.xlsx(
  paste0(rootpath, "/data/derived/legacy_study_C_VCG_results_table.xlsx")
)
#*******************************************************************************
#*extract the ones where there were discrepancies in the statistics and check
#*whether the conclusion of the study director would lead to the same
#*outcome
study_C_discrepancies_extracted <- study_C_trr %>%
  #remove the ones where there wasn't any statistical discrepancy
  filter(!consistent)
#*******************************************************************************
#check the inconsistently significant onec
study_C_incon_sig <- study_C_discrepancies_extracted %>%
  filter(
    grepl("CCG: (?!\\*)", where_was_the_discrepancy, perl = TRUE)
  )

#*check the ones where the findings were not test_substance_related.
#*These are the trouble makers
study_C_incon_sig_TR <- study_C_incon_sig %>%
  filter(test_substance_related)
#*******************************************************************************
#check the inconsistently non-significant ones
study_C_incon_non_sig <- study_C_discrepancies_extracted %>%
  filter(
    grepl("VCG: (?!\\*)", where_was_the_discrepancy, perl = TRUE)
  )

#*check the ones which are now treatment related in comparison to the original
#*study. These may alter the toxicological outcome
study_C_incon_non_sig_not_test_substance_related <- study_C_incon_non_sig %>%
  filter(test_substance_related_in_original_study)

#*get all treatment related noteworthy findings, irrelevant whether they were
#*consistent with the original study or not
study_C_sig_TR <- study_C_trr %>% filter(test_substance_related)


study_C_discrepancies_extracted %>%
  filter(!test_substance_related_in_original_study,
         grepl("CCG: (?!\\*)", where_was_the_discrepancy, perl = TRUE)
         ) %>% nrow()
#*******************************************************************************
#*******************************************************************************
#*get a summary of how often things were classified as not test-substance
#*related----
#*******************************************************************************
study_A_CCG_decision_summary <- decision_summary(study_A_CCG_trr)
study_A_VCG_decision_summary <- decision_summary(study_A_trr)

study_B_CCG_decision_summary <- decision_summary(study_B_CCG_trr)
study_B_VCG_decision_summary <- decision_summary(study_B_trr)

study_C_CCG_decision_summary <- decision_summary(study_C_CCG_trr)
study_C_VCG_decision_summary <- decision_summary(study_C_trr)

mean(c(
  study_A_CCG_decision_summary %>% filter(value == "percentage for non-relevace decision") %>% pull(percentage),
  study_B_CCG_decision_summary %>% filter(value == "percentage for non-relevace decision") %>% pull(percentage),
  study_C_CCG_decision_summary %>% filter(value == "percentage for non-relevace decision") %>% pull(percentage)
))

mean(c(
  study_A_VCG_decision_summary %>% filter(value == "percentage for non-relevace decision") %>% pull(percentage),
  study_B_VCG_decision_summary %>% filter(value == "percentage for non-relevace decision") %>% pull(percentage),
  study_C_VCG_decision_summary %>% filter(value == "percentage for non-relevace decision") %>% pull(percentage)
))

n_all_non_tr_findings <- rbind(
  study_A_CCG_decision_summary, study_A_VCG_decision_summary,
  study_B_CCG_decision_summary, study_B_VCG_decision_summary, 
  study_C_CCG_decision_summary, study_C_VCG_decision_summary
  ) %>%
  group_by(value) %>%
  summarize(summed = sum(n)) %>%
  filter(value == "not test-substance related") %>%
  pull(summed)

rbind(
  study_A_CCG_decision_summary, study_A_VCG_decision_summary,
  study_B_CCG_decision_summary, study_B_VCG_decision_summary, 
  study_C_CCG_decision_summary, study_C_VCG_decision_summary
) %>%
  group_by(value) %>%
  summarize(
    mean_percentage = round(sum(n) / n_all_non_tr_findings * 100)
    ) %>%
  arrange(desc(mean_percentage))


study_A_trr %>% mutate(across(everything(), as.character)) %>%
  bind_rows(study_B_trr %>% mutate(across(everything(), as.character))) %>%
  bind_rows(study_C_trr %>% mutate(across(everything(), as.character))) %>%
  bind_rows(study_A_CCG_trr %>% mutate(across(everything(), as.character))) %>%
  bind_rows(study_B_CCG_trr %>% mutate(across(everything(), as.character))) %>%
  bind_rows(study_C_CCG_trr %>% mutate(across(everything(), as.character))) %>%
  filter(
    reason_1 == "different deviations in different sexes" |
      reason_2 == "different deviations in different sexes" |
      reason_3 == "different deviations in different sexes" | 
      reason_4 == "different deviations in different sexes" |
      reason_5 == "different deviations in different sexes"
  ) %>% View()
