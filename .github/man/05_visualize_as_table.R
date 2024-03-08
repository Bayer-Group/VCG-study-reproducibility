#*05_visualize_as_table.R

#*This script generates a gt table showing significant results of the
#*legacy study and the concatenated results of the VCG resampling experiment.
#*The goal is to decide whether treatment-relatedness occurs by observing not
#*only statistical significance but also other points which might speak for
#*treatment relatedness (such as dose dependency or whether CG is outisde of
#*5th/95th percentile of HCD).
#*Addittionally, the tables are exported as XLSX so that they can be edited by
#*subject matter experts and study directors.
##### program starts here
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*load libraries----
#*****************************************************************************
require(tidyverse)
require(gt)
require(webshot2)
require(openxlsx)
require(caret)
#*******************************************************************************
#*******************************************************************************
#Read RData for fast recovery of results. Uncomment if needed----
#*******************************************************************************
# CCG_legacy_study_A <- readRDS(
#   paste0(rootpath, "/data/Derived/CCG_legacy_study_A.RData")
# )
# VCG_legacy_study_A <- readRDS(
#   paste0(rootpath, "/data/Derived/VCG_legacy_study_A.RData")
# )
# VCG_legacy_study_A_imp_median <- readRDS(
#   paste0(rootpath, "/data/Derived/VCG_legacy_study_A_imp_median.RData")
# )
# VCG_legacy_study_A_imp_rs <- readRDS(
#   paste0(rootpath, "/data/Derived/VCG_legacy_study_A_imp_rs.RData")
# )
# VCG_legacy_study_A_imp_pmm <- readRDS(
#   paste0(rootpath, "/data/Derived/VCG_legacy_study_A_imp_pmm.RData")
# )
# 
# CCG_legacy_study_B <- readRDS(
#   paste0(rootpath, "/data/Derived/CCG_legacy_study_B.RData")
# )
# VCG_legacy_study_B <- readRDS(
#   paste0(rootpath, "/data/Derived/VCG_legacy_study_B.RData")
# )
# CCG_legacy_study_B_recovery <- readRDS(
#   paste0(rootpath, "/data/Derived/CCG_legacy_study_B_recovery.RData")
# )
# VCG_legacy_study_B_recovery <- readRDS(
#   paste0(rootpath, "/data/Derived/VCG_legacy_study_B_recovery.RData")
# )
# VCG_legacy_study_B_imp_median <- readRDS(
#   paste0(rootpath, "/data/Derived/VCG_legacy_study_B_imp_median.RData")
# )
# VCG_legacy_study_B_imp_rs <- readRDS(
#   paste0(rootpath, "/data/Derived/VCG_legacy_study_B_imp_rs.RData")
# )
# VCG_legacy_study_B_imp_pmm <- readRDS(
#   paste0(rootpath, "/data/Derived/VCG_legacy_study_B_imp_pmm.RData")
# )
# 
# CCG_legacy_study_C <- readRDS(
#   paste0(rootpath, "/data/Derived/CCG_legacy_study_C.RData")
# )
# VCG_legacy_study_C <- readRDS(
#   paste0(rootpath, "/data/Derived/VCG_legacy_study_C.RData")
# )
# VCG_legacy_study_C_imp_median <- readRDS(
#   paste0(rootpath, "/data/Derived/VCG_legacy_study_C_imp_median.RData")
# )
# VCG_legacy_study_C_imp_rs <- readRDS(
#   paste0(rootpath, "/data/Derived/VCG_legacy_study_C_imp_rs.RData")
# )
# VCG_legacy_study_C_imp_pmm <- readRDS(
#   paste0(rootpath, "/data/Derived/VCG_legacy_study_C_imp_pmm.RData")
# )

#beginning of function
#*******************************************************************************
#*******************************************************************************
study_results_table <- function(
  legacy_study_results = CCG_legacy_study_C,
  resampling_results = VCG_legacy_study_C
)
  {
  #*****************************************************************************
  #*Extract the numeric values form the legacy study result and the resampling
  #*result
  resampling_results_flatten <- resampling_results[[1]] %>%
    #create a column storing the LBTESTCD, LBSPEC, and unit of measurement
    mutate(
      code_unit_spec = paste0(
        LBTESTCD, " [", LBORRESU, "] in ", LBSPEC
      )
    ) %>%
    #*remove body weight gain parameters for now as they were only experimental
    filter(!grepl("BG_D", LBTESTCD))
  
  legacy_study_results_flatten <- legacy_study_results[[1]] %>%
    #create a column storing the LBTESTCD, LBSPEC, and unit of measurement
    mutate(
      code_unit_spec = paste0(
        LBTESTCD, " [", LBORRESU, "] in ", LBSPEC
      )
    ) %>%
    select(names(resampling_results_flatten))
  #*****************************************************************************
  #*****************************************************************************
  #*Compute reproducibility percentage----
  #*****************************************************************************
  #*calculate in how many cases the statistical result of the CCG study is in
  #*accordance with the statistical result of the VCG
  #*****************************************************************************
  #Get results of legacy study with concurrent control
  CCG_stat_res <- legacy_study_results_flatten %>%
    #drop parameters which were omitted by VCGs
    inner_join(
      resampling_results_flatten %>% select(LBTESTCD, LBSPEC, LBORRESU),
      by = c("LBTESTCD", "LBSPEC", "LBORRESU")
      ) %>%
    select(
      LBTESTCD, LBSPEC, LBORRESU, starts_with(c("CG_", "LD_", "MD_", "HD_"))
      ) %>%
    mutate(across(
      starts_with(c("CG_", "LD_", "MD_", "HD_")),
      function(x){x = if_else(grepl("\\*", x), TRUE, FALSE)}
    ))
  #*****************************************************************************
  #*Further specify the direction of the significant difference: increase
  CCG_stat_res_pos <- legacy_study_results_flatten %>%
    #drop parameters which were omitted by VCGs
    inner_join(
      resampling_results_flatten %>% select(LBTESTCD, LBSPEC, LBORRESU),
      by = c("LBTESTCD", "LBSPEC", "LBORRESU")
    ) %>%
    select(
      LBTESTCD, LBSPEC, LBORRESU, starts_with(c("CG_", "LD_", "MD_", "HD_"))
    ) %>%
    mutate(across(
      starts_with(c("CG_", "LD_", "MD_", "HD_")),
      function(x){x = if_else(grepl("\\+", x), TRUE, FALSE)}
    )) %>%
    select(-LBTESTCD, -LBORRESU, -LBSPEC)
  #add a suffix to the increases
  names(CCG_stat_res_pos) <- paste0(names(CCG_stat_res_pos), "_pos")
  
  #*Further specify the direction of the significant difference: decrease
  CCG_stat_res_neg <- legacy_study_results_flatten %>%
    #drop parameters which were omitted by VCGs
    inner_join(
      resampling_results_flatten %>% select(LBTESTCD, LBSPEC, LBORRESU),
      by = c("LBTESTCD", "LBSPEC", "LBORRESU")
    ) %>%
    select(
      LBTESTCD, LBSPEC, LBORRESU, starts_with(c("CG_", "LD_", "MD_", "HD_"))
    ) %>%
    mutate(across(
      starts_with(c("CG_", "LD_", "MD_", "HD_")),
      function(x){x = if_else(grepl("\\-", x), TRUE, FALSE)}
    )) %>%
    select(-LBTESTCD, -LBORRESU, -LBSPEC)
  #add a suffix to the decreases
  names(CCG_stat_res_neg) <- paste0(names(CCG_stat_res_neg), "_neg")
  #*****************************************************************************
  #*****************************************************************************
  #*Do the same for VCGs
  #Get results of legacy study with virtual control instead
  VCG_stat_res <- resampling_results_flatten %>%
    ungroup() %>%
    #drop parameters which are not present in CCG results table
    inner_join(
      legacy_study_results_flatten %>% select(LBTESTCD, LBSPEC, LBORRESU),
      by = c("LBTESTCD", "LBSPEC", "LBORRESU")
    ) %>%
    select(
      LBTESTCD, LBSPEC, LBORRESU, starts_with(c("CG_", "LD_", "MD_", "HD_"))
    ) %>%
    mutate(across(
      starts_with(c("CG_", "LD_", "MD_", "HD_")),
      function(x){x = if_else(grepl("\\*", x), TRUE, FALSE)}
    ))
  #*****************************************************************************
  #*Further specify the direction of the significant difference: increase
  VCG_stat_res_pos <- resampling_results_flatten %>%
    ungroup() %>%
    #drop parameters which are not present in CCG results table
    inner_join(
      legacy_study_results_flatten %>% select(LBTESTCD, LBSPEC, LBORRESU),
      by = c("LBTESTCD", "LBSPEC", "LBORRESU")
    ) %>%
    select(
      LBTESTCD, LBSPEC, LBORRESU, starts_with(c("CG_", "LD_", "MD_", "HD_"))
    ) %>%
    mutate(across(
      starts_with(c("CG_", "LD_", "MD_", "HD_")),
      function(x){x = if_else(grepl("\\+", x), TRUE, FALSE)}
    )) %>%
    select(-LBTESTCD, -LBORRESU, -LBSPEC)
  #add a suffix to the increases
  names(VCG_stat_res_pos) <- paste0(names(VCG_stat_res_pos), "_pos")
  
  #*Further specify the direction of the significant difference: decrease
  VCG_stat_res_neg <- resampling_results_flatten %>%
    ungroup() %>%
    #drop parameters which are not present in CCG results table
    inner_join(
      legacy_study_results_flatten %>% select(LBTESTCD, LBSPEC, LBORRESU),
      by = c("LBTESTCD", "LBSPEC", "LBORRESU")
    ) %>%
    select(
      LBTESTCD, LBSPEC, LBORRESU, starts_with(c("CG_", "LD_", "MD_", "HD_"))
    ) %>%
    mutate(across(
      starts_with(c("CG_", "LD_", "MD_", "HD_")),
      function(x){x = if_else(grepl("\\-", x), TRUE, FALSE)}
    )) %>%
    select(-LBTESTCD, -LBORRESU, -LBSPEC)
  #add a suffix to the decreases
  names(VCG_stat_res_neg) <- paste0(names(VCG_stat_res_neg), "_neg")
  #*****************************************************************************
  #*****************************************************************************
  #Check how many times the significance was in accordance between CCG and VCG
  VCG_res_matching_with_CCG <- apply(
    cbind(
      CCG_stat_res,
      CCG_stat_res_pos,
      CCG_stat_res_neg,
      VCG_stat_res,
      VCG_stat_res_pos,
      VCG_stat_res_neg
      ), 1, function(x){
    identical(
      x[4:ncol(cbind(CCG_stat_res, CCG_stat_res_pos, CCG_stat_res_neg))],
      x[
        (ncol(cbind(CCG_stat_res, CCG_stat_res_pos, CCG_stat_res_neg)) + 4):
          (2 * ncol(cbind(CCG_stat_res, CCG_stat_res_pos, CCG_stat_res_neg)))
        ]
    )
  })
  #Calculate the percentage from it
  matching_percentage <- round(
    sum(VCG_res_matching_with_CCG) / length(VCG_res_matching_with_CCG) * 100
  )
  #*****************************************************************************
  #*Specify, where exactly the discrepancies occurred (in which columns)----
  #*****************************************************************************
  CCG_which_res <- legacy_study_results_flatten %>%
    #drop parameters which were omitted by VCGs
    inner_join(
      resampling_results_flatten %>% select(LBTESTCD, LBSPEC, LBORRESU),
      by = c("LBTESTCD", "LBSPEC", "LBORRESU")
    ) %>%
    select(
      LBTESTCD, LBSPEC, LBORRESU, starts_with(c("CG_", "LD_", "MD_", "HD_"))
    ) %>%
    mutate(across(
      starts_with(c("CG_", "LD_", "MD_", "HD_")),
      #extract all the "*" and the value in the brackets "(+)" or "(-)"
      ~sapply(
        str_extract_all(., "\\*+|\\(\\+\\)|\\(\\-\\)"), paste, collapse = ""
        ),
      .names = "{.col}"
    ))
  
  #Get results of legacy study with virtual control instead
  VCG_which_res <- resampling_results_flatten %>%
    ungroup() %>%
    #drop parameters which were omitted by CCGs
    inner_join(
      legacy_study_results_flatten %>% select(LBTESTCD, LBSPEC, LBORRESU),
      by = c("LBTESTCD", "LBSPEC", "LBORRESU")
    ) %>%
    select(
      LBTESTCD, LBSPEC, LBORRESU, starts_with(c("CG_", "LD_", "MD_", "HD_"))
    ) %>%
    mutate(across(
      starts_with(c("CG_", "LD_", "MD_", "HD_")),
      #extract all the "*" and the value in the brackets "(+)" or "(-)"
      ~sapply(
        str_extract_all(., "\\*+|\\(\\+\\)|\\(\\-\\)"), paste, collapse = ""
        ),
      .names = "{.col}"
    ))
  
  
  #bind columns
  bound_CCG_and_VCG <- cbind(CCG_which_res, VCG_which_res)
  #add a suffix to the CCG results and VCG results
  names(bound_CCG_and_VCG) <- c(
    paste0(names(CCG_stat_res), "_CCG"),
    paste0(names(VCG_stat_res), "_VCG")
  )
  
  discrepancy_check <- bound_CCG_and_VCG %>%
    select(-LBTESTCD_VCG, -LBSPEC_VCG, -LBORRESU_VCG) %>%
    rename(
      "LBTESTCD" = "LBTESTCD_CCG",
      "LBSPEC" = "LBSPEC_CCG",
      "LBORRESU" = "LBORRESU_CCG"
      )
  #Create a function to check for identical values
  discrepancy_check_function <- function(row) {
    #empty vector to store discrepancies
    discrepancy_string <- c()
    
    #loop over each row to find whether discrepancies occurred
    for (i in 4:(ncol(CCG_which_res))) {
      #*compare whether column 4 correspond to column 11, column 5 to column 12,
      #*and so on
      col_name <- str_remove_all(names(CCG_which_res)[i], "_CCG")
      col_val_CCG <- row[i]
      col_val_VCG <- row[i + ncol(CCG_which_res) - 3]
      #*If both columns are true or false, skip. Otherwise print a message
      if (
        #check if both are significant and go into positive direction
        (
          grepl("\\*", col_val_CCG) & grepl("\\*", col_val_VCG)
          &
          grepl("\\+", col_val_CCG) & grepl("\\+", col_val_VCG)
        ) | #or if both are significant and go into negative direction
        (
          grepl("\\*", col_val_CCG) & grepl("\\*", col_val_VCG)
          &
          grepl("\\-", col_val_CCG) & grepl("\\-", col_val_VCG)
        ) | #or if both are not significant
        (!grepl("\\*", col_val_CCG) & !grepl("\\*", col_val_VCG))
        ) {
        next
      } else {
        discrepancy_string <- c(
          discrepancy_string,
          paste(col_name, ": CCG:", col_val_CCG, ", VCG:", col_val_VCG))
      }
    }
    #return it in a format which can be directly executed in EXCEL
    values_for_XLSX <- paste0(
      "=VERKETTEN(\"",
      paste(
        discrepancy_string,
        collapse = "\";ZEICHEN(10);\""
        ),
      "\")")
    
    return(values_for_XLSX)
  }
  #Apply the function to each row
  discrepancy_strings <- apply(discrepancy_check, 1, discrepancy_check_function)
  
  discrepancy_check$where_was_the_discrepancy <- discrepancy_strings
  #*****************************************************************************
  #*****************************************************************************
  #*Summarize over all results which we could not reproduce to obtain something
  #*which resembles a confusion matrix
  #*****************************************************************************
  #Create a function to check for consistency between the two data types
  consistency_check <- function(CCG, VCG) {
    #Convert A and B into vectors of logical values
    CCG_sig <- grepl("\\*", CCG)
    VCG_sig <- grepl("\\*", VCG)
    
    #check direction in which the significant difference goes
    CCG_dir <- grepl("\\(\\+\\)", CCG)
    VCG_dir <- grepl("\\(\\+\\)", VCG)
    
    #Compute the error counts
    Con_sig <- sum(CCG_sig & VCG_sig)
    Con_non_sig <- sum(!CCG_sig & !VCG_sig)
    Incon_sig   <- sum(!CCG_sig & VCG_sig)
    Incon_non_sig  <- sum(CCG_sig & !VCG_sig)
    Inv_sig <- sum(CCG_sig & VCG_sig & (CCG_dir != VCG_dir))
    
    # Return the error counts as a data frame
    return(
      tibble(
        con_sig = Con_sig,
        con_non_sig = Con_non_sig,
        incon_sig = Incon_sig,
        incon_non_sig = Incon_non_sig,
        inv_sig = Inv_sig
        )
      )
  }
  
  #Apply consistency_check function from rowwise and then column wise
  consistency_counts <- do.call(
    rbind,
    mapply(
      consistency_check,
      CCG_which_res[4:ncol(CCG_which_res)] %>% select(-starts_with("CG_")),
      VCG_which_res[4:ncol(VCG_which_res)] %>% select(-starts_with("CG_")),
      SIMPLIFY = FALSE
      )
    ) %>%
    summarize(across(everything(), sum)) %>%
    mutate(
      n_stats_results = sum(
        con_sig, con_non_sig, incon_sig, incon_non_sig, inv_sig
        ),
      n_sig_stats_results = sum(con_sig, incon_non_sig, inv_sig),
      n_non_sig_stats_results = sum(con_non_sig, incon_sig),
      con_sig_performance = paste0(
        con_sig, " out of ", sum(con_sig, incon_non_sig, inv_sig),
        " (", round(con_sig / sum(con_sig, incon_non_sig, inv_sig) * 100), "%)"
      ),
      con_non_sig_performance = paste0(
        con_non_sig, " out of ", sum(con_non_sig, incon_sig),
        " (", round(con_non_sig / sum(con_non_sig, incon_sig) * 100), "%)"
      ),
      incon_sig_performance = paste0(
        incon_sig, " out of ", sum(con_non_sig, incon_sig),
        " (", round(incon_sig / sum(con_non_sig, incon_sig) * 100), "%)"
      ),
      incon_non_sig_performance = paste0(
        incon_non_sig, " out of ", sum(con_sig, incon_non_sig, inv_sig),
        " (", round(
          incon_non_sig / sum(con_sig, incon_non_sig, inv_sig) * 100
          ), "%)"
      ),
      inv_sig_performance = paste0(
        inv_sig, " out of ", sum(con_sig, incon_non_sig, inv_sig),
        " (", round(inv_sig / sum(con_sig, incon_non_sig, inv_sig) * 100), "%)"
      )
    )
  #*****************************************************************************
  #*****************************************************************************
  #add parameter-wise consistency to the table----
  #*****************************************************************************
  parameter_wise_consistency_check <- function(CCG, VCG) {
    #Convert A and B into vectors of logical values
    CCG_sig <- CCG %>%
      as_tibble() %>%
      mutate(across(everything(), function(x) {grepl("\\*", x)}))
    VCG_sig <- VCG %>%
      as_tibble() %>%
      mutate(across(everything(), function(x) {grepl("\\*", x)}))
    
    #check direction in which the significant difference goes
    CCG_dir <- CCG %>%
      as_tibble() %>%
      mutate(across(everything(), function(x) {grepl("\\(\\+\\)", x)}))
    VCG_dir <- VCG %>%
      as_tibble() %>%
      mutate(across(everything(), function(x) {grepl("\\(\\+\\)", x)}))
    
    #*Combine the rows together. This will be used for logic checking
    combined_sig <- CCG_sig %>%
      #add row id which is used to combine CCG and VCG
      mutate(row_id = row_number()) %>%
      pivot_longer(-row_id, names_to = "col", values_to = "CCG") %>%
      #join VCG results
      left_join(
        VCG_sig %>%
          mutate(row_id = row_number()) %>%
          pivot_longer(-row_id, names_to = "col", values_to = "VCG"),
        by = c("row_id", "col"))
    
    combined_dir <- CCG_dir %>%
      #add row id which is used to combine CCG and VCG
      mutate(row_id = row_number()) %>%
      pivot_longer(-row_id, names_to = "col", values_to = "CCG") %>%
      #join VCG results
      left_join(
        VCG_dir %>%
          mutate(row_id = row_number()) %>%
          pivot_longer(-row_id, names_to = "col", values_to = "VCG"),
        by = c("row_id", "col"))
    
    #Count the number of total parameters
    n_parameters <- VCG_sig %>% nrow()
    #Compute the error counts
    #*count number of parameters where at least one group was significantly
    #*different to the control
    orig_sig <- apply(CCG_sig, 1, any) %>% sum()
    #*count the number of parameters where none of the groups show significant
    #*differences to the original control
    orig_non_sig <- apply(!CCG_sig, 1, all) %>% sum()
    #*count number of parameters where CCGs and VCGs were consistent in all
    #*significance classes
    Con <- apply(CCG_sig == VCG_sig, 1, all) %>% sum()
    #*Count numbers of consistently significant parameters, i.e., check if (1)
    #*there's at least one significantly different parameter in the CCG results
    #*and then check if (2) the significant results show the same structure as
    #*in the VCG results.
    Con_sig <- combined_sig %>%
      #if any CCG result is positive, check whether CCG and VCG are identical
      group_by(row_id) %>%
      summarise(flag = any(CCG) & all(CCG == VCG, na.rm = TRUE)) %>%
      #count the frequency of this occurrence
      pull(flag) %>%
      sum()
    #*Count numbers of consistently non-significant parameters, i.e., where
    #*all dose groups are not significantly different to the respective
    #*controls (CCG or VCG)
    Con_non_sig <- apply(CCG_sig | VCG_sig, 1, function(x) {all(!x)}) %>% sum()
    #*Inconsistently significant is defined as all CCG values being non-
    #*significant but there is at least one significance in VCG
    Incon_sig <- combined_sig %>%
      group_by(row_id) %>%
      summarise(flag = all(!CCG) & any(VCG)) %>%
      pull(flag) %>%
      sum()
    #*Inconsistently non significant is defined as there is at least one
    #*significant value in CCG group but this one is missing in VCG
    Incon_non_sig  <- combined_sig %>%
      group_by(row_id) %>%
      summarise(flag = any(CCG) & any(CCG != VCG)) %>%
      pull(flag) %>%
      sum()
    #*Inverse significant is defined as there is at least one significant row
    #*in CCG which is also reproduced by VCGs, but the direction is wrong.
    Inv_sig <- combined_sig %>%
      left_join(
        combined_dir, by = c("row_id", "col"), suffix = c("_sig", "_dir")
        ) %>%
      group_by(row_id) %>%
      summarise(flag = any(CCG_sig & VCG_sig & (CCG_dir != VCG_dir))) %>%
      pull(flag) %>%
      sum()

    # Return the error counts as a data frame
    return(
      tibble(
        n_parameters,
        orig_sig,
        orig_non_sig, 
        con = Con,
        con_sig = Con_sig,
        con_non_sig = Con_non_sig,
        incon_sig = Incon_sig,
        incon_non_sig = Incon_non_sig,
        inv_sig = Inv_sig
      )
    )
  }
  
  parameter_wise_discrepancies <- parameter_wise_consistency_check(
    CCG = CCG_which_res[4:ncol(CCG_which_res)] %>% select(-starts_with("CG_")),
    VCG = VCG_which_res[4:ncol(CCG_which_res)] %>% select(-starts_with("CG_"))
  ) %>%
    mutate(
      n_sig_stats_results = sum(con_sig, incon_non_sig, inv_sig),
      n_non_sig_stats_results = sum(con_non_sig, incon_sig),
      con_sig_performance = paste0(
        con_sig, " out of ", orig_sig,
        " (", round(con_sig / orig_sig * 100), "%)"
      ),
      con_non_sig_performance = paste0(
        con_non_sig, " out of ", orig_non_sig,
        " (", round(con_non_sig / orig_non_sig * 100), "%)"
      ),
      incon_sig_performance = paste0(
        incon_sig, " out of ", orig_non_sig,
        " (", round(incon_sig / orig_non_sig * 100), "%)"
      ),
      incon_non_sig_performance = paste0(
        incon_non_sig, " out of ", orig_sig,
        " (", round(
          incon_non_sig / orig_sig * 100
        ), "%)"
      ),
      inv_sig_performance = paste0(
        inv_sig, " out of ", orig_sig,
        " (", round(inv_sig / orig_sig * 100), "%)"
      )
    )
  t(parameter_wise_discrepancies)
  #*****************************************************************************
  #*Append whether the statistical results between legacy and resampling are
  #*consistent towards each other
  resampling_with_consistency <- resampling_results_flatten
  resampling_with_consistency$consistent <- VCG_res_matching_with_CCG
  
  #*add the consistent results only to those rows where a sig-test was made with
  #*both, CCG and VCG
  legacy_with_consistency <- legacy_study_results_flatten %>%
    left_join(
      resampling_with_consistency %>%
        select(LBTESTCD, LBORRESU, LBSPEC, consistent),
      by = c("LBTESTCD", "LBORRESU", "LBSPEC")
    )
  
  #*****************************************************************************
  #*****************************************************************************
  #*Style the tables----
  #*****************************************************************************
  #*Make the table headers more readable and appealing and reorder the columns
  styled_table <- function(res_table){
    res_table %>%
      select(
        code_unit_spec,
        consistent,
        starts_with(c("CG_", "LD_", "MD_", "HD_")),
        starts_with("out_of_LON_"),
        starts_with("dose_dependency_"),
        starts_with("LLN")
      )
  }
  
  #*get column names which will be automatically renamed if they are present
  #*in the styled table data frame
  column_names <- c(
    "code_unit_spec",
    "consistent",
    "LLN/ULN_M",
    "CG_M", "out_of_LON_CG_M",
    "LD_M", "out_of_LON_LD_M",
    "MD_M", "out_of_LON_MD_M",
    "HD_M", "out_of_LON_HD_M",
    "dose_dependency_M",
    "LLN/ULN_F",
    "CG_F", "out_of_LON_CG_F",
    "LD_F", "out_of_LON_LD_F",
    "MD_F", "out_of_LON_MD_F",
    "HD_F", "out_of_LON_HD_F",
    "dose_dependency_F",
    "where_was_the_discrepancy"
    )
  #how the column names will be renamed afterwards
  renaming <- c(
    "Parameter short name [unit] in site of measurment",
    "Statistical result consistent with both CCG and VCG",
    rep(c(
      "LLN/ULN",
      "Control", "Out of LON",
      "Low dose", "Out of LON",
      "Mid dose", "Out of LON",
      "High dose", "Out of LON",
      "Dose dependency?"
      ), 2),
    "Effect in only one sex?"
  )
  
  #Create the pair list with the column names and labels
  labeling_list <- setNames(renaming, column_names)
  
  #*Create a function to conditionally find an asterisk and color the cell
  #*accordingly
  asterisk_finder <- function(x, stng){
    cells_body(
      columns = !!sym(x), 
      rows =  grepl(stng, !!sym(x))
    )
  }
  #*****************************************************************************
  #Get legacy study with concurrent control groups
  styled_legacy <- styled_table(legacy_with_consistency)
  
  #Check whether all columns are present in the legacy study
  column_names_in_CCG <- column_names %in% colnames(legacy_with_consistency)
  #*****************************************************************************
  #Get legacy study with virtual control groups
  styled_VCG <- styled_table(
    resampling_with_consistency %>% ungroup() %>% select(-LBTESTCD, -LBSPEC)
      )
  
  #*also add the exact discrepancies between CCG and VCG as a separate column
  #*to the styled VCG which is going to be printed as XLSX
  styled_VCG_for_XLSX <- merge(
    styled_VCG,
    discrepancy_check %>%
      #create a column storing the LBTESTCD, LBSPEC, and unit of measurement
      mutate(
        code_unit_spec = paste0(
          LBTESTCD, " [", LBORRESU, "] in ", LBSPEC
        )
      ) %>% select(code_unit_spec, where_was_the_discrepancy),
    by = c("code_unit_spec")
  ) %>%
    #*rearrange the table so that the ones that aren't reproduced are on top of
    #*the table followed by code_unit_spec in alphabetical order
    #*however, keep body weight (BW), body weight gain (BG), food consumption
    #*(FC) and water consumption (WC) values together and place them last.
    arrange(
      #*this tweak creates numerical values and stores everything at the end
      #*of the table if it starts with these prefixes
      as.integer(grepl("^(FC_|WC_|BW_|BG_)", code_unit_spec)),
      #ignore the TRUE/FALSE separation if parameters start with these values
      ifelse(grepl("^(FC_|WC_|BW_|BG_)", code_unit_spec), 1, consistent),
      code_unit_spec
      )
  
  #Check whether all columns are present in the legacy study
  column_names_in_VCG <- column_names %in% colnames(styled_VCG)
  #*****************************************************************************
  #*Sort styled tables by consistency and alphabetical order
  styled_legacy <- styled_legacy %>%
    arrange(
      .,
      as.integer(grepl("^(FC_|WC_|BW_|BG_)", code_unit_spec)),
      #ignore the TRUE/FALSE separation if parameters start with these values
      ifelse(grepl("^(FC_|WC_|BW_|BG_)", code_unit_spec), 1, consistent),
      code_unit_spec
      )
  
  styled_VCG <- styled_VCG %>%
    arrange(
      .,
      as.integer(grepl("^(FC_|WC_|BW_|BG_)", code_unit_spec)),
      #ignore the TRUE/FALSE separation if parameters start with these values
      ifelse(grepl("^(FC_|WC_|BW_|BG_)", code_unit_spec), 1, consistent),
      code_unit_spec
      )
  #*****************************************************************************
  #*Reorder columns. The tables are going to be exported as EXCEL sheets
  reorder_cols <- function(stat_res_table){
    my_order <- column_names
    #Check which columns are present in data frame and drop NA columns
    present_in_df <- match(column_names, colnames(stat_res_table))
    present_in_df <- present_in_df[!is.na(present_in_df)]
    stat_res_table <- stat_res_table[, present_in_df]
    return(stat_res_table)
  }
  
  styled_legacy <- reorder_cols(styled_legacy)
  styled_VCG <- reorder_cols(styled_VCG)
  styled_VCG_for_XLSX <- reorder_cols(styled_VCG_for_XLSX)
  #*****************************************************************************
  #*****************************************************************************
  #*Render results in a GT Table----
  #*******************************************************************************
  table_creator <- function(res_table, tabletitle, tablesubtitle){
    #***************************************************************************
    #***************************************************************************
    #*Create function to draw borders and add spanners----
    #***************************************************************************
    #*Create a function which will create tab spanners across all columns which
    #*share the same suffix. In this case, "_M", or "_F"
    create_spanners <- function(gt_table, suffix, spannertitle) {
      columns <- column_names[column_names %in% colnames(res_table)]
      columns_with_suffix <- columns[grepl(paste0("_", suffix, "$"), columns)]
      if (length(columns_with_suffix) > 0) {
        gt_table <- gt_table %>% 
          tab_spanner(
            label = md(spannertitle),
            columns = columns_with_suffix
          )
      }
      gt_table
    }
    #***************************************************************************
    #Create a function to add a border to columns with prefix "out_of_LON_"
    HCD_border <- function(gt_table, prefix) {
      columns_to_use <- column_names[column_names %in% colnames(res_table)]
      prefix_pattern <- paste0("^", prefix)
      columns_with_prefix <- columns_to_use[
        grepl(prefix_pattern, columns_to_use)
        ]
      if (length(columns_with_prefix) > 0) {
        for (col_name in columns_with_prefix) {
          gt_table <- gt_table %>% 
            tab_style(
              style = cell_borders(
                sides = c("right"),
                weight = px(1)),
              locations = cells_body(columns = col_name)
            )
        }
      }
      gt_table
    }
    #***************************************************************************
    #Do the same for the column labels, not only the column body
    HCD_border_column_labels <- function(gt_table, prefix) {
      columns_to_use <- column_names[column_names %in% colnames(res_table)]
      prefix_pattern <- paste0("^", prefix)
      columns_with_prefix <- columns_to_use[
        grepl(prefix_pattern, columns_to_use)
        ]
      if (length(columns_with_prefix) > 0) {
        for (col_name in columns_with_prefix) {
          gt_table <- gt_table %>% 
            tab_style(
              style = cell_borders(
                sides = c("right"),
                weight = px(1)),
              locations = cells_column_labels(col_name)
            )
        }
      }
      gt_table
    }
    #***************************************************************************
    #***************************************************************************
    #*render the table----
    #***************************************************************************
    res_table %>%
      select(all_of(column_names[column_names_in_CCG])) %>%
      gt() %>%
      #Relabel the columns
      cols_label(
        !!!labeling_list[column_names_in_CCG]
      ) %>%
      tab_header(
        title = md(tabletitle),
        subtitle = md(tablesubtitle)
      ) %>%
      #highlight statistically significant increases in red
      tab_style(
        style = cell_fill(color = "#fdba99"),
        locations = lapply(names(res_table), asterisk_finder, "\\*.*\\(\\+\\)")
      ) %>%
      #highlight statistically significant decreases in cyan
      tab_style(
        style = cell_fill(color = "#a6cee3"),
        locations = lapply(names(res_table), asterisk_finder, "\\*.*\\(\\-\\)")
      ) %>%
      #highlight "1" in dose dependency questions in red
      tab_style(
        style = cell_fill(color = "#fdba99"),
        locations = lapply(
          c(
            "dose_dependency_F",
            "dose_dependency_M"
          ),
          asterisk_finder,
          "1"
        )
      ) %>%
      #highlight "0" in dose dependency questions in blue
      tab_style(
        style = cell_fill(color = "#a6cee3"),
        locations = lapply(
          c("dose_dependency_F",
            "dose_dependency_M"
          ),
          asterisk_finder,
          "0"
        )
      ) %>%
      #Create a big label above all females
      create_spanners(., "F", "**Females**") %>%
      #Create a big label above all males
      create_spanners(., "M", "**Males**") %>%
      #Add a border separating females from males
      tab_style(
        style = list(
          cell_borders(
            sides = "everything",
            color = "#000000",
            weight = px(2)
          ),
          cell_fill(color = "#f0f0f0")
        ),
        locations = cells_column_spanners(
          spanners = c(
            md("**Females**"), md("**Males**")
          )
        )
      ) %>%
      #draw the border through the end of the "males" column
      tab_style(
        style = list(
          cell_borders(
            sides = "right",
            color = "#000000",
            weight = px(2)
          ),
          cell_fill(color = "#f0f0f0")
        ),
        locations = cells_column_labels(dose_dependency_M)
      ) %>%
      #Do the same for the cells body
      tab_style(
        style = cell_borders(
          sides = c("right"),
          weight = px(2)),
        locations = cells_body(columns = dose_dependency_M)
      ) %>%
      #Add a border left of control
      tab_style(
        style = 
          cell_borders(
            sides = "left",
            color = "#000000",
            weight = px(2)
          ),
        locations = cells_column_labels(
          c(`LLN/ULN_M`, `LLN/ULN_F`)
        )
      ) %>%
      #Add a border left of control values
      tab_style(
        style = 
          cell_borders(
            sides = "left",
            color = "#000000",
            weight = px(2)
          ),
        locations = cells_body(
          columns = c(`LLN/ULN_M`, `LLN/ULN_F`)
        )
      ) %>%
      #Add a border to separate the limits of normal tab
      HCD_border(., "LLN") %>%
      HCD_border_column_labels(., "LLN") %>%
      #Add a border separating each dose group
      HCD_border(., "out_of_LON_") %>%
      HCD_border_column_labels(., "out_of_LON_") %>%
      #Color the parameter names in a light yellow
      tab_style(
        style = list(
          cell_text(weight = "bold")
        ),
        locations = cells_body(columns = code_unit_spec)
      ) %>%
      tab_style(
        style = list(
          cell_text(weight = "bold")
        ),
        locations = cells_column_labels(columns = code_unit_spec)
      ) %>%
      cols_align(align = "center", columns = code_unit_spec) %>%
      #Color the header in grey
      tab_style(
        style = cell_fill(color = "#f0f0f0"),
        locations = cells_column_labels(everything())
      ) %>%
      fmt_number()
  }
  #Return the tables out of function
  tabled_legacy_results <- table_creator(
    res_table = styled_legacy,
    tabletitle = "**Legacy Study original results**",
    tablesubtitle = "**Original results**"
  )
  
  tabled_VCG_results <- table_creator(
    res_table = styled_VCG,
    tabletitle = "**Results after CCG replaced by VCGs**",
    tablesubtitle = paste0(
      "**",
      "Statistical results of ", sum(VCG_res_matching_with_CCG),
      " out of ",
      length(VCG_res_matching_with_CCG),
      " parameters (", matching_percentage,
      " %) reproduced**"
    )
  )
  
  return(list(
    tabled_legacy_results,
    tabled_VCG_results,
    styled_legacy,
    styled_VCG_for_XLSX,
    consistency_counts,
    parameter_wise_discrepancies
    ))
}

#*******************************************************************************
#*******************************************************************************
#*Get results for legacy study A----
#*******************************************************************************
tables_legacy_study_A <- study_results_table(
  legacy_study_results = CCG_legacy_study_A,
  resampling_results = VCG_legacy_study_A
  )

tables_legacy_study_A_imp_median <- study_results_table(
  legacy_study_results = CCG_legacy_study_A,
  resampling_results = VCG_legacy_study_A_imp_median
)


tables_legacy_study_A_imp_rs <- study_results_table(
  legacy_study_results = CCG_legacy_study_A,
  resampling_results = VCG_legacy_study_A_imp_rs
)

tables_legacy_study_A_imp_pmm <- study_results_table(
  legacy_study_results = CCG_legacy_study_A,
  resampling_results = VCG_legacy_study_A_imp_pmm
)

#*******************************************************************************
print("Legacy study A")
print(t(tables_legacy_study_A[[5]]))
print(t(tables_legacy_study_A[[6]]))
tables_legacy_study_A[[1]]

#save the table
gtsave(tables_legacy_study_A[[1]], "legacy_study_A_CCG_results_table.html")
gtsave(tables_legacy_study_A[[1]], "legacy_study_A_CCG_results_table.png", vwidth = 3000)
#save the results as Excel sheet
write.xlsx(tables_legacy_study_A[[3]], "legacy_study_A_CCG_results_table.xlsx")
#*******************************************************************************
tables_legacy_study_A[[2]]
#save the table
gtsave(tables_legacy_study_A[[2]], "legacy_study_A_VCG_results_table.html")
gtsave(tables_legacy_study_A[[2]], "legacy_study_A_VCG_results_table.png", vwidth = 3000)
#save the results as Excel sheet
write.xlsx(tables_legacy_study_A[[4]], "legacy_study_A_VCG_results_table.xlsx")
#*******************************************************************************
tables_legacy_study_A_imp_median[[2]]
#save the table
gtsave(tables_legacy_study_A_imp_median[[2]], "legacy_study_A_VCG_imp_median_results_table.html")
gtsave(tables_legacy_study_A_imp_median[[2]], "legacy_study_A_VCG_imp_median_results_table.png", vwidth = 3000)
#save the results as Excel sheet
write.xlsx(tables_legacy_study_A_imp_median[[4]], "legacy_study_A_VCG_imp_median_results_table.xlsx")
#*******************************************************************************
tables_legacy_study_A_imp_rs[[2]]
#save the table
gtsave(tables_legacy_study_A_imp_rs[[2]], "legacy_study_A_VCG_imp_rs_results_table.html")
gtsave(tables_legacy_study_A_imp_rs[[2]], "legacy_study_A_VCG_imp_rs_results_table.png", vwidth = 3000)
#save the results as Excel sheet
write.xlsx(tables_legacy_study_A_imp_rs[[4]], "legacy_study_A_VCG_imp_rs_results_table.xlsx")
#*******************************************************************************
tables_legacy_study_A_imp_pmm[[2]]
#save the table
gtsave(tables_legacy_study_A_imp_pmm[[2]], "legacy_study_A_VCG_imp_pmm_results_table.html")
gtsave(tables_legacy_study_A_imp_pmm[[2]], "legacy_study_A_VCG_imp_pmm_results_table.png", vwidth = 3000)
#save the results as Excel sheet
write.xlsx(tables_legacy_study_A_imp_pmm[[4]], "legacy_study_A_VCG_imp_pmm_results_table.xlsx")
#*******************************************************************************
#*******************************************************************************
#*Get results for legacy study B----
#*******************************************************************************
tables_legacy_study_B <- study_results_table(
  legacy_study_results = CCG_legacy_study_B,
  resampling_results = VCG_legacy_study_B
)

tables_legacy_study_B_recovery <- study_results_table(
  legacy_study_results = CCG_legacy_study_B_recovery,
  resampling_results = VCG_legacy_study_B_recovery
)

tables_legacy_study_B_imp_median <- study_results_table(
  legacy_study_results = CCG_legacy_study_B,
  resampling_results = VCG_legacy_study_B_imp_median
)

tables_legacy_study_B_imp_rs <- study_results_table(
  legacy_study_results = CCG_legacy_study_B,
  resampling_results = VCG_legacy_study_B_imp_rs
)

tables_legacy_study_B_imp_pmm <- study_results_table(
  legacy_study_results = CCG_legacy_study_B,
  resampling_results = VCG_legacy_study_B_imp_pmm
)

#*******************************************************************************
print("Legacy study B")
print(t(tables_legacy_study_B[[5]]))
print(t(tables_legacy_study_B[[6]]))

tables_legacy_study_B[[1]]
#save the table
gtsave(tables_legacy_study_B[[1]], "legacy_study_B_CCG_results_table.html")
gtsave(tables_legacy_study_B[[1]], "legacy_study_B_CCG_results_table.png", vwidth = 3000)
#save the results as Excel sheet
write.xlsx(tables_legacy_study_B[[3]], "legacy_study_B_CCG_results_table.xlsx")
#*******************************************************************************
tables_legacy_study_B[[2]]
#save the table
gtsave(tables_legacy_study_B[[2]], "legacy_study_B_VCG_results_table.html")
gtsave(tables_legacy_study_B[[2]], "legacy_study_B_VCG_results_table.png", vwidth = 3000)
#save the results as Excel sheet
write.xlsx(tables_legacy_study_B[[4]], "legacy_study_B_VCG_results_table.xlsx")
#*******************************************************************************
print("Legacy study B recovery")
print(tables_legacy_study_B_recovery[[5]])
tables_legacy_study_B_recovery[[1]]
#save the table
gtsave(tables_legacy_study_B_recovery[[1]], "legacy_study_B_recovery_CCG_results_table.html")
gtsave(tables_legacy_study_B_recovery[[1]], "legacy_study_B_recovery_CCG_results_table.png", vwidth = 3000)
#save the results as Excel sheet
write.xlsx(tables_legacy_study_B_recovery[[3]], "legacy_study_B_recovery_CCG_results_table.xlsx")
#*******************************************************************************
tables_legacy_study_B_recovery[[2]]
#save the table
gtsave(tables_legacy_study_B_recovery[[2]], "legacy_study_B_recovery_VCG_results_table.html")
gtsave(tables_legacy_study_B_recovery[[2]], "legacy_study_B_recovery_VCG_results_table.png", vwidth = 3000)
#save the results as Excel sheet
write.xlsx(tables_legacy_study_B_recovery[[4]], "legacy_study_B_recovery_VCG_results_table.xlsx")
#*******************************************************************************
tables_legacy_study_B_imp_median[[2]]
#save the table
gtsave(tables_legacy_study_B_imp_median[[2]], "legacy_study_B_VCG_imp_median_results_table.html")
gtsave(tables_legacy_study_B_imp_median[[2]], "legacy_study_B_VCG_imp_median_results_table.png", vwidth = 3000)
#save the results as Excel sheet
write.xlsx(tables_legacy_study_B_imp_median[[4]], "legacy_study_B_VCG_imp_median_results_table.xlsx")
#*******************************************************************************
tables_legacy_study_B_imp_rs[[2]]
#save the table
gtsave(tables_legacy_study_B_imp_rs[[2]], "legacy_study_B_VCG_imp_rs_results_table.html")
gtsave(tables_legacy_study_B_imp_rs[[2]], "legacy_study_B_VCG_imp_rs_results_table.png", vwidth = 3000)
#save the results as Excel sheet
write.xlsx(tables_legacy_study_B_imp_rs[[4]], "legacy_study_B_VCG_imp_rs_results_table.xlsx")
#*******************************************************************************
tables_legacy_study_B_imp_pmm[[2]]
#save the table
gtsave(tables_legacy_study_B_imp_pmm[[2]], "legacy_study_B_VCG_imp_pmm_results_table.html")
gtsave(tables_legacy_study_B_imp_pmm[[2]], "legacy_study_B_VCG_imp_pmm_results_table.png", vwidth = 3000)
#save the results as Excel sheet
write.xlsx(tables_legacy_study_B_imp_pmm[[4]], "legacy_study_B_VCG_imp_pmm_results_table.xlsx")
#*******************************************************************************
#*******************************************************************************
#*Get results for legacy study C----
#*******************************************************************************
tables_legacy_study_C <- study_results_table(
  legacy_study_results = CCG_legacy_study_C,
  resampling_results = VCG_legacy_study_C
)

tables_legacy_study_C_imp_median <- study_results_table(
  legacy_study_results = CCG_legacy_study_C,
  resampling_results = VCG_legacy_study_C_imp_median
)

tables_legacy_study_C_imp_rs <- study_results_table(
  legacy_study_results = CCG_legacy_study_C,
  resampling_results = VCG_legacy_study_C_imp_rs
)

tables_legacy_study_C_imp_pmm <- study_results_table(
  legacy_study_results = CCG_legacy_study_C,
  resampling_results = VCG_legacy_study_C_imp_pmm
)
#*******************************************************************************
print("Legacy study C")
print(t(tables_legacy_study_C[[5]]))
print(t(tables_legacy_study_C[[6]]))
tables_legacy_study_C[[1]]
#save the table
gtsave(tables_legacy_study_C[[1]], "legacy_study_C_CCG_results_table.html")
gtsave(tables_legacy_study_C[[1]], "legacy_study_C_CCG_results_table.png", vwidth = 3000)
#save the results as Excel sheet
write.xlsx(tables_legacy_study_C[[3]], "legacy_study_C_CCG_results_table.xlsx")
#*******************************************************************************
tables_legacy_study_C[[2]]
#save the table
gtsave(tables_legacy_study_C[[2]], "legacy_study_C_VCG_results_table.html")
gtsave(tables_legacy_study_C[[2]], "legacy_study_C_VCG_results_table.png", vwidth = 3000)
#save the results as Excel sheet
write.xlsx(tables_legacy_study_C[[4]], "legacy_study_C_VCG_results_table.xlsx")
# #*******************************************************************************
tables_legacy_study_C_imp_median[[2]]
#save the table
gtsave(tables_legacy_study_C_imp_median[[2]], "legacy_study_C_VCG_imp_median_results_table.html")
gtsave(tables_legacy_study_C_imp_median[[2]], "legacy_study_C_VCG_imp_median_results_table.png", vwidth = 3000)
#save the results as Excel sheet
write.xlsx(tables_legacy_study_C_imp_median[[4]], "legacy_study_C_VCG_imp_median_results_table.xlsx")
#*******************************************************************************
tables_legacy_study_C_imp_rs[[2]]
#save the table
gtsave(tables_legacy_study_C_imp_rs[[2]], "legacy_study_C_VCG_imp_rs_results_table.html")
gtsave(tables_legacy_study_C_imp_rs[[2]], "legacy_study_C_VCG_imp_rs_results_table.png", vwidth = 3000)
#save the results as Excel sheet
write.xlsx(tables_legacy_study_C_imp_rs[[4]], "legacy_study_C_VCG_imp_rs_results_table.xlsx")
#*******************************************************************************
tables_legacy_study_C_imp_pmm[[2]]
#save the table
gtsave(tables_legacy_study_C_imp_pmm[[2]], "legacy_study_C_VCG_imp_pmm_results_table.html")
gtsave(tables_legacy_study_C_imp_pmm[[2]], "legacy_study_C_VCG_imp_pmm_results_table.png", vwidth = 3000)
#save the results as Excel sheet
write.xlsx(tables_legacy_study_C_imp_pmm[[4]], "legacy_study_C_VCG_imp_pmm_results_table.xlsx")
#*******************************************************************************
#*******************************************************************************
#*Check results where you replaced legacy study C with CCG animals from A----
#*******************************************************************************
tables_legacy_study_C_with_A_CCGs <- study_results_table(
  legacy_study_results = CCG_legacy_study_C,
  resampling_results = CCG_A_legacy_study_C
)
