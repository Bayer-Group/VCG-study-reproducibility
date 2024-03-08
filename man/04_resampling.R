#04_resampling.R

#*This script is based on "02_make_statistical_test.R". However, instead of the
#*concurrent control group (CCG) of the legacy studies, this set randomly
#*samples control group values (without replacement) of the HCD. These values
#*are then used to replace the CCGs.
#*It can be selected whether all CCGs are removed or only a fraction of them.
#*Furthermore, the number of iterations can be selected.
#*This script needs the location of the CSV files containing SEND formatted data
#*stored in Data/Original/legacy_study and Data/Original/quantitative_parameters.





#*Create a function which does the following:
#*  - Replace the original control group values with a random sample (without
#*    replacement) of values from HCD and recalculate the Dunnett test
#*  - Repeat previous step n times
#*  - Summarize results in a list giving the percentage of how often the
#*    original result was reproduced (per dose group)
#*    as well as the mean values of each sample of the VCGs
#*******************************************************************************
#Check if packages are installed------------------------------------------------
#*******************************************************************************
#*******************************************************************************
require(tidyverse)
#*******************************************************************************
#*Get the function which calculates respective significance test
source(paste0(rootpath, "/man/02_make_statistical_test.R"))
#*Get function which plots the INITBW distribution before and after matching by
#*ININTBW
source(paste0(rootpath, "/man/10_BW_selection_visualization.R"))
#*Beginning of function---------------------------------------------------------
#*******************************************************************************
#*******************************************************************************
resampling_VCG <- function(
    legacy_study = legacy_study_A,
    HCD = HCD,
    iterations = 100,
    replacement_aim = 100,
    match_by_initbw = TRUE,
    sensitivity_analyses = FALSE,
    plottitle = "Legacy study A",
    #todo****************************
    sampling = "random",
    refined_resampling = FALSE,
    simulate_missing_parameters = FALSE
    #todo end************************
)
{
  #*get summary statistics of HCD for each parameter by each sex
  HCD_refdata <- HCD %>%
    group_by(SEX, LBTESTCD, LBSPEC, LBORRESU) %>%
    summarise(
      `Percentile2.5th` = quantile(LBORRES, .025),
      `Percentile97.5th` = quantile(LBORRES, .975)
    ) %>%
    suppressMessages()
  #*****************************************************************************
  #*Remove all animals who not match initial body weight of dose groups in
  #*legacy study (if selected)----
  HCD_initbw_matched <- if(match_by_initbw){
    
    minmax_initbw <- legacy_study %>%
      filter(
        trial_set_description != "CG",
        LBTESTCD == "BW_D001"
      ) %>%
      #extract highest and lowest values of dose groups per sex
      group_by(SEX) %>%
      summarize(
        min_initbw = min(LBORRES),
        max_initbw = max(LBORRES)
      )
    
    legacy_minmax <- split(minmax_initbw, minmax_initbw$SEX)
    
    legacy_minmax_list <- lapply(
      legacy_minmax,
      function(x){
        x %>% select(-SEX) %>% as_vector()
      }
    )
    #***************************************************************************
    #*Filter HCD so that animals are only selected who match INITBW of legacy 
    filtered_HCD_USUBJIDs <- split(
      HCD %>%
        filter(LBTESTCD == "BW_D001") %>%
        select(USUBJID, LBORRES),
      HCD %>%
        filter(LBTESTCD == "BW_D001") %>%
        pull(SEX)
    )
    
    #Select subject IDs from the animals within INITBW
    filtered_by_initbw <- mapply(function(x,y){
      x %>% 
        filter(between(LBORRES, y[1], y[2])) %>%
        pull(USUBJID)
    },
    x = filtered_HCD_USUBJIDs, y = legacy_minmax_list, SIMPLIFY = F
    ) %>%
      unlist()
    
    #filter now only those animals who match in their INITBW
    HCD_filtered_by_initbw <- HCD %>% filter(USUBJID %in% filtered_by_initbw)
  #*****************************************************************************
  #Plot the INITBW distribution of HCD before and after matching by INITBW
  bw_selection_plot(HCD, HCD_filtered_by_initbw, legacy_study, plottitle)
  #*****************************************************************************
  HCD_filtered_by_initbw
    
  }else{
    HCD
  }
  #*****************************************************************************
  #*****************************************************************************
  #*Select only those LBTESTCD-values which are present in the HCD.
  #*Return a warning message listing parameters which were present in the 
  #*legacy studies but not in the HCD.----
  #*****************************************************************************
  HCD_LBTESTCD_LBSPEC_LBORRESU <- unique(
    paste(
      HCD_initbw_matched$LBTESTCD,
      HCD_initbw_matched$LBSPEC,
      HCD_initbw_matched$LBORRESU,
      sep = "_"
    )
  )
  
  missing_LBTESTCD <- legacy_study %>%
    mutate(
      LBTESTCD_LBSPEC_LBORRESU = paste(LBTESTCD, LBSPEC, LBORRESU, sep = "_")
    ) %>%
    filter(!LBTESTCD_LBSPEC_LBORRESU %in% HCD_LBTESTCD_LBSPEC_LBORRESU) %>%
    pull(LBTESTCD_LBSPEC_LBORRESU) %>%
    unique() %>% sort()
  if(!is.null(missing_LBTESTCD)){
    warning(
      "The following parameters are not present in the HCD and will therefore be omitted:\n - ",
      paste(sort(missing_LBTESTCD), collapse = "\n - "))
  }
  
  legacy_study_only_LBTESTCD_from_HCD_with_LBSPECinfo <- legacy_study %>%
    mutate(
      LBTESTCD_LBSPEC_LBORRESU = paste(LBTESTCD, LBSPEC, LBORRESU, sep = "_")
    ) %>%
    filter(!LBTESTCD_LBSPEC_LBORRESU %in% missing_LBTESTCD)
  
  #*Do this the same way around. Omit all parameters in HCD data which were not
  #*calculated in the legacy study
  HCD_selected_LBTESTCDs <- HCD_initbw_matched %>%
    mutate(
      LBTESTCD_LBSPEC_LBORRESU = paste(LBTESTCD, LBSPEC, LBORRESU, sep = "_")
    ) %>%
    filter(
      LBTESTCD_LBSPEC_LBORRESU %in%
        unique(
          legacy_study_only_LBTESTCD_from_HCD_with_LBSPECinfo %>%
            pull(LBTESTCD_LBSPEC_LBORRESU)
        )
    ) %>%
    select(-LBTESTCD_LBSPEC_LBORRESU)
  
  #remove concatenated column
  legacy_study_only_LBTESTCD_from_HCD <-
    legacy_study_only_LBTESTCD_from_HCD_with_LBSPECinfo %>%
    select(-LBTESTCD_LBSPEC_LBORRESU)
  #*****************************************************************************
  #*****************************************************************************
  #Keep sentinel animals--------------------------------------------------------
  #*****************************************************************************
  #*The VCG sample population is selected to maintain the location parameters
  #*of the CCG as close as possible
  #*Remove individuals from the CCG, based on selection.
  
  #How many animals do you want to keep?
  n_keep <- 
    (
      (
        legacy_study_only_LBTESTCD_from_HCD %>%
          filter(trial_set_description == "CG") %>%
          pull(USUBJID) %>%
          n_distinct()
      ) -
        (
          legacy_study_only_LBTESTCD_from_HCD %>%
            filter(trial_set_description == "CG") %>%
            pull(USUBJID) %>%
            n_distinct()
        ) *
        replacement_aim / 100 %>%
        round()
    ) /
    #divide the whole number by the number of sexes present in the legacy study
    n_distinct(legacy_study_only_LBTESTCD_from_HCD$SEX)
  
  #Get unique subject IDs from control group of the case study ordered by
  #initial body weight (BW_D1)
  CCG_by_sex <- split(
    legacy_study_only_LBTESTCD_from_HCD %>%
      filter(
        trial_set_description == "CG",
        LBTESTCD == "BW_D001"
      ) %>%
      select(LBORRES, USUBJID),
    legacy_study_only_LBTESTCD_from_HCD %>%
      filter(
        trial_set_description == "CG",
        LBTESTCD == "BW_D001"
      ) %>%
      pull(SEX)
  )
  
  sentinel_animals_USUBJIDs <- lapply(CCG_by_sex, function(x){
    if(n_keep == 0) {
      NULL
    } else if(n_keep %% 2 != 0) {
      #*get one animal from the middle if the amount of sentinel animals is odd
      x %>%
        arrange(desc(LBORRES)) %>%
        slice(
          1:floor(n_keep / 2),
          (ceiling(nrow(x)) / 2),
          (nrow(x) - floor(n_keep / 2) + 1):nrow(x)
        ) %>%
        pull(USUBJID)
    } else if(n_keep %% 2 == 0) {
      #if the amount of sentinel animals is even, get only the extremes
      x %>%
        arrange(desc(LBORRES)) %>%
        slice(
          1:floor(n_keep / 2),
          (nrow(x) - floor(n_keep / 2) + 1):nrow(x)
        ) %>%
        pull(USUBJID)
    }
  }) %>%
    unlist()
  
  #select sentinel animals based on USUBJIDs
  sentinel_animals <- legacy_study_only_LBTESTCD_from_HCD %>%
    filter(USUBJID %in% sentinel_animals_USUBJIDs)
  #*****************************************************************************
  #*****************************************************************************
  #*Replace CCG with VCGs and calculate significance test----
  #*****************************************************************************
  #*make a loop where you sample from the selected subset and recalculate the
  #*respective significance test
  # collected_vcg_samples <- list()
  resampling_significance_test_results <- list()
  
  #*sample n animals from the sampling group with n being the number of animals per group
  #*make a stratified sampling so each group is equally represented in the sample
  #*****************************************************************************
  #*begin of loop
  #*****************************************************************************
  for (i in 1:iterations) {
    #Print number of iterations
    print(paste0("iteration ", i, " of ", iterations))
    #Get INITBW of HCD by sex. This is used to extract USUBJIDs of VCG animals
    HCD_by_sex <- split(
      HCD_selected_LBTESTCDs %>%
        filter(LBTESTCD == "BW_D001") %>%
        select(USUBJID),
      HCD_selected_LBTESTCDs %>%
        filter(LBTESTCD == "BW_D001") %>%
        pull(SEX)
    )
    
    #Sample n animals selected by USUBJID on initial body weight
    vcg_sample_USUBJID <- mapply(function(x,y) {
      lapply(x, function(row){slice_sample(x, n = nrow(y) - n_keep)})
    }, x = HCD_by_sex, y = CCG_by_sex) %>%
      unlist()
    
    
    
    #filter VCGs based on selected USUBJID
    vcg_sample <- HCD_selected_LBTESTCDs %>%
      filter(USUBJID %in% vcg_sample_USUBJID)
    
    #Combine VCG sample and sentinel animal values
    vcg_sample_with_sentinel <- vcg_sample %>%
      rbind(sentinel_animals, vcg_sample)
    
    
    #Collect the values for all VCG samples
    # collected_vcg_samples[[i]] <- vcg_sample_with_sentinel
    #***************************************************************************
    #*Calculate the significance tests with the VCG instead of CCG
    #***************************************************************************
    #*Replace CCG values of legacy studies with VCGs
    legacy_study_with_VCG <- legacy_study_only_LBTESTCD_from_HCD %>%
      filter(trial_set_description != "CG") %>%
      rbind(vcg_sample_with_sentinel) %>%
      #security step: remove duplicates
      unique()
    
    #Calculate significance test
    vcg_significance_res <- tox_statistical_test(
      studydata = legacy_study_with_VCG,
      refdata = HCD_refdata,
      power_calculations = sensitivity_analyses
    ) %>%
      suppressMessages()
    
    #Collect the results in list
    resampling_significance_test_results[[i]] <- vcg_significance_res
  }
  #*****************************************************************************
  #end of loop
  #*****************************************************************************
  #*****************************************************************************
  #*Collect all results from the resampling approach----------------------------
  #*****************************************************************************
  #*Extract results of significance tests from VCGs
  vcg_sigtest_results <- lapply(
    resampling_significance_test_results,
    function(x){
      x[[1]]
    }
  )
  vcg_results_flatten <- vcg_sigtest_results %>%
    imap_dfr(~ .x %>% as_tibble(), .id = "iteration")
  #*****************************************************************************
  #*****************************************************************************
  #compute majority vote for significance tests----
  #*****************************************************************************
  #*The majority vote is done in the following approach:
  #*1st vote: do we have more significant results then non-significant ones?
  #*(in case of a tie, vote for significant results).
  #*2nd vote: if we have a majority of iterations leading to a statistically
  #*significante result, we check, which classification is given the most
  #*(i.e., "*" for p-val <= 0.05, "**" for pval <= 0.01,
  #*and "***" for p-val <= 0.001).
  
  #Extract asterisks of all results
  vcg_res_asterisks <- vcg_results_flatten %>%
    select(
      LBTESTCD, LBSPEC, LBORRESU, starts_with(c("CG", "LD", "MD", "HD"))
    ) %>%
    mutate(
      across(
        starts_with(c("CG", "LD", "MD", "HD")),
        #Get only the asterisks
        function(x){gsub("[^\\*]", "", x)}
      )
    )
  
  #check whether majority of results vote for significance
  vcg_vote_sig_nonsig <- vcg_res_asterisks %>%
    group_by(LBTESTCD, LBSPEC, LBORRESU) %>%
    mutate(
      across(
        starts_with(c("CG", "LD", "MD", "HD")),
        function(x){if_else(str_count(x, "\\*") > 0, "1", "0")}
      )
    ) %>%
    #upon a tie, classify as statistically significant
    summarize(
      across(
        starts_with(c("CG", "LD", "MD", "HD")),
        function(x){if_else(sum(x == "1") >= sum(x == "0"), TRUE, FALSE)}
      )
    ) %>%
    #pivot into long structure
    pivot_longer(
      cols = starts_with(c("CG", "LD", "MD", "HD")),
      names_to = "colnames",
      values_to = "is_significant"
    ) %>%
    #*create a unique ID for further filtering out of LBTESTCD, LBSPEC, LBORRESU,
    #*and colnames
    mutate(LLLc = paste(LBTESTCD, LBSPEC, LBORRESU, colnames, sep = ".")) %>%
    suppressMessages()
  
  
  #Calculate percentage of cases where no significance is detected
  vcg_vote_no_sig <- vcg_res_asterisks %>%
    group_by(LBTESTCD, LBSPEC, LBORRESU) %>%
    mutate(
      across(
        starts_with(c("CG", "LD", "MD", "HD")),
        function(x){str_count(x, "\\*")}
      )
    ) %>%
    #pivot into long structure to count majority votes
    pivot_longer(
      cols = starts_with(c("CG", "LD", "MD", "HD")),
      names_to = "colnames",
      values_to = "asterisk_count"
    ) %>%
    #*create a unique ID for further filtering out of LBTESTCD, LBSPEC,
    #*LBORRESU, and colnames
    mutate(LLLc = paste(LBTESTCD, LBSPEC, LBORRESU, colnames, sep = ".")) %>%
    #*extract those rows where no significance is the majority of votes
    filter(
      LLLc %in% (
        vcg_vote_sig_nonsig %>% filter(!is_significant) %>% pull(LLLc)
      )
    ) %>%
    #drop LLLc column
    select(-LLLc) %>%
    #group by everything
    group_by(LBTESTCD, LBSPEC, LBORRESU, colnames, asterisk_count) %>%
    #add the frequency of occurred significances
    summarise(n = n(), .groups = "drop") %>%
    #drop the ones where significances were detected
    filter(asterisk_count == 0) %>%
    #write non significant results with percentage leading to this res
    mutate(
      result = paste0("n.s. {", round(n / iterations * 100), " %}")
    ) %>%
    select(-asterisk_count, -n)
  
  #calculate absolute percentage of votes if there was a significant result
  vcg_percentage_sig <- vcg_res_asterisks %>%
    group_by(LBTESTCD, LBSPEC, LBORRESU) %>%
    mutate(
      across(
        starts_with(c("CG", "LD", "MD", "HD")),
        function(x){str_count(x, "\\*")}
      )
    ) %>%
    #pivot into long structure to count majority votes
    pivot_longer(
      cols = starts_with(c("CG", "LD", "MD", "HD")),
      names_to = "colnames",
      values_to = "asterisk_count"
    ) %>%
    #*create a unique ID for further filtering out of LBTESTCD, LBSPEC,
    #*LBORRESU, and colnames
    mutate(LLLc = paste(LBTESTCD, LBSPEC, LBORRESU, colnames, sep = ".")) %>%
    #*extract those rows where no significance is the majority of votes
    filter(
      LLLc %in% (
        vcg_vote_sig_nonsig %>% filter(is_significant) %>% pull(LLLc)
      )
    ) %>%
    #drop LLLc column
    select(-LLLc) %>%
    #drop the ones where significances were not detected
    filter(asterisk_count != 0) %>%
    #group by everything
    group_by(LBTESTCD, LBSPEC, LBORRESU, colnames) %>%
    #add the frequency of occurred significances
    summarise(n = n(), .groups = "drop") %>%
    #write non significant results with percentage leading to this res
    mutate(
      percentage = round(n / iterations * 100)
    ) %>%
    select(-n)
  
  
  #check for the most frequent number of asterisks
  vcg_vote_most_frequent_sig <- vcg_res_asterisks %>%
    group_by(LBTESTCD, LBSPEC, LBORRESU) %>%
    mutate(
      across(
        starts_with(c("CG", "LD", "MD", "HD")),
        function(x){str_count(x, "\\*")}
      )
    ) %>%
    #pivot into long structure to count majority votes
    pivot_longer(
      cols = starts_with(c("CG", "LD", "MD", "HD")),
      names_to = "colnames",
      values_to = "asterisk_count"
    ) %>%
    #*create a unique ID for further filtering out of LBTESTCD, LBSPEC,
    #*LBORRESU, and colnames
    mutate(LLLc = paste(LBTESTCD, LBSPEC, LBORRESU, colnames, sep = ".")) %>%
    #*extract those rows where no significance is the majority of votes
    filter(
      LLLc %in% (
        vcg_vote_sig_nonsig %>% filter(is_significant) %>% pull(LLLc)
      )
    ) %>%
    #drop LLLc column
    select(-LLLc) %>%
    #group by everything
    group_by(LBTESTCD, LBSPEC, LBORRESU, colnames, asterisk_count) %>%
    #add the frequency of occurred significances
    summarise(n = n(), .groups = "drop") %>%
    #*remove the ones which have no significance as they are already covered
    #*in the previous data frame
    filter(asterisk_count != 0) %>%
    #*arrange in descending order, so that the one with the highest number of
    #*asterisks is shown first
    arrange(
      LBTESTCD, LBSPEC, LBORRESU, colnames, desc(n), desc(asterisk_count)
    ) %>%
    group_by(LBTESTCD, LBSPEC, LBORRESU, colnames) %>%
    #extract the one with the highest frequency
    filter(row_number() == 1) %>%
    #append absolute percentage of all votes leading to a significant result
    left_join(
      vcg_percentage_sig, by = c("LBTESTCD", "LBSPEC", "LBORRESU", "colnames")
      ) %>%
    #replace "n" with the percentage
    mutate(
      result = paste0(strrep("*", asterisk_count), " {", percentage, " %}")
    ) %>%
    select(-asterisk_count, -n, -percentage)
  
  #*add the ones where no statistical test was performed in general
  vcg_vote_no_test <- vcg_vote_sig_nonsig %>%
    filter(is.na(is_significant)) %>%
    #make results column with NA as value
    mutate(result = NA) %>%
    select(-is_significant, -LLLc)
  
  #combine all results
  vcg_sig_res_majority_votes <- rbind.data.frame(
    vcg_vote_no_sig,
    vcg_vote_most_frequent_sig,
    vcg_vote_no_test
  ) %>%
    #separate dose group and sex
    separate(colnames, into = c("trial_set_description", "SEX"), sep = "_")
  #*****************************************************************************
  #*****************************************************************************
  #*Calculate the median of everything numeric----
  #*****************************************************************************
  #get the median value of the mean value from each iteration
  median_vals_of_vcg_res_means <- vcg_results_flatten %>%
    select(
      LBTESTCD, LBSPEC, LBORRESU, starts_with(c("CG", "LD", "MD", "HD"))) %>%
    mutate(
      across(
        starts_with(c("CG", "LD", "MD", "HD")),
        #Get only the mean values. Drop everything after "\u00B1"
        function(x){str_extract(x, "[^\u00B1]+") %>% parse_number(x)}
      )
    ) %>%
    group_by(LBTESTCD, LBSPEC, LBORRESU) %>%
    summarise(
      #*Attention! Here, not the median over the entire set is calculated but
      #*the median over all means of each iteration!
      across(
        starts_with(c("CG", "LD", "MD", "HD")),
        function(x) signif(median(x, na.rm = T),3))
    ) %>%
    suppressMessages()
  #*****************************************************************************
  #*Calculate median standard deviations
  #*Get the average standard deviation values for each LBTESTCD
  #*pivot data into long structure
  median_vals_of_vcg_res_sd <- vcg_results_flatten  %>%
    select(LBTESTCD, LBSPEC, LBORRESU, starts_with(c("CG", "LD", "MD", "HD"))
    ) %>%
    mutate(
      across(
        starts_with(c("CG", "LD", "MD", "HD")),
        #Get only the standard deviation values. Drop everything before "\u00B1"
        function(x){str_extract(x, "\u00B1 (.+)") %>% parse_number(x)}
      )
    ) %>%
    group_by(LBTESTCD, LBSPEC, LBORRESU) %>%
    summarise(
      #*Attention! Here, not the median over the entire set is calculated but
      #*the median over all sds of each iteration!
      across(
        starts_with(c("CG", "LD", "MD", "HD")),
        function(x) signif(median(x, na.rm = T), 2))
    ) %>%
    #the warning comes from trying to extract empty values. we can ignore this
    suppressWarnings() %>%
    suppressMessages()
  
  #*****************************************************************************
  #*Extract information on how many dose group animals are within the HCD range
  #*****************************************************************************
  #*Since the information is the same in every iteration (due to the fact that
  #*HCD does not change nor the dose group location parameters), the results of
  #*the very first iteration are extracted and then attached to the final
  #*results table
  outside_of_HCD_doses <- resampling_significance_test_results[[1]][[1]] %>%
    #remove VCGs. They do change in each iteration.
    select(
      -out_of_LON_CG_F,
      -out_of_LON_CG_M
    ) %>%
    select(
      LBTESTCD,
      LBSPEC,
      LBORRESU,
      starts_with("out")
    )
  #*****************************************************************************
  #*Extract information on how many VCG animals are within the HCD range
  #*****************************************************************************
  outside_of_HCD_VCGs <- 
    #get the fist data frame of the nested list
    lapply(resampling_significance_test_results, function(x){x[[1]]}) %>%
    #flatten list
    imap_dfr(~ .x %>% as_tibble(), .id = "iteration") %>%
    select(
      iteration,
      LBTESTCD,
      LBSPEC,
      LBORRESU,
      out_of_LON_CG_F,
      out_of_LON_CG_M
    ) %>%
    #extract the first number from the "out_of_LON" columns
    mutate(
      #extract the first character and convert to numeric            
      animal_number_F = as.integer(substr(out_of_LON_CG_F, 1, 1)),
      animal_number_M = as.integer(substr(out_of_LON_CG_M, 1, 1)),
      #extract everything after the first character
      suffix_F = substring(out_of_LON_CG_F, 2),
      suffix_M = substring(out_of_LON_CG_M, 2)
      ) %>%
    #remove original columns
    select(-out_of_LON_CG_F, -out_of_LON_CG_M) %>%
    #compute median for numerical results
    group_by(LBTESTCD, LBSPEC, LBORRESU) %>%
    summarise(
      median_animal_number_F = round(median(animal_number_F)),
      median_animal_number_M = round(median(animal_number_M)),
      suffix_F = first(suffix_F),
      suffix_M = first(suffix_M),
      .groups = "drop"
    ) %>%
    #combine the median numbers with the suffix
    mutate(
      out_of_LON_CG_F = paste0(median_animal_number_F, suffix_F),
      out_of_LON_CG_M = paste0(median_animal_number_M, suffix_M)
    ) %>%
    select(LBTESTCD, LBSPEC, LBORRESU, out_of_LON_CG_F, out_of_LON_CG_M)
  #*****************************************************************************
  #*Combine median VCG results and dose group results
  #***************************************************************************** 
  outside_of_HCD <- outside_of_HCD_VCGs %>%
    left_join(outside_of_HCD_doses, by = c("LBTESTCD", "LBSPEC", "LBORRESU"))
  #*****************************************************************************
  #*Extract the HCD range, i.e., the LON (limit of normal)
  #*****************************************************************************
  #*again, the limits of normal remain the same in each iteration, therefore,
  #*only the first limits of normal are extracted
  LON <- resampling_significance_test_results[[1]][[1]] %>%
    select(
      LBTESTCD,
      LBSPEC,
      LBORRESU,
      starts_with("LLN")
    )
  #*****************************************************************************
  #*Check for dose dependency (i.e., does the difference between means rise with
  #*rising dose?)
  #extract VCG mean of each parameter and for each sex for that
  vcg_means_longer <- median_vals_of_vcg_res_means %>%
    select(LBTESTCD, LBSPEC, LBORRESU, CG_F, CG_M) %>%
    pivot_longer(
      -c(LBTESTCD, LBSPEC, LBORRESU),
      names_to = c("trial_set_description", "SEX"),
      names_pattern = "(CG)_(M|F)"
    ) %>%
    rename("VCG_mean" = "value") %>%
    select(-trial_set_description)
  
  #*get the average p-values which are going to be used to delete rows where no
  #*significance is present.
  asterisks <- vcg_sig_res_majority_votes %>%
    # pivot_longer(
    #   -c(LBTESTCD, LBSPEC, LBORRESU),
    #   names_to = c("trial_set_description", "SEX"),
    #   names_pattern = "(CG|.D)_(M|F)"
    # ) %>%
    mutate(significance = if_else(grepl("\\*", result), "*", NA)) #%>%
  # select(-value)
  
  #check for dose dependency
  VCG_assigned_dose_dependency <- median_vals_of_vcg_res_means %>%
    pivot_longer(
      -c(LBTESTCD, LBSPEC, LBORRESU),
      names_to = c("trial_set_description", "SEX"),
      names_pattern = "(CG|.D)_(M|F)"
    ) %>%
    left_join(
      vcg_means_longer,
      by = c("LBTESTCD", "LBSPEC", "LBORRESU", "SEX")
    ) %>%
    left_join(
      asterisks,
      by = c("LBTESTCD", "LBSPEC", "LBORRESU", "trial_set_description", "SEX")
    ) %>%
    mutate(
      #calculate absolute difference between control and each dose group
      difference = abs(value - VCG_mean)
    ) %>%
    group_by(LBTESTCD, LBSPEC, LBORRESU, SEX) %>%
    #sort by dose group
    mutate(
      trial_set_description = factor(
        trial_set_description,
        levels = c("CG", "LD", "MD", "HD")
      )
    ) %>%
    arrange(LBTESTCD, LBSPEC, LBORRESU, SEX, trial_set_description) %>%
    #*check whether the difference in means in rising doses is increasing
    mutate(
      dose_dependency_per_group = difference > lag(difference),
      dose_dependency = if_else(
        any(!dose_dependency_per_group, na.rm = TRUE), "0", "1"
      ),
      #remove dose dependencies where there is no significant difference
      dose_dependency = if_else(
        any(!is.na(significance)), unique(dose_dependency), ""
      )
    ) %>%
    select(LBTESTCD, LBSPEC, LBORRESU, SEX, dose_dependency) %>%
    unique() %>%
    #pivot the structure
    pivot_wider(
      id_cols = c(LBTESTCD, LBSPEC, LBORRESU),
      names_from = c(SEX),
      values_from = dose_dependency
    ) %>%
    rename(
      "dose_dependency_F" = "F",
      "dose_dependency_M" = "M"
    ) %>%
    #add correlations and sexes in different directions in future
    mutate(
      correlations_F = "tbd",
      correlations_M = "tbd",
      sexes_HD_in_diff_directions = "tbd"
    )
  #*****************************************************************************
  #*Concatenate all numeric results into one data frame and assign asterisks
  #*if the corresponding p-value falls below 0.05
  vcg_asterisks <- median_vals_of_vcg_res_means %>%
    pivot_longer(
      -c(LBTESTCD, LBSPEC, LBORRESU),
      names_to = c("trial_set_description", "SEX"),
      names_pattern = "(CG|.D)_(M|F)"
    ) %>%
    #join longer format of average standard deviations
    left_join(
      median_vals_of_vcg_res_sd %>% pivot_longer(
        -c(LBTESTCD, LBSPEC, LBORRESU),
        names_to = c("trial_set_description", "SEX"),
        names_pattern = "(CG|.D)_(M|F)"
      ) %>%
        rename("sd" = "value"),
      by = c("LBTESTCD", "LBSPEC", "LBORRESU", "trial_set_description", "SEX")
    ) %>%
    #*add control group value to check in which direction the significant
    #*differences show
    left_join(
      vcg_means_longer,
      by = c("LBTESTCD", "LBSPEC", "LBORRESU", "SEX")
    ) %>%
    #add the majority votes of the significances
    left_join(
      vcg_sig_res_majority_votes,
      by = c("LBTESTCD", "LBSPEC", "LBORRESU", "trial_set_description", "SEX")
    ) %>%
    #check in which direction the significant difference goes
    mutate(
      direction = case_when(
        grepl("\\*", result) & (value > VCG_mean) ~ "(+)",
        grepl("\\*", result) & (value < VCG_mean) ~ "(-)"
      ),
      #turn NA into empty values
      direction = if_else(is.na(direction), "", direction),
      #*turn mean value into string and add standard deviation as well as the
      #*significance asterisks, the direction, and the percentage of the
      #*majority votes (if not control group)
      param_mean_string = case_when(
        is.nan(value) ~ "",
        !is.nan(value) & trial_set_description == "CG" ~ paste(
          as.character(value), "\u00B1", as.character(sd), sep = " "
        ),
        TRUE ~ paste(
          as.character(value), "\u00B1", as.character(sd), result, direction,
          sep = " "
        )
      )
    ) %>%
    # #remove everything from control which has to do with significance tests
    # #Get only the standard deviation values. Drop everything before "\u00B1"
    # function(x){sub(".\u00B1(.+)", "\\1", x) %>% parse_number()}
    select(
      LBTESTCD, LBSPEC, LBORRESU, trial_set_description, SEX, param_mean_string
    ) %>%
    #pivot back into wide structure
    pivot_wider(
      id_cols = c(LBTESTCD, LBSPEC, LBORRESU),
      names_from = c(trial_set_description, SEX),
      values_from = param_mean_string
    )
  #*****************************************************************************
  #*Concatenate all results into one data frame
  VCG_logic_results <- vcg_asterisks %>%
    left_join(
      VCG_assigned_dose_dependency,
      by = c("LBTESTCD", "LBSPEC", "LBORRESU")
    ) %>%
    #add the limits of normal to the results table
    left_join(LON, by = c("LBTESTCD", "LBSPEC", "LBORRESU")) %>%
    #Add the count on how many animals were ouside of HCD reference range
    left_join(
      outside_of_HCD,
      by = c("LBTESTCD", "LBSPEC", "LBORRESU")
    )
  #*****************************************************************************
  #*****************************************************************************
  #*Collect the effect sizes of all iterations----
  #*****************************************************************************
  effsize_resampling_results <- lapply(
    resampling_significance_test_results,
    function(x){
      x[[2]]
    }
  )
  effsize_resampling_results_flatten <- effsize_resampling_results %>%
    imap_dfr(~ .x %>% as_tibble(), .id = "iteration")
  #*****************************************************************************
  #*****************************************************************************
  #*Collect power of all iterations----
  #*****************************************************************************
  #*Collect the effect sizes of all iterations
  power_resampling_results <- lapply(
    resampling_significance_test_results,
    function(x){
      x[[3]]
    }
  )
  power_resampling_results_flatten <- power_resampling_results %>%
    imap_dfr(~ .x %>% as_tibble(), .id = "iteration")
  #*****************************************************************************
  
  #*return results from function
  return(list(
    VCG_logic_results,
    effsize_resampling_results_flatten,
    power_resampling_results_flatten
  ))
}
#End of function
#*******************************************************************************
#*******************************************************************************
#*Generate legacy study original results
#*******************************************************************************
#*Use the "study_results_table" function for this
print("CCG_legacy_study_A")
Sys.time()
CCG_legacy_study_A <- tox_statistical_test(
  studydata = legacy_study_A,
  refdata = HCD_study_A,
  power_calculations = TRUE
)
Sys.time()
#*******************************************************************************
#*******************************************************************************
#*Apply function to legacy study where CCGs are replaced by VCGs via resampling
#*******************************************************************************
print("VCG_legacy_study_A")
VCG_legacy_study_A <- resampling_VCG(
  legacy_study = legacy_study_A,
  HCD = HCD,
  iterations = 100,
  replacement_aim = 100,
  sensitivity_analyses = TRUE
)

print("VCG_legacy_study_A_imp_median")
VCG_legacy_study_A_imp_median <- resampling_VCG(
  legacy_study = legacy_study_A,
  HCD = HCD_imputed_median,
  iterations = 100,
  replacement_aim = 100
)
print("VCG_legacy_study_A_imp_rs")
VCG_legacy_study_A_imp_rs <- resampling_VCG(
  legacy_study = legacy_study_A,
  HCD = HCD_imputed_rs,
  iterations = 100,
  replacement_aim = 100
)
print("VCG_legacy_study_A_imp_pmm")
VCG_legacy_study_A_imp_pmm <- resampling_VCG(
  legacy_study = legacy_study_A,
  HCD = HCD_imputed_pmm,
  iterations = 100,
  replacement_aim = 100
)
#*******************************************************************************
#*******************************************************************************
#*Compare studies where CCGs from study A and C are exchanged. How large is
#*the reproducibility percantage there?
#*******************************************************************************
print("CCG_legacy_study_A_replaced_C")
CCG_C_legacy_study_A <- tox_statistical_test(
  studydata = legacy_study_A %>%
    filter(trial_set_description != "CG") %>%
    bind_rows(legacy_study_C %>% filter(trial_set_description == "CG")),
  refdata = HCD_study_A,
  power_calculations = FALSE
) %>%
  suppressMessages()

print("CCG_legacy_study_C_replaced_A")
CCG_A_legacy_study_C <- tox_statistical_test(
  studydata = legacy_study_C %>%
    filter(trial_set_description != "CG") %>%
    bind_rows(legacy_study_A %>% filter(trial_set_description == "CG")),
  refdata = HCD_study_C,
  power_calculations = FALSE
) %>%
  suppressMessages()
#*******************************************************************************
#*******************************************************************************
#*Generate legacy study original results
#*******************************************************************************
#*Use the "study_results_table" function for this
print("CCG_legacy_study_B")
CCG_legacy_study_B <- tox_statistical_test(
  studydata = legacy_study_B,
  refdata = HCD_study_B,
  power_calculations = TRUE
) %>%
  suppressMessages()

#study B recovery groups
print("CCG_legacy_study_B_recovery")
CCG_legacy_study_B_recovery <- tox_statistical_test(
  studydata = legacy_study_B_recovery,
  refdata = HCD_study_B,
  power_calculations = TRUE
)

#*******************************************************************************
#*******************************************************************************
#*Apply function to legacy study where CCGs are replaced by VCGs via resampling
#*******************************************************************************
print("VCG_legacy_study_B")
VCG_legacy_study_B <- resampling_VCG(
  legacy_study = legacy_study_B,
  HCD = HCD,
  iterations = 100,
  replacement_aim = 100,
  plottitle = "Legacy study B",
  sensitivity_analyses = TRUE
)
# 
saveRDS(
  CCG_legacy_study_B,
  paste0(rootpath, "/data/Derived/CCG_legacy_study_B.RData")
)
saveRDS(
  VCG_legacy_study_B,
  paste0(rootpath, "/data/Derived/VCG_legacy_study_B.RData")
)

#study B recovery groups
print("VCG_legacy_study_B_recovery")
VCG_legacy_study_B_recovery <- resampling_VCG(
  legacy_study = legacy_study_B_recovery,
  HCD = HCD_recovery,
  iterations = 100,
  replacement_aim = 100,
  match_by_initbw = FALSE,
  plottitle = "Legacy study B recovery",
  sensitivity_analyses = TRUE
)

saveRDS(
  CCG_legacy_study_B_recovery,
  paste0(rootpath, "/data/Derived/CCG_legacy_study_B_recovery.RData")
)
saveRDS(
  VCG_legacy_study_B_recovery,
  paste0(rootpath, "/data/Derived/VCG_legacy_study_B_recovery.RData")
)

print("VCG_legacy_study_B_imp_median")
VCG_legacy_study_B_imp_median <- resampling_VCG(
  legacy_study = legacy_study_B,
  HCD = HCD_imputed_median,
  iterations = 100,
  replacement_aim = 100,
  plottitle = "Legacy study B imp median"
)
print("VCG_legacy_study_B_imp_rs")
VCG_legacy_study_B_imp_rs <- resampling_VCG(
  legacy_study = legacy_study_B,
  HCD = HCD_imputed_rs,
  iterations = 100,
  replacement_aim = 100,
  plottitle = "Legacy study B imp rs"
)
print("VCG_legacy_study_B_imp_pmm")
VCG_legacy_study_B_imp_pmm <- resampling_VCG(
  legacy_study = legacy_study_B,
  HCD = HCD_imputed_pmm,
  iterations = 100,
  replacement_aim = 100,
  plottitle = "Legacy study B imp pmm"
)

#*******************************************************************************
#*******************************************************************************
#*Generate legacy study original results
#*******************************************************************************
#*Use the "study_results_table" function for this
print("CCG_legacy_study_C")
CCG_legacy_study_C <- tox_statistical_test(
  studydata = legacy_study_C,
  refdata = HCD_study_C,
  power_calculations = TRUE

) %>%
  suppressMessages()
#*******************************************************************************
#*******************************************************************************
#*Apply function to legacy study where CCGs are replaced by VCGs via resampling
#*******************************************************************************
print("VCG_legacy_study_C")
VCG_legacy_study_C <- resampling_VCG(
  legacy_study = legacy_study_C,
  HCD = HCD,
  iterations = 100,
  replacement_aim = 100,
  plottitle = "Legacy study C",
  sensitivity_analyses = TRUE
)

saveRDS(
  CCG_legacy_study_C,
  paste0(rootpath, "/data/Derived/CCG_legacy_study_C.RData")
)
saveRDS(
  VCG_legacy_study_C,
  paste0(rootpath, "/data/Derived/VCG_legacy_study_C.RData")
)
print("VCG_legacy_study_C_imp_median")
VCG_legacy_study_C_imp_median <- resampling_VCG(
  legacy_study = legacy_study_C,
  HCD = HCD_imputed_median,
  iterations = 100,
  replacement_aim = 100,
  plottitle = "Legacy study C imp median"
)
print("VCG_legacy_study_C_imp_rs")
VCG_legacy_study_C_imp_rs <- resampling_VCG(
  legacy_study = legacy_study_C,
  HCD = HCD_imputed_rs,
  iterations = 100,
  replacement_aim = 100,
  plottitle = "Legacy study C imp rs"
)
print("VCG_legacy_study_C_imp_pmm")
VCG_legacy_study_C_imp_pmm <- resampling_VCG(
  legacy_study = legacy_study_C,
  HCD = HCD_imputed_pmm,
  iterations = 100,
  replacement_aim = 100,
  plottitle = "Legacy study C imp pmm"
)
#*******************************************************************************
#*******************************************************************************
#*Generate legacy study original results
#*******************************************************************************
#*Use the "study_results_table" function for this
print("CCG_legacy_study_NHP")
CCG_legacy_study_NHP <- tox_statistical_test(
  studydata = legacy_study_NHP,
  refdata = HCD_study_C
)
#*******************************************************************************
#*******************************************************************************
#*Save results as R-files for fast recovery of results. Uncomment if needed----
#*******************************************************************************
saveRDS(
  CCG_legacy_study_A,
  paste0(rootpath, "/data/Derived/CCG_legacy_study_A.RData")
  )
saveRDS(
  VCG_legacy_study_A,
  paste0(rootpath, "/data/Derived/VCG_legacy_study_A.RData")
)
saveRDS(
  VCG_legacy_study_A_imp_median,
  paste0(rootpath, "/data/Derived/VCG_legacy_study_A_imp_median.RData")
)
saveRDS(
  VCG_legacy_study_A_imp_rs,
  paste0(rootpath, "/data/Derived/VCG_legacy_study_A_imp_rs.RData")
)
saveRDS(
  VCG_legacy_study_A_imp_pmm,
  paste0(rootpath, "/data/Derived/VCG_legacy_study_A_imp_pmm.RData")
)

saveRDS(
  CCG_legacy_study_B,
  paste0(rootpath, "/data/Derived/CCG_legacy_study_B.RData")
)
saveRDS(
  VCG_legacy_study_B,
  paste0(rootpath, "/data/Derived/VCG_legacy_study_B.RData")
)
saveRDS(
  CCG_legacy_study_B_recovery,
  paste0(rootpath, "/data/Derived/CCG_legacy_study_B_recovery.RData")
)
saveRDS(
  VCG_legacy_study_B_recovery,
  paste0(rootpath, "/data/Derived/VCG_legacy_study_B_recovery.RData")
)
saveRDS(
  VCG_legacy_study_B_imp_median,
  paste0(rootpath, "/data/Derived/VCG_legacy_study_B_imp_median.RData")
)
saveRDS(
  VCG_legacy_study_B_imp_rs,
  paste0(rootpath, "/data/Derived/VCG_legacy_study_B_imp_rs.RData")
)
saveRDS(
  VCG_legacy_study_B_imp_pmm,
  paste0(rootpath, "/data/Derived/VCG_legacy_study_B_imp_pmm.RData")
)

saveRDS(
  CCG_legacy_study_C,
  paste0(rootpath, "/data/Derived/CCG_legacy_study_C.RData")
)
saveRDS(
  VCG_legacy_study_C,
  paste0(rootpath, "/data/Derived/VCG_legacy_study_C.RData")
)
saveRDS(
  VCG_legacy_study_C_imp_median,
  paste0(rootpath, "/data/Derived/VCG_legacy_study_C_imp_median.RData")
)
saveRDS(
  VCG_legacy_study_C_imp_rs,
  paste0(rootpath, "/data/Derived/VCG_legacy_study_C_imp_rs.RData")
)
saveRDS(
  VCG_legacy_study_C_imp_pmm,
  paste0(rootpath, "/data/Derived/VCG_legacy_study_C_imp_pmm.RData")
)