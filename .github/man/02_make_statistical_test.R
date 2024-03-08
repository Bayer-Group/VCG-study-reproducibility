#* 02_make_statistical_test.R

#*Create a function which does the following:
#*  - Calculate a statistical test with respect to the endpoint of interest
#*    for each sex separately using the individual values
#*  - Summarize results in a table and return the table

#Check if packages are installed
require(tidyverse)
require(PMCMRplus)
require(rstatix)
require(effsize)
#*******************************************************************************
#*Source the function for sensitivity analysis
source(paste0(rootpath, "/man/09_sensitivity_analysis.R"))
#*Read the table listing the needed statistical test for evaluating the results
#*with respect to each endpoint. This table need to consist of two columns:
#*The first contains the LBTESTCD (i.e., the coding for the endpoint).
#*The second column needs to contain the name of the respective statistical
#*test.
#*In addition to that, the Cliff's delta effect size and a post-hoc power
#*calculation is performed in order to gain an understanding how large the
#*effect needs to be in order to have a significant effect.
#*Currently, the following tests are included:
#*  - Dunnett's test
#*  - Dunnett's test with heterogenicity correction
#*  - Bonferroni adjusted Mann Whitney U-test
#*  - Dunnett's test after logarithmic transformation
#*  - NA:: no statistical test
#*  - Additionally, calculate Cliff's delta effect size and
#*    post-hoc power analysis.
#*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#*Warning: in the statistical list, there's a test, calles "het_dunn". With
#*this, the "Dunn test" is not meant but the "exact heterogenous Dunnett's"
#*test, i.e. the Dunnett's T3 test!
#*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

significance_test_table <- read.csv(
  paste0(rootpath, "/data/Original/significance_test_list.csv"),
  fileEncoding="UTF-8-BOM"
  )

tox_statistical_test <- function(
  studydata,
  refdata,
  power_calculations = FALSE
  )
{
  #*****************************************************************************
  #*****************************************************************************
  #*Create a function which formats numeric values to have three digits and adds
  #*tailing zeros, if needed
  #*****************************************************************************
  signif3 <- function(x) {
    # In case of empty values, return NA
    if(is.na(x) | x == "") {
      return(NA_character_)
    }
    
    #*Round to 3 significant figures and transform "," to "." and
    #*remove whitespace, if necessary
    digits3 <- x %>%
      gsub(",", ".", .) %>%
      str_trim() %>%
      as.numeric(.) %>%
      signif(., 3)
    
    # Count the number of digits in the present data
    num_digits <- ifelse(digits3 == 0, 1, floor(log10(abs(digits3))) + 1)
    
    # Check the vectorized conditions using sapply
    formatted_numbers <- sapply(1:length(x), function(i) {
      if (num_digits[i] == 1) {
        return(sprintf("%.2f", digits3[i]))
      } else if (num_digits[i] == 2) {
        return(sprintf("%.1f", digits3[i]))
      } else {
        return(as.character(digits3[i]))
      }
    })
    
    return(formatted_numbers)
  }
  
  format_number_vectorized <- Vectorize(signif3)
  #*****************************************************************************
  #*****************************************************************************
  #*Check how many individuals (per does group) are outside of HCD----
  #*****************************************************************************
  #*HCD needs to be provided in the "refdata" variable
  outside_of_HCD_LON_range <- merge(
    studydata,
    refdata %>%
      select(
        SEX, LBTESTCD, LBSPEC, LBORRESU, `Percentile2.5th`, `Percentile97.5th`
        ),
    by = c("SEX", "LBTESTCD", "LBSPEC", "LBORRESU")
    ) %>%
    mutate(
      #turn data into numerics, remove whitespace, and replace "," with "."
      across(
        c(`Percentile2.5th`, `Percentile97.5th`),
        function(x){as.numeric(gsub(",", ".", str_trim(x)))}
        ),
      outside_LON_range = case_when(
        !between(LBORRES, `Percentile2.5th`, `Percentile97.5th`) ~ "1",
        between(LBORRES, `Percentile2.5th`, `Percentile97.5th`) ~ "0"
      )
    ) %>%
    group_by(
      trial_set_description, SEX, LBTESTCD, LBSPEC, LBORRESU
      ) %>%
    #summarize the results per sex and per parameter
    count(outside_LON_range) %>%
    rename("n_outside_LON_range" = n) %>%
    filter(outside_LON_range == "1")
  
  #*Create a data frame containing the reference ranges, i.e., the lower limit
  #*of normal (LLN) and the upper limit of normal (ULN)
  LONs <- refdata %>%
    select(
      SEX, LBTESTCD, LBSPEC, LBORRESU, `Percentile2.5th`, `Percentile97.5th`
    ) %>%
  mutate(
    #define the lower limit of normal and the upper limit of normal
    `LLN/ULN` = paste0(
      map_chr(`Percentile2.5th`, ~signif3(.x)),
      "/",
      map_chr(`Percentile97.5th`, ~signif3(.x))
    )
  ) %>%
    select(-`Percentile2.5th`, -`Percentile97.5th`)
  
  #merge individual parameters of study with reference data
  studydata_with_reference_long <- merge(
    studydata,
    refdata %>%
      select(
        SEX, LBTESTCD, LBSPEC, LBORRESU, `Percentile2.5th`, `Percentile97.5th`
      ),
    by = c("SEX", "LBTESTCD", "LBSPEC", "LBORRESU")
  ) %>%
    group_by(trial_set_description, SEX, LBTESTCD, LBSPEC, LBORRESU) %>%
    #summarize the results per sex and per parameter
    summarise(population = n()) %>%
    #add limit of normal
    left_join(LONs, by = c("SEX", "LBTESTCD", "LBSPEC", "LBORRESU")) %>%
    #delete duplicates which appear for some reason
    unique() %>%
    #add a number on how many animals are outside of the limits of normal
    left_join(
      outside_of_HCD_LON_range %>%
        select(
          trial_set_description,
          SEX,
          LBTESTCD,
          LBSPEC,
          LBORRESU,
          n_outside_LON_range
          ),
      by = c(
        "trial_set_description",
        "SEX",
        "LBTESTCD",
        "LBSPEC",
        "LBORRESU"
      )
    ) %>%
    mutate(
      n_outside_LON_range = coalesce(n_outside_LON_range, 0L),
      n_outside_LON_range = paste0(n_outside_LON_range, " out of ", population)
    ) %>%
    select(-population) %>%
    ungroup() %>%
    suppressMessages()
  
  #Pivot study data with reference into a wide forma
  studydata_with_reference <- studydata_with_reference_long %>%
    select(-`LLN/ULN`) %>%
    pivot_wider(
      names_from = c(trial_set_description, SEX),
      values_from = n_outside_LON_range,
      names_prefix = "out_of_LON_"
    ) %>% 
    left_join(
      studydata_with_reference_long %>%
        select(-n_outside_LON_range, -trial_set_description) %>%
        unique() %>%
        pivot_wider(
          names_from = c(SEX),
          values_from = `LLN/ULN`,
          names_prefix = "LLN/ULN_"
        ),
      by = c("LBTESTCD", "LBSPEC", "LBORRESU")
      )
  #*****************************************************************************
  #*Create a warning message if a LBTESTCD is measured which is not present in 
  #*the significance test list (ignore body weight parameters)
  studydata_LBTESTCD <- studydata %>%
    select(LBTESTCD, LBSPEC, LBORRESU) %>%
    unique()
  
  LBTESTCD_not_in_sig_test_list <- studydata_LBTESTCD %>%
    #get only values which are not present in significance test table
    anti_join(
      significance_test_table %>% select(LBTESTCD, LBSPEC),
      by = c("LBTESTCD", "LBSPEC")
      ) %>%
    filter(
      !grepl("BW_", LBTESTCD),
      !grepl("WEIGHT_", LBTESTCD),
      !grepl("FC_", LBTESTCD),
      !grepl("WC_", LBTESTCD),
      !grepl("BG_", LBTESTCD)
    )
  if(nrow(LBTESTCD_not_in_sig_test_list) > 0){
    warning("The following LBTESTCD parameters are not captured in the significance test list and will be omitted:\n - ",
            paste(
              paste(
                LBTESTCD_not_in_sig_test_list$LBTESTCD,
                LBTESTCD_not_in_sig_test_list$LBSPEC,
                sep = ","
                ),
              collapse = "\n - "
              )
    )
  }
  #*****************************************************************************  
  #*Create a warning message if some parameters are missing a control group
  LBTESTCD_missing_CG <- studydata %>%
    mutate(
      LBTESTCD_LBSPEC_LBORRESU_SEX = paste(
        LBTESTCD, LBSPEC, LBORRESU, SEX, sep = "_"
      )
    ) %>%
    group_by(LBTESTCD_LBSPEC_LBORRESU_SEX) %>%
    filter(
      !LBTESTCD_LBSPEC_LBORRESU_SEX %in% (
        studydata %>%
          mutate(
            LBTESTCD_LBSPEC_LBORRESU_SEX = paste(
              LBTESTCD, LBSPEC, LBORRESU, SEX, sep = "_"
              )
            ) %>%
          filter(trial_set_description == "CG") %>%
          pull(LBTESTCD_LBSPEC_LBORRESU_SEX) %>%
          unique()
      )
    ) %>%
    pull(LBTESTCD_LBSPEC_LBORRESU_SEX) %>%
    unique()
  
  if(length(LBTESTCD_missing_CG) != 0){
    warning("The following LBTESTCD parameters are missing a control group and will be omitted:\n - ",
            paste(LBTESTCD_missing_CG, collapse = "\n - ")
    )
  }
  #*****************************************************************************
  #*Create a warning message if only one animal value is present in a group.
  #*By this, no statistical calculation can be made.
  LBTESTCD_only_one_CG <- studydata %>%
    mutate(
      LBTESTCD_LBSPEC_LBORRESU_SEX = paste(
        LBTESTCD, LBSPEC, LBORRESU, SEX, sep = "_"
      )
      ) %>%
    filter(
      !LBTESTCD_LBSPEC_LBORRESU_SEX %in% LBTESTCD_missing_CG,
      trial_set_description == "CG"
      ) %>%
    group_by(LBTESTCD_LBSPEC_LBORRESU_SEX, trial_set_description) %>%
    summarise(population = n()) %>%
    filter(population <= 1) %>%
    pull(LBTESTCD_LBSPEC_LBORRESU_SEX) %>%
    unique() %>%
    suppressMessages()
  
  if(length(LBTESTCD_only_one_CG) != 0){
    warning("The following LBTESTCD parameters have only one control group value and will be omitted:\n - ",
            paste(LBTESTCD_only_one_CG, collapse = "\n - ")
    )
  } 

  #*Keep only LBTESTCD which were used in the legacy study
  LBTESTCD_in_legacy_study <- significance_test_table %>%
    #add body weight, body weight gain, and organ weight to this list
    rbind(
      tibble(
        studydata_LBTESTCD %>%
          select(-LBORRESU) %>% #drop units here
          filter(
            grepl("BW_", LBTESTCD) |
            grepl("WEIGHT_", LBTESTCD) |
            grepl("BG_", LBTESTCD)
            ),
        significance_test = "dunnett"
      )
    ) %>%
    #add food and water consumption to this list
    rbind(
      tibble(
        studydata_LBTESTCD %>%
          select(-LBORRESU) %>% #drop units here
          filter(LBSPEC == "FW"),
        significance_test = "u_test"
      )
    ) %>%
    inner_join(
      studydata %>% select(LBTESTCD, LBSPEC) %>% unique(),
      by = c("LBTESTCD", "LBSPEC")
    )
  
  #add significance test to the studydata which needs to be used
  studydata_with_significance_test <- merge(
    studydata,
    LBTESTCD_in_legacy_study, by = c("LBTESTCD", "LBSPEC")
    ) %>%
    mutate(
      LBTESTCD_LBSPEC_LBORRESU_SEX = paste(
        LBTESTCD, LBSPEC, LBORRESU, SEX, sep = "_"
        )
      ) %>%
    filter(
      !LBTESTCD_LBSPEC_LBORRESU_SEX %in% LBTESTCD_missing_CG,
      !LBTESTCD_LBSPEC_LBORRESU_SEX %in% LBTESTCD_only_one_CG
        ) %>%
    select(-LBTESTCD_LBSPEC_LBORRESU_SEX)
  
  #*Turn data frame into list divided into respective parameters and sex
  #*define this data frame globally, if a sensitivity analysis needs to be done
  if(power_calculations){
    studydata_list <<- split(
      studydata_with_significance_test %>%
        select(trial_set_description, LBORRES, significance_test),
      paste(
        studydata_with_significance_test$LBTESTCD,
        studydata_with_significance_test$LBSPEC,
        studydata_with_significance_test$LBORRESU,
        studydata_with_significance_test$SEX,
        sep = "_"
      )
    )
  }else{
    studydata_list <- split(
      studydata_with_significance_test %>%
        select(trial_set_description, LBORRES, significance_test),
      paste(
        studydata_with_significance_test$LBTESTCD,
        studydata_with_significance_test$LBSPEC,
        studydata_with_significance_test$LBORRESU,
        studydata_with_significance_test$SEX,
        sep = "_"
      )
    )
  }

  #*****************************************************************************
  #*****************************************************************************
  #*Excecute significance test based on significance test list----
  #*****************************************************************************
  statistical_test_results_no_bg <- lapply(
    studydata_list, function(parameter){
      if(
        is.na(unique(parameter$significance_test))
        ){
        #no test----
        # print("no test")
      }else if(
        #*A prerequisite of the Dunnett's test is that at least 3 groups are
        #*present. If this is not the case (i.e., if only 2 groups are present)
        #*a t-test is calculated instead.
        unique(parameter$significance_test) == "dunnett" &
        n_distinct(parameter$trial_set_description) == 2
      ){
        # print("t-test")
        parameter %>%
          droplevels() %>%
          t_test(
            data = .,
            LBORRES ~ trial_set_description,
            ref.group = "CG",
            p.adjust.method = "none",
            var.equal = TRUE
            ) %>%
          select(group2, p) %>%
          as_tibble() %>%
          rename(
            "trial_set_description" = "group2",
            "p_value" = "p"
            ) %>%
          mutate(
            significance = case_when(
              p_value > .05 ~ NA_character_,
              p_value <= .05 & p_value > .01 ~ "*",
              p_value <= .01 & p_value > .001 ~ "**",
              p_value <= .001 ~ "***"
            )
          )
      }else if(
        unique(parameter$significance_test) == "dunnett" &
        n_distinct(parameter$trial_set_description) > 2
        ){
        #*Dunnett's exact homogenous test----
        # print("dunnett")
        dunnettTest(
          parameter$LBORRES,
          parameter$trial_set_description
        )$p.value %>%
        as_tibble(rownames = "trial_set_description") %>%
        rename("p_value" = "CG") %>%
        mutate(
          significance = case_when(
            p_value > .05 ~ NA_character_,
            p_value <= .05 & p_value > .01 ~ "*",
            p_value <= .01 & p_value > .001 ~ "**",
            p_value <= .001 ~ "***"
          )
        )
      }
      else if(
        #*A prerequisite of the Dunnett's test is that at least 3 groups are
        #*present. If this is not the case (i.e., if only 2 groups are present)
        #*a t-test is calculated instead.
        unique(parameter$significance_test) == "het_dunn" &
        n_distinct(parameter$trial_set_description) == 2
      ){
        # print("Welch-adapted t-test")
        parameter %>%
          droplevels() %>%
          t_test(
            data = .,
            LBORRES ~ trial_set_description,
            ref.group = "CG",
            p.adjust.method = "none",
            var.equal = FALSE
          ) %>%
          select(group2, p) %>%
          as_tibble() %>%
          rename(
            "trial_set_description" = "group2",
            "p_value" = "p"
          ) %>%
          mutate(
            significance = case_when(
              p_value > .05 ~ NA_character_,
              p_value <= .05 & p_value > .01 ~ "*",
              p_value <= .01 & p_value > .001 ~ "**",
              p_value <= .001 ~ "***"
            )
          )
        }else if(
          unique(parameter$significance_test) == "het_dunn" &
          n_distinct(parameter$trial_set_description) > 2
          ){
        #Dunnett's exact heterogenous test----
        # print("hetdunn")
        dunnettT3Test(
          parameter$LBORRES,
          parameter$trial_set_description
        )$p.value %>%
          as_tibble(rownames = "trial_set_description") %>%
          select("trial_set_description", "CG") %>%
          rename("p_value" = "CG") %>%
          mutate(
            significance = case_when(
              p_value > .05 ~ NA_character_,
              p_value <= .05 & p_value > .01 ~ "*",
              p_value <= .01 & p_value > .001 ~ "**",
              p_value <= .001 ~ "***"
            )
          )
      }
      else if(
        #*A prerequisite of the Dunnett's test is that at least 3 groups are
        #*present. If this is not the case (i.e., if only 2 groups are present)
        #*a t-test is calculated instead.
        unique(parameter$significance_test) == "log_trans_dunnett" &
        n_distinct(parameter$trial_set_description) == 2
      ){
        # print("log trans t-test")
        parameter %>%
          droplevels() %>%
          mutate(LBORRES = log10(LBORRES)) %>%
          t_test(
            data = .,
            LBORRES ~ trial_set_description,
            ref.group = "CG",
            p.adjust.method = "none",
            var.equal = FALSE
          ) %>%
          select(group2, p) %>%
          as_tibble() %>%
          rename(
            "trial_set_description" = "group2",
            "p_value" = "p"
          ) %>%
          mutate(
            significance = case_when(
              p_value > .05 ~ NA_character_,
              p_value <= .05 & p_value > .01 ~ "*",
              p_value <= .01 & p_value > .001 ~ "**",
              p_value <= .001 ~ "***"
            )
          )
      }else if(
        unique(parameter$significance_test) == "log_trans_dunnett" &
        n_distinct(parameter$trial_set_description) > 2
        ){
        #Dunnett's exact homogenous test after logarhitmic transformation----
        # print("logtrandunn")
        dunnettTest(
          log10(parameter$LBORRES),
          parameter$trial_set_description
        )$p.value %>%
          as_tibble(rownames = "trial_set_description") %>%
          rename("p_value" = "CG") %>%
          mutate(
            significance = case_when(
              p_value > .05 ~ NA_character_,
              p_value <= .05 & p_value > .01 ~ "*",
              p_value <= .01 & p_value > .001 ~ "**",
              p_value <= .001 ~ "***"
            )
          )
      }else if(unique(parameter$significance_test) == "u_test"){
        #*An important prerequesite of the Bonferroni-adjusted Mann Whitney U
        #*test is that at least 6 data points per group are present. If this
        #*is not given, no test is performed and a warning message is returned.
        if(
          any(
            parameter %>%
            group_by(trial_set_description) %>%
            summarize(population = n()) %>%
            pull(population) <= 6
          )
        ){
          # print("no u-test")
          warning(
            paste0(
            "At least one group in the parameter", deparse(substitute((parameter))), " has less than 6 values. No U-test is calculated"
            )
          )
          tibble(
            trial_set_description = parameter %>%
              pull(trial_set_description) %>%
              unique(),
            p_value = NaN,
            significance = "(!)"
          )
        }else{
          #Bonferroni adjusted Mann Whitney U-test---
          # print("Bonferroni u-test")
          parameter %>%
            droplevels() %>%
            wilcox_test(LBORRES ~ trial_set_description, ref.group = "CG") %>%
            adjust_pvalue(method = "bonferroni") %>%
            add_significance("p.adj") %>%
            select(group2, p.adj) %>%
            rename("trial_set_description" = "group2", "p_value" = "p.adj") %>%
            mutate(
              significance = case_when(
                p_value > .05 ~ NA_character_,
                p_value <= .05 & p_value > .01 ~ "*",
                p_value <= .01 & p_value > .001 ~ "**",
                p_value <= .001 ~ "***"
              )
            )
        }
      }
    }
  )
  #*****************************************************************************
  statistical_test_results <- c(
    statistical_test_results_no_bg
  )
  
  #sort names alphabetically
  statistical_test_results <- statistical_test_results[
    order(names(statistical_test_results))
    ]
  #*****************************************************************************
  #*****************************************************************************
  #*Calculate effect size----
  #*****************************************************************************
  #*In addition, calculate effect size (Cliff's delta as an non-parameteric
  #*effect size)
  
  #*In some cases, the number of groups may vary (e.g., due to mortality). To
  #*account for that, the number of groups per parameter is calculated first,
  #*and then the effect size is calculated between the control and each dose
  #*group separately.
  cliff_delta_calculator <- function(endpoint, dose_group){
    #get dose groups present in parameter
    dose_groups <- endpoint %>%
      filter(trial_set_description != "CG") %>%
      droplevels()
    
    #get control group
    control_group <- endpoint %>%
      filter(trial_set_description == "CG") %>%
      droplevels()
    
    #Calculate effect size between control and dose group.
    estimate <- cliff.delta(
      endpoint %>%
        droplevels()%>%
        filter(trial_set_description == "CG") %>%
        pull(LBORRES),
      endpoint %>%
        droplevels() %>%
        filter(trial_set_description == dose_group) %>%
        pull(LBORRES)
    )$estimate
    
    return(tibble(trial_set_description = dose_group, effsize = estimate))
  }
  #*****************************************************************************
  #*****************************************************************************
  #*Compute post-hoc power----
  #*****************************************************************************
  if(power_calculations){
    test_power <- sensitivity_analysis()
    #*flatten results of post-hoc power
    test_power_flatten <- test_power %>%
      imap_dfr(~ .x %>% as_tibble(), .id = "LBTESTCD_LBSPEC_LBORRESU_SEX") %>%
      rename(power = value) %>%
      #split LBTESTCD, LBSPEC, and SEX into three different columns
      #*note that this splitting looks for underscores in the values which doesn't
      #*work for the BW/OM parameters, as they have one underscore too much.
      #*In this case, we remove the "BW_" prefix and add it in the end again.
      mutate(
        #store prefixes in separate column
        prefix = if_else(
          grepl(
            "^(BW_|WEIGHT_|OWBW_|FC_|WC_|BG_)", LBTESTCD_LBSPEC_LBORRESU_SEX
          ), 
          gsub("^(.*?_).*", "\\1", LBTESTCD_LBSPEC_LBORRESU_SEX), ""),
        new_LBTESTCD_LBSPEC_LBORRESU_SEX = if_else(
          prefix != "", 
          gsub("^(.*?_)(.*)", "\\2", LBTESTCD_LBSPEC_LBORRESU_SEX), 
          LBTESTCD_LBSPEC_LBORRESU_SEX
        )
      ) %>%
      #separate columns
      separate(
        new_LBTESTCD_LBSPEC_LBORRESU_SEX,
        c("LBTESTCD", "LBSPEC", "LBORRESU", "SEX"),
        sep = "_",
        remove = FALSE
      ) %>%
      #add prefixes again
      mutate(
        LBTESTCD = if_else(prefix != "", paste0(prefix, LBTESTCD), LBTESTCD)
      ) %>%
      #remove columns
      select(
        -prefix, -LBTESTCD_LBSPEC_LBORRESU_SEX, -new_LBTESTCD_LBSPEC_LBORRESU_SEX
      )
  }else{test_power_flatten <- NULL}
  #*****************************************************************************
  #Apply function on the parameters of the selected study
  study_effect_size <- lapply(
    studydata_list,
    function(parameter){
      dose_group_names <- parameter %>%
        pull(trial_set_description) %>%
        droplevels() %>%
        unique()
      
      effsizes <- lapply(
        dose_group_names,
        cliff_delta_calculator,
        endpoint = parameter
      )
      
      results <- do.call(rbind, effsizes)
      
      return(results)
    }
  ) %>%
    #*Suppress warning message "The samples are fully disjoint, using
    #*approximate Confidence Interval estimation". This message appears if there
    #*is no overlap between one group and another.
    suppressWarnings()
 #******************************************************************************
 #******************************************************************************
 #*Unlist the statistical results----
 #******************************************************************************
 #*Create a data frame where each row represents a tested parameter and each
 #*column is the respective dose group, shown with its mean value and whether
 #*there's a significant difference
  flatten_test_results <- statistical_test_results %>%
    imap_dfr(~ .x %>% as_tibble(), .id = "LBTESTCD_LBSPEC_LBORRESU_SEX")
  
  #*Extract the parameter where no significance test were done and append them
  #*to the results data frame
  stat_res_names_not_done <- tibble(
    LBTESTCD_LBSPEC_LBORRESU_SEX = rep(
      names(statistical_test_results), each = 3
      ),
    p_value = NaN,
    trial_set_description = rep(
      c("LD", "MD", "HD"), length(statistical_test_results)
      ),
    significance = NA_character_
    ) %>%
    filter(
      !LBTESTCD_LBSPEC_LBORRESU_SEX %in%
        flatten_test_results$LBTESTCD_LBSPEC_LBORRESU_SEX
      )
  
  #Append empty values to flatten test results
  test_results_tibble <- rbind(
    flatten_test_results,
    stat_res_names_not_done
    ) %>%
    arrange(LBTESTCD_LBSPEC_LBORRESU_SEX) %>%
    #split LBTESTCD, LBSPEC, LBORRESU, and SEX into three different columns
    #*note that this splitting looks for underscores in the values which doesn't
    #*work for the BW/OM parameters, as they have one underscore too much.
    #*In this case, we remove the "BW_" prefix and add it in the end again.
    mutate(
      #store prefixes in separate column
      prefix = if_else(
        grepl("^(BW_|WEIGHT_|OWBW_|FC_|WC_|BG_)", LBTESTCD_LBSPEC_LBORRESU_SEX), 
        gsub("^(.*?_).*", "\\1", LBTESTCD_LBSPEC_LBORRESU_SEX), ""),
      new_LBTESTCD_LBSPEC_LBORRESU_SEX = if_else(
        prefix != "", 
        gsub("^(.*?_)(.*)", "\\2", LBTESTCD_LBSPEC_LBORRESU_SEX), 
        LBTESTCD_LBSPEC_LBORRESU_SEX
        )
      ) %>%
    #separate columns
    separate(
      new_LBTESTCD_LBSPEC_LBORRESU_SEX,
      c("LBTESTCD", "LBSPEC", "LBORRESU", "SEX"),
      sep = "_",
      remove = FALSE
      ) %>%
    #add prefixes again
    mutate(
      LBTESTCD = if_else(prefix != "", paste0(prefix, LBTESTCD), LBTESTCD)
      ) %>%
    #remove columns
    select(
      -prefix,
      -LBTESTCD_LBSPEC_LBORRESU_SEX,
      -new_LBTESTCD_LBSPEC_LBORRESU_SEX
      )
  #*****************************************************************************
  #*****************************************************************************
  #*Extract the mean value of each quantitative parameter----
  #*****************************************************************************
  means_of_parameters <- studydata %>%
    group_by(LBTESTCD, LBSPEC, LBORRESU, SEX, trial_set_description) %>%
    summarise(
      param_mean = mean(LBORRES),
      param_sd = sd(LBORRES),
      param_population = n()
      ) %>%
    #attach mean and sd from the slopes of the body weight gain linear models
    ungroup() %>%
    mutate(
      param_mean_str = signif(param_mean, 3),
      param_sd_str = signif(param_sd, 2)
      ) %>%
    #join the significance stars to this data frame
    left_join(
      test_results_tibble,
      by = c("SEX", "LBSPEC", "LBTESTCD", "LBORRESU", "trial_set_description")
      ) %>%
    suppressMessages()
  
  #*Add a temporary column containing CCG mean values to determine direction
  #*of significant changes. This will be deleted later
  CCG_means <- means_of_parameters %>%
    filter(trial_set_description == "CG") %>%
    select(LBTESTCD, LBSPEC, LBORRESU, SEX, param_mean) %>%
    rename("CCG_mean" = "param_mean")
  
  means_of_parameters_CG_sep <- merge(
    means_of_parameters,
    CCG_means,
    by = c("LBTESTCD", "LBSPEC", "LBORRESU", "SEX"),
    all.x = TRUE
    ) %>%
    mutate(
      #add a (+), if significant increase and a (-) if significant decrease
      direction = case_when(
        grepl("\\*", significance) & (param_mean > CCG_mean) ~ "(+)",
        grepl("\\*", significance) & (param_mean < CCG_mean) ~ "(-)"
      )
    ) %>%
    #append the stars (if any) to param_mean
    mutate(param_mean_str = case_when(
      !is.na(significance) ~ paste0(
        param_mean_str, " \u00B1 ", param_sd_str, significance, direction
        ),
      TRUE ~ paste0(param_mean_str, " \u00B1 ", param_sd_str)
    )
    )
  #*****************************************************************************
  #*Pivot the data frame (numeric results for further logic test)
  significance_test_numeric_results <- means_of_parameters %>%
    select(
      LBTESTCD,
      LBSPEC,
      LBORRESU,
      trial_set_description,
      SEX,
      param_mean,
      param_population
      ) %>%
    pivot_wider(
      names_from = c(trial_set_description, SEX),
      values_from = c(param_mean, param_population)
    ) %>%
    rename_with(~str_remove(., "param_mean_"))
  
  #*****************************************************************************
  #*Check for dose dependency (i.e., does the difference between means rise with
  #*rising dose?)
  rising_difference <- means_of_parameters_CG_sep %>%
    mutate(
      #calculate absolute distance between control and each dose group
      difference = abs(param_mean - CCG_mean)
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
      #check if the dependency is TRUE in all cases per group
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
    )
  #*****************************************************************************
  #*****************************************************************************
  #*Add logic rules to the data frame----
  #*****************************************************************************
  #*These are rules which were often used in toxicological decision making
  #*when a statistical significance was classified as not treatment related.
  #*Namely:
  #* 1) Do we see a dose dependency? (this is defined here as an increasing
  #*    absolute difference between control and increasing dose)
  #* 2) How many individuals are outside of HCD 2*sd range? (todo)
  #* 3) Do we observe effects in correlating endpoints? (todo)
  #* 4) Do M and F go into different directions?
  #* 5) Is the effect only visible in one sex?
  results_with_logic_rules <- significance_test_numeric_results %>%
    #add dose dependency results
    left_join(rising_difference, by = c("LBTESTCD", "LBSPEC", "LBORRESU")) %>%
    mutate(
      correlations_F = "tbd",
      correlations_M = "tbd",
      sexes_HD_in_diff_directions = "tbd") %>%
    #*Here, we check how many values of individual dose group are within the
    #*2*sd range of the supplied historical control data (refdata)
    left_join(
      studydata_with_reference,
      by = c("LBTESTCD", "LBSPEC", "LBORRESU"))
  #*****************************************************************************
  #*Pivot the data frame (character results)
  significance_test_results_as_strings <- means_of_parameters_CG_sep %>%
    select(
      LBTESTCD,
      LBSPEC,
      LBORRESU,
      trial_set_description,
      SEX,
      param_mean_str
      ) %>%
    pivot_wider(
      names_from = c(trial_set_description, SEX),
      values_from = param_mean_str
    ) 
  #*****************************************************************************
  #Check whether significances are present in both sexes
  #add a logic test to see whether significance is visible only in one sex
  #this includes also significances of e.g. LD in M and MD and HD in F!
  both_sexes <- significance_test_results_as_strings %>%
    mutate(
      across(
        ends_with("D_M"), ~grepl("\\*", .), .names = "has_significance_{.col}"
      ),
      across(
        ends_with("D_F"), ~grepl("\\*", .), .names = "has_significance_{.col}"
      ),
    )
  #Check if significances are present per sex
  both_sexes$significance_in_males <- apply(
    both_sexes %>%
      select(starts_with("has_")) %>%
      select(ends_with("_M")),
    1,
    any,
    na.rm = TRUE
    )
  both_sexes$significance_in_females <- apply(
    both_sexes %>%
      select(starts_with("has_")) %>%
      select(ends_with("_F")),
    1,
    any,
    na.rm = TRUE
  )
  #*****************************************************************************
  #*Join results together
  significance_test_results <- significance_test_results_as_strings %>%
    #join the results from previous logic rule tests
    inner_join(
      results_with_logic_rules %>%
        select(
          LBTESTCD,
          LBSPEC,
          LBORRESU,
          starts_with("LLN/ULN"),
          dose_dependency_F,
          dose_dependency_M,
          correlations_F,
          correlations_M,
          sexes_HD_in_diff_directions,
          starts_with("param_population"),
          starts_with("out_of_LON")
        ),
      by = c("LBTESTCD", "LBSPEC", "LBORRESU")
      ) %>%
    #add other results (pivoted) from means_of_parameters_CG_sep
    left_join(
      means_of_parameters_CG_sep %>%
        select(
          trial_set_description, LBTESTCD, LBSPEC, LBORRESU, SEX, p_value
          ) %>%
        pivot_wider(
          names_from = c(trial_set_description, SEX),
          values_from = p_value,
          names_prefix = "pval_"
        ) %>%
        #remove p value information. comment out this line, if needed
        select(!starts_with("pval_"), !starts_with("param_population_")),
      by = c("LBTESTCD", "LBSPEC", "LBORRESU")
    )
  #*****************************************************************************
  #*Transform effect sizes into a tibble containing effect size estimates per
  #*dose group and sex and parameter along with significance test results
  flatten_effect_sizes <- study_effect_size %>%
    imap_dfr(~ .x %>% as_tibble(), .id = "LBTESTCD_LBSPEC_LBORRESU_SEX") %>%
    #*join significance test results by right join (i.e., only values will be
    #*shown which were also calculated in the significance tests)
    right_join(
      test_results_tibble %>%
        mutate(
          LBTESTCD_LBSPEC_LBORRESU_SEX = paste(
            LBTESTCD, LBSPEC, LBORRESU, SEX, sep = "_"
            ),
          significance = case_when(
            grepl("\\*", significance) ~ "1",
            TRUE ~ "0")) %>%
        select(LBTESTCD_LBSPEC_LBORRESU_SEX, trial_set_description, significance) %>%
        unique(),
      by = c("LBTESTCD_LBSPEC_LBORRESU_SEX", "trial_set_description")
    ) %>%
    #split LBTESTCD, LBSPEC, and SEX into three different columns
    #*note that this splitting looks for underscores in the values which doesn't
    #*work for the BW/OM parameters, as they have one underscore too much.
    #*In this case, we remove the "BW_" prefix and add it in the end again.
    mutate(
      #store prefixes in separate column
      prefix = if_else(
        grepl("^(BW_|WEIGHT_|OWBW_|FC_|WC_|BG_)", LBTESTCD_LBSPEC_LBORRESU_SEX), 
        gsub("^(.*?_).*", "\\1", LBTESTCD_LBSPEC_LBORRESU_SEX), ""),
      new_LBTESTCD_LBSPEC_LBORRESU_SEX = if_else(
        prefix != "", 
        gsub("^(.*?_)(.*)", "\\2", LBTESTCD_LBSPEC_LBORRESU_SEX), 
        LBTESTCD_LBSPEC_LBORRESU_SEX
      )
    ) %>%
    #separate columns
    separate(
      new_LBTESTCD_LBSPEC_LBORRESU_SEX,
      c("LBTESTCD", "LBSPEC", "LBORRESU", "SEX"),
      sep = "_",
      remove = FALSE
    ) %>%
    #add prefixes again
    mutate(
      LBTESTCD = if_else(prefix != "", paste0(prefix, LBTESTCD), LBTESTCD)
    ) %>%
    #remove columns
    select(
      -prefix, -LBTESTCD_LBSPEC_LBORRESU_SEX, -new_LBTESTCD_LBSPEC_LBORRESU_SEX
      )
  
  return(list(
    significance_test_results,
    flatten_effect_sizes,
    test_power_flatten
    ))
}
