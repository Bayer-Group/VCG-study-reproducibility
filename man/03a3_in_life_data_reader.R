#*03a3_in_life_data_reader.R
#*
#*This script reads the in-life data table results extracted from PDF study
#*reports. It extracts the summary statistics and the significant differences
#*which will be compared to the statistical results of the R-script
#*02_make_statistical_test.R. If the results do not overlap, the significance
#*tests will be adjusted accordingly
#*******************************************************************************
#*******************************************************************************
#*read scripts----
#*******************************************************************************
#*#Get the function which calculates respective significance test
source(paste0(rootpath, "/man/02_make_statistical_test.R"))
#*******************************************************************************
#*******************************************************************************
#*load libraries----
#*******************************************************************************
require(tidyverse)
require(openxlsx)
#*******************************************************************************
#*******************************************************************************
#*Read files for cleaning----
#*******************************************************************************
#read mapping table to harmonize data names
inlife_data_mapping <- read.xlsx(
  paste0(
    rootpath,
    "/data/Original/legacy_studies/inlife_data/inlife_data_mapping.xlsx"
  )
) %>%
  rename(
    "LBTESTCD" = "LBTESTCD_old",
    "LBSPEC" = "LBSPEC_old",
    "LBORRESU" = "LBORRESU_old",
  )
#*******************************************************************************
#*******************************************************************************
#*create a function to read and transform in life data----
#*******************************************************************************
read_and_pivot_inlife <- function(inlife_data){
  #read in life data
  inlife_data <- read.xlsx(
    paste0(
      rootpath,
      "/data/Original/legacy_studies/inlife_data/",
      inlife_data
    )
  ) %>%
    mutate(across(everything(), as.character))
  
  #select only mean values
  inlife_mean <- inlife_data %>%
    filter(numbers == "Means") %>%
    select(-numbers)
  
  #pivot the data to long format
  inlife_long <- inlife_mean %>%
    pivot_longer(
      c(-trial_set_description, -SEX),
      names_to = c("LBTESTCD", "LBSPEC", "LBDY", "LBORRESU"),
      names_sep = "_"
    ) %>%
    #*remove empty mean values.Don't be on alert, this is normal because
    #*measurements from males and females are usually taken on different days
    filter(!is.na(value)) %>%
    #rename the units back to original structure
    mutate(
      LBORRESU = gsub("X", "/", LBORRESU),
      LBORRESU = gsub("00", " ", LBORRESU),
      #turn day of visit back to numeric
      LBDY = parse_number(LBDY)
    )
  
  #*harmonize the names to match the ones of the legacy study
  inlife_cleaned <- inlife_long %>%
    left_join(
      inlife_data_mapping,
      by = c("LBTESTCD", "LBSPEC", "LBORRESU")
    ) %>%
    select(-LBTESTCD, -LBSPEC, -LBORRESU) %>%
    rename(
      "LBTESTCD" = "LBTESTCD_new",
      "LBTEST" = "LBTEST_new",
      "LBSPEC" = "LBSPEC_new",
      "LBORRESU" = "LBORRESU_new"
    )
  
  #create a column showing all statistically significant results
  inlife_with_significance <- inlife_cleaned %>%
    mutate(
      significance = gsub("[0-9\\.E\\-]", "", value),
      significance = case_when(
        significance == "" ~ FALSE,
        TRUE ~ TRUE
      )
    ) %>%
    #add domain name and BWDY as "LBTESTCD"
    mutate(
      LBTEST = case_when(
        LBSPEC == "BW" ~ paste0(
          "Body Weight Day ", str_pad(LBDY, 3, pad = "0")
          ),
        TRUE ~ LBTEST
      ),
      LBTESTCD = case_when(
        LBSPEC == "BW" ~ paste0("BW_D", str_pad(LBDY, 3, pad = "0")),
        TRUE ~ LBTESTCD
      )
    )
}

#*******************************************************************************
#*******************************************************************************
#*read the in life data files----
#*******************************************************************************
inlife_study_A <- read_and_pivot_inlife("Study_A_inlife_data.xlsx")
inlife_study_B <- read_and_pivot_inlife("Study_B_inlife_data.xlsx")
inlife_study_C <- read_and_pivot_inlife("Study_C_inlife_data.xlsx")

#for cleaning purposes: uncomment if needed
# concatenated_inlife <- inlife_study_A %>%
#   select(LBTESTCD, LBSPEC, LBORRESU) %>%
#   rbind(inlife_study_B %>% select(LBTESTCD, LBSPEC, LBORRESU)) %>%
#   rbind(inlife_study_B %>% select(LBTESTCD, LBSPEC, LBORRESU)) %>%
#   unique()
# 
# write.xlsx(concatenated_inlife, "inlife_data_mapping.xlsx")

#*******************************************************************************
#*******************************************************************************
#*Compare to original results----
#*******************************************************************************
#*Study A
#*******************************************************************************
#*Use the "study_results_table" function for this
CCG_legacy_study_A <- tox_statistical_test(
  studydata = legacy_study_A,
  refdata = HCD_study_A
)[[1]]

#pivot into longer format
CCG_legacy_study_A_longer <- CCG_legacy_study_A %>%
  select(LBTESTCD, LBSPEC, LBORRESU, matches("^[A-Z]{2}_M$|^[A-Z]{2}_F$")) %>%
  pivot_longer(
    -c(LBTESTCD, LBSPEC, LBORRESU),
    names_to = c("trial_set_description", "SEX"),
    names_sep = "_"
  ) %>%
  mutate(significance = if_else(grepl("\\*", value), TRUE, FALSE))

print(
  nrow(
    CCG_legacy_study_C_longer %>%
         filter(
           LBSPEC %in% c("SERUM", "PLASMA", "WHOLE BLOOD", "URINE"),
           trial_set_description != "CG"
           )
    )
  )

#check which studies match with the original results
matched_studies_A <- inlife_study_A %>%
  inner_join(
    CCG_legacy_study_A_longer,
    by = c("trial_set_description", "LBTESTCD", "LBSPEC", "LBORRESU", "SEX")
    ) %>%
  rename(
    "studyreport_significance" = "significance.x",
    "reproduced_significance" = "significance.y"
  ) %>%
  mutate(
    overlapping_significances = if_else(
      studyreport_significance == reproduced_significance, T, F)
      ) %>%
  filter(overlapping_significances == F)
#*******************************************************************************
#*Study B
#*******************************************************************************
#*Use the "study_results_table" function for this
CCG_legacy_study_B <- tox_statistical_test(
  studydata = legacy_study_B,
  refdata = HCD_study_A
)[[1]]

#pivot into longer format
CCG_legacy_study_B_longer <- CCG_legacy_study_B %>%
  select(LBTESTCD, LBSPEC, matches("^[A-Z]{2}_M$|^[A-Z]{2}_F$")) %>%
  pivot_longer(
    -c(LBTESTCD, LBSPEC),
    names_to = c("trial_set_description", "SEX"),
    names_sep = "_"
  ) %>%
  mutate(significance = if_else(grepl("\\*", value), TRUE, FALSE))

#check which studies match with the original results
matched_studies_B <- inlife_study_B %>%
  inner_join(
    CCG_legacy_study_B_longer,
    by = c("trial_set_description", "LBTESTCD", "LBSPEC", "SEX")
  ) %>%
  rename(
    "studyreport_significance" = "significance.x",
    "reproduced_significance" = "significance.y"
  ) %>%
  mutate(
    overlapping_significances = if_else(
      studyreport_significance == reproduced_significance, T, F)
  ) %>%
  filter(overlapping_significances == F)
#*******************************************************************************
#*Study C
#*******************************************************************************
#*Use the "study_results_table" function for this
CCG_legacy_study_C <- tox_statistical_test(
  studydata = legacy_study_C,
  refdata = HCD_study_A
)[[1]]

#pivot into longer format
CCG_legacy_study_C_longer <- CCG_legacy_study_C %>%
  select(LBTESTCD, LBSPEC, matches("^[A-Z]{2}_M$|^[A-Z]{2}_F$")) %>%
  pivot_longer(
    -c(LBTESTCD, LBSPEC),
    names_to = c("trial_set_description", "SEX"),
    names_sep = "_"
  ) %>%
  mutate(significance = if_else(grepl("\\*", value), TRUE, FALSE))

matched_studies_C <- inlife_study_C %>%
  inner_join(
    CCG_legacy_study_C_longer,
    by = c("trial_set_description", "LBTESTCD", "LBSPEC", "SEX")
  ) %>%
  rename(
    "studyreport_significance" = "significance.x",
    "reproduced_significance" = "significance.y"
  ) %>%
  mutate(
    overlapping_significances = if_else(
      studyreport_significance == reproduced_significance, T, F)
  ) %>%
  filter(overlapping_significances == F)

rbind(
  matched_studies_A, matched_studies_B, matched_studies_C
) %>%
  select(LBTESTCD, SEX, trial_set_description) %>% View()
