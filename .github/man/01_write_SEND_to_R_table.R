#01_write_SEND_to_R_table.R

#*This script selects the quantitative results of the legacy study inputs for
#*males and females respectively.
#*These results will be exported as tables to the data/derived folder.
#*******************************************************************************
#*******************************************************************************
#*Load libraries----
#*******************************************************************************
library(data.table)
library(tidyverse)
library(openxlsx)
#*******************************************************************************
#*******************************************************************************
#*Read scripts----
#*******************************************************************************
#*Read the script which calculates significance tests
source(paste0(rootpath, "/man/02_make_statistical_test.R"))
#*Read the script that reads reference data (HCD historical control data)
source(paste0(rootpath, "/man/03a4_refdata_reader.R"))
#*Create function which reads LB, BW, and OM parameters from the legacy_study
#*folder and joins them to one data frame----
#*******************************************************************************
SEND_reader <- function(
    study,
    lbdylow = 21
    ,
    lbdyhigh = 35
    ,
    omdylow = 21
    ,
    omdyhigh = 35
    ,
    bwdylow = 1
    ,
    bwdyhigh = 35
    ,
    fwdylow = 1
    ,
    fwdyhigh = 35
    ){
  #*****************************************************************************
  #*****************************************************************************
  #*Load the harmonization mapping files----
  #*****************************************************************************
  lb_mapping <- read.xlsx(paste0(rootpath, "/data/LB_mapping.xlsx")) %>%
    rename("LBTEST" = "LBTEST_old")
  
  om_mapping <- read.xlsx(paste0(rootpath, "/data/OM_mapping.xlsx")) %>%
    rename(
      "LBTEST" = "OMTEST",
      "LBSPEC_new" = "OMSPEC_new",
      "LBTESTCD_new" = "OMTESTCD_new",
      "LBTESTCD" = "OMTESTCD_old"
    )
  om_mapping <- om_mapping[!duplicated(om_mapping$LBTESTCD),]
  
  #*****************************************************************************
  #*Read unit mapping table to harmonize parameter names and units
  unit_mapping <- read.xlsx(paste0(rootpath, "/data/LB_unit_mapping.xlsx")) %>%
    rename(
      "LBTEST" = "LBTEST_old",
      "LBTESTCD" = "LBTESTCD_old",
      "LBSPEC" = "LBSPEC_old",
      "LBORRESU" = "LBORRESU_old"
    ) %>%
    select(-meanval)
  
  #*****************************************************************************
  #*****************************************************************************
  #*Read SEND study files as CSV----
  #*****************************************************************************
  #read ts data from studies
  study_ts <- fread(paste0(
    rootpath, "/data/Original/", study, "/ts.csv"
  )) %>%
    #extract year of the study from studydata
    mutate(START_YEAR_SPREFID = paste(START_YEAR, SPREFID, sep = "_")) %>%
    select(
      START_YEAR,
      SPREFID,
      START_YEAR_SPREFID,
      TRTV,
      AGE
    )
  
  #read lb data from legacy study
  study_lb <- fread(paste0(
    rootpath, "/data/Original/", study, "/lb.csv"
    )) %>%
    #remove LB measurements before day 21 and after day 35
    filter(between(LBDY, lbdylow, lbdyhigh)) %>%
    ungroup()
  
  #Harmonize LB parameters
  lb_cleaned <- study_lb %>%
    mutate(LBTEST = toupper(LBTEST)) %>%
    select(-LBTESTCD) %>%
    left_join(lb_mapping, by = "LBTEST") %>%
    mutate(LBTEST = LBTEST_new) %>%
    select(-LBTEST_new) %>%
    filter(!is.na(LBTEST)) %>%
    #*Additional cleaning step: there are urine parameters who are named the
    #*same even though they are different parameters. We need to change them
    #*here
    mutate(
      LBTESTCD = if_else(
        grepl("Creatinine per sampling period UVOL", LBMETHOD),
        "CREATUVOC",
        LBTESTCD
      ),
      LBTEST = if_else(
        grepl("Creatinine per sampling period UVOL", LBMETHOD),
        "PROTEIN PER URINE VOLUME",
        LBTEST
      )
    )
  #****************************************************************************
  #read bw data from legacy study
  study_bw <- fread(paste0(
    rootpath, "/data/Original/", study, "/bw.csv"
    )) %>%
    #*remove LB measurements manually to select time window. Always include day
    #*1 is this one is used to select initial body weight
    filter(between(BWDY, bwdylow, bwdyhigh) | BWDY == 1) %>%
    #remove animal data by IDs not present in the LB domain
    filter(USUBJID %in% unique(study_lb$USUBJID)) %>%
    #rename body weight day to labday (LBDY)
    rename(
      "LBDY" = "BWDY",
      "LBORRES" = "BWORRES",
      "LBORRESU" = "BWORRESU"
           ) %>%
    #add domain name and BWDY as "LBTESTCD"
    mutate(
      LBTEST = paste0("Body Weight Day ", str_pad(LBDY, 3, pad = "0")),
      LBTESTCD = paste0("BW_D", str_pad(LBDY, 3, pad = "0")),
      LBSPEC = "BW"
      )
  #****************************************************************************
  #derive body weight gain (bg) out of body weight
  study_bg <- study_bw %>%
    #arrange by LBDY (day of measurment)
    arrange(USUBJID, LBDY) %>%
    #calculate weight gain of each individual animal
    group_by(USUBJID) %>%
    mutate(
      BGORRES = LBORRES - lag(LBORRES),
      #rename parameters
      LBSPEC = "BG",
      LBTEST = gsub("Weight Day", "Weight Gain Day", LBTEST),
      LBTESTCD = gsub("BW_D", "BG_D", LBTESTCD)
      ) %>%
    #remove day 1 values and minimal values (if they are negative)
    filter(LBDY != 1, LBDY != min(study_bw$LBDY)) %>%
    #remove old body weight data and rename BGORRES to LBORRES
    select(-LBORRES) %>%
    rename("LBORRES" = "BGORRES")
  #****************************************************************************
  #read om data from legacy study
  study_om <- fread(paste0(
    rootpath, "/data/Original/", study, "/om.csv"
  )) %>%
    #remove excluded values
    filter(is.na(OMEXCLFL)) %>%
    #*remove LB measurements before day 21 and after day 35
    #*note: this filter may lead to heavy data loss. Longer term studies which
    #*might carry blood parameters in that time window do not necessarily have
    #*organ measurements on that same day.
    filter(between(OMDY, omdylow, omdyhigh)) %>%
    #remove animal data by IDs not present in the LB domain
    filter(USUBJID %in% unique(study_lb$USUBJID)) %>%
    #if measurements were taken twice during that time, take the later one
    group_by(USUBJID, OMSPEC) %>%
    unique() %>%
    top_n(1, abs(OMDY)) %>%
    #rename body weight day to labday (LBDY)
    rename(
      "LBEXCLFL" = "OMEXCLFL",
      "LBSTAT" = "OMSTAT",
      "LBREASND" = "OMREASND",
      "LBTESTCD" = "OMTESTCD",
      "LBTEST" = "OMTEST",
      "LBSPEC" = "OMSPEC",
      "LBORRES" = "OMSTRES",
      "LBORRESU" = "OMSTRESU",
      "LBDY" = "OMDY"
    ) %>%
    #Combine LBSPEC and LBTESTCD as the weight per organ is the interesting one
    mutate(
      LBTESTCD = paste(LBTESTCD, LBSPEC, sep = "_")
      ) %>%
    ungroup()
  
  #Harmonize LB parameters
  om_cleaned <- study_om %>%
    mutate(LBTEST = toupper(LBTEST)) %>%
    select(-LBTEST) %>%
    left_join(om_mapping, by = "LBTESTCD") %>%
    mutate(
      LBTESTCD = LBTESTCD_new,
      LBSPEC = LBSPEC_new
      ) %>%
    select(-LBTESTCD_new, -LBSPEC_new) %>%
    filter(!is.na(LBTEST))
  #****************************************************************************
  #read fw data from legacy study
  study_fw <- fread(paste0(
    rootpath, "/data/Original/", study, "/fw.csv"
  )) %>%
    #remove excluded values
    filter(FWSTAT == "Y") %>%
    #*remove FW measurements before day 1 and after day 35
    filter(between(FWDY, fwdylow, fwdyhigh)) %>%
    #remove animal data by IDs not present in the LB domain
    filter(USUBJID %in% unique(study_lb$USUBJID)) %>%
    #calculate week of measurement
    mutate(FWWK = (FWDY - 1) %/% 7 + 1) %>%
    #rename body weight day to labday (LBDY)
    rename(
      "LBSTAT" = "FWSTAT",
      "LBREASND" = "FWREASND",
      "LBTEST" = "FWTEST",
      "LBTESTCD" = "FWTESTCD",
      "LBWK" = "FWWK",
      "LBORRES" = "FWORRES",
      "LBORRESU" = "FWORRESU",
      "LBDY" = "FWDY",
      "LBENDY" = "FWENDY"
    ) %>%
    #add domain name and FWDY as "LBTESTCD"
    mutate(
      LBTEST = paste0(LBTEST, " Week ", str_pad(LBWK, 2, pad = "0")),
      LBTESTCD = paste0(LBTESTCD, "_W", str_pad(LBWK, 2, pad = "0")),
      LBSPEC = "FW"
    )
  #*****************************************************************************
  #Join the domains into one
  concatenated_study <- lb_cleaned %>%
    bind_rows(
      study_bw %>%
        select(
          SPREFID,
          SEX,
          trial_set_description,
          LBDY,
          USUBJID,
          LBORRES,
          LBORRESU,
          LBTEST,
          LBTESTCD,
          LBSPEC
          )
      ) %>%
    bind_rows(
      study_bg %>%
        select(
          SPREFID,
          SEX,
          trial_set_description,
          LBDY,
          USUBJID,
          LBORRES,
          LBORRESU,
          LBTEST,
          LBTESTCD,
          LBSPEC
        )
    ) %>%
    bind_rows(
      om_cleaned %>%
        select(
          SPREFID,
          SEX,
          trial_set_description,
          LBDY,
          USUBJID,
          LBEXCLFL,
          LBSTAT,
          LBREASND,
          LBTESTCD,
          LBTEST,
          LBSPEC,
          LBORRES,
          LBORRESU,
          LBDY
        )
      ) %>%
    bind_rows(
      study_fw %>%
        select(
          SPREFID,
          SEX,
          trial_set_description,
          LBDY,
          LBWK,
          USUBJID,
          LBSTAT,
          LBREASND,
          LBTESTCD,
          LBTEST,
          LBSPEC,
          LBORRES,
          LBORRESU,
          LBDY
        )
    ) %>%
    left_join(study_ts, by = "SPREFID") %>%
    #*Transform the dose groups in the "trial set description" into a character
    #*with 4 levels: "CG", "LD", "MD", "HD"
    mutate(
      trial_set_description = factor(
        trial_set_description,
        levels = c("CG", "LD", "MD", "HD")
      )
      ) %>%
    ungroup() %>%
    #***************************************************************************
    #harmonize the units
    #*turn first all empty values to NA
    mutate(across(where(is.character), ~na_if(., ""))) %>%
    #*Harmonize the parameter names and the units.
    #*Use the unit-mapping table for that
    left_join(
      unit_mapping,
      by = c("LBTEST", "LBTESTCD", "LBSPEC", "LBORRESU")
    ) %>%
    #*if the LBTESTCD is the body weight, or in food and water consumption
    #*translate it directly to LBTESTCD_new
    mutate(
      LBTESTCD_new = if_else(
        grepl("BW_D", LBTESTCD) |
          grepl("WC_W", LBTESTCD) |
          grepl("FC_W", LBTESTCD) |
          grepl("BG_", LBTESTCD),
        LBTESTCD,
        LBTESTCD_new
        ),
      LBTEST_new = if_else(
        grepl("BW_D", LBTESTCD) |
          grepl("WC_W", LBTESTCD) |
          grepl("FC_W", LBTESTCD) |
          grepl("BG_", LBTESTCD),
        LBTEST,
        LBTEST_new
      ),
      LBSPEC_new = case_when(
        grepl("BW_D", LBTESTCD) ~ "BW",
        grepl("BG_D", LBTESTCD) ~ "BG",
        grepl("WC_W", LBTESTCD) | grepl("FC_W", LBTESTCD) ~ "FW",
        TRUE ~ LBSPEC_new
      ),
      LBORRESU_new = case_when(
        grepl("BW_D", LBTESTCD) ~ "g",
        grepl("BG_D", LBTESTCD) ~ "g",
        grepl("WC_W", LBTESTCD) | grepl("FC_W", LBTESTCD) ~ "g/day",
        TRUE ~ LBORRESU_new
      ),
      ) %>%
    select(-LBTEST, -LBTESTCD, -LBSPEC, -LBORRESU) %>%
    rename(
      "LBTEST" = "LBTEST_new",
      "LBTESTCD" = "LBTESTCD_new",
      "LBSPEC" = "LBSPEC_new",
      "LBORRESU" = "LBORRESU_new"
    ) %>%
    #remove values with empty lab values or lab names
    filter(
      !is.na(LBTESTCD),
      !is.na(LBORRES)
      )
    
 return(concatenated_study) 
}
#*******************************************************************************
#*******************************************************************************
#*Read the study data from the legacy studies and the HCD
#*******************************************************************************
#*Read data from study A
legacy_study_A <- SEND_reader(study = "legacy_studies/study_A")

#*Reference data used in study A
HCD_study_A <- refdata_reader(
  study = "study_A",
  param = "refdata_quantitative"
)
#*******************************************************************************
#*Read data from study B
legacy_study_B <- SEND_reader(
  study = "legacy_studies/study_B",
  #cap body weight observation to day 31 (the day when LB parameters were read)
  bwdyhigh = 31,
  fwdyhigh = 31
)

#*Reference data used in study B
HCD_study_B <- refdata_reader(
  study = "study_B",
  param = "refdata_quantitative"
)
#*******************************************************************************
#*Read recovery group data from study B
legacy_study_B_recovery <- SEND_reader(
  study = "legacy_studies/study_B",
  #cap body weight observation to day 31 (the day when LB parameters were read)
  lbdylow = 36,
  lbdyhigh = 49,
  omdylow = 36,
  omdyhigh = 49,
  bwdylow = 32,
  bwdyhigh = 49
)
#*******************************************************************************
#*Read data from study C
legacy_study_C <- SEND_reader(study = "legacy_studies/study_C")

#*Reference data used in study C
HCD_study_C <- refdata_reader(
  study = "study_C",
  param = "refdata_quantitative"
)