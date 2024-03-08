#*01a1_write_qualitative_SEND_data_to_R_table.R
#*
#*2023 08 17
#*
#*This script reads and concatenates CSV data of the qualitative parameters
#*(i.e., microscopic findings (MI), macroscopic findings (MA), and clinical
#*observations (CO)) into one table where the incidences can be calculated from.
#*
#*Input:  CSV data in SEND format.
#*Steps:  - Read data.
#*        - Harmonize column names.
#*        - Concatenate data into one table.
#*        - Export table.
#*Output: Concatenated table with qualitative parameter data and values.
#*******************************************************************************
#*******************************************************************************
#*load libraries----
#*******************************************************************************
require(data.table)
require(tidyverse)
#*******************************************************************************
#*******************************************************************************
#*beginning of function
#*******************************************************************************
qualitative_SEND_reader <- function(
    study,
    madylow = 1
    ,
    madyhigh = 100
    ,
    midylow = 1
    ,
    midyhigh = 100
    ,
    cldylow = 1
    ,
    cldyhigh = 29
    ){
  #*****************************************************************************
  #*****************************************************************************
  #*renamer function----
  #*****************************************************************************
  #*Create a function which renames all prefixes (CL, MA, MI) into QP
  #*(qualitative parameters)
  renamer <- function(df){
    #read and replace names of data frame
    gsub("^MA|^MI|^CL", "QP", names(df))
  }
  #*****************************************************************************
  #*****************************************************************************
  #read MI findings----
  #*****************************************************************************
  #*microscopic findings (MI)
  mi <- fread(paste0(
    rootpath, "/data/Original/", study, "/mi.csv"
  )) %>%
    #filter for study day
    filter(between(MIDY, midylow, midyhigh)) %>%
    #turn severity findings into factor
    mutate(
      MISEV = as.character(MISEV),
      MISEV = factor(MISEV, levels = c(NA, "1", "2", "3", "4", "5")),
      #turn severity grade into text
      MISEVTXT = case_when(
        is.na(MISEV) ~ "no finding",
        MISEV == "1" ~ "minimal",
        MISEV == "2" ~ "mild",
        MISEV == "3" ~ "moderate",
        MISEV == "4" ~ "marked",
        MISEV == "5" ~ "severe"
      ),
      MISEVTXT = factor(
        MISEVTXT,
        levels = c(
          "no finding", "minimal", "mild", "moderate", "marked", "severe"
        )
      ),
      #turn dose groups into factor
      trial_set_description = factor(
        trial_set_description, levels = c("CG", "LD", "MD", "HD")
      ),
      #add name of domain
      DOMAIN = "MI"
    )
  #rename prefixes
  names(mi) <- renamer(mi)
  #*****************************************************************************
  #*****************************************************************************
  #read MA findings----
  #*****************************************************************************
  #*macroscopic findings (MA)
  ma <- fread(paste0(
    rootpath, "/data/Original/", study, "/ma.csv"
  )) %>%
    #filter for study day
    filter(between(MADY, madylow, madyhigh)) %>%
    #turn severity findings into factor
    mutate(
      MASEV = as.character(MASEV),
      MASEV = factor(MASEV, levels = c(NA, "1", "2", "3", "4", "5")),
      #turn severity grade into text
      MASEVTXT = case_when(
        is.na(MASEV) ~ "no finding",
        MASEV == "1" ~ "minimal",
        MASEV == "2" ~ "mild",
        MASEV == "3" ~ "moderate",
        MASEV == "4" ~ "marked",
        MASEV == "5" ~ "severe"
      ),
      MASEVTXT = factor(
        MASEVTXT,
        levels = c(
          "no finding", "minimal", "mild", "moderate", "marked", "severe"
        )
      ),
      #turn dose groups into factor
      trial_set_description = factor(
        trial_set_description, levels = c("CG", "LD", "MD", "HD")
      ),
      #add name of domain
      DOMAIN = "MA"
    )
  #rename prefixes
  names(ma) <- renamer(ma)
  #*****************************************************************************
  #*****************************************************************************
  #read CL findings----
  #*****************************************************************************
  #*clinical observations (CL)
  cl <- fread(paste0(
    rootpath, "/data/Original/", study, "/cl.csv"
  )) %>%
    #filter for study day
    filter(between(CLDY, cldylow, cldyhigh)) %>%
    #turn severity findings into factor
    mutate(
      CLSEV = as.character(CLSEV),
      CLSEV = factor(CLSEV, levels = c(NA, "1", "2", "3", "4", "5")),
      #turn severity grade into text
      CLSEVTXT = case_when(
        is.na(CLSEV) ~ "no finding",
        CLSEV == "1" ~ "minimal",
        CLSEV == "2" ~ "mild",
        CLSEV == "3" ~ "moderate",
        CLSEV == "4" ~ "marked",
        CLSEV == "5" ~ "severe"
      ),
      CLSEVTXT = factor(
        CLSEVTXT,
        levels = c(
          "no finding", "minimal", "mild", "moderate", "marked", "severe"
        )
      ),
      #turn dose groups into factor
      trial_set_description = factor(
        trial_set_description, levels = c("CG", "LD", "MD", "HD")
      ),
      #add name of domain
      DOMAIN = "CL"
    ) %>%
    #rename CLCAT to CLSPEC to that it can be combined with other domains
    rename("CLSPEC" = "CLCAT")
  #rename prefixes
  names(cl) <- renamer(cl)
  #*****************************************************************************
  #*****************************************************************************
  #concatenate all findings into one data frame----
  #*****************************************************************************
  qualitative_concatenated <- ma %>% bind_rows(mi) %>% bind_rows(cl)
}
#*******************************************************************************
#*******************************************************************************
#*Call function to create qualitative data for HCD and legacy studies----
#*******************************************************************************
#legacy study A
legacy_study_A_qualitative <- qualitative_SEND_reader("legacy_studies/study_A")
#legacy study B
#specify study days, so that main group and recovery group is not mixed up
legacy_study_B_qualitative <- qualitative_SEND_reader(
  study =  "legacy_studies/study_B",
  cldyhigh = 29
  )
legacy_study_B_recovery_qualitative <- qualitative_SEND_reader(
  study = "legacy_studies/study_B",
  madylow = 36, madyhigh = 49,
  midylow = 36, midyhigh = 49,
  cldylow = 30, cldyhigh = 49
)
#legacy study C
legacy_study_C_qualitative <- qualitative_SEND_reader(
  study = "legacy_studies/study_C",
  #include the study data from female high dose which was taken on day 6
  madylow = 1,
  midylow = 1
  )
#historical control data (HCD)
HCD_qualitative <- qualitative_SEND_reader("HCD_RAT")
