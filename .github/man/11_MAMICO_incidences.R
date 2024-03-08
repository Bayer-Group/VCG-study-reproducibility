#*11_MAMICO_incidences.R
#*
#*This is a standalone script which means it is not executed by 00_master.R
#*
#*This script takes the results of the legacy studies and compares the
#*percentage of incidences qualitative findings (MAcroscopic findings,
#*MIcroscopic findings, Clinical Observations) of  HCD.
#*For instance, if digging and cleaning movements are observed in animals in
#*the legacy study, you can compare how many animals in the historical control
#*data (HCD) showed the same behavior.
#*
#*Input:  HCD and legacy study data. Especially, MA, MI, CO data
#*Steps:  - Extract the findings of interest.
#*        - Calculate the percentage of findings per population.
#*        - Do the same for HCD.
#*        - Visualize this using the 06a2_MI_plot.R plot.
#*Output: Percentage of finding with background incidences.
#*******************************************************************************
#*******************************************************************************
#*Load libraries----
#*******************************************************************************
require(tidyverse)
#*******************************************************************************
#*******************************************************************************
#*Read qualitative parameters----
#*******************************************************************************
source(paste0(rootpath, "/man/01a1_write_qualitative_SEND_data_to_R_table.R"))
#*******************************************************************************
#*******************************************************************************
#*beginning of function
#*******************************************************************************
MAMICO_incidences_calculator <- function(
    legacy_study,
    HCD,  
    QPSPEC_input = NULL,
    QPORRES_input = NULL,
    DOMAIN_input = NULL
    ){
  #*****************************************************************************
  #*filter data by pre-select parameters
  legacy_study_filtered <- legacy_study %>%
    filter(
      grepl(QPSPEC_input, QPSPEC, perl = TRUE),
      grepl(QPORRES_input, QPORRES, perl = TRUE),
      DOMAIN == DOMAIN_input
    )
  
  HCD_filtered <- HCD %>%
    filter(
      grepl(QPSPEC_input, QPSPEC, perl = TRUE),
      grepl(QPORRES_input, QPORRES, perl = TRUE),
      DOMAIN == DOMAIN_input
    )
  
  #filter HCD and CCG also for domain
  legacy_study_unfiltered <- legacy_study %>% filter(DOMAIN == DOMAIN_input)
  
  HCD_unfiltered <- HCD %>% filter(DOMAIN == DOMAIN_input)
  
  test <- HCD %>%
    select(SPREFID, USUBJID, SEX, trial_set_description, DOMAIN) %>%
    unique()
  
  test_CL <- test %>% filter(DOMAIN == "CL")
  test_MI <- test %>% filter(DOMAIN == "MI")
  
  test2 <- test_CL %>% filter(!USUBJID %in% test_MI$USUBJID)
  n_distinct(test2$SPREFID)
  
  #*****************************************************************************
  #*Split data into control and dose groups
  data_unfiltered_control <- legacy_study_unfiltered %>%
    filter(trial_set_description == "CG")
  
  data_unfiltered_dose <- legacy_study_unfiltered %>%
    filter(trial_set_description != "CG")
  
  data_filtered_control <- legacy_study_filtered %>%
    filter(trial_set_description == "CG")
  
  data_filtered_dose <- legacy_study_filtered %>%
    filter(trial_set_description != "CG")
  #*****************************************************************************
  #calculate sample size per sex
  sample_size_calculator <- function(data){
    data %>%
      select(USUBJID, SEX, trial_set_description) %>%
      unique() %>%
      group_by(SEX, trial_set_description) %>%
      summarize(population = n(), .groups = "drop")
  }
  #sample size of all animals per sex to calculate incidence
  CG_size_unfiltered <- sample_size_calculator(data_unfiltered_control)
  dose_size_unfiltered <- sample_size_calculator(data_unfiltered_dose)
  HCD_size_unfiltered <- sample_size_calculator(HCD_unfiltered)
  #sample size of animals having exactly this finding
  CG_size_filtered <- sample_size_calculator(data_filtered_control)
  dose_size_filtered <- sample_size_calculator(data_filtered_dose)
  HCD_size_filtered <- sample_size_calculator(HCD_filtered)
  #*****************************************************************************
  #*calculate % of incidences in concurrent controls and in HCD
  #join number of incidences and total number of animals
  CCG_incidences <- merge(
    CG_size_filtered,
    CG_size_unfiltered,
    by = c("SEX", "trial_set_description"),
    all = TRUE
    ) %>%
    #turn all NA into 0
    mutate(
      across(starts_with("population"),function(x) {if_else(is.na(x), 0, x)}),
      incidences = round(population.x / population.y * 100),
      incidence_text = paste0(
        population.x, " out of ", population.y, " (", incidences, " %)"
        ),
      trial_set_description = "CCG"
      )
  
  HCD_incidences <- merge(
    HCD_size_filtered,
    HCD_size_unfiltered,
    by = c("SEX", "trial_set_description"),
    all = TRUE
  ) %>%
    #turn all NA into 0
    mutate(
      across(starts_with("population"),function(x) {if_else(is.na(x), 0, x)}),
      incidences = round(population.x / population.y * 100),
      incidence_text = paste0(
        population.x, " out of ", population.y, " (", incidences, " %)"
      ),
      trial_set_description = "HCD"
    )
  
  dose_incidences <- merge(
    dose_size_filtered,
    dose_size_unfiltered,
    by = c("SEX", "trial_set_description"),
    all = TRUE
  ) %>%
    #turn all NA into 0
    mutate(
      across(starts_with("population"),function(x) {if_else(is.na(x), 0, x)}),
      incidences = round(population.x / population.y * 100),
      incidence_text = paste0(
        population.x, " out of ", population.y, " (", incidences, " %)"
      )
    )
  
  #*****************************************************************************
  #*summarize results in a table----
  #*****************************************************************************
  incidences <- dose_incidences %>%
    bind_rows(CCG_incidences) %>%
    bind_rows(HCD_incidences) %>%
    #order dose groups by making a factor out of character vector
    mutate(
      trial_set_description = factor(
        trial_set_description, levels = c("HCD", "CCG", "LD", "MD", "HD")
        )
      ) %>%
    arrange(SEX, trial_set_description)
  
  return(incidences)
}
#end of function
#*******************************************************************************
#*******************************************************************************
#*Call function to calculate background incidences----
#*******************************************************************************
#*******************************************************************************
#*legacy study B----
#*******************************************************************************
#study B follicular cell hypertrophy in the thyroid gland----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_B_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*gland)(?=.*thyroid)",
  QPORRES_input = "(?i)(?=.*follicul)(?=.*cell)(?=.*hypertroph)",
  DOMAIN_input = "MI"
)

#study B macrophages in the thymus----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_B_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*thymus)",
  QPORRES_input = "(?i)(?=.*macrophages)",
  DOMAIN_input = "MI"
)

#study B degeneration regeneration of Harderian gland----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_B_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*gland)(?=.*harderian)",
  QPORRES_input = "(?i)(?=.*degeneration)(?=.*regeneration)",
  DOMAIN_input = "MI"
)

#study B extramedullary hematopoiesis in the spleen----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_B_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*spleen)",
  QPORRES_input = "(?i)(?=.*extramedullary)(?=.*hematopoiesis)",
  DOMAIN_input = "MI"
)
#*******************************************************************************
#*legacy study C----
#*******************************************************************************
#check for names. un-comment if needed
legacy_study_C_qualitative %>%
  filter(
    # DOMAIN == "MA",
    # !grepl("UNREMARKABLE", QPORRES, ignore.case = T),
    grepl("congest", QPORRES, ignore.case = T),
    # grepl("fat", QPSPEC, ignore.case = T)
  ) %>%
  select(-SPREFID, -USUBJID) %>%
  unique()

#study C Dilatation/distension in small intestine----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*small)(?=.*intest)",
  QPORRES_input = "(?i)(?=.*(?:dilat.*|disten.*))",
  DOMAIN_input = "MA"
)
#study C Dilatation/distension or change in contents in large intestine----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*large)(?=.*intest)",
  QPORRES_input = "(?i)(?=.*(?:dilat.*|disten.*|chang.*cont.*))",
  DOMAIN_input = "MA"
)
#study C Dilatation/distension in stomach----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*stom)",
  QPORRES_input = "(?i)(?=.*(?:dilat.*|disten.*))",
  DOMAIN_input = "MA"
)
#study C Change in contents in small intestine----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*small)(?=.*intest)",
  QPORRES_input = "(?i)(?=.*chang)(?=.*cont)",
  DOMAIN_input = "MA"
)
#study C Change in contents in large intestine----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*larg)(?=.*intest)",
  QPORRES_input = "(?i)(?=.*chang)(?=.*cont)",
  DOMAIN_input = "MA"
)
#study C Change in contents in stomach----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.stom)",
  QPORRES_input = "(?i)(?=.*chang)(?=.*cont)",
  DOMAIN_input = "MA"
)
#study C Hypertrophy of mesenteric/cardiac arteries in vascular system----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*mesentery)",
  QPORRES_input = "(?i)(?=.*hypertrophy)(?=.*arter)",
  DOMAIN_input = "MI"
)

#study C inflammatory infiltrates in the mesentery----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*mesentery)",
  QPORRES_input = "(?i)(?=.*Infiltrate)(?=.*inflammatory)",
  DOMAIN_input = "MI"
)

#study C paneth cell hypertrophy in the small intestine----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*small)(?=.*intestin)",
  QPORRES_input = "(?i)(?=.*Paneth)(?=.*Hypertrophy)",
  DOMAIN_input = "MI"
)
#study C inflammatory/degenrative changes in the GI tract----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?:intest.*|stom.*)",
  QPORRES_input = "(?i)(?:inflam.*|degen.*)",
  DOMAIN_input = "MI"
)

#study C hypertrophy/hyperplasia in zona glomerulosa in adrenal cortex----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*adrenal)(?=.*gland)",
  QPORRES_input = "(?i)(?:Hypertr.*|hyperpl.*)(?=.*glomerulosa)",
  DOMAIN_input = "MI"
)

#study C hypertrophy/hyperplasia in zona fasciculata in adrenal cortex----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*adrenal)(?=.*gland)",
  QPORRES_input = "(?i)(?:Hypertr.*|hyperpl.*)(?=.*fascicul)",
  DOMAIN_input = "MI"
)

#study C congestion/hemorraghe in the adrenal cortex----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*adrenal)(?=.*gland)",
  QPORRES_input = "(?i)(?:congest.*|hem.or.ag*)",
  DOMAIN_input = "MI"
)

#study C cytoplasmic alteration in the adrenal medulla----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*adrenal)(?=.*gland)",
  QPORRES_input = "(?i)(?=.*Cytoplasm)",
  DOMAIN_input = "MI"
)
#study C decreased cellularity in bone marrow----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*marrow)",
  QPORRES_input = "(?i)(?=.*cellul)(?=.*decr)",
  DOMAIN_input = "MI"
)
#study C congestion in bone marrow----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*marrow)",
  QPORRES_input = "(?i)(?=.*cong)",
  DOMAIN_input = "MI"
)
#study C accumulation of adipocytes in bone marrow----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*marrow)",
  QPORRES_input = "(?i)(?=.*adipocyt)(?=.*accum)",
  DOMAIN_input = "MI"
)
#study C decreased cellularity in the thymus----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*thym.*s)",
  QPORRES_input = "(?i)(?=.*decr)(?=.*cellul)",
  DOMAIN_input = "MI"
)

#study C increase of tingible body macrophages in the thymus----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*thym.*s)",
  QPORRES_input = "(?i)(?=.*incr)(?=.*ting)(?=.*body)(?=.*macrop)",
  DOMAIN_input = "MI"
)

#study C apoptosis/necrosis in the thymus----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*thym.*s)",
  QPORRES_input = "(?i)(?=.*(?:apop.*|ne.ros.*))",
  DOMAIN_input = "MI"
)

#study C decreased cellularity in the spleen----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*spleen)",
  QPORRES_input = "(?i)(?=.*decr)(?=.*cellul)",
  DOMAIN_input = "MI"
)

#study C decreased cellularity in the lymph nodes----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*lymph)",
  QPORRES_input = "(?i)(?=.*decr)(?=.*cellul)",
  DOMAIN_input = "MI"
)

#study C decreased cellularity in the Peyer's patches----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*peyer)",
  QPORRES_input = "(?i)(?=.*decr)(?=.*cellul)",
  DOMAIN_input = "MI"
)

#study C glycogen depletion in the liver----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*liver)",
  QPORRES_input = "(?i)(?=.*glyc)(?=.*depl.t)",
  DOMAIN_input = "MI"
)

#study C atrophy of acinar cells in salivary glands----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*saliv)",
  QPORRES_input = "(?i)(?=.*atrop)(?=.*acina)",
  DOMAIN_input = "MI"
)

#study C atrophy of acinar cells in the pancreas----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*pan.r)",
  QPORRES_input = "(?i)(?=.*atrop)(?=.*acina)",
  DOMAIN_input = "MI"
)

#study C vacuolation of the Harderian gland----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*harder)",
  QPORRES_input = "(?i)(?=.*vac..lat)",
  DOMAIN_input = "MI"
)

#study C atrophy of follicular epithelial cell in thyroid gland----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*th.r)",
  QPORRES_input = "(?i)(?=.*atrop)(?=.*foll)(?=.*epit.el)",
  DOMAIN_input = "MI"
)

#study C decreased secretion in prostate----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*prost)",
  QPORRES_input = "(?i)(?=.*decr)(?=.*secr)",
  DOMAIN_input = "MI"
)

#study C decreased secretion in seminal vesicles----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*semin)(?=.*vesic)",
  QPORRES_input = "(?i)(?=.*decr)(?=.*secr)",
  DOMAIN_input = "MI"
)

#study C atrophy of epithelial cells in the vagina----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*vagin)",
  QPORRES_input = "(?i)(?=.*atrop)(?=.*epit)",
  DOMAIN_input = "MI"
)

#study C atrophy of adipocates in the adipose tissue----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*)",
  QPORRES_input = "(?i)(?=.*atrop)(?=.*adipo)",
  DOMAIN_input = "MI"
)

#study C clinical observation: piloerection----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*clinic)(?=.*observ)",
  QPORRES_input = "(?i)(?=.*piloere)",
  DOMAIN_input = "CL"
)

#study C clinical observation: digging and cleaning movements----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*clinic)(?=.*observ)",
  QPORRES_input = "(?i)(?=.*digg)(?=.*clean)",
  DOMAIN_input = "CL"
)

#study C clinical observation: increased salivation----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*clinic)(?=.*observ)",
  QPORRES_input = "(?i)(?=.*inc)(?=.*saliv)",
  DOMAIN_input = "CL"
)

#study C clinical observation: resistance during handling----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*clinic)(?=.*observ)",
  QPORRES_input = "(?i)(?=.*resist)(?=.*handling)",
  DOMAIN_input = "CL"
)

#study C clinical observation: changed feces consistency----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*clinic)(?=.*observ)",
  QPORRES_input = "(?i)(?=.*fec)(?=.*consist)",
  DOMAIN_input = "CL"
)

#study C clinical observation: high stepping gait----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*)",
  QPORRES_input = "(?i)(?=.*stepp)(?=.*gait)",
  DOMAIN_input = "CL"
)

#study C clinical observation: reduced motility----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*)",
  QPORRES_input = "(?i)(?=.*reduc)(?=.*motility)",
  DOMAIN_input = "CL"
)

#study C clinical observation: erect tail----
MAMICO_incidences_calculator(
  legacy_study = legacy_study_C_qualitative,
  HCD = HCD_qualitative, 
  QPSPEC_input = "(?i)(?=.*)",
  QPORRES_input = "(?i)(?=.*erect)(?=.*tail)",
  DOMAIN_input = "CL"
)

