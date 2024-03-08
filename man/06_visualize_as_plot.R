#*06_visualize_as_plot.R

#*This script generates a ggplot showing significant results of the
#*legacy study and the concatenated results of the VCG resampling experiment.
#*The goal is to see whether the CCG values are strongly outside of the VCG range
#*resulting poor reproducibility.
##### program starts here
#*********************************************************************************;
#*********************************************************************************;
#*********************************************************************************;
#*****************************************************************************
#*****************************************************************************
#*load libraries----
#*****************************************************************************
require(tidyverse)
require(effsize)
#*********************************************************************************;
study_results_plot <- function(
  legacy_study_results = CCG_legacy_study_C,
  resampling_results = VCG_legacy_study_C,
  keep_consistent_results = TRUE,
  color_by = "power"
){
  #*****************************************************************************
  #*Get the effect size results from each iteration
  effsizes_legacy_study <- legacy_study_results[[2]]
  effsizes_resampling <- resampling_results[[2]]
  #*****************************************************************************
  #*Get the power calculation results (if done) and add them to the results tab
  power_of_legacy <- legacy_study_results[[3]]
  power_of_resampling <- resampling_results[[3]]
  #*****************************************************************************
  #*Add power to effect size table
  legacy <- merge(
    effsizes_legacy_study,
    power_of_legacy,
    by = c("LBTESTCD", "LBSPEC", "LBORRESU", "SEX"),
    all.x = TRUE
    )
  #*Add power to effect size table
  resampling <- merge(
    effsizes_resampling,
    power_of_resampling,
    by = c("LBTESTCD", "LBSPEC", "LBORRESU", "SEX", "iteration"),
    all.x = TRUE
  )
  #*****************************************************************************
  #Turn dose groups into factors to order them in the plot
  table_styler <- function(unstyled_table){
    styled_table <- unstyled_table %>%
      mutate(
        trial_set_description = case_when(
          trial_set_description == "LD" ~ "Low Dose",
          trial_set_description == "MD" ~ "Mid Dose",
          trial_set_description == "HD" ~ "High Dose"
        ),
        trial_set_description = factor(
          trial_set_description,
          levels = c("Low Dose", "Mid Dose", "High Dose")
        ),
        SEX = case_when(
          SEX == "M" ~ "Males",
          SEX == "F" ~ "Females"
        ),
        #*add a color conditionally to the data frames based on whether you
        #*want to color by power or by significance
        colorings = case_when(
          color_by == "power" ~ power,
          color_by == "p_value" ~ as.numeric(significance)
        )
      )
    
    #*change to discrete values separately as this isn't allowed in a case_when
    #*function
    corrected_styled_table <- if(color_by == "p_value"){
      styled_table %>% mutate(colorings = as.character(colorings))
    }else{styled_table}
    
    return(corrected_styled_table)
  }
  #apply function to tables
  legacy_study_styled <- table_styler(legacy)
  resampling_styled <- table_styler(resampling)
  
  
  #*Select whether you only want to observe significant changes and those which
  #*were inconsistent with the results of the legacy study
  #get all significant LBTESTCD from parameters in VCG
  significant_VCG_res_LBTESTCD <- resampling_styled %>%
    filter(significance == "1") %>%
    pull(LBTESTCD) %>%
    unique()
  
  #*select CCG results which were either significant or inconsistent to VCG
  #*results in their significance test
  legacy_study_styled <- if(keep_consistent_results == T){
    legacy_study_styled %>%
      #drop empty values
      filter(!is.na(trial_set_description))
  }else{
    legacy_study_styled %>%
      filter(
        significance == "1" |
          LBTESTCD %in% significant_VCG_res_LBTESTCD,
        #drop empty values
        !is.na(trial_set_description)
      )
  }
  
  #filter by inconsistent results and remove empty values (bugfix)
  resampling_styled <- if(keep_consistent_results == T){
    resampling_styled %>%
      #drop empty values
      filter(!is.na(trial_set_description))
  }else{
    resampling_styled %>%
      filter(
        LBTESTCD %in% resampling_styled$LBTESTCD,
        #drop empty values
        !is.na(trial_set_description)
        )
  }
  #*****************************************************************************
  #*****************************************************************************
  #Define title and colorings----
  #*****************************************************************************
  #*add the title based on whether the plot is colored by the significance
  #*or by the power
  plottitle <- ifelse(
    color_by == "power",
    "Effect sizes colored by statistical power",
    "Effect sizes colored by significance test"
    )
  #*****************************************************************************
  #*****************************************************************************
  #Create a plot showing effect sizes of each parameter by dose group----
  #*****************************************************************************
  p <- ggplot(
    data = resampling_styled,
    aes(
      x = effsize,
      y = LBTESTCD,
      color = colorings
    )
  ) +
    #add the VCG data as lines
    geom_point(
      shape="|",
      alpha=.75,
      size=2
    ) +
    #add the cenral line indicating that the dose group is 0
    geom_vline(
      xintercept = 0, col="#a5a5a5"
    ) +
    #add the CCG as a rhombus
    geom_point(
      data = legacy_study_styled,
      shape = 23,
      size = 2
    ) +
    #Create a theme to the plot
    theme_bw() +
    #add title to the plot
    ggtitle(plottitle) +
    #label x-axis
    xlab("Effect size (Cliff's delta)") +
    #label y-axis
    ylab("Parameter") +
    #*reverse order of y-axis so that it's in alphabetical order
    scale_y_discrete(limits = rev) +
    #Split the plot into sex and dose groups
    facet_grid(SEX ~ trial_set_description)
  
  if(color_by == "p_value"){
    p <- p +
    scale_color_manual(
      values = c("0" = "#1f78b4", "1" = "#e31a1c"),
      name = "Significance", labels = c("n.s.", "p <= .05")
      )
  }else{
    p <- p + scale_color_gradient2(
      low =  "#e31a1c",
      mid =  "#10384f",
      high = "#1f78b4",
      midpoint = .8,
      name = "Power"
      )
  }
  p
}

#*******************************************************************************
#*******************************************************************************
#*Write and save images----
#*******************************************************************************
#*Legacy study A
#*******************************************************************************
effsize_legacy_study_A <- study_results_plot(
    legacy_study_results = CCG_legacy_study_A,
    resampling_results = VCG_legacy_study_A,
    keep_consistent_results = F,
    color_by = "p_value"
)
ggsave(
  plot = effsize_legacy_study_A,
  filename = paste0(path_res, "/effsize_study_A_no_imp.png"),
  dpi = 600,
  width = 20,
  height = 30,
  scale = 2.5,
  units = "cm"
  )

power_legacy_study_A <- study_results_plot(
  legacy_study_results = CCG_legacy_study_A,
  resampling_results = VCG_legacy_study_A,
  keep_consistent_results = T,
  color_by = "power"
)

ggsave(
  plot = power_legacy_study_A,
  filename = paste0(path_res, "/power_study_A_no_imp.png"),
  dpi = 600,
  width = 20,
  height = 30,
  scale = 2.5,
  units = "cm"
)


effsize_VCG_study_A_median_imp <- study_results_plot(
  legacy_study_results = CCG_legacy_study_A,
  resampling_results = VCG_legacy_study_A_imp_median,
  keep_consistent_results = F
)
ggsave(
  plot = effsize_VCG_study_A_median_imp,
  filename = paste0(path_res, "/effsize_study_A_median_imp.png"),
  dpi = 600,
  width = 20,
  height = 30,
  scale = 2.5,
  units = "cm"
)

effsize_VCG_study_A_rs_imp <- study_results_plot(
  legacy_study_results = CCG_legacy_study_A,
  resampling_results = VCG_legacy_study_A_imp_rs,
  keep_consistent_results = F
)
ggsave(
  plot = effsize_VCG_study_A_rs_imp,
  filename = paste0(path_res, "/effsize_study_A_rs_imp.png"),
  dpi = 600,
  width = 20,
  height = 30,
  scale = 2.5,
  units = "cm"
)

effsize_VCG_study_A_pmm_imp <- study_results_plot(
  legacy_study_results = CCG_legacy_study_A,
  resampling_results = VCG_legacy_study_A_imp_pmm,
  keep_consistent_results = F
)
ggsave(
  plot = effsize_VCG_study_A_pmm_imp,
  filename = paste0(path_res, "/effsize_study_A_pmm_imp.png"),
  dpi = 600,
  width = 20,
  height = 30,
  scale = 2.5,
  units = "cm"
)
#*******************************************************************************
#*Legacy study B
#*******************************************************************************
effsize_legacy_study_B <- study_results_plot(
  legacy_study_results = CCG_legacy_study_B,
  resampling_results = VCG_legacy_study_B,
  keep_consistent_results = F
)
ggsave(
  plot = effsize_legacy_study_B,
  filename = paste0(path_res, "/effsize_study_B_no_imp.png"),
  dpi = 600,
  width = 20,
  height = 30,
  scale = 2.5,
  units = "cm"
)

power_legacy_study_B <- study_results_plot(
  legacy_study_results = CCG_legacy_study_B,
  resampling_results = VCG_legacy_study_B,
  keep_consistent_results = T,
  color_by = "power"
)

ggsave(
  plot = power_legacy_study_B,
  filename = paste0(path_res, "/power_study_B_no_imp.png"),
  dpi = 600,
  width = 20,
  height = 30,
  scale = 2.5,
  units = "cm"
)

power_legacy_study_B_recovery <- study_results_plot(
  legacy_study_results = CCG_legacy_study_B_recovery,
  resampling_results = VCG_legacy_study_B_recovery,
  keep_consistent_results = T,
  color_by = "power"
)

ggsave(
  plot = power_legacy_study_B_recovery,
  filename = paste0(path_res, "/power_study_B_recovery_no_imp.png"),
  dpi = 600,
  width = 20,
  height = 30,
  scale = 2.5,
  units = "cm"
)

effsize_VCG_study_B_median_imp <- study_results_plot(
  legacy_study_results = CCG_legacy_study_B,
  resampling_results = VCG_legacy_study_B_imp_median,
  keep_consistent_results = F
)
ggsave(
  plot = effsize_VCG_study_B_median_imp,
  filename = paste0(path_res, "/effsize_study_B_median_imp.png"),
  dpi = 600,
  width = 20,
  height = 30,
  scale = 2.5,
  units = "cm"
)

effsize_VCG_study_B_rs_imp <- study_results_plot(
  legacy_study_results = CCG_legacy_study_B,
  resampling_results = VCG_legacy_study_B_imp_rs,
  keep_consistent_results = F
)
ggsave(
  plot = effsize_VCG_study_B_rs_imp,
  filename = paste0(path_res, "/effsize_study_B_rs_imp.png"),
  dpi = 600,
  width = 20,
  height = 30,
  scale = 2.5,
  units = "cm"
)

effsize_VCG_study_B_pmm_imp <- study_results_plot(
  legacy_study_results = CCG_legacy_study_B,
  resampling_results = VCG_legacy_study_B_imp_pmm,
  keep_consistent_results = F
)
ggsave(
  plot = effsize_VCG_study_B_pmm_imp,
  filename = paste0(path_res, "/effsize_study_B_pmm_imp.png"),
  dpi = 600,
  width = 20,
  height = 30,
  scale = 2.5,
  units = "cm"
)

#*******************************************************************************
#*Legacy study C
#*******************************************************************************
effsize_legacy_study_C <- study_results_plot(
  legacy_study_results = CCG_legacy_study_C,
  resampling_results = VCG_legacy_study_C,
  keep_consistent_results = F
)
ggsave(
  plot = effsize_legacy_study_C,
  filename = paste0(path_res, "/effsize_study_C_no_imp.png"),
  dpi = 600,
  width = 20,
  height = 30,
  scale = 2.5,
  units = "cm"
)

power_legacy_study_C <- study_results_plot(
  legacy_study_results = CCG_legacy_study_C,
  resampling_results = VCG_legacy_study_C,
  keep_consistent_results = T,
  color_by = "power"
)

ggsave(
  plot = power_legacy_study_C,
  filename = paste0(path_res, "/power_study_C_no_imp.png"),
  dpi = 600,
  width = 20,
  height = 30,
  scale = 2.5,
  units = "cm"
)

effsize_VCG_study_C_median_imp <- study_results_plot(
  legacy_study_results = CCG_legacy_study_C,
  resampling_results = VCG_legacy_study_C_imp_median,
  keep_consistent_results = F
)
ggsave(
  plot = effsize_VCG_study_C_median_imp,
  filename = paste0(path_res, "/effsize_study_C_median_imp.png"),
  dpi = 600,
  width = 20,
  height = 30,
  scale = 2.5,
  units = "cm"
)

effsize_VCG_study_C_rs_imp <- study_results_plot(
  legacy_study_results = CCG_legacy_study_C,
  resampling_results = VCG_legacy_study_C_imp_rs,
  keep_consistent_results = F
)
ggsave(
  plot = effsize_VCG_study_C_rs_imp,
  filename = paste0(path_res, "/effsize_study_C_rs_imp.png"),
  dpi = 600,
  width = 20,
  height = 30,
  scale = 2.5,
  units = "cm"
)

effsize_VCG_study_C_pmm_imp <- study_results_plot(
  legacy_study_results = CCG_legacy_study_C,
  resampling_results = VCG_legacy_study_C_imp_pmm,
  keep_consistent_results = F
)
ggsave(
  plot = effsize_VCG_study_C_pmm_imp,
  filename = paste0(path_res, "/effsize_study_C_pmm_imp.png"),
  dpi = 600,
  width = 20,
  height = 30,
  scale = 2.5,
  units = "cm"
)