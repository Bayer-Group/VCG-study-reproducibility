#*06a2_MI_plot.R
#*
#*This script aims to visualize the percentage of microscopic findings per
#*dose group of a legacy study. Once with the original concurrent controls
#*(CCG) and once with all historical control data (HCD). The goal of this plot
#*is to gain an understanding in "how the number of organ findings change" if
#*using virtual control groups (VCGs) instead of CCGs.
#*
#*Input:  legacy study data and historical control data as CSV
#*
#*Steps:  - match HCD to the initial body weight (INITBW) of the dose group
#*         animals of the legacy study.
#*        - visualize the original findings (for liver and kidney as the
#*         predominant organs in toxicological assessments)
#*        - do the same visualization but with HCDs instead of CCGs.
#*        - we use percentages instead of absolute number of findings to
#*         account for the imbalance of HCD compared to dose group animals.
#*        - males and females are regarded separately.
#*Output: ggplot plot exported as JPEG and PDF

#*******************************************************************************
#*******************************************************************************
#*load libraries----
#*******************************************************************************
require(data.table)
require(tidyverse)
require(gridExtra)
#*******************************************************************************
#*******************************************************************************
#*create visualizer function----
#*******************************************************************************
mi_visualizer <- function(study, HCD_path, target_organ, titletext_study){
  #*****************************************************************************
  #*****************************************************************************
  #*read data----
  #*****************************************************************************
  #read legacy study MI findings
  mi_legacy <- fread(paste0(
    rootpath, "/data/Original/legacy_studies/", study, "/mi.csv"
  )) %>%
    #filter for target organ
    filter(MISPEC == target_organ) %>%
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
      #change the name of SEX into full names, so it looks nicer on the plots
      SEX = if_else(SEX == "F", "Females", "Males")
      )
  
  #*split data frames into dose groups and control group. They will be combined
  #*later in the plot
  mi_CCG <- mi_legacy %>% filter(trial_set_description == "CG")
  mi_dose <- mi_legacy %>% filter(trial_set_description != "CG")
  
  #read HCD MI findings
  mi_HCD <- fread(paste0(
    rootpath, "/data/Original/", HCD_path, "/mi.csv"
  ))  %>%
    #filter for target organ
    filter(MISPEC == target_organ) %>%
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
      #change the name of SEX into full names, so it looks nicer on the plots
      SEX = if_else(SEX == "F", "Females", "Males")
    )
  #****************************************************************************
  #calculate sample size per sex and per dose group
  sample_size_calculator <- function(data){
    data %>%
      group_by(SEX, trial_set_description) %>%
      summarize(population = n(), .groups = "drop")
  }
  
  CCG_size <- sample_size_calculator(mi_CCG)
  dose_size <- sample_size_calculator(mi_dose)
  HCD_size <- sample_size_calculator(mi_HCD)
  #****************************************************************************
  #*calculate percentage of findings (i.e., out of n animals, how many showed
  #*a finding with the severity grade k?)
  percentage_calculator <- function(data, size_data){
    data %>%
      group_by(SEX, trial_set_description, MISEVTXT) %>%
      reframe(number_of_findings = n()) %>%
      left_join(size_data, by = c("SEX", "trial_set_description")) %>%
      mutate(percentage = round(number_of_findings / population * 100))
  }

  CCG_percentage <- percentage_calculator(mi_CCG, CCG_size)
  dose_percentage <- percentage_calculator(mi_dose, dose_size)
  HCD_percentage <- percentage_calculator(mi_HCD, HCD_size)
  
  #add color scaling and names to the plot
  color_labels <- c(
    "no finding", "minimal", "mild", "moderate", "marked", "severe"
    )
  color_values <- c(
    "#c2decc", "#7bb9cf", "#5482b4", "#334a8f", "#172a63", "#051033"
    )
  
  names(color_values) <- color_labels
  #*****************************************************************************
  #*****************************************************************************
  #*Create function to visualize MI findings----
  #*****************************************************************************
  mi_plot <- function(plotdata, control_data, titletext_group){
    ggplot(
      data = plotdata %>% rbind.data.frame(control_data), #%>% filter(SEX == sex),
      aes(x = trial_set_description, y = percentage, fill = MISEVTXT)
    ) +
      geom_col(position = position_dodge2(width = .9, preserve = "single")) +
      #adding the population in each dose group as a text above the bars
      geom_text(
        aes(label = paste0("n = ", population), y = 100)
                ) +
      #custom colors
      scale_fill_manual(
        values = color_values,
        name = "Severity grade",
        # labels = color_labels
        # values = c(
        #   "#c2decc", "#7bb9cf", "#5482b4", "#334a8f", "#172a63", "#051033"
        #   )
      ) +
      #make scale of y axis from 0 to 100
      scale_y_continuous(
        limits = c(0, 100),
        #add a percent to each axis tick
        labels = function(x){paste0(x, " %")}
        ) +
      #relabel y axis
      ylab("Percentage of findings") +
      xlab("Trial set groups") +
      #split the plot into two grids separated by sex
      facet_wrap(~SEX, dir = "h") +
      #title of the plot. Add legacy study name to it
      ggtitle(paste0(
        "Percentage of microscopic findings in \"",
        titletext_study,
        "\" using ",
        titletext_group
        )) +
      theme_minimal() +
      theme(
        #add border to each facet strip
        strip.background = element_rect(
          fill = "#ffffff", colour = "#000000", linewidth = 1
          ),
        #increase the margin between the facet strips
        panel.spacing = unit(2, "lines")
      )
  }
  #*****************************************************************************
  #*****************************************************************************
  #*Create the plots for CCG and VCG----
  #*****************************************************************************
  CCG_MI_plot <- mi_plot(
    plotdata = dose_percentage,
    control_data = CCG_percentage,
    titletext_group = "concurrent control group (CCG)"
  )
  
  VCG_MI_plot <- mi_plot(
    plotdata = dose_percentage,
    control_data = HCD_percentage,
    titletext_group = "historical control data (HCD)"
  )

  #*****************************************************************************
  #*****************************************************************************
  #*Return the plots----
  #*****************************************************************************
  return(list(CCG_MI_plot, VCG_MI_plot))
  }

#*******************************************************************************
#*******************************************************************************
#*Call the function and compare to legacy studies----
#*******************************************************************************
#*Do this for liver----
#*******************************************************************************
#*Legacy study A
#*******************************************************************************
#*Create plot
mi_legacy_study_A <- mi_visualizer(
  study = "study_A",
  HCD_path = "HCD_RAT",
  target_organ = "LIVER",
  titletext_study = "Legacy study A"
  )
#*******************************************************************************
#*Export plot as JPEG
# open a new jpeg file
jpeg(
  paste0(path_res, "/legacy_study_A_MI_findings_liver.jpeg"),
  width = 1200, height = 600
  )

#Combine the plots
grid.arrange(
  mi_legacy_study_A[[1]], mi_legacy_study_A[[2]], ncol = 2
)

#Close the jpeg file
dev.off()
#*******************************************************************************
#*Legacy study B
#*******************************************************************************
#*Create plot
mi_legacy_study_B <- mi_visualizer(
  study = "study_B",
  HCD_path = "HCD_RAT",
  target_organ = "LIVER",
  titletext_study = "Legacy study B"
)
#*******************************************************************************
#*Export plot as JPEG
# open a new jpeg file
jpeg(
  paste0(path_res, "/legacy_study_B_MI_findings_liver.jpeg"),
  width = 1200, height = 600
)

#Combine the plots
grid.arrange(
  mi_legacy_study_B[[1]], mi_legacy_study_B[[2]], ncol = 2
)

#Close the jpeg file
dev.off()
#*******************************************************************************
#*Legacy study C
#*******************************************************************************
#*Create plot
mi_legacy_study_C <- mi_visualizer(
  study = "study_C",
  HCD_path = "HCD_RAT",
  target_organ = "LIVER",
  titletext_study = "Legacy study C"
)
#*******************************************************************************
#*Export plot as JPEG
# open a new jpeg file
jpeg(
  paste0(path_res, "/legacy_study_C_MI_findings_liver.jpeg"),
  width = 1200, height = 600
)

#Combine the plots
grid.arrange(
  mi_legacy_study_C[[1]], mi_legacy_study_C[[2]], ncol = 2
)

#Close the jpeg file
dev.off()
#*******************************************************************************
#*Do this for kidneys----
#*******************************************************************************
#*Legacy study A
#*******************************************************************************
#*Create plot
mi_legacy_study_A <- mi_visualizer(
  study = "study_A",
  HCD_path = "HCD_RAT",
  target_organ = "KIDNEY",
  titletext_study = "Legacy study A"
)
#*******************************************************************************
#*Export plot as JPEG
# open a new jpeg file
jpeg(
  paste0(path_res, "/legacy_study_A_MI_findings_kidney.jpeg"),
  width = 1200, height = 600
)

#Combine the plots
grid.arrange(
  mi_legacy_study_A[[1]], mi_legacy_study_A[[2]], ncol = 2
)

#Close the jpeg file
dev.off()
#*******************************************************************************
#*Legacy study B
#*******************************************************************************
#*Create plot
mi_legacy_study_B <- mi_visualizer(
  study = "study_B",
  HCD_path = "HCD_RAT",
  target_organ = "KIDNEY",
  titletext_study = "Legacy study B"
)
#*******************************************************************************
#*Export plot as JPEG
# open a new jpeg file
jpeg(
  paste0(path_res, "/legacy_study_B_MI_findings_kidney.jpeg"),
  width = 1200, height = 600
)

#Combine the plots
grid.arrange(
  mi_legacy_study_B[[1]], mi_legacy_study_B[[2]], ncol = 2
)

#Close the jpeg file
dev.off()
#*******************************************************************************
#*Legacy study C
#*******************************************************************************
#*Create plot
mi_legacy_study_C <- mi_visualizer(
  study = "study_C",
  HCD_path = "HCD_RAT",
  target_organ = "KIDNEY",
  titletext_study = "Legacy study C"
)
#*******************************************************************************
#*Export plot as JPEG
# open a new jpeg file
jpeg(
  paste0(path_res, "/legacy_study_C_MI_findings_kidney.jpeg"),
  width = 1200, height = 600
)

#Combine the plots
grid.arrange(
  mi_legacy_study_C[[1]], mi_legacy_study_C[[2]], ncol = 2
)

#Close the jpeg file
dev.off()
