#*10_BW_selction_visualizaton.R
#*
#*This is a standalone script which means it is not executed by 00_master.R
#*
#*This script creates a histogram showing the distribution of the initial body
#*weights of the HCD and the INITBW values of the dose groups.
#*
#*Input:  HCD_filtered_by_initbw, HCD, and legacy study from the script
#*        04_resampling.R
#*Steps:  - Read the data.
#*        - Plot a histogram (before filtering by INITBW) and add dashed lines
#*         illustrating where the mean, the SD and the ranges of the dose group
#*         INITW are.
#*        - Plot the histogram (after filtering by INITBW) with the same ranges
#*Output: Histogram plots
#*******************************************************************************
#*******************************************************************************
#*Load libraries----
#*******************************************************************************
require(tidyverse)
require(cowplot)
#*******************************************************************************
#*******************************************************************************
#*Beginning of function----
#*******************************************************************************
bw_selection_plot <- function(
    HCD_unfiltered,
    HCD_filtered,
    legacy_study_data,
    study_name
    ){
  #*****************************************************************************
  #*****************************************************************************
  #*Extract mean, sd, and range of legacy study
  #*****************************************************************************
  legacy_summary <- legacy_study_data %>%
    #extract the initial body weight
    filter(LBTESTCD == "BW_D001", trial_set_description != "CG") %>%
    #do this by sex
    group_by(SEX) %>%
    summarise(
      mean_value = mean(LBORRES, na.rm = TRUE),
      lower_2sd = mean(LBORRES, na.rm = TRUE) - 2 * sd(LBORRES, na.rm = TRUE),
      upper_2sd = mean(LBORRES, na.rm = TRUE) + 2 * sd(LBORRES, na.rm = TRUE),
      min_value = min(LBORRES, na.rm = TRUE),
      max_value = max(LBORRES, na.rm = TRUE)
    ) %>%
    #round everything
    mutate(
      across(c(everything(), -SEX), round),
      #write animal sex as whole words for better visualization
      SEX = if_else(SEX == "M", "Males", "Females"),
      SEX = factor(SEX, levels = c("Males", "Females"))
    )
  #*****************************************************************************
  #*****************************************************************************
  #join the HCD_unfiltered and HCD_filtered with the summarized legacy study----
  #*****************************************************************************
  #HCD unfiltered
  HCD_unfiltered_with_summary <- HCD_unfiltered %>%
    #write animal sex as whole words for better visualization
    mutate(
      SEX = if_else(SEX == "M", "Males", "Females"),
      SEX = factor(SEX, levels = c("Males", "Females"))
    ) %>%
    #extract the initial body weight
    filter(LBTESTCD == "BW_D001") %>%
    left_join(legacy_summary, by = "SEX")
  #*****************************************************************************
  #*Add counts of the data
  HCD_unfiltered_count <- HCD_unfiltered %>%
    #write animal sex as whole words for better visualization
    mutate(
      SEX = if_else(SEX == "M", "Males", "Females"),
      SEX = factor(SEX, levels = c("Males", "Females"))
    ) %>%
    select(USUBJID, SEX) %>%
    unique() %>%
    group_by(SEX) %>%
    summarize(population = n())
  #*****************************************************************************
  #*****************************************************************************
  #HCD pre-filtered by INITBW
  HCD_filtered_with_summary <- HCD_filtered %>%
    #write animal sex as whole words for better visualization
    mutate(
      SEX = if_else(SEX == "M", "Males", "Females"),
      SEX = factor(SEX, levels = c("Males", "Females"))
    ) %>%
    #extract the initial body weight
    filter(LBTESTCD == "BW_D001") %>%
    left_join(legacy_summary, by = "SEX")
  #*****************************************************************************
  #*Add counts of the data
  HCD_filtered_count <- HCD_filtered %>%
    #write animal sex as whole words for better visualization
    mutate(
      SEX = if_else(SEX == "M", "Males", "Females"),
      SEX = factor(SEX, levels = c("Males", "Females"))
      ) %>%
    select(USUBJID, SEX) %>%
    unique() %>%
    group_by(SEX) %>%
    summarize(population = n())
  #*****************************************************************************
  #*Create a count to select minimum and max value of both data frames so that
  #*the x-axes have the same limits in the end
  combined_xlim <- c(
    min(c(
      min(HCD_unfiltered_with_summary$LBORRES),
      min(HCD_filtered_with_summary$LBORRES)
      )) - 10,
    max(c(
      max(HCD_unfiltered_with_summary$LBORRES),
      max(HCD_filtered_with_summary$LBORRES)
      )) + 10
    )
  #*****************************************************************************
  #*****************************************************************************
  #*Generate histograms----
  #*****************************************************************************
  #*create a function to generate histograms and execute this function for the
  #*histograms before and after filtering
  histogram_plot <- function(which_HCD, HCD_count, plottitle){
    #For the unfiltered HCD
    which_HCD %>%
      #set LBORRES as x-axis values
      ggplot(aes(x = LBORRES)) +
      #define histogram and customize binwidth
      geom_histogram(binwidth = 5, fill = "#a5a5a5", color = "#000000") + 
      #add mean value of legacy dose groups as a vertical straight line
      geom_vline(
        aes(xintercept = mean_value), linewidth = .5
      ) +
      #add lower and upper 2*sd value of legacy dose groups as a dahed lines
      geom_vline(
        aes(xintercept = lower_2sd), linetype = "dashed", linewidth = .5
      ) +
      geom_vline(
        aes(xintercept = upper_2sd), linetype = "dashed", linewidth = .5
      ) +
      #add lower and upper range values of legacy dose groups as dotted lines
      geom_vline(
        aes(xintercept = min_value), linetype = "dotted", linewidth = .5
      ) +
      geom_vline(
        aes(xintercept = max_value), linetype = "dotted", linewidth = .5
      ) +
      #Create invisible lines in order to create a legend
      geom_line(
        aes(linetype = "Mean", y = 0), linewidth = 0
      ) +
      geom_line(
        aes(linetype = "Mean \u00B1 2*SD", y = 0), linewidth = 0
      ) +
      geom_line(
        aes(linetype = "Min-Max ranges", y = 0), linewidth = 0
      ) +
      #annotate the line types so that they appear in the legend
      scale_linetype_manual(
        paste0(study_name, "\nlocation parameters"),
        values = c(
          "Mean" = "solid", "Mean \u00B1 2*SD" = "dashed", "Min-Max ranges" = "dotted"
        )
      ) +
      #split histogram by sex
      facet_wrap(~SEX, as.table = FALSE) +
      #label the x and y axis and a title
      labs(
        title = plottitle,
        x = "Body weight [g]",
        y = "Count"
      ) +
      #minimal theme
      theme_minimal() +
      #set the limits of the x-axis
      xlim(combined_xlim) +
      #make the facet wrap text bigger
      theme(strip.text = element_text(face = "bold", size = 12)) +
      #make the legend elements bigger
      theme(legend.text = element_text(size = 10)) +
      #add counts of the data to the plot
      geom_text(
        data = HCD_count, 
        aes(label = paste("Population:", population), y = Inf, x = Inf), 
        hjust = 1.1, vjust = 1.1
      ) +
      #*manually override the legend. Otherwise saving the plot destroys the
      #*legend lines
      guides(
        linetype = guide_legend(
          override.aes = list(
            linetype = c("solid", "dashed", "dotted"),
            linewidth = .5
          )
        )
      )
    
  }
  #*****************************************************************************
  #*****************************************************************************
  #*Call histogram function----
  #*****************************************************************************
  #*Create histogram of HCD *BEFORE* filtering by INITBW of legacy dose groups
  histogram_before_filtering <- histogram_plot(
    HCD_unfiltered_with_summary,
    HCD_unfiltered_count,
    "A"# "HCD before filtering by INITBW of legacy dose groups"
    )
  #*****************************************************************************
  #*Create histogram of HCD *AFTER* filtering by INITBW of legacy dose groups
  histogram_after_filtering <- histogram_plot(
    HCD_filtered_with_summary,
    HCD_filtered_count,
    "B"# "HCD after filtering by INITBW of legacy dose groups"
  )
  
  #Arrange the plots into a grid
  result_plot <- plot_grid(
    histogram_before_filtering, histogram_after_filtering, ncol = 1
    )
  #*****************************************************************************
  #*****************************************************************************
  #write plot as JPEG----
  #*****************************************************************************
  ggsave2(
    #*change the legacy study name so it is consistent to the other plots in
    #*this project
    paste0(
      path_res, "/", study_name %>% gsub("L", "l", .) %>% gsub(" ", "_", .),
      "_INITBW_filtering_process.jpeg"
      ),
    result_plot,
    width = 20,
    height = 15,
    units = "cm",
    dpi = 300
  )
}
