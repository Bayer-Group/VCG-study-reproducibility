#*06a3_BW_plot.R
#*
#*This script aims to visualize the body weight development of animals with as
#*a line plot with respect to the study day.
#*Each dose group is shown as a mean value along with their standard deviations
#*as error bars. In case if a dose group has a statistically significant
#*difference to the control on a specific day, this group (on this day) is
#*marked with an asterisk.
#*Same plot is done for body weight gain (BG) instead of BW. But now, instead
#*of body weight, the gain is plotted there.
#*
#*Two plots are shown side by side: on the left, the plot with the concurrent
#*control groups of the legacy study. On the right, the same legacy study, but
#*with the mean values of all iterations from the VCGs.
#*
#*Input:  CCG_legacy_study_X and VCG_legacy_X obtained from 04_resampling.R.
#*        From there, body weight values are extracted.
#*Steps:  - Extract the body weight values.
#*        - Split by sex, and CCG/VCG.
#*        - Visualize everything.
#*Output: JPEG image of body weight development of the animals.
#*******************************************************************************
#*******************************************************************************
#*load libraries----
#*******************************************************************************
require(data.table)
require(tidyverse)
require(gridExtra)
#*******************************************************************************
#*******************************************************************************
#*Create function to visualize body weight curves----
#*******************************************************************************
BW_plot <- function(
    CCG_study,
    VCG_study,
    parameter,
    BWDY_thershold = 28
    ){
  #*****************************************************************************
  #*****************************************************************************
  #*Read BW data from studies----
  #*****************************************************************************
  data_extractor <- function(resdata, CCGVCG) {
    resdata %>%
    #security step: ungroup data frames if they were grouped beofre together
    ungroup() %>%
    #*extract all parameters carrying the character "BW_D" (body weight_day)
    #*or "BG_D" (body weight gain day). Specify by function argument.
    filter(grepl(parameter, LBTESTCD)) %>%
      select(
        LBTESTCD, LBSPEC, LBORRESU, starts_with(c("CG", "LD", "MD", "HD"))
      ) %>%
      #pivot into long structure separated by SEX and DOSE_GROUP
      pivot_longer(
        cols = starts_with(c("CG", "LD", "MD", "HD")),
        names_to = "DOSE_SEX",
        values_to = "values"
      ) %>%
      mutate(
        #split the values column into mean and SD
        #get the mean values
        mean_values = parse_number(str_extract(values, "[^\u00B1]+")),
        #get sd values
        sd_values = parse_number(str_extract(values, "\u00B1 (.+)")),
        #get boolean value to check if there was a significance or not
        significance = if_else(grepl("\\*", values), TRUE, FALSE),
        #extract study daty from LBTESTCD
        BWDY = parse_number(LBTESTCD),
        #split DOSE_SEX into dose group and sex, which are spearated by "_"
        #get dose group
        trial_set_description = str_extract(DOSE_SEX, "[^_]+"),
        #turn dose group into factor
        trial_set_description = factor(
          trial_set_description, levels = c("CG", "LD", "MD", "HD")
        ),
        #get animal sex
        SEX = str_sub(DOSE_SEX, -1),
        #write animal sex as whole words for better visualization
        SEX = if_else(SEX == "M", "males", "females"),
        #note whether this set belongs to CCG or to VCG
        #make sex into a factor so that males are shown left of females
        SEX = factor(SEX, levels = c("males", "females")),
        CCGVCG = CCGVCG
      ) %>%
      #remove data from days past the maximum selected threshold
      filter(BWDY <= BWDY_thershold) %>%
      #drop epmty values
      filter(!is.na(mean_values)) %>%
      #drop the character values column and LBTESTCD
      select(-values, -LBTESTCD, -LBSPEC, -DOSE_SEX) %>%
      #rename LBORRESU to BWORRESU
      rename(BWORRESU = LBORRESU) %>%
      #*a warning message may occur if the sd values are empty. This may happen
      #*if there weren't enough animals to compute an sd value. The message
      #*is dropped.
      suppressWarnings()
  }
  #*****************************************************************************
  #*extract BW data from CCG and VCG and combine into a list
  BW_data <- rbind.data.frame(
    data_extractor(CCG_study[[1]], "CCG"),
    data_extractor(VCG_study[[1]], "VCG")
  )
  #split data frame by sex
  BW_data_list <- split(BW_data, BW_data$SEX)
  #*****************************************************************************
  #*****************************************************************************
  #*Create a plot title and an y-axis title----
  #*****************************************************************************
  plottitle <- ifelse(
    parameter == "BW_D", ", body weight", ", body weight gain"
  )
  
  yaxistitle <- ifelse(
    parameter == "BW_D", "Body weight [", "Body weight gain ["
  )
  #*****************************************************************************
  #*****************************************************************************
  #*Create the plot----
  #*****************************************************************************
  BW_growth_plot <- lapply(BW_data_list, function(which_subdata){
    ggplot(
      #select data to plot
      which_subdata,
      #make study day as x axis and mean value as y axis. separate by dose_group
      aes(
        x = BWDY,
        y = mean_values,
        color = trial_set_description,
        shape = trial_set_description
        )
      ) +
      #make a line plot out of this
      geom_line(linewidth = .4) +
      #add point to the lines so that they have different shapes too
      geom_point(size = 1) +
      #add sd as error bars
      geom_errorbar(
        aes(ymin = mean_values - sd_values, ymax = mean_values + sd_values),
        width = .3,
        size = .3
        ) +
      #upon significant difference, mark the point as an "*"
      geom_text(
        data = subset(which_subdata, significance),
        aes(label = "*"),
        size = 5,
        color = "#000000",
        fontface = "bold"
        ) +
      #separate by sex so the plots are side by side
      facet_grid(CCGVCG ~ .) +
      #*create a dynamic title which changes depending on sex and whether you
      #*observe CCG results or VCG results
      labs(
        title = paste0(unique(which_subdata$SEX), plottitle),
        x = "Study day",
        y = paste0(yaxistitle, unique(which_subdata$BWORRESU), "]"),
        color = "Dose group",
        shape = "Dose group"
        ) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5))
  })
  
  #Arrange the plots into a grid
  result_plot <- grid.arrange(grobs = BW_growth_plot, ncol = 2)
  
  return(result_plot)
}
#*******************************************************************************
#*******************************************************************************
#*write the plots----
#*******************************************************************************
#*Legacy study A
#*body weight
ggsave2(
  paste0(path_res, "/legacy_study_A_BW.jpeg"),
  BW_plot(CCG_legacy_study_A, VCG_legacy_study_A, "BW_D"),
  width = 20,
  height = 10,
  units = "cm",
  dpi = 300
)
#*body weight gain
ggsave2(
  paste0(path_res, "/legacy_study_A_BG.jpeg"),
  BW_plot(CCG_legacy_study_A, VCG_legacy_study_A, "BG_D"),
  width = 20,
  height = 10,
  units = "cm",
  dpi = 300
)
#*******************************************************************************
#*Legacy study B
#*body weight
ggsave2(
  paste0(path_res, "/legacy_study_B_BW.jpeg"),
  BW_plot(CCG_legacy_study_B, VCG_legacy_study_B, "BW_D"),
  width = 20,
  height = 10,
  units = "cm",
  dpi = 300
)
#*body weight gain
ggsave2(
  paste0(path_res, "/legacy_study_B_BG.jpeg"),
  BW_plot(CCG_legacy_study_B, VCG_legacy_study_B, "BG_D"),
  width = 20,
  height = 10,
  units = "cm",
  dpi = 300
)

#combine BW and BG into one plot and label
ggsave2(
  paste0(path_res, "/legacy_study_B_BW_BG.jpeg"),
  plot_grid(
    BW_plot(CCG_legacy_study_B, VCG_legacy_study_B, "BW_D"),
    BW_plot(CCG_legacy_study_B, VCG_legacy_study_B, "BG_D"),
    ncol = 1,
    labels = c("A", "B")
  ),
  width = 20,
  height = 20,
  units = "cm",
  dpi = 300
)

#*******************************************************************************
#*Legacy study C
#*body weight
ggsave2(
  paste0(path_res, "/legacy_study_C_BW.jpeg"),
  BW_plot(CCG_legacy_study_C, VCG_legacy_study_C, "BW_D"),
  width = 20,
  height = 10,
  units = "cm",
  dpi = 300
)
#*body weight gain
ggsave2(
  paste0(path_res, "/legacy_study_C_BG.jpeg"),
  BW_plot(CCG_legacy_study_C, VCG_legacy_study_C, "BG_D"),
  width = 20,
  height = 10,
  units = "cm",
  dpi = 300
)
