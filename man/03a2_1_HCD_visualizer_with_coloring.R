#*03a2_1_HCD_visualizer_with_coloring.R
#*
#*This script is a version of 03a2_HCD_visualizer.R. But here, instead of
#*batch-plotting everything as a box plot and/or as a histogram, you can
#*select the x-axis of the boxes and color it accordingly. For instance, you
#*can color the boxes so that they represent certain vehicles used in a study.
#*
#*This script creates histograms and box plots (denounced histbox) for each
#*parameter of interest. This script is used for quality control of HCD.
#*As an input, historical control data is fed into the function and then
#*all parameters are visualized below each other creating a long HTML list.
#*
#*As an additional quality control, studies are highlighted if their quartiles
#*are outside of the "grand median" value (i.e., the study seems to be different
#*towards all other studies).

##### program starts here
#*******************************************************************************
#*******************************************************************************
#*load libraries----
#*******************************************************************************
require(tidyverse)
require(gridExtra)
#*******************************************************************************
#*******************************************************************************
#*Begin of function
HCD_visualizer_with_coloring <- function(data_to_plot, plottitle, color_by = NULL){
  #*****************************************************************************
  #*****************************************************************************
  #*Prepare data----
  #*****************************************************************************
  #Arrange data by year
  data_to_plot <- data_to_plot %>% arrange(START_YEAR_SPREFID)
  
  #calculate median value of each parameter
  data_with_param_median <- data_to_plot %>%
    group_by(SEX, LBTESTCD, LBORRESU, LBSPEC) %>%
    mutate(
      param_median = median(LBORRES)
    ) %>%
    ungroup() %>%
    #*calculate percentiles of each study to check whether the median value
    #*over all parameters is within or outside of the quartiles of a study
    group_by(SEX, LBTESTCD, LBORRESU, LBSPEC, SPREFID) %>%
    mutate(
      lower_quartile = quantile(LBORRES)[["25%"]],
      upper_quartile = quantile(LBORRES)[["75%"]],
      out_of_range = case_when(
        lower_quartile > param_median | upper_quartile < param_median ~ TRUE,
        TRUE ~ FALSE
      )
    ) %>%
    ungroup() %>%
    group_by(SEX, LBTESTCD, LBORRESU, LBSPEC) %>%
    mutate(
      #*Add a sign which is going to be attached to the axis title signifying
      #*that there is some study with abnormally high or low values
      abnormality_text = case_when(
        any(out_of_range) == TRUE ~ "(!)",
        TRUE == FALSE ~ ""
      ))
  
  #Split the data by SEX, LBTESTCD, LBSPEC, and unit to make separate plots
  data_to_plot_list <- split(
    data_with_param_median,
    list(
      data_with_param_median$LBTESTCD, 
      data_with_param_median$LBORRESU,
      data_with_param_median$LBSPEC
    ),
    drop = T
  )
  
  #Create a nested list out of the list and split by SEX
  nested_list <- lapply(
    data_to_plot_list,
    function(x){split(x, x$SEX, drop = T)}
  )
  
  # Function to generate the plots
  plot_data <- function(data_list, param_name){
    
    if("M" %in% names(data_list)){
      df_males <- data_list$M
      df_males$sex <- "Male"
    } else {
      df_males <- data.frame()
    }
    
    if("F" %in% names(data_list)){
      df_females <- data_list$`F`
      df_females$sex <- "Female"
    } else {
      df_females <- data.frame()
    }
    
    df <- rbind.data.frame(df_males, df_females)
    df$param <- param_name
    
    df
  }
  
  nested_list_df <- imap_dfr(nested_list, plot_data)
  
  # Open the PDF device
  pdf(paste0(plottitle, ".pdf"), width = 28, height = 8)
  # Loop through the parameters and generate plots for each parameter
  for(param_name in unique(nested_list_df$param)){
    
    df_param <- nested_list_df %>% filter(param == param_name)
    
    #calculate population in parameter for each sex
    population_m <- df_param %>% filter(SEX == "M") %>% nrow()
    population_f <- df_param %>% filter(SEX == "F") %>% nrow()
    
    #Determine min and max 'value' for current parameter
    min_val <- min(df_param$LBORRES, na.rm = TRUE)
    max_val <- max(df_param$LBORRES, na.rm = TRUE)
    
    #Determine binwidth automatically for each parameter
    breaks <- pretty(
      range(df_param$LBORRES), n = nclass.FD(df_param$LBORRES), min.n = 1
      )
    bwidth <- breaks[2]-breaks[1]
    
    #Create automatic filling if the color_fill is not empty
    fill_mapping <- if(!is.null(color_by)) aes_string(fill = color_by) else NULL
    #***************************************************************************
    #***************************************************************************
    #*Create the plots----
    #***************************************************************************
    p_male_box <- ggplot(
      df_param %>% filter(sex == "Male"),
      aes(
        x = START_YEAR_SPREFID,
        y = LBORRES
        )
      ) +
      fill_mapping +
      geom_boxplot() +
      ggtitle("Boxplot (Male)") +
      scale_y_continuous(limits = c(min_val, max_val)) +
      theme(axis.text.x = element_text(angle = 90)) +
      theme_minimal()
    #***************************************************************************
    p_female_box <- ggplot(
      df_param %>% filter(sex == "Female"),
      aes(
        x = START_YEAR_SPREFID,
        y = LBORRES
        )
      ) +
      fill_mapping +
      geom_boxplot() +
      ggtitle("Boxplot (Female)") +
      scale_y_continuous(limits = c(min_val, max_val)) +
      theme(axis.text.x = element_text(angle = 180)) +
      theme_minimal()
    
    # Arrange the plots into a grid, and add a title for the parameter
    gridExtra::grid.arrange(
      # p_male_hist, p_male_box, p_female_hist, p_female_box,
      p_male_box, p_female_box,
      top = grid::textGrob(
        param_name, gp = grid::gpar(fontsize = 16, fontface = "bold")
        ),
      nrow = 1)
    
  } %>%
    suppressWarnings()

  #***************************************************************************
  # Close the PDF device
  dev.off()
}
