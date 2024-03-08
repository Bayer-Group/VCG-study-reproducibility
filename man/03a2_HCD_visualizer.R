#*03a2_HCD_visualizer.R
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
HCD_visualizer <- function(
    data_to_plot = HCD
    , 
    plottitle = "HCD_no_imputation"
    ){
  #*****************************************************************************
  #*****************************************************************************
  #*Prepare data----
  #*****************************************************************************
  #Arrange data by year
  data_to_plot <- data_to_plot %>% arrange(START_YEAR_SPREFID)
  #*****************************************************************************
  #*Detect extreme studies----
  #*****************************************************************************
  #*As an "extreme study, I define every study, whose median value (of each
  #*respective parameter) is outside of the 1.5*IQR of the overall data.
  
  #*Calculate the interquartile range across all studies
  data_with_param_median <- data_to_plot %>%
    group_by(SEX, LBTESTCD, LBORRESU, LBSPEC) %>%
    mutate(
      #calculate interquartile range (IQR)
      grand_IQR = IQR(LBORRES),
      #calculate median
      grand_median = median(LBORRES),
      #Calculate the 1.5*IQR ranges which are used to mark extreme studies
      grand_lower_fence = grand_median - 1.5 * grand_IQR,
      grand_upper_fence = grand_median + 1.5 * grand_IQR
      ) %>%
    ungroup() %>%
    #check whether the median of each study is within the IQR over all studies
    group_by(SEX, LBTESTCD, LBORRESU, LBSPEC, SPREFID) %>%
    mutate(
      #calculate median of each study
      study_median = median(LBORRES),
      #Flag if the study median is outside of the grand fences
      extreme_flag = if_else(
        between(study_median, grand_lower_fence, grand_upper_fence),
        FALSE,
        TRUE
      )
      )
  #*****************************************************************************
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
  pdf(paste0(plottitle, ".pdf"), width = 14, height = 8)
  
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
    breaks <- pretty(range(df_param$LBORRES), n = nclass.FD(df_param$LBORRES), min.n = 1)
    bwidth <- breaks[2]-breaks[1]
    #***************************************************************************
    #*plots for males----
    #***************************************************************************
    #Create histogram for males
    p_male_hist <- ggplot(
      df_param %>% filter(sex == "Male"),
      #do the plot sideways, so make y the variable
      aes(y = LBORRES)
      ) + 
      geom_histogram(
        aes(x = ..density..),
        binwidth = bwidth,
        color = "#000000",
        fill = "#ffffff"
        ) +
      #add a density plot
      geom_density(alpha = .2, fill = "#ff6666", aes(x = ..density..)) +
      #add a title
      ggtitle("Histogram (Male)") +
      #*Scale from lowest to highest value (in both sexes so that the axes)
      #*are the same
      scale_y_continuous(limits = c(min_val, max_val)) +
      scale_x_reverse() +
      theme_minimal() +
      #Add the population count to the data
      annotation_custom(
        grob = textGrob(
          paste0("n = ", population_m),
          hjust = 0, vjust = 1, gp = gpar(col = "#000000")
          ), 
        xmin = -Inf, xmax = -Inf, ymin = Inf, ymax = Inf
      )
    #***************************************************************************
    #Create box plot for males
    p_male_box <- ggplot(
      df_param %>% filter(sex == "Male"),
      aes(x = START_YEAR_SPREFID, y = LBORRES)
      ) +
      #make a box plot and fill in extreme values
      geom_boxplot(aes(fill = extreme_flag)) +
      #Add dashed lines for the grand median and fences
      geom_hline(
        yintercept = df_param %>%
          filter(sex == "Male") %>% pull(grand_median) %>% unique(),
        linetype = "dashed"
        ) +
      geom_hline(
        yintercept = df_param %>%
          filter(sex == "Male") %>% pull(grand_lower_fence) %>% unique(),
        linetype = "dashed"
        ) +
      geom_hline(
        yintercept = df_param %>%
          filter(sex == "Male") %>% pull(grand_upper_fence) %>% unique(),
        linetype = "dashed"
        ) +
      #Add annotations for the lines
      annotate(
        "text", x = max(df_param$START_YEAR_SPREFID),
        y = df_param %>%
          filter(sex == "Male") %>% pull(grand_median) %>% unique(),
        label = "Median", vjust = -1, hjust = 1
        ) +
      annotate(
        "text", x = max(df_param$START_YEAR_SPREFID),
        y = df_param %>%
          filter(sex == "Male") %>% pull(grand_lower_fence) %>% unique(),
        label = "Lower fence", vjust = -1, hjust = 1
        ) +
      annotate(
        "text", x = max(df_param$START_YEAR_SPREFID),
        y = df_param %>%
          filter(sex == "Male") %>% pull(grand_upper_fence) %>% unique(),
        label = "Upper fence", vjust = -1, hjust = 1
        ) +
      #color the boxes if extreme
      scale_fill_manual(values = c("FALSE" = "#ffffff", "TRUE" = "#fb9a99")) +
      ggtitle("Boxplot (Male)") +
      #*Scale from lowest to highest value (in both sexes so that the axes)
      #*are the same
      scale_y_continuous(limits = c(min_val, max_val)) +
      theme_minimal() +
      #rotate axis text
      theme(axis.text.x = element_text(angle = 90))
    #***************************************************************************
    #*plots for females----
    #***************************************************************************
    #Create histogram for females
    p_female_hist <- ggplot(
      df_param %>% filter(sex == "Female"),
      #do the plot sideways, so make y the variable
      aes(y = LBORRES)
    ) + 
      geom_histogram(
        aes(x = ..density..),
        binwidth = bwidth,
        color = "#000000",
        fill = "#ffffff"
      ) +
      #add a density plot
      geom_density(alpha = .2, fill = "#ff6666", aes(x = ..density..)) +
      #add a title
      ggtitle("Histogram (Female)") +
      #*Scale from lowest to highest value (in both sexes so that the axes)
      #*are the same
      scale_y_continuous(limits = c(min_val, max_val)) +
      scale_x_reverse() +
      theme_minimal() +
      #Add the population count to the data
      annotation_custom(
        grob = textGrob(
          paste0("n = ", population_m),
          hjust = 0, vjust = 1, gp = gpar(col = "#000000")
        ), 
        xmin = -Inf, xmax = -Inf, ymin = Inf, ymax = Inf
      )
    #***************************************************************************
    #Create box plot for females
    p_female_box <- ggplot(
      df_param %>% filter(sex == "Female"),
      aes(x = START_YEAR_SPREFID, y = LBORRES)
    ) +
      #make a box plot and fill in extreme values
      geom_boxplot(aes(fill = extreme_flag)) +
      #Add dashed lines for the grand median and fences
      geom_hline(
        yintercept = df_param %>%
          filter(sex == "Female") %>% pull(grand_median) %>% unique(),
        linetype = "dashed"
      ) +
      geom_hline(
        yintercept = df_param %>%
          filter(sex == "Female") %>% pull(grand_lower_fence) %>% unique(),
        linetype = "dashed"
      ) +
      geom_hline(
        yintercept = df_param %>%
          filter(sex == "Female") %>% pull(grand_upper_fence) %>% unique(),
        linetype = "dashed"
      ) +
      #Add annotations for the lines
      annotate(
        "text", x = max(df_param$START_YEAR_SPREFID),
        y = df_param %>%
          filter(sex == "Female") %>% pull(grand_median) %>% unique(),
        label = "Median", vjust = -1, hjust = 1
      ) +
      annotate(
        "text", x = max(df_param$START_YEAR_SPREFID),
        y = df_param %>%
          filter(sex == "Female") %>% pull(grand_lower_fence) %>% unique(),
        label = "Lower fence", vjust = -1, hjust = 1
      ) +
      annotate(
        "text", x = max(df_param$START_YEAR_SPREFID),
        y = df_param %>%
          filter(sex == "Female") %>% pull(grand_upper_fence) %>% unique(),
        label = "Upper fence", vjust = -1, hjust = 1
      ) +
      #color the boxes if extreme
      scale_fill_manual(values = c("FALSE" = "#ffffff", "TRUE" = "#fb9a99")) +
      ggtitle("Boxplot (Female)") +
      #*Scale from lowest to highest value (in both sexes so that the axes)
      #*are the same
      scale_y_continuous(limits = c(min_val, max_val)) +
      theme_minimal() +
      #rotate axis text
      theme(axis.text.x = element_text(angle = 90))
    #***************************************************************************
    # Arrange the plots into a grid, and add a title for the parameter
    gridExtra::grid.arrange(
      p_male_hist, p_male_box, p_female_hist, p_female_box,
      top = grid::textGrob(
        param_name, gp = grid::gpar(fontsize = 16, fontface = "bold")
        ),
      nrow = 1
      )
  } %>%
    suppressWarnings()
  
  # Close the PDF device
  dev.off()
}