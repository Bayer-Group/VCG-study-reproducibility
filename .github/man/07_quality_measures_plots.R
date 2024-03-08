#*07_quality_measures_plot.R

#*This script generates a ggplot showing significant results of the
#*legacy study and the concatenated results of the VCG resampling experiment.
#*The goal is to see whether the CCG values are strongly outside of the VCG range
#*resulting poor reproducibility
##### program starts here
#*****************************************************************************
#*****************************************************************************
#*load libraries----
#*****************************************************************************
require(tidyverse)
require(plotly)
#*****************************************************************************
#*Get the effect size results from each iteration
legacy_study <- legacy_study_C
HCD_unfiltered <- HCD

#filter for initial body weight. un-comment if needed
minmax_initbw <- legacy_study %>%
  filter(
    trial_set_description != "CG",
    LBTESTCD == "BW_D001"
  ) %>%
  #extract highest and lowest values of dose groups per sex
  group_by(SEX) %>%
  summarize(
    min_initbw = min(LBORRES),
    max_initbw = max(LBORRES)
  )

legacy_minmax <- split(minmax_initbw, minmax_initbw$SEX)

legacy_minmax_list <- lapply(
  legacy_minmax,
  function(x){
    x %>% select(-SEX) %>% as_vector()
  }
)
#***************************************************************************
#*Filter HCD so that animals are only selected who match INITBW of legacy
filtered_HCD_USUBJIDs <- split(
  HCD_unfiltered %>%
    filter(LBTESTCD == "BW_D001") %>%
    select(USUBJID, LBORRES),
  HCD_unfiltered %>%
    filter(LBTESTCD == "BW_D001") %>%
    pull(SEX)
)

#Select subject IDs from the animals within INITBW
filtered_by_initbw <- mapply(function(x,y){
  x %>%
    filter(between(LBORRES, y[1], y[2])) %>%
    pull(USUBJID)
},
x = filtered_HCD_USUBJIDs, y = legacy_minmax_list, SIMPLIFY = F
) %>%
  unlist()

#filter now only those animals who match in their INITBW
HCD_filtered <- HCD_unfiltered %>% filter(USUBJID %in% filtered_by_initbw)

#*****************************************************************************

HCD_filtered <- HCD_filtered %>% arrange(START_YEAR_SPREFID)

hist <- plot_ly(opacity = 0.5, y = ~LBORRES) %>%
  add_histogram(
    # ybins = list(size = xbinsize, start = xbinstart),
    data = HCD_filtered %>% filter(SEX == "F", LBTESTCD == "WC_W02"),
    marker = list(color = "#a5a5a5", line = list(color = "#000000", width = 2)),
    name = "HCD",
    # legendgroup = "anesthetic2",
    # showlegend = F#ifelse(sex == "M", T, F)
  ) %>%
  add_histogram(
    ybins = list(size = .1),
    data = legacy_study %>%
      filter(
        trial_set_description == "CG",
        SEX == "F",
        LBTESTCD == "WC_W02"
      ),
    marker = list(color = "#6a3d9a", line = list(color = "#000000", width = 2)),
    name = "legacy study",
    # legendgroup = "anesthetic1",
    # showlegend = F#ifelse(sex == "M", T, F)
  ) %>%
  layout(
    font = list(size = 18),
    yaxis = list(
      title = ~paste0(unique(LBTESTCD), " [", unique(LBORRESU),"]"),
      showline = T,
      showgrid = T,
      linewidth = 2,
      ticks = "outside",
      tickwidth = 2,
      ticklen = 10,
      mirror = T
      # range = c(axisrange_l, axisrange_u),
      # dtick = xticks
    ),
    xaxis = list(
      title = "Frequency",
      showline = T,
      showgrid = T,
      linewidth = 2,
      ticks = "outside",
      tickwidth = 2,
      ticklen = 10,
      mirror = T,
      autorange = "reversed"
    ),
    barmode = "stack"
  )

#box plot
box <- plot_ly(
  x = ~START_YEAR_SPREFID,
  y = ~LBORRES,
  opacity = 0.8) %>%
  #add values from CO2 subset
  add_boxplot(
    data = HCD_filtered %>% filter(LBTESTCD == "WC_W02", SEX == "F"),
    marker = list(color = "#a5a5a5",
                  size = 7, opacity = 0.5,
                  symbol = "triangle-up",
                  line = list(color = "#000000", width = 2)),
    line = list(color = "#000000", width = 2),
    fillcolor = "#a5a5a5",
    boxpoints = "all", jitter = "0.7", pointpos = "0",
    name = "HCD",
    # legendgroup = "anesthetic2",
    showlegend = F
  ) %>%
  add_boxplot(
    data = legacy_study %>%
      filter(
        LBTESTCD == "WC_W02",
        SEX == "F",
        trial_set_description == "CG"
      ),
    marker = list(color = "#6a3d9a",
                  size = 7, opacity = 0.5,
                  symbol = "pentagon",
                  line = list(color = "#000000", width = 2)),
    line = list(color = "#000000", width = 2),
    fillcolor = "#6a3d9a",
    boxpoints = "all", jitter = "0.7", pointpos = "0",
    name = "legacy study",
    # legendgroup = "anesthetic1",
    showlegend = T
  ) %>%
  #Add the count of each group
  add_annotations(
    text = ~paste0(
      "<i>n</i><sub>legacy study</sub>: ", legacy_study %>% filter(SEX == "F", trial_set_description == "CG", LBTESTCD == "WC_W02") %>% nrow(),
      "\n<i>n</i><sub>HCD</sub>: ", HCD_filtered %>% filter(SEX == "F", LBTESTCD == "WC_W02") %>% nrow()),
    size = 18,
    x = 0.95,
    y = 0.95,
    xref = "paper",
    yref = "paper",
    xanchor = "right",
    yanchor = "top",
    align = "right",
    showarrow = FALSE
  ) %>%
  layout(
    font = list(size = 18),
    xaxis = list(
      title = "Year of the study begin",
      tickmode = "linear",
      showline = T,
      showgrid = T,
      linewidth = 2,
      ticks = "outside",
      tickwidth = 2,
      ticklen = 10,
      mirror = T
    ),
    yaxis = list(
      title = "",#paste0(axistitle),
      showline = T,
      showgrid = T,
      showticklabels = F,
      linewidth = 2,
      ticks = "outside",
      tickwidth = 2,
      ticklen = 10,
      mirror = T
      # dtick = xticks,
      # range = c(axisrange_l, axisrange_u)
    )
  )


#align both plots into one subplot
histbox <-
subplot(
  hist,
  box,
  nrows = 1 ,shareY = F, shareX = F, titleX = T, titleY = T
)
    