#03a1_missing_data_calculator.R

#*Based on the historical control data (HCD) extracted in
#*"03_write_HCD_to_R_table.R", this script calculates the missing data per
#*measured endpoint and exports a table containing a list of all endpoints and
#*how well they are present in the data.

#*******************************************************************************
#*******************************************************************************
#*Load libraries----
#*******************************************************************************
require(data.table)
require(tidyverse)
require(gt)
require(webshot2)
#*******************************************************************************
#*******************************************************************************
#*Import HCD----
#*******************************************************************************
HCD <- HCD

#*Extract endpoints measured in each animal and pivot the data so that every
#*endpoint is a new column
HCD_endpoints <- HCD %>%
  select(USUBJID, LBORRES, LBTESTCD, LBSPEC) %>%
  pivot_wider(
    id_cols = USUBJID,
    names_from = c(LBTESTCD, LBSPEC),
    values_from = LBORRES,
    values_fn = function(x) paste(x, collapse=",")
  ) %>%
  #Transform back to original structure and keep the NAs
  pivot_longer(
    cols = c(everything(), -USUBJID),
    names_to = "LBTESTCD",
    values_to = "LBORRES"
  ) %>%
  #remove everything after second underscore if presend
  mutate(LBTESTCD = str_extract(LBTESTCD, "[^_]*_[^_]*")) %>%
  group_by(LBTESTCD) %>%
  summarize(
    na_count = sum(is.na(LBORRES)),
    total_population = HCD %>% pull(USUBJID) %>% n_distinct(),
    available = round(100 - (na_count / total_population * 100))
  ) %>%
  ungroup() %>%
  arrange(desc(available))
#*******************************************************************************
#*Do the same for qualitative HCD parameters----
HCD_qualitative <- HCD_qualitative

#*Extract endpoints measured in each animal and pivot the data so that every
#*endpoint is a new column
HCD_qual_endpoints <- HCD_qualitative %>%
  select(USUBJID, DOMAIN) %>%
  unique() %>%
  group_by(DOMAIN) %>%
  summarize(population = n(), .groups = "drop")
#*******************************************************************************
#*******************************************************************************
#*******************************************************************************
#*Create a table listing missing data----
#*******************************************************************************
missing_data_table <- HCD_endpoints %>%
  gt() %>%
  #Relable the columns
  cols_label(
    LBTESTCD = "Parameter short name",
    na_count = "Missing values",
    total_population = "Total population",
    available = "Measured data [%]"
  ) %>%
  tab_header(
    title = md("**Mising values**"),
    subtitle = md("Quantitative Parameters of Legacy Study")
  ) %>%
  #highlight data below 80 % percentage in red
  tab_style(
    style = cell_fill(color = "#fdba99"),
    locations = cells_body(columns = available, rows = available < 80)
  )

#*******************************************************************************
#*******************************************************************************
#*Export missing data restults as HTML, JPEG, and CSV----
#*******************************************************************************
fwrite(HCD_endpoints, paste0(rootpath, "/data/Derived/missing_data_table.csv"))
gtsave(missing_data_table, paste0(path_res, "/missing_data_table.html"))
gtsave(missing_data_table, paste0(path_res, "/missing_data_table.png"))
