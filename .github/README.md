# VCG-study-reproducibility
Reproduce statistical results of preclinical legacy studies. Replace concurrent controls of these studies with virtual control groups (VCG) and check whether the treatment relatedness can be reproduced.

## Overview
This repository consists of R scripts for replacing a legacy study's concurrent control group (CCG) with virtual control groups (VCGs) derived from historical control data (HCD).
This is followed by statistical reanalysis of the significance tests used in the legacy study to assess the performance of VCGs in their ability to reproduce the original results.
Furthermore, this code creates Excel sheets carrying the following information:
- Are the dose-groups within or outside the 5th and 95th percentile of the HCD?
- Is there a dose-dependency, i.e. is there an increase/a decrease of the mean value with increasing test-substance dose?
The mean values of the legacy study, once with the CCG and once with VCGs, are exported and can be evaluated by study directors and subject matter experts (SMEs) to assess whether the statistical significant findings are of actual biological relevance.
This helps in replicating the study-interpretation process by study directors and ultimately, supports in finding a study conclusion.

## Dependencies
In order to run the R-scripts the following R-packages were used:
- [`data.table`](https://CRAN.R-project.org/package=data.table)
- [`tidyverse`](https://doi.org/10.21105/joss.01686)
- [`openxlsx`](https://CRAN.R-project.org/package=openxlsx) 
- [`PMCMRplus`](https://CRAN.R-project.org/package=PMCMRplus) 
- [`rstatix`](https://CRAN.R-project.org/package=rstatix) 
- [`effsize`](https://CRAN.R-project.org/package=effsize) 
- [`gt`](https://CRAN.R-project.org/package=gt) 
- [`webshot2`](https://CRAN.R-project.org/package=webshot2) 
- [`mice`](https://doi.org/10.18637/jss.v045.i03) 
- [`janitor`](https://CRAN.R-project.org/package=janitor) 
- [`grid`](https://cran.r-project.org/package=grid) 
- [`cowplot`](https://CRAN.R-project.org/package=cowplot) 
- [`gridExtra`](https://CRAN.R-project.org/package=gridExtra) 
- [`pwr`](https://CRAN.R-project.org/package=pwr) 
- [`parallel`](https://www.R-project.org/) 

## The data
All data are stored under the `~/data/` file. However, legacy study information was removed due to proprietary and confidentiality reasons. To make this code run, one must add a folder under `data/Original` named `/legacy_studies/your_study_name`, place CSV files with all study information as CSV files in CDISC SEND format, and map the study name in the code respectively, in `man/01_write_SEND_to_R_table.R`.
To further understand how the files need to look like, please refer to "data/Original/HCD_RAT" where the data format ist stored accordingly.

#### HCD and mapping tables are provided in the script and consist of the following:
- `Original/HCD_RAT` contains all historical control data stored in CDISC SEND format as CSV files.
- `Original/LBTEST_dictionary.csv` carries translations for laboratory parameter short names into long names.
- `LB_mapping.xlsx`, `LB_unit_mapping.xlsx`, and `OM_mapping.xlsx` are harmonizing steps to unify different laboratory parameters, e.g., units like "U/L", "U/l", "u/l" are harmonized as "U/L".

## Scripts
The code sinsists of 23 R-scripts and is controlled over the `00_master.R` script. In the following, a short description of every script and its function is given:

### `00_master.R`
This program is the execution file of the calculation and visualization steps of the remaining R-scripts. The functions of the R-scripts are executed here with respect to the entered values.
### `01_write_SEND_to_R_table.R`
This script selects the quantitative results of the legacy study inputs for males and females respectively. These results will be exported as tables to the data/derived folder.
### `01a1_write_qualitative_SEND_data_to_R_table.R`
This script reads and concatenates CSV data of the qualitative parameters (i.e., microscopic findings (MI), macroscopic findings (MA), and clinical observations (CO)) into one table where the background incidences can be calculated from.
### `02_make_statistical_test.R`
This script creates a function which does the following:
- Calculate a statistical test with respect to the endpoint of interest for each sex separately using the individual values.
- Summarize results in a table and return the table.
### `03_write_HCD_to_R_table.R`
This script extracts the historical control data stored in ~/data/Original/quantitative_parameters and prepares them into a table in the same fashion as in "01_write_legacy_study_original_results.R". The same function is used for that.
### `03a1_missing_data_calculator.R`
This is a standalone script, i.e., it is not called by the master script.
Based on the historical control data (HCD) extracted in "03_write_HCD_to_R_table.R", this script calculates the missing data per measured endpoint and exports a table containing a list of all endpoints and how well they are present in the data.
### `03a2_1_HCD_visualizer_with_coloring.R`
This script is a version of 03a2_HCD_visualizer.R. But here, instead of batch-plotting everything as a box plot and/or as a histogram, you can select the x-axis of the boxes and color it accordingly. For instance, you can color the boxes so that they represent certain vehicles used in a study.
### `03a2_HCD_visualizer.R`
This script creates histograms and box plots (denounced histbox) for each parameter of interest. This script is used for quality control of HCD. As an input, historical control data is fed into the function and then all parameters are visualized below each other creating a long HTML list.
As an additional quality control, studies are highlighted if their quartiles are outside of the "grand median" value (i.e., the study seems to be different towards all other studies).
### `03a3_in_life_data_reader.R`
This script reads the in-life data table results extracted from PDF study reports. It extracts the summary statistics and the significant differences which will be compared to the statistical results of the R-script 02_make_statistical_test.R. If the results do not overlap, the significance tests will be adjusted accordingly.
### `03a4_refdata_reader.R`
This script reads the reference data which is entered as DOCX documents and transforms them in to tables easily readable by R. These reference values are then used to compare whether animal values are within or outside of the HCD 5th/95th percentile reference range.
### `04_resampling.R`
This script is based on "02_make_statistical_test.R". However, instead of the concurrent control group (CCG) of the legacy studies, this set randomly samples control group values (without replacement) of the HCD. These values are then used to replace the CCGs. It can be selected whether all CCGs are removed or only a fraction of them. Furthermore, the number of iterations can be selected. This script needs the location of the CSV files containing SEND formatted data stored in Data/Original/legacy_study and Data/Original/quantitative_parameters.
This function does the following:
- Replace the original control group values with a random sample (without replacement) of values from HCD and recalculate the Dunnett test
- Repeat previous step n times
- Summarize results in a list giving the percentage of how often the original result was reproduced (per dose group) as well as the mean values of each sample of the VCGs.
### `04a1_simulate_missing_parameters.R`
This function aims to simulate parameters if no sufficient data is present in the VCG data set. This is done in the following fashion:
Input:  historical control data (HCD) and if present, sentinel animals.
Steps:  - imute missing data by mean value of non-missing data
        - impute missing data by randomly sampling non-missing data
        - impute missing data by predictive mean matching (pmm)      
Output: Imputed data
### `05_visualize_as_table.R`
This script generates a gt table showing significant results of the legacy study and the concatenated results of the VCG resampling experiment. The goal is to decide whether treatment-relatedness occurs by observing not only statistical significance but also other points which might speak for treatment relatedness (such as dose dependency or whether CG is outisde of 5th/95th percentile-area of HCD).
Addittionally, the tables are exported as XLSX so that they can be edited by subject matter experts and study directors.
### `06_visualize_as_plot.R`
This script generates a ggplot showing significant results of the legacy study and the concatenated results of the VCG resampling experiment. The goal is to see whether the CCG values are strongly outside of the VCG range resulting poor reproducibility.
### `06a2_MI_plot.R`
This script aims to visualize the percentage of microscopic findings per dose group of a legacy study. Once with the original concurrent controls (CCG) and once with all historical control data (HCD). The goal of this plot is to gain an understanding in "how the number of organ findings change" if using virtual control groups (VCGs) instead of CCGs.
Input:  legacy study data and historical control data as CSV
Steps:  - match HCD to the initial body weight (INITBW) of the dose group animals of the legacy study.
        - visualize the original findings (for liver and kidney as the predominant organs in toxicological assessments).
        - do the same visualization but with HCDs instead of CCGs.
        - we use percentages instead of absolute number of findings to account for the imbalance of HCD compared to dose group animals.
        - males and females are regarded separately.
Output: ggplot plot exported as JPEG and PDF.
### `06a3_BW_plot.R`
This is a standalone script which means it is not executed by 00_master.R
This script aims to visualize the body weight development of animals with as a line plot with respect to the study day. Each dose group is shown as a mean value along with their standard deviations as error bars. In case if a dose group has a statistically significant difference to the control on a specific day, this group (on this day) is marked with an asterisk. Same plot is done for body weight gain (BG) instead of BW. But now, instead of body weight, the gain is plotted there.
Two plots are shown side by side: on the left, the plot with the concurrent control groups of the legacy study. On the right, the same legacy study, but with the mean values of all iterations from the VCGs.
### `07_quality_measures_plot.R`
This is a standalone script which means it is not executed by 00_master.R
This script generates a ggplot showing significant results of the legacy study and the concatenated results of the VCG resampling experiment. The goal is to see whether the CCG values are strongly outside of the VCG range resulting poor reproducibility.
### `08_compare_test_substance_relatedness.R`
This is a standalone script which means it is not executed by 00_master.R
This script collects the XLSX documents created in 05_visualize_as_table.R and summarizes the results of them. The XLSX files serve to support a study director in decision making process. Basically, the study directors can use the XLSX file to classify whether a significant finding is treatment related or not and - if not - provide reasons for their decisions. In this script, we summarize the decisions and the discrepancies between CCG results and VCG results.
### `09_sensitivity_analysis.R`
This is a standalone script which means it is not executed by 00_master.R
This script is executed in 02_make_statistical_test.R. It assesses the power of the statistical results in order to check how far away the observed effect needs to be in order to obtain 80 % of power. The point of interest here is to check whether the results of the CCG had enough power from the beginning and whether we'd be able to detect the result with VCGs head on.
### `10_BW_selction_visualizaton.R`
This is a standalone script which means it is not executed by 00_master.R
This script creates a histogram showing the distribution of the initial body weights of the HCD and the INITBW values of the dose groups.
### `11_MAMICO_incidences.R`
This is a standalone script which means it is not executed by 00_master.R
This script takes the results of the legacy studies and compares the percentage of incidences of qualitative findings (MAcroscopic findings, MIcroscopic findings, Clinical Observations) of HCD. For instance, if digging and cleaning movements are observed in animals in the legacy study, you can compare how many animals in the historical control data (HCD) showed the same behavior.
### `12_effect_size_interpretation.R`
This is a standalone script which means it is not executed by 00_master.R
This script interprets the results of Cliff's Delta effect sizes. The following questions I wanted to answer:
1) Is the CCG between 5 % and 95 % percentile of VCG range? Y/N,
2) If so, does the majority vote for the same thing? Y/N
3) How big is the smallest effect leading to a significant difference?
4) How big is the biggest effect leading to a non significant difference?
#*For this script, the effect size results and the voting results of the VCG script are used.
### `13_power_interpretations`
This is a standalone script which means it is not executed by 00_master.R
This script extracts the power calculations of the data and checks which of the parameter values were actually well powered.

## Figures
The resulting figures are stored in the `~/Results`-subfolder.

## Authors
Alexander Gurjanov - [Alexander Gurjanov](mailto:alexander.gurjanov@bayer.com)
Carlos Vieira e Vieira - [Carlos Vieira e Vieira](mailto:carlos.vieiraevieira@bayer.com)
