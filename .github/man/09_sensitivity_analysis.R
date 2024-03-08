#*09_sensitivity_analysis.R.
#*
#*2023 07 27
#*
#*This script assesses the power of the statistical results in order to check
#*how far away the observed effect needs to be in order to obtain 80 % of power.
#*The point of interest here is to check whether the results of the CCG had
#*enough power from the beginning and whether we'd be able to detect the
#*result with VCGs head on.
#*
#*Input:  Statistical results of 04_resampling.R
#*Steps:  - Extract significance result for each one.
#*        - Compute Cohen's D for each difference between CG and dose group.
#*        - Calculate post-hoc power.
#*        - Plot the results where the effect size is on the x-axis, the
#*          power is the color and the statistical results are different shapes.
#*Output: Table with the power and graph showing how far away the results need
#*        to be in order to obtain enough power.
#*******************************************************************************
#*******************************************************************************
#*Divide work to many clusters----
#*******************************************************************************

#*******************************************************************************
#*******************************************************************************
#beginning of function----
#*******************************************************************************
sensitivity_analysis <- function()
{
  #*****************************************************************************
  #*****************************************************************************
  #*Divide work to many clusters----
  #*****************************************************************************
  n_cores <- detectCores()
  #create cluster
  cluster <- makeCluster(n_cores - 1)
  #*****************************************************************************
  #*****************************************************************************
  #*Load libraries----
  #*****************************************************************************
  #export libraries to cluster
  clusterEvalQ(cluster,{library(tidyverse)})

  #*Load data into cluster
  clusterExport(cluster, "studydata_list")
  #*****************************************************************************
  #*****************************************************************************
  #*Compute post-hoc power analysis----
  #*****************************************************************************
  powers <- parLapply(cluster, studydata_list, function(parameter){
    #get the statistical test which was used in the data
    if(
      is.na(unique(parameter$significance_test))
    ){
      #no test----
      # print("no test")
    }

    else if(
      #*A prerequisite of the Dunnett's test is that at least 3 groups are
      #*present. If this is not the case (i.e., if only 2 groups are present)
      #*a t-test is calculated instead.
      unique(parameter$significance_test) == "dunnett" &
      n_distinct(parameter$trial_set_description) == 2
    ){
      # print("t-test")
      pwr.t.test(
        #*get the population of the smallest group size (if the sizes are un-
        #*balanced, the more conservative approach is chosen)
        n = parameter %>%
          group_by(trial_set_description) %>%
          summarize(population = n()) %>%
          pull(population) %>%
          min(),
        #get Cohen's D as an effect size
        d = cohen.d(
          parameter %>%
            filter(trial_set_description == "CG") %>%
            pull(LBORRES),
          parameter %>%
            filter(trial_set_description == "HD") %>%
            pull(LBORRES)
        )$estimate
      )$power

    }else if(
      unique(parameter$significance_test) == "dunnett" &
      n_distinct(parameter$trial_set_description) > 2
    ){
      #*Dunnett's exact homogenous test----
      # print("dunnett")
      powerMCTests(
        mu = c(
          parameter %>%
            group_by(trial_set_description) %>%
            summarise(mean_vals = mean(LBORRES)) %>%
            pull(mean_vals)
        ),
        n = c(
          parameter %>%
            group_by(trial_set_description) %>%
            summarise(population = n()) %>%
            pull(population)
        ),
        #*adjust the standard deviation. A normal distribution for errors is
        #*assumed. The max error range is selected for a conservative approach
        parms = list(
          mean = 0,
          sd = parameter %>%
            group_by(trial_set_description) %>%
            summarise(sd_vals = sd(LBORRES)) %>%
            pull(sd_vals) %>%
            max()
        ),
        test = "dunnettTest",
        replicates = 500
      )$anypair

    }else if(
      #*A prerequisite of the Dunnett's test is that at least 3 groups are
      #*present. If this is not the case (i.e., if only 2 groups are present)
      #*a t-test is calculated instead.
      unique(parameter$significance_test) == "het_dunn" &
      n_distinct(parameter$trial_set_description) == 2
    ){
      #*Note: t-test doesn't support unequal variances. The power may be
      #*overestimated here!
      # print("t-test")
      pwr.t.test(
        #*get the population of the smallest group size (if the sizes are un-
        #*balanced, the more conservative approach is chosen)
        n = parameter %>%
          group_by(trial_set_description) %>%
          summarize(population = n()) %>%
          pull(population) %>%
          min(),
        #get Cohen's D as an effect size
        d = cohen.d(
          parameter %>%
            filter(trial_set_description == "CG") %>%
            pull(LBORRES),
          parameter %>%
            filter(trial_set_description == "HD") %>%
            pull(LBORRES)
        )$estimate
      )$power

    }else if(
      unique(parameter$significance_test) == "het_dunn" &
      n_distinct(parameter$trial_set_description) > 2
    ){
      #Dunnett's exact heterogenous test----
      # print("hetdunn")
      powerMCTests(
        mu = c(
          parameter %>%
            group_by(trial_set_description) %>%
            summarise(mean_vals = mean(LBORRES)) %>%
            pull(mean_vals)
        ),
        n = c(
          parameter %>%
            group_by(trial_set_description) %>%
            summarise(population = n()) %>%
            pull(population)
        ),
        #*adjust the standard deviation. A normal distribution for errors is
        #*assumed. The max error range is selected for a conservative approach
        parms = list(
          mean = 0,
          sd = parameter %>%
            group_by(trial_set_description) %>%
            summarise(sd_vals = sd(LBORRES)) %>%
            pull(sd_vals) %>%
            max()
        ),
        test = "dunnettT3Test",
        replicates = 500
      )$anypair

    }else if(
      #*A prerequisite of the Dunnett's test is that at least 3 groups are
      #*present. If this is not the case (i.e., if only 2 groups are present)
      #*a t-test is calculated instead.
      unique(parameter$significance_test) == "log_trans_dunnett" &
      n_distinct(parameter$trial_set_description) == 2
    ){
      # print("log trans t-test")
      pwr.t.test(
        #*get the population of the smallest group size (if the sizes are un-
        #*balanced, the more conservative approach is chosen)
        n = parameter  %>%
          group_by(trial_set_description) %>%
          summarize(population = n()) %>%
          pull(population) %>%
          min(),
        #get Cohen's D as an effect size
        d = cohen.d(
          parameter  %>%
            filter(trial_set_description == "CG") %>%
            pull(LBORRES) %>%
            log10(),
          parameter  %>%
            filter(trial_set_description == "HD") %>%
            pull(LBORRES) %>%
            log10()
        )$estimate
      )$power

    }else if(
      unique(parameter$significance_test) == "log_trans_dunnett" &
      n_distinct(parameter$trial_set_description) > 2
    ){
      #Dunnett's exact homogeneous test after logarithmic transformation----
      # print("logtrandunn")
      powerMCTests(
        mu = c(
          parameter %>%
            group_by(trial_set_description) %>%
            summarise(mean_vals = mean(LBORRES)) %>%
            pull(mean_vals) %>%
            log10()
        ),
        n = c(
          parameter %>%
            group_by(trial_set_description) %>%
            summarise(population = n()) %>%
            pull(population)
        ),
        #*adjust the standard deviation. A normal distribution for errors is
        #*assumed. The max error range is selected for a conservative approach
        parms = list(
          mean = 0,
          sd = parameter %>%
            group_by(trial_set_description) %>%
            summarise(sd_vals = sd(LBORRES)) %>%
            pull(sd_vals) %>%
            max()
        ),
        test = "dunnettTest",
        replicates = 500
      )$anypair

    }else if(unique(parameter$significance_test) == "u_test"){
      #*An important prerequesite of the Bonferroni-adjusted Mann Whitney U
      #*test is that at least 6 data points per group are present.
      if(
        any(
          parameter %>%
          group_by(trial_set_description) %>%
          summarize(population = n()) %>%
          pull(population) < 6
        )
      ){
        # print("no u-test")
      }else if(n_distinct(parameter$trial_set_description) == 2){
        #U-test
        # print("u-test")
        NaN
      }else{
        #Bonferroni adjusted Mann Whitney U-test---
        # print("u-test")
        powerMCTests(
          mu = c(
            parameter  %>%
              group_by(trial_set_description) %>%
              summarise(mean_vals = mean(LBORRES)) %>%
              pull(mean_vals)
          ),
          n = c(
            parameter  %>%
              group_by(trial_set_description) %>%
              summarise(population = n()) %>%
              pull(population)
          ),
          #*adjust the standard deviation. A normal distribution for errors is
          #*assumed. The max error range is selected for a conservative approach
          parms = list(
            mean = 0,
            sd = parameter %>%
              group_by(trial_set_description) %>%
              summarise(sd_vals = sd(LBORRES)) %>%
              pull(sd_vals) %>%
              max()
          ),
          test = "pairwise.wilcox.test",
          p.adjust.method = "bonferroni",
          replicates = 500
        )$anypair
      }
    }
  })
  #stop the cluster
  stopCluster(cluster)
  return(powers)
}

