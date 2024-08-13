################################### HEADER #####################################
### SCRIPT: Protection/VE                                                    ###
### AUTHOR: Savannah Hammerton                                               ###
### PURPOSE: This script will generate the protection values and estimate    ###
###          VE from those values.                                           ###
################################# END HEADER ###################################

# First steps-------------------------------------------------------------------
## Load packages----------------------------------------------------------------
library(tidyverse)
library(here)
library(boot)

## Load data--------------------------------------------------------------------
# Processed homologous dataset 
homologous_dat <- readRDS(here::here("data/raw-data/homologous_dat.rds"))

# Pulled CDC reported seasonal VE 
seasonal_ve <- readRDS(here::here("data/raw-data/seasonal_ve.rds"))
## Script setup-----------------------------------------------------------------
here::here() #set paths
source(here::here("code/objects.R")) #load objects 
source(here::here("code/functions.R")) #load functions

# Analysis----------------------------------------------------------------------

## Main analysis----------------------------------------------------------------
# Need to nest the data so that it doesn't use group; that way I can get the 
# group contrasts
# First, taking out anyone under 18 years old, as we're not using kids in main
# analysis
seasons <-
  homologous_dat |> 
  # only keeping those 18+ years old 
  dplyr::filter(age >= 18) |> 
  protect_df() |> 
  tidyr::nest(data = -c(season, strain_type))

# Set seed 
set.seed(42)

# Calculate VEs and IVEs using `get_boots()` function from `functions.R` script
ive_boots_list <- 
  purrr::map(seasons$data, ~get_boots(.x, nstat=9), .progress = TRUE)

# Wrangle VE estimates and CIs and into usable final data set
boot_ve <- 
  dplyr::select(seasons, !data) |> 
  dplyr::mutate(boot = ive_boots_list) |> 
  tidyr::unnest(boot) |> 
  dplyr::filter(!grepl("ive", statistic)) |> 
  tidyr::separate(statistic, into = c("statistic", "group"),
                  sep = "_") |> 
  dplyr::mutate(group = toupper(group),
                data = "Main analysis") |> 
  dplyr::rename("boot_ve" = boot_est,
                "boot_ve_lwr" = boot_est_lwr,
                "boot_ve_upr" = boot_est_upr)
# change variable names to be more specific to what we've calculated 

saveRDS(boot_ve, file = here::here("results/output/boot_ve.rds"))

# Wrangle IVE estimates and CIs and into usable final data set
ive_boot <- 
  dplyr::select(seasons, !data) |> 
  dplyr::mutate(boot = ive_boots_list) |> 
  tidyr::unnest(boot) |> 
  dplyr::filter(grepl("ive", statistic)) |> 
  tidyr::separate(statistic, into = c("type", "comparison1", "comparison2"),
                  sep = "_") |> 
  tidyr::unite(col = "comparison", comparison1, comparison2) |> 
  dplyr::mutate(type = ifelse(type == "rive", "Relative", "Absolute"),
                comparison = toupper(comparison),
                data = "Main analysis") 

# Save bootstrapped IVEs
saveRDS(ive_boot, file = here::here("results/output/boot_ive.rds"))

## Individual-based VE----------------------------------------------------------
# This VE uses the same population as the main analysis, but is based on the 
# comparison of each person's pre- and post- protection, rather than comparing
# the pre- and post-vaccine "group" protections


# Set seed 
set.seed(42)

# Calculate VEs and IVEs using `get_boots()` function from `functions.R` script
ive_boots_list_ind <- 
  purrr::map(seasons$data, ~get_boots(.x, nstat=9, statistic = est_IVE_ind), 
             .progress = TRUE)

# Wrangle VE estimates and CIs and into usable final data set
boot_ve_ind <- 
  dplyr::select(seasons, !data) |> 
  dplyr::mutate(boot = ive_boots_list_ind) |> 
  tidyr::unnest(boot) |> 
  dplyr::filter(!grepl("ive", statistic)) |> 
  tidyr::separate(statistic, into = c("statistic", "group"),
                  sep = "_") |> 
  dplyr::mutate(group = toupper(group),
                data = "Individual VE") |> 
  dplyr::rename("boot_ve" = boot_est,
                "boot_ve_lwr" = boot_est_lwr,
                "boot_ve_upr" = boot_est_upr)
# change variable names to be more specific to what we've calculated 

saveRDS(boot_ve_ind, file = here::here("results/output/boot_ve_ind.rds"))

# Wrangle IVE estimates and CIs and into usable final data set
ive_boot_ind <- 
  dplyr::select(seasons, !data) |> 
  dplyr::mutate(boot = ive_boots_list_ind) |> 
  tidyr::unnest(boot) |> 
  dplyr::filter(grepl("ive", statistic)) |> 
  tidyr::separate(statistic, into = c("type", "comparison1", "comparison2"),
                  sep = "_") |> 
  tidyr::unite(col = "comparison", comparison1, comparison2) |> 
  dplyr::mutate(type = ifelse(type == "rive", "Relative", "Absolute"),
                comparison = toupper(comparison),
                data = "Individual VE") 

# Save bootstrapped IVEs
saveRDS(ive_boot_ind, file = here::here("results/output/boot_ive_ind.rds")) 

## LoD VE-----------------------------------------------------------------------
# This is the sensitivity analysis excluding anyone with a titer at the limit
# of detection (5) at either the pre-vaccine or post-vaccine titer 
# Filter out LoD data, get protection estimates, and nest for VE/IVE
# calculations
nolod_dat <- 
  homologous_dat |> 
  dplyr::filter(pretiter != 5 & postiter != 5,
                age >= 18, 
                # there are not enough data in one of the groups in 2020, 
                #so we'll skip that 
                season != 2020) |> 
  protect_df() |> 
  tidyr::nest(data = -c(season, strain_type)) 


# Set seed 
set.seed(42)

# Calculate VEs and IVEs using `get_boots()` function from `functions.R` script
ive_boots_list_nolod <- 
  purrr::map(nolod_dat$data, ~get_boots(.x, nstat=9), .progress = TRUE)

# Wrangle VE estimates and CIs and into usable final data set
boot_ve_nolod <- 
  dplyr::select(nolod_dat, !data) |> 
  dplyr::mutate(boot = ive_boots_list_nolod) |> 
  tidyr::unnest(boot) |> 
  dplyr::filter(!grepl("ive", statistic)) |> 
  tidyr::separate(statistic, into = c("statistic", "group"),
                  sep = "_") |> 
  dplyr::mutate(group = toupper(group),
                data = "No LoD") |> 
  dplyr::rename("boot_ve" = boot_est,
                "boot_ve_lwr" = boot_est_lwr,
                "boot_ve_upr" = boot_est_upr)
# change variable names to be more specific to what we've calculated 

saveRDS(boot_ve_nolod, file = here::here("results/output/boot_ve_nolod.rds"))

# Wrangle IVE estimates and CIs and into usable final data set
ive_boot_nolod <- 
  dplyr::select(nolod_dat, !data) |> 
  dplyr::mutate(boot = ive_boots_list_nolod) |> 
  tidyr::unnest(boot) |> 
  dplyr::filter(grepl("ive", statistic)) |> 
  tidyr::separate(statistic, into = c("type", "comparison1", "comparison2"),
                  sep = "_") |> 
  tidyr::unite(col = "comparison", comparison1, comparison2) |> 
  dplyr::mutate(type = ifelse(type == "rive", "Relative", "Absolute"),
                comparison = toupper(comparison),
                data = "No LoD") 

# Save bootstrapped IVEs
saveRDS(ive_boot_nolod, file = here::here("results/output/boot_ive_nolod.rds")) 


## Including minors-------------------------------------------------------------

# Get protection estimates, and nest for VE/IVE calculations; not excluding 
# minors 
kid_dat <- 
  homologous_dat |> 
  protect_df() |> 
  tidyr::nest(data = -c(season, strain_type)) 

# Set seed 
set.seed(42)

# Calculate VEs and IVEs using `get_boots()` function from `functions.R` script
ive_boots_list_kid <- 
  purrr::map(kid_dat$data, ~get_boots(.x, nstat=9), .progress = TRUE)

# Wrangle VE estimates and CIs and into usable final data set
boot_ve_kid <- 
  dplyr::select(kid_dat, !data) |> 
  dplyr::mutate(boot = ive_boots_list_kid) |> 
  tidyr::unnest(boot) |> 
  dplyr::filter(!grepl("ive", statistic)) |> 
  tidyr::separate(statistic, into = c("statistic", "group"),
                  sep = "_") |> 
  dplyr::mutate(group = toupper(group),
                data = "All ages") |> 
  dplyr::rename("boot_ve" = boot_est,
                "boot_ve_lwr" = boot_est_lwr,
                "boot_ve_upr" = boot_est_upr)
# change variable names to be more specific to what we've calculated 

saveRDS(boot_ve_kid, file = here::here("results/output/boot_ve_kid.rds"))

# Wrangle IVE estimates and CIs and into usable final data set
ive_boot_kid <- 
  dplyr::select(kid_dat, !data) |> 
  dplyr::mutate(boot = ive_boots_list_kid) |> 
  tidyr::unnest(boot) |> 
  dplyr::filter(grepl("ive", statistic)) |> 
  tidyr::separate(statistic, into = c("type", "comparison1", "comparison2"),
                  sep = "_") |> 
  tidyr::unite(col = "comparison", comparison1, comparison2) |> 
  dplyr::mutate(type = ifelse(type == "rive", "Relative", "Absolute"),
                comparison = toupper(comparison),
                data = "All ages") 

# Save bootstrapped IVEs
saveRDS(ive_boot_kid, file = here::here("results/output/boot_ive_kid.rds")) 

## VE for comparing to CDC estimates--------------------------------------------
# Setup data frame with relevant grouping for age/dose groups, seasons, and 
# strains
cdc_groups <- 
  protect_df(homologous_dat) |> 
  dplyr::filter(season < 2020) |> 
  dplyr::mutate(age_group = dplyr::case_when(
    age <= 17 ~ age_groups[[1]],
    age < 50 & age > 17 ~ age_groups[[2]],
    age > 50 & age < 65 ~ age_groups[[3]],
    age >= 65 ~ age_groups[[4]]
  )) |> 
  # not including anyone under 18 because many seasons did not have include 
  # those under 18
  dplyr::filter(age_group != age_groups[[1]]) |> 
  tidyr::nest(data = -c(age_group, season, strain_type))

# Set seed 
set.seed(42)

# Bootstrap VE estimates and CIs and wrangle into usable final data set
boot_ve_cdc <- 
  purrr::map(cdc_groups$data, bootstrap, .progress = TRUE) |> 
  purrr::map(bootstrap.ci) |> 
  dplyr::bind_rows() |> 
  cbind(cdc_groups) |> 
  # cbind back to groups dataset to get seasons/strain/group info matched to
  # ve estimate
  dplyr::select(!data) |> 
  # remove data list variable from the dataset
  dplyr::rename("boot_ve" = boot_est,
                "boot_ve_lwr" = boot_est_lwr,
                "boot_ve_upr" = boot_est_upr)
# change variable names to be more specific to what we've calculated 

# Cleaning up data to match the pulled CDC data 
ve_compare <-
  boot_ve_cdc |> 
  dplyr::mutate(
    dominating_strain = dplyr::case_when(
      season == 2013 ~ "H1N1",
      season == 2014 ~ "H3N2",
      season == 2015 ~ "H1N1",
      season == 2016 ~ "H3N2",
      season == 2017 ~ "H3N2",
      season == 2018 ~ "H1N1 & H3N2",
      season == 2019 ~ "H1N1"
    ),
    season = dplyr::case_when(
      season == 2013 ~ "2013-2014", 
      season == 2014 ~ "2014-2015", 
      season == 2015 ~ "2015-2016", 
      season == 2016 ~ "2016-2017", 
      season == 2017 ~ "2017-2018", 
      season == 2018 ~ "2018-2019", 
      season == 2019 ~ "2019-2020"
    ),
    Type = 
      factor(ifelse(stringr::str_detect(dominating_strain, 
                                        as.character(strain_type)), 
                    1, 0),
             levels = c(0,1),
             labels = c("Not Dominant", "Dominant"))
  ) |> 
  dplyr::left_join(seasonal_ve, 
                   by = c("season", "age_group", "strain_type"))

# Save VE comparisons dataset 
saveRDS(ve_compare, here::here("results/output/ve_compare.rds"))

