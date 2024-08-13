################################### HEADER #####################################
### SCRIPT: Protection/VE                                                    ###
### AUTHOR: Savannah Hammerton                                               ###
### PURPOSE: This script will generate the protection values and conduct     ###
###          any data cleaning/prep needed for other scripts                 ###
################################# END HEADER ###################################


## Load data--------------------------------------------------------------------
# Processed homologous dataset 
homologous_dat <- readRDS(here::here("data/raw-data/homologous_dat.rds"))

## Script setup-----------------------------------------------------------------
here::here() #set paths
ggplot2::theme_set(ggplot2::theme_minimal()) #set plot theme
source(here::here("code/objects.R")) #load objects 
source(here::here("code/functions.R")) #load functions


## Protection-------------------------------------------------------------------
# Calculate predicted protections and 95% CIs
protection <-
  homologous_dat |> 
  dplyr::mutate(pre_protect = predict_protect(pretiter, 3.09, 1.42),
                post_protect = predict_protect(postiter, 3.09, 1.42),
                inverse_pre = 1 - pre_protect,
                inverse_post = 1 - post_protect)

saveRDS(protection, here::here("results/output/protection.rds"))
