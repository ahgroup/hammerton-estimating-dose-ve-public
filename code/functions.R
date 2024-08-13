################################### HEADER #####################################
### SCRIPT: Functions                                                        ###
### AUTHOR: Savannah Hammerton                                               ###
### PURPOSE: This script will hold all the functions used across various     ###
###          scripts in this project. Running `source(here::here(            ###
###          "code/functions.R"))` will place all the functions in the       ###
###          script into the global environment and get them ready for use.  ###
################################# END HEADER ###################################

library(boot)

# Calculations------------------------------------------------------------------
# Function to predict VE 
# This function will calculated an estimated protection probability given an
# input titer value and alpha and beta. Alpha and beta, while pre-determined,
# will vary based on whether we are calculating protection from the main curve,
# lower credible bound, or upper credible bound
predict_protect <- function(titer, alpha, beta) {
  1 - (1/(1 + exp(beta*(log(titer) - alpha))))
}


# Function to add protection probabilities to dataframe
# This function will calculate predicted protections and 95% CIs for an input 
# dataframe. It will also calculate risk values, which are just the inverse of 
# the protection values (one risk value for each protection value)
protect_df <- function(dat) {
output <-
  dat |> 
  dplyr::mutate(pre_protect = predict_protect(pretiter, 3.09, 1.42),
                post_protect = predict_protect(postiter, 3.09, 1.42),
                risk_pre = 1 - pre_protect,
                risk_post = 1 - post_protect)
}




# Function to calculate bootstrap distribution, estimate, confidence interval
# Default arguments for VE statistics, 9999 resamples
# Can change these defaults as needed for different calculations 
bootstrap <- function(data, statistic = est_VE, R = 9999, sim = "ordinary") {
  # Calculate bootstrap distribution, needed for getting CI in next step
  dist <- boot::boot(data = data, statistic = statistic, R = R, sim = sim)
  
  return(dist)
}

# Function to get actual estimate and CI using `dist` object produced from 
# boostrap function 
bootstrap.ci <- function(dist, ...) {
  # Calculate boostrapped estimate and 95% ci, BCa method
  ci <- boot::boot.ci(dist, conf = 0.95, type = "bca", ...)
  
  # Pull out relevant info from ci object and name variables - keeping names 
  # vague here so I can change them to names that match the statistic as I
  # implement the function
  out <- c("boot_est" = ci$t0, 
           "boot_est_lwr" = ci$bca[,4], 
           "boot_est_upr" = ci$bca[,5])
  
  return(out)
}



# Function to calculate relative and absolute increase in VE across groups
# The VE estimate used here is more similar to an efficacy trial estimate/other 
# published estimates. It essentially takes a ratio of the mean risk in the
# post/pre vaccine groups
est_IVE <- function(dat, idx) {
  # setting up the second argument for the boot function (indices for 
  # resample) and defining new data object
  di <- dat[idx, ]
  # get ve for oahd group
  oahd <- dplyr::filter(di, group == "OAHD")
  ve_oahd <- 
    (1 - sum(oahd$risk_post)/sum(oahd$risk_pre))
  # get ve for oasd group
  oasd <- dplyr::filter(di, group == "OASD")
  ve_oasd <- 
    (1 - sum(oasd$risk_post)/sum(oasd$risk_pre))
  # get ve for yasd group
  yasd <- dplyr::filter(di, group == "YASD")
  ve_yasd <-
    (1 - sum(yasd$risk_post)/sum(yasd$risk_pre))
  
  # Get relative and absolute increase in VE from previously calculated VEs 
  rive_oahd_oasd <- (ve_oahd - ve_oasd) / ve_oasd
  aive_oahd_oasd <- ve_oahd - ve_oasd
  rive_oahd_yasd <- (ve_oahd - ve_yasd) / ve_yasd
  aive_oahd_yasd <- ve_oahd - ve_yasd
  rive_oasd_yasd <- (ve_oasd - ve_yasd) / ve_yasd
  aive_oasd_yasd <- ve_oasd - ve_yasd
  
  return(c(ve_oahd, ve_oasd, ve_yasd,
           rive_oahd_oasd, aive_oahd_oasd,
           rive_oahd_yasd, aive_oahd_yasd,
           rive_oasd_yasd, aive_oasd_yasd))
}

# Function to calculate relative and absolute increase in VE across groups
# The VE estimate used here compares each individual's post-vaccine risk to  
# their pre-vaccine risk. While this isn't as easily comparable to traditional  
# VE estimates, it is a benefit of these methods that we can compare each 
# individual's "outcome" to their "counterfactual" since we have their pre- 
# and post- risks
est_IVE_ind <- function(dat, idx) {
  # setting up the second argument for the boot function (indices for 
  # resample) and defining new data object
  di <- dat[idx, ]
  # get ve for oahd group
  oahd <- dplyr::filter(di, group == "OAHD")
  ve_oahd <- 
    mean(1 - (oahd$risk_post/oahd$risk_pre))
  # get ve for oasd group
  oasd <- dplyr::filter(di, group == "OASD")
  ve_oasd <- 
    mean(1 - (oasd$risk_post/oasd$risk_pre))
  # get ve for yasd group
  yasd <- dplyr::filter(di, group == "YASD")
  ve_yasd <-
    mean(1 - (yasd$risk_post/yasd$risk_pre))
  
  # Get relative and absolute increase in VE from previously calculated VEs 
  rive_oahd_oasd <- (ve_oahd - ve_oasd) / ve_oasd
  aive_oahd_oasd <- ve_oahd - ve_oasd
  rive_oahd_yasd <- (ve_oahd - ve_yasd) / ve_yasd
  aive_oahd_yasd <- ve_oahd - ve_yasd
  rive_oasd_yasd <- (ve_oasd - ve_yasd) / ve_yasd
  aive_oasd_yasd <- ve_oasd - ve_yasd
  
  return(c(ve_oahd, ve_oasd, ve_yasd,
           rive_oahd_oasd, aive_oahd_oasd,
           rive_oahd_yasd, aive_oahd_yasd,
           rive_oasd_yasd, aive_oasd_yasd))
}


# Final function that incorporates all previous functions and returns
# estimates and CIs for all VEs and DVEs
get_boots <- function(dat, statistic = est_IVE, nstat, IVE = TRUE){
  
  test.dist <-  bootstrap(dat, statistic = statistic) 
  
  test.ci <- purrr::map_dfr(1:nstat, ~bootstrap.ci(test.dist, index=.x)) 
  
  if(isTRUE(IVE)) {
    test.ci$statistic = c("ve_oahd", "ve_oasd", "ve_yasd",
                          "rive_oahd_oasd", "aive_oahd_oasd",
                          "rive_oahd_yasd", "aive_oahd_yasd",
                          "rive_oasd_yasd", "aive_oasd_yasd")
  }
  
  return(test.ci)
}


# Function to calculate VE - group specific ratios 
# This specific function will only be used when we are not getting contrasts;
# namely, for the age group VE estimates to compare to CDC reported VE 
est_VE <- function(dat, idx) {
  # setting up the second argument for the boot function (indices for 
  # resample) and defining new data object
  di <- dat[idx, ]
  # define ve with new data object
  ve <- (1 - sum(di$risk_post)/sum(di$risk_pre))
}



