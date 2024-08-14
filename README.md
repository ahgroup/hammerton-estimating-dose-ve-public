# Overview


This is the repository for the analysis of ""Estimating standard-dose and high-dose Fluzone vaccine efficacies for influenza A based on HAI titers". All analysis reproduction should be possible from this repo.


# Repository Description 

  * All initial, "raw" data can be found in the `data/raw-data` folder 
  * All code used to conduct this analysis and generate results can be found in the `code` folder 
  * All generated results can be found in the `results` folder. The `figures` and `tables` folders contain all figures and tables presented in the main text and supplement, while the `output` folder contains .rds files with data sets generated through the analysis 
  * The files used to generate the main text and supplement of the manuscript can be found in the `products/manuscript` and `products/manuscript/supplement` folders, respectively
  * More information as needed is provided in the README files within each folder


# Analysis Replication 

This project uses `renv` to help ensure reproducibility. Upon downloading the repository, running `renv::restore()` will download all used packages and dependencies in the versions used in conducting this analysis. 

To completely replicate our analysis, run the code scripts in the following order: 

  1) processing-code.R
  2) statistical-analysis.R
  3) fig-gen.R
  4) table-gen.R
  
The functions.R and objects.R scripts will each be called at the beginning of every script listed above, so they do not need to run explicitly. However, the functions.R script in particular contains the code responsible for much of the analysis, such as the VE and DVE calculations and bootstrapped CIs, so alterations to that part of the analysis would take place in that script. 

The CDC reported estimates were initially directly webscraped from the CDC website. However, several of the seasons' pages have since been archived, altering the URLs used to access them. We have therefore placed the data set with CDC reported VE in the `raw-data` folder. All values can still be checked by accessing each season's reported VE and accompanying information at <https://www.cdc.gov/flu/vaccines-work/past-seasons-estimates.html>. 

