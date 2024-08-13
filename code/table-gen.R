################################### HEADER #####################################
### SCRIPT: Table Generation                                                 ###
### AUTHOR: Savannah Hammerton                                               ###
### PURPOSE: This script will generate all the tables for the manuscript.    ###
################################# END HEADER ###################################

# First steps-------------------------------------------------------------------
## Load packages----------------------------------------------------------------
library(tidyverse)
library(here)
library(patchwork)
library(ggside)
library(ggbeeswarm)
library(flextable)
library(paletteer)

## Load data--------------------------------------------------------------------
# Processed homologous dataset 
homologous_dat <- readRDS(here::here("data/raw-data/homologous_dat.rds"))

# Calculated VEs, linking back to seasonal info from initial data set 
ve <- readRDS(here::here("results/output/boot_ve.rds")) |> 
  dplyr::left_join(
    dplyr::select(homologous_dat, c(season, strain_type, strains_short)),
    by = c("season", "strain_type"),
    relationship = "many-to-many"
  ) |> 
  dplyr::distinct()

# Calculated DVEs, linking back to seasonal info from initial data set 
ive <- readRDS(here::here("results/output/boot_ive.rds")) |> 
  dplyr::left_join(
    dplyr::select(homologous_dat, c(season, strain_type, strains_short)),
    by = c("season", "strain_type"),
    relationship = "many-to-many"
  ) |> 
  dplyr::distinct()

# Reported CDC seasonal VE 
ve_compare <- readRDS(here::here("results/output/ve_compare.rds"))

# Pulling out basic seasonal info
strains <- 
  ve |> 
  dplyr::ungroup() |> 
  dplyr::select(season, strain_type, strains_short) |> 
  dplyr::distinct()


## Script setup-----------------------------------------------------------------
here::here() #set paths
ggplot2::theme_set(ggplot2::theme_minimal()) #set plot theme
source(here::here("code/objects.R")) #load objects 
source(here::here("code/functions.R")) #load functions

# Tables------------------------------------------------------------------------
## Table 1----------------------------------------------------------------------
desc_table <-
  homologous_dat |> 
  dplyr::filter(age >= 18) |> 
  dplyr::group_by(group, season, strain_type) |> 
  dplyr::mutate(n = n()) |> 
  dplyr::select(n, season, strain_type, strains_fullname, group, age) |> 
  dplyr::summarise(strain = strains_fullname[1], 
                   age_med = median(age),
                   age_1st_quart = quantile(age, 0.25),
                   age_3rd_quart = quantile(age, 0.75),
                   n = n[1]) |> 
  tidyr::pivot_wider(names_from = strain_type, values_from = strain) |> 
  dplyr::mutate(Strains = paste0(H1N1, "\n", H3N2),
                season = as.character(season),
                Age = paste0(round(age_med, digits = 0), " (",
                             round(age_1st_quart, digits = 0), ", ",
                             round(age_3rd_quart, digits = 0), ")")) |> 
  dplyr::select(!c(H1N1, H3N2, age_med, age_1st_quart, age_3rd_quart)) |> 
  dplyr::arrange(season) |>  
  dplyr::relocate(Strains, .after = season) |> 
  dplyr::relocate(group, .after = season) |> 
  dplyr::rename("Season" = season,
                "Group" = group) |> 
  flextable::flextable() |> 
  flextable::merge_at(i = 1:3, j = 1, part = "body") |> 
  flextable::merge_at(i = 4:6, j = 1, part = "body") |> 
  flextable::merge_at(i = 7:9, j = 1, part = "body") |>
  flextable::merge_at(i = 10:12, j = 1, part = "body") |> 
  flextable::merge_at(i = 13:15, j = 1, part = "body") |> 
  flextable::merge_at(i = 16:18, j = 1, part = "body") |> 
  flextable::merge_at(i = 19:21, j = 1, part = "body") |> 
  flextable::merge_at(i = 22:24, j = 1, part = "body") |> 
  flextable::merge_at(i = 25:27, j = 1, part = "body") |>
  flextable::hline(i = 3, part = "body") |> 
  flextable::hline(i = 6, part = "body") |> 
  flextable::hline(i = 9, part = "body") |> 
  flextable::hline(i = 12, part = "body") |> 
  flextable::hline(i = 15, part = "body") |> 
  flextable::hline(i = 18, part = "body") |> 
  flextable::hline(i = 21, part = "body") |> 
  flextable::hline(i = 24, part = "body") |> 
  flextable::hline(i = 27, part = "body") |> 
  flextable::line_spacing(space = 1.1, part = "all") |> 
  flextable::footnote(i = 1, j = 5, part = "header", ref_symbols = "1",
                      value = flextable::as_paragraph("Median (IQR)"))


desc_table

saveRDS(desc_table, here::here("results/tables/tb1.rds"))

# Second table one, wide format for slides/presentations
desc_table_wide <-
  homologous_dat |> 
  dplyr::filter(age >= 18) |> 
  dplyr::group_by(group, season, strain_type) |> 
  dplyr::mutate(n = n()) |> 
  dplyr::select(n, season, strain_type, strains_fullname, group, age) |> 
  dplyr::summarise(strain = strains_fullname[1], 
                   age_med = median(age),
                   age_1st_quart = quantile(age, 0.25),
                   age_3rd_quart = quantile(age, 0.75),
                   n = n[1]) |> 
  tidyr::pivot_wider(names_from = strain_type, values_from = strain) |> 
  dplyr::mutate(Strains = paste0(H1N1, "\n", H3N2),
                season = as.character(season),
                Age = paste0(round(age_med, digits = 0), " (",
                             round(age_1st_quart, digits = 0), ", ",
                             round(age_3rd_quart, digits = 0), ")"),
                value = paste0(n, "\n", Age)) |> 
  dplyr::select(!c(H1N1, H3N2, age_med, age_1st_quart, age_3rd_quart,
                   Age, n)) |> 
  dplyr::arrange(season) |>  
  dplyr::relocate(Strains, .after = season) |> 
  tidyr::pivot_wider(names_from = group, values_from = value) |> 
  dplyr::rename("Season" = season) |> 
  flextable::flextable() |> 
  flextable::autofit() |> 
  flextable::footnote(i = 1, j = 3:5, part = "header", ref_symbols = "1",
                      value = flextable::as_paragraph("Number in group \n Median (IQR)"))


desc_table_wide

flextable::save_as_image(desc_table_wide, 
                         here::here("results/tables/wide_tb1.png")
                         )

## VEs--------------------------------------------------------------------------
ve_res <- 
  ve |> 
  dplyr::mutate(value = paste0(round(boot_ve*100, 2), " \n (", 
                               round(boot_ve_lwr*100, 2), ", ",
                               round(boot_ve_upr*100, 2), ")"),
                Season = factor(season),
                Strain = paste0(strain_type, "\n", strains_short)) |> 
  dplyr::select(!c(contains("boot"),data, statistic, season, strain_type, 
                                      strains_short)) |> 
  tidyr::pivot_wider(names_from = group, values_from = value) |> 
  flextable::flextable() |> 
  flextable::merge_at(i = 1:2, j = 1, part = "body") |> 
  flextable::merge_at(i = 3:4, j = 1, part = "body") |> 
  flextable::merge_at(i = 5:6, j = 1, part = "body") |> 
  flextable::merge_at(i = 7:8, j = 1, part = "body") |> 
  flextable::merge_at(i = 9:10, j = 1, part = "body") |> 
  flextable::merge_at(i = 11:12, j = 1, part = "body") |> 
  flextable::merge_at(i = 13:14, j = 1, part = "body") |> 
  flextable::merge_at(i = 15:16, j = 1, part = "body") |> 
  flextable::merge_at(i = 17:18, j = 1, part = "body") |> 
  flextable::hline(i = 2, part = "body") |> 
  flextable::hline(i = 4, part = "body") |> 
  flextable::hline(i = 6, part = "body") |> 
  flextable::hline(i = 8, part = "body") |> 
  flextable::hline(i = 10, part = "body") |> 
  flextable::hline(i = 12, part = "body") |> 
  flextable::hline(i = 14, part = "body") |> 
  flextable::hline(i = 16, part = "body") |> 
  flextable::hline(i = 18, part = "body") 

ve_res

saveRDS(ve_res, here::here("results/tables/ve.rds"))

## DVEs-------------------------------------------------------------------------
dve_res <- 
  ive |> 
  dplyr::mutate(value = paste0(round(boot_est*100, 2), " \n (", 
                               round(boot_est_lwr*100, 2), ", ",
                               round(boot_est_upr*100, 2), ")"),
                Season = factor(season),
                Strain = paste0(strain_type, "\n", strains_short)) |> 
  dplyr::select(!c(contains("boot"),data, season, strain_type, 
                   strains_short)) |> 
  dplyr::rename("Type" = type) |> 
  dplyr::relocate(Type, .after = Strain) |> 
  tidyr::pivot_wider(names_from = comparison, values_from = value) |> 
  flextable::flextable() |> 
  flextable::merge_at(i = 1:4, j = 1, part = "body") |>
  flextable::merge_at(i = 5:8, j = 1, part = "body") |>
  flextable::merge_at(i = 9:12, j = 1, part = "body") |>
  flextable::merge_at(i = 13:16, j = 1, part = "body") |>
  flextable::merge_at(i = 17:20, j = 1, part = "body") |>
  flextable::merge_at(i = 21:24, j = 1, part = "body") |>
  flextable::merge_at(i = 25:28, j = 1, part = "body") |>
  flextable::merge_at(i = 29:32, j = 1, part = "body") |>
  flextable::merge_at(i = 33:36, j = 1, part = "body") |> 
  flextable::hline(i = 4, part = "body") |>
  flextable::hline(i = 8, part = "body") |>
  flextable::hline(i = 12, part = "body") |>
  flextable::hline(i = 16, part = "body") |>
  flextable::hline(i = 20, part = "body") |>
  flextable::hline(i = 24, part = "body") |>
  flextable::hline(i = 28, part = "body") |>
  flextable::hline(i = 32, part = "body") |>
  flextable::hline(i = 36, part = "body")
  
dve_res

saveRDS(dve_res, here::here("results/tables/dve.rds"))

