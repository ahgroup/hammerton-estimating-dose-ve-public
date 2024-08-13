################################### HEADER #####################################
### SCRIPT: Figure Generation                                                ###
### AUTHOR: Savannah Hammerton                                               ###
### PURPOSE: This script will generate all the figures for the manuscript.   ###
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
library(ggrepel)

## Load data--------------------------------------------------------------------
# Processed homologous dataset 
homologous_dat <- readRDS(here::here("data/raw-data/homologous_dat.rds"))
# Predicted protection
protection <- readRDS(here::here("results/output/protection.rds"))
# Bootstrapped VEs
ve <- readRDS(here::here("results/output/boot_ve.rds"))
# Bootstrapped IVEs
ive <- readRDS(here::here("results/output/boot_ive.rds"))
# Bootstrapped VEs - individual VE
ve_ind <- readRDS(here::here("results/output/boot_ve_ind.rds"))
# Bootstrapped IVEs - individual VE 
ive_ind <- readRDS(here::here("results/output/boot_ive_ind.rds"))
# No LoD
ive_nolod <- readRDS(here::here("results/output/boot_ive_nolod.rds"))
# With minors 
ive_kids <- readRDS(here::here("results/output/boot_ive_kid.rds")) 
# Reported CDC seasonal VE 
ve_compare <- readRDS(here::here("results/output/ve_compare.rds"))
# Extracted Coudeville 2010 curve data 
extracted_curve <- 
  readr::read_csv(here::here("data/raw-data/coudeville_curve.txt"), skip = 1, 
                  col_names = c("postiter", "protect"))


## Script setup-----------------------------------------------------------------
here::here() #set paths
ggplot2::theme_set(ggplot2::theme_minimal()) #set plot theme
source(here::here("code/objects.R")) #load objects 
source(here::here("code/functions.R")) #load functions

# Figures-----------------------------------------------------------------------

## Curve replication------------------------------------------------------------

# Generate values between 1 and 350 for x axis/HAI titer
# Not starting at 0 to avoid problems taking the log of the titer
coude_hi <- data.frame(postiter = seq(from = 1, to = 350, length.out = 1000))

# Function to generate protection probabilities based on input data and 
# alpha/beta parameter estimates 
curve_fx <- function(dat, alpha, beta) {
  
  dat$protect <- 1 - (1/(1 + exp(beta*(log(dat$postiter) - alpha))))
  
  return(dat)
  
} 

# Applying curve application function to generated HAI titers with 3.09 and 
# 1.42 alpha and beta coefficients 
dat <-
  curve_fx(dat = coude_hi, alpha = 3.09, beta = 1.42)

# Figure with extracted and estimated curves overlayed 
curve <- 
  extracted_curve |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = postiter, y = protect) +
  ggplot2::geom_line(linetype = 1, linewidth = 1.25,
                     ggplot2::aes(color = "Extracted Curve")) +
  ggplot2::geom_line(data = dat, linetype = 2, linewidth = 1.25,
                     ggplot2::aes(color = "Estimated Curve")) +
  ggplot2::scale_color_manual(values = c("Extracted Curve" = "black",
                                         "Estimated Curve" = "red"),
                              name = "Curve Version") +
  ggplot2::xlab("HAI Titer") +
  ggplot2::ylab("Protection") +
  ggplot2::theme(legend.position = "bottom")
curve

ggplot2::ggsave(here::here("results/figures/curve.png"),
                curve,
                height = 4, width = 6.5)

## Ages-------------------------------------------------------------------------
under18 <- 
  homologous_dat |> 
  dplyr::group_by(season) |> 
  dplyr::filter(strain_type == "H1N1") |> 
  dplyr::mutate(minors = ifelse(age < 18, 1, 0)) |> 
  dplyr::summarise(under18 = sum(minors)) 


kids <- homologous_dat |> 
  dplyr::filter(group == "YASD", strain_type == "H1N1") |> 
  dplyr::left_join(under18, by = "season") |> 
  ggplot2::ggplot()+
  ggplot2::aes(x = season, y = age, label = under18) +
  ggplot2::geom_count() +
  ggplot2::geom_label(ggplot2::aes(y = 8)) +
  ggplot2::geom_hline(ggplot2::aes(yintercept = 18), color = "red") +
  ggplot2::scale_x_continuous(breaks = 2013:2021) +
  ggplot2::scale_y_continuous(breaks = seq(8, 66, 4)) +
  ggplot2::labs(x= "Season",y="Age", title = "Ages of YASD group")

kids

ggplot2::ggsave(here::here("results/figures/kids.png"),
                height = 5, width = 4.25)

## Titer distributions----------------------------------------------------------
### H1N1------------------------------------------------------------------------
h1n1_pretiters <- 
  homologous_dat |> 
  dplyr::filter(strain_type=="H1N1",
                age >= 18) |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = season, y = log(pretiter), 
               group = interaction(season, group),
               fill = group) +
  ggplot2::scale_x_continuous(breaks = 2013:2021,
                              labels = labelsh1n1) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::geom_boxplot(alpha = 0.85) +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::labs(fill = "Group",
                x = "Season",
                y = "Log HAI Titer") +
  ggplot2::scale_fill_manual(values = group_colors) +
  ggplot2::labs(x = NULL) +
  ggplot2::ggtitle("H1N1 Pre-Vaccination")

h1n1_pretiters

ggplot2::ggsave(plot = h1n1_pretiters, 
                filename = here::here("results/figures/h1n1_pretiters.png"),
                width = 5, height = 3.5)

h1n1_postiters <- 
  homologous_dat |> 
  dplyr::filter(strain_type=="H1N1",
                age >= 18) |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = season, y = log(postiter), 
               group = interaction(season, group),
               fill = group) +
  ggplot2::scale_x_continuous(breaks = 2013:2021,
                              labels = labelsh1n1) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::geom_boxplot(alpha = 0.85) +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::labs(fill = "Group",
                x = "Season",
                y = "Log HAI Titer") +
  ggplot2::scale_fill_manual(values = group_colors) +
  ggplot2::labs(x = NULL) +
  ggplot2::ggtitle("H1N1 Post-Vaccination")


h1n1_postiters

ggplot2::ggsave(plot = h1n1_postiters, 
                filename = here::here("results/figures/h1n1_postiters.png"),
                width = 5, height = 3.5)


# Combined pre and post
h1n1_post_patch <- 
  h1n1_postiters +
  ggplot2::theme(axis.text.x = ggplot2::element_blank())

h1n1_combined_titers <- 
  (h1n1_pretiters + h1n1_post_patch) + 
  patchwork::plot_layout(guides = "collect", nrow = 1) &
  ggplot2::theme(legend.position = "bottom",
                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 0.7)) &
  ggplot2::ylim(1.6,8.6)

h1n1_combined_titers

ggplot2::ggsave(h1n1_combined_titers,
                filename = here::here("results/figures/h1n1_combined_titers.png"),
                width = 8.5, height = 3.5)

### H3N2------------------------------------------------------------------------
h3n2_pretiters <- 
  homologous_dat |> 
  dplyr::filter(strain_type=="H3N2",
                age >= 18) |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = season, y = log(pretiter), group = interaction(season, group),
               fill = group) +
  ggplot2::scale_x_continuous(breaks = 2013:2021,
                              labels = labelsh3n2) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::geom_boxplot(alpha = 0.85) +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::labs(fill = "Group",
                x = "Season",
                y = "Log HAI Titer") +
  ggplot2::scale_fill_manual(values = group_colors) +
  ggplot2::labs(x = NULL) +
  ggplot2::ggtitle("H3N2 Pre-Vaccination")

h3n2_pretiters

ggplot2::ggsave(plot = h3n2_pretiters, 
                filename = here::here("results/figures/h3n2_pretiters.png"),
                width = 5, height = 3.5)

h3n2_postiters <- 
  homologous_dat |> 
  dplyr::filter(strain_type=="H3N2",
                age >= 18) |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = season, y = log(postiter), 
               group = interaction(season, group),
               fill = group) +
  ggplot2::scale_x_continuous(breaks = 2013:2021,
                              labels = labelsh3n2) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::geom_boxplot(alpha = 0.85) +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::labs(fill = "Group",
                x = "Season",
                y = "Log HAI Titer") +
  ggplot2::scale_fill_manual(values = group_colors) +
  ggplot2::labs(x = NULL) +
  ggplot2::ggtitle("H3N2 Post-Vaccination")
h3n2_postiters

ggplot2::ggsave(plot = h3n2_postiters, 
                filename = here::here("results/figures/h3n2_postiters.png"),
                width = 5, height = 3.5)

# Combined pre and post
h3n2_post_patch <- 
  h3n2_postiters +
  ggplot2::theme(axis.text.x = ggplot2::element_blank())

h3n2_combined_titers <- 
  (h3n2_pretiters + h3n2_post_patch) + 
  patchwork::plot_layout(guides = "collect", nrow = 1) &
  ggplot2::theme(legend.position = "bottom",
                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 0.7)) &
  ggplot2::ylim(1.6,8.6)

h3n2_combined_titers

ggplot2::ggsave(h3n2_combined_titers,
                filename = here::here("results/figures/h3n2_combined_titers.png"),
                width = 8.5, height = 3.5)


### All combined----------------------------------------------------------------
titer_patch <- 
  (h1n1_pretiters + h1n1_postiters) / (h3n2_pretiters + h3n2_postiters) +
  patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom",
                 axis.text.x = ggplot2::element_text(
                   angle = 45, hjust = 1
                 )) &
  ggplot2::ylim(1.6,8.6)

titer_patch

ggplot2::ggsave(titer_patch,
                filename = here::here("results/figures/combined_titers.png"),
                width = 8, height = 5)

## LoDs-------------------------------------------------------------------------
h1n1_lod <- 
  homologous_dat |> 
  dplyr::filter(strain_type == "H1N1",
                age >= 18) |> 
  tidyr::pivot_longer(c(pretiter, postiter), 
                      names_to = "stat", values_to = "titer") |> 
  dplyr::mutate(LoD = factor(ifelse(titer == 5, 1, 0),
                             labels = c("Not LoD", "LoD"))) |> 
  ggplot2::ggplot()+
  ggplot2::aes(x = season, y = log(titer), color = LoD)+
  ggplot2::geom_count() +
  ggplot2::scale_x_continuous(breaks = 2013:2021,
                              labels = labelsh1n1) +
  ggplot2::facet_wrap(~fct_rev(stat) + group) +
  ggplot2::labs(x = "Season", y = "Log(Titer)", title = "H1N1 Titers") +
  ggplot2::scale_size_area(breaks = seq(25, 125, 25),
                           limits = c(0, 200)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))


h1n1_lod

ggplot2::ggsave(here::here("results/figures/h1n1_lod.png"),
                height = 5, width = 11)


h3n2_lod <- 
  homologous_dat |> 
  dplyr::filter(strain_type == "H3N2",
                age >= 18) |> 
  tidyr::pivot_longer(c(pretiter, postiter), 
                      names_to = "stat", values_to = "titer") |> 
  dplyr::mutate(LoD = factor(ifelse(titer == 5, 1, 0),
                             labels = c("Not LoD", "LoD"))) |> 
  ggplot2::ggplot()+
  ggplot2::aes(x = season, y = log(titer), color = LoD)+
  ggplot2::geom_count() +
  ggplot2::scale_x_continuous(breaks = 2013:2021,
                              labels = labelsh3n2) +
  ggplot2::facet_wrap(~fct_rev(stat) + group) +
  ggplot2::labs(x = "Season", y = "Log(Titer)", title = "H3N2 Titers") +
  ggplot2::scale_size_area(breaks = seq(25, 125, 25),
                           limits = c(0, 200)) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))


h3n2_lod

ggplot2::ggsave(here::here("results/figures/h3n2_lod.png"),
                height = 5, width = 11)

lod <- 
  (h1n1_lod / h3n2_lod ) /
  patchwork::plot_layout(guides = "collect",
                         axis_titles = "collect") +
  patchwork::plot_annotation(tag_levels = 'A') &
  ggplot2::theme(legend.position = "bottom",
                 axis.text.x = ggplot2::element_text(size = 7),
                 axis.title.x = ggplot2::element_text(size = 8))


lod


ggplot2::ggsave(here::here("results/figures/lod.png"), plot = lod,
                height = 9, width = 6)


## Probability distributions----------------------------------------------------
### H1N1------------------------------------------------------------------------
h1n1_preprotect <- 
  protection |> 
  dplyr::filter(strain_type=="H1N1",
                age >= 18) |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = season, y = pre_protect, 
               group = interaction(season, group),
               fill = group) +
  ggplot2::scale_x_continuous(breaks = 2013:2021,
                              labels = labelsh1n1) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::geom_boxplot(alpha = 0.85) +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::labs(fill = "Group",
                x = "Season",
                y = "Probability of Protection") +
  ggplot2::scale_fill_manual(values = group_colors) +
  ggplot2::labs(x = NULL) +
  ggplot2::ggtitle("H1N1 Pre-Vaccination")

h1n1_preprotect

ggplot2::ggsave(plot = h1n1_preprotect, 
                filename = here::here("results/figures/h1n1_preprotect.png"),
                width = 5, height = 3.5)

h1n1_postiters <- 
  homologous_dat |> 
  dplyr::filter(strain_type=="H1N1",
                age >= 18) |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = season, y = log(postiter), 
               group = interaction(season, group),
               fill = group) +
  ggplot2::scale_x_continuous(breaks = 2013:2021,
                              labels = labelsh1n1) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::geom_boxplot(alpha = 0.85) +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::labs(fill = "Group",
                x = "Season",
                y = "Log HAI Titer") +
  ggplot2::scale_fill_manual(values = group_colors) +
  ggplot2::labs(x = NULL) +
  ggplot2::ggtitle("H1N1 Post-Vaccination")


h1n1_postiters

ggplot2::ggsave(plot = h1n1_postiters, 
                filename = here::here("results/figures/h1n1_postiters.png"),
                width = 5, height = 3.5)

h1n1_postprotect <- 
  protection |> 
  dplyr::filter(strain_type=="H1N1") |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = season, y = post_protect, 
               group = interaction(season, group),
               fill = group) +
  ggplot2::scale_x_continuous(breaks = 2013:2021,
                              labels = labelsh1n1) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::geom_boxplot(alpha = 0.85) +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::labs(fill = "Group",
                x = "Season",
                y = "Probability of Protection") +
  ggplot2::scale_fill_manual(values = group_colors) +
  ggplot2::labs(x = NULL) +
  ggplot2::ggtitle("H1N1 Post-Vaccination")

h1n1_postprotect

ggplot2::ggsave(plot = h1n1_postprotect, 
                filename = here::here("results/figures/h1n1_postprotect.png"),
                width = 5, height = 3.5)

## VEs--------------------------------------------------------------------------
### H1N1------------------------------------------------------------------------
h1n1_ve <- 
  ggplot2::ggplot(data = dplyr::filter(ve, strain_type == "H1N1",
                                       data == "Main analysis"),
  ) +
  ggplot2::aes(x = season, y = (boot_ve*100), 
               group = factor(group, labels =  c("YASD", "OAHD", "OASD"),
                              levels = c("YASD", "OAHD", "OASD")), 
               color = factor(group, labels =  c("YASD", "OAHD", "OASD"),
                              levels = c("YASD", "OAHD", "OASD")),
               fill =  factor(group, labels =  c("YASD", "OAHD", "OASD"),
                              levels = c("YASD", "OAHD", "OASD")),
               shape = factor(group, labels =  c("YASD", "OAHD", "OASD"),
                              levels = c("YASD", "OAHD", "OASD"))) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_ve_lwr*100, ymax = boot_ve_upr*100),
    linewidth = 0.75, alpha = 0.6, 
    position = ggplot2::position_dodge(width = 0.75)) +
  ggplot2::geom_point(size = 2, alpha = 0.85,
                      position = ggplot2::position_dodge(width = 0.75)) +
  ggplot2::labs(color = "Group",
                shape = "Group",
                fill = "Group",
                x = "Season",
                y = "VE and 95% CI",
                title = "H1N1") +
  ggplot2::scale_x_continuous(breaks = 2013:2021,
                              labels = labelsh1n1) +
  ggplot2::scale_color_manual(values = group_colors) +
  ggplot2::scale_shape_manual(values = group_shapes) +
  
  ggplot2::scale_fill_manual(values = group_colors) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::ylim(c(0,100)) 


h1n1_ve

ggplot2::ggsave(plot = h1n1_ve, 
                filename = here::here("results/figures/h1n1_ve.png"),
                width = 5, height = 3.5)

### H3N2------------------------------------------------------------------------
h3n2_ve <-
  ggplot2::ggplot(data = dplyr::filter(ve, strain_type == "H3N2",
                                       data == "Main analysis")) +
  ggplot2::aes(x = season, y = (boot_ve*100), 
               group = factor(group, labels =  c("YASD", "OAHD", "OASD"),
                              levels = c("YASD", "OAHD", "OASD")), 
               color = factor(group, labels =  c("YASD", "OAHD", "OASD"),
                              levels = c("YASD", "OAHD", "OASD")),
               shape =  factor(group, labels =  c("YASD", "OAHD", "OASD"),
                               levels = c("YASD", "OAHD", "OASD")),
               fill =  factor(group, labels =  c("YASD", "OAHD", "OASD"),
                              levels = c("YASD", "OAHD", "OASD"))) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_ve_lwr*100, ymax = boot_ve_upr*100),
    linewidth = 0.75, alpha = 0.6, 
    position = ggplot2::position_dodge(width = 0.75)) +
  ggplot2::geom_point(size = 2, alpha = 0.85,
                      position = ggplot2::position_dodge(width = 0.75)) +
  ggplot2::labs(color = "Group",
                shape = "Group",
                fill = "Group",
                x = "Season",
                y = "VE and 95% CI",
                title = "H3N2") +
  ggplot2::scale_x_continuous(breaks = 2013:2021,
                              labels = labelsh3n2) +
  ggplot2::scale_color_manual(values = group_colors) +
  ggplot2::scale_shape_manual(values = group_shapes) +
  ggplot2::scale_fill_manual(values = group_colors) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::ylim(c(0,100))

h3n2_ve

ggplot2::ggsave(plot = h3n2_ve, 
                filename = here::here("results/figures/h3n2_ve.png"),
                width = 5, height = 3.5)

### Combined--------------------------------------------------------------------
h3n2_patch <- 
  h3n2_ve +
  ggplot2::theme(axis.text.y = ggplot2::element_blank()) +
  ggplot2::labs(y = NULL)

combined_ve <- 
  (h1n1_ve + h3n2_patch) + patchwork::plot_layout( guides = "collect") & 
  ggplot2::theme(legend.position = "bottom",
                 axis.text.x = ggplot2::element_text(angle = 45, hjust = 0.7)) &
  ggplot2::labs(x = NULL)

combined_ve

ggplot2::ggsave(plot = combined_ve,
                filename = here::here("results/figures/combined_ve.png"),
                width = 8, height = 3.5)

## Individual vs Group VE comparison--------------------------------------------
ve_comp_h1n1 <-
  ve |> 
  dplyr::filter(data == "Main analysis", strain_type == "H1N1") |> 
  dplyr::bind_rows(dplyr::filter(ve_ind, strain_type == "H1N1")) |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = season, y = (boot_ve*100), color = data) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_ve_lwr*100, ymax = boot_ve_upr*100),
    linewidth = 0.75, alpha = 0.6, 
    position = ggplot2::position_dodge(width = 0.75)) +
  ggplot2::geom_point(size = 2, alpha = 0.85,
                      position = ggplot2::position_dodge(width = 0.75)) +
  ggplot2::facet_grid(group~strain_type) +
  ggplot2::scale_x_continuous(breaks = 2013:2021,
                              labels = labelsh1n1) +
  ggplot2::labs(x = "Season",
                y = "VE and 95% CI",
                color = "VE Source") 

ve_comp_h1n1

ve_comp_h3n2 <-
  ve |> 
  dplyr::filter(data == "Main analysis", strain_type == "H3N2") |> 
  dplyr::bind_rows(dplyr::filter(ve_ind, strain_type == "H3N2")) |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = season, y = (boot_ve*100), color = data) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_ve_lwr*100, ymax = boot_ve_upr*100),
    linewidth = 0.75, alpha = 0.6, 
    position = ggplot2::position_dodge(width = 0.75)) +
  ggplot2::geom_point(size = 2, alpha = 0.85,
                      position = ggplot2::position_dodge(width = 0.75)) +
  ggplot2::facet_grid(group~strain_type) +
  ggplot2::scale_x_continuous(breaks = 2013:2021,
                              labels = labelsh3n2) +
  ggplot2::labs(x = "Season",
                y = "VE and 95% CI",
                color = "VE Source") 

ve_comp_h3n2  

ve_versions <- 
  ve_comp_h1n1 + ve_comp_h3n2 + patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom", 
                 axis.text = ggplot2::element_text(angle = 30, hjust = 0.75)) 

ve_versions

ggplot2::ggsave(here::here("results/figures/ve_versions.png"),
                width = 8, height = 4)
## ADVE-------------------------------------------------------------------------
### OAHD_OASD-------------------------------------------------------------------
oahd_oasd <- 
  ggplot2::ggplot(data = dplyr::filter(ive, 
                                       comparison == "OAHD_OASD",
                                       type == "Absolute",
                                       data == "Main analysis")) +
  ggplot2::aes(x = season, y = boot_est*100, color = strain_type,
               shape = strain_type) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_est_lwr*100, ymax = boot_est_upr*100),
    alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(alpha = 0.95, position = ggplot2::position_dodge(1)) +
  ggplot2::scale_x_continuous(limits = c(2012.5, 2021.5),
                              breaks = 2013:2021,
                              labels = labelsboth) +
  ggplot2::coord_cartesian(ylim=c(-100, 150)) +
  ggplot2::scale_color_manual(values = ve_comp_colors, name = "Type") +
  ggplot2::labs(x = "Season \n(H1N1 Strain / H3N2 Strain)",
                y = NULL,
                title = "OAHD vs. OASD",
                shape = "Type") +
  ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::annotate("label", x = 2014, y = -70, 
                    label = "Better VE with SD \n vaccine in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::annotate("label", x = 2014, y = 120, 
                    label = "Better VE with HD \n vaccines in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::theme(legend.position = "bottom")

oahd_oasd

ggplot2::ggsave(here::here("results/figures/oahd_oasd.png"), oahd_oasd,
                height = 5, width = 4.5)

# For slides 
oahd_oasd_slides <- 
  ggplot2::ggplot(data = dplyr::filter(ive, 
                                       comparison == "OAHD_OASD",
                                       type == "Absolute",
                                       data == "Main analysis")) +
  ggplot2::aes(x = season, y = boot_est*100, color = strain_type) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_est_lwr*100, ymax = boot_est_upr*100),
    alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(alpha = 0.95, position = ggplot2::position_dodge(1)) +
  ggplot2::scale_x_continuous(limits = c(2012.5, 2021.5),
                              breaks = 2013:2021,
                              labels = labelsboth) +
  ggplot2::coord_cartesian(ylim=c(-100, 150)) +
  ggplot2::scale_color_manual(values = ve_comp_colors, name = "Type") +
  ggplot2::labs(x = "Season \n(H1N1 Strain / H3N2 Strain)",
                y = NULL,
                title = "OAHD vs. OASD") +
  ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::annotate("label", x = 2014, y = -70, 
                    label = "Better VE with SD \n vaccine in older adults",
                    fill = "#F3F3F370", color = "#BA0C2F", size = 2.5) +
  ggplot2::annotate("label", x = 2014, y = 120, 
                    label = "Better VE with HD \n vaccines in older adults",
                    fill = "#F3F3F370", color = "#BA0C2F", size = 2.5) +
  ggplot2::theme(legend.position = "right")
oahd_oasd_slides

ggplot2::ggsave(here::here("results/figures/oahd_oasd_slides.png"),
                oahd_oasd_slides,
                height = 5, width = 6.5)

### OAHD_YASD-------------------------------------------------------------------
oahd_yasd <- 
  ggplot2::ggplot(data = dplyr::filter(ive, 
                                       comparison == "OAHD_YASD",
                                       type == "Absolute",
                                       data == "Main analysis")) +
  ggplot2::aes(x = season, y = boot_est*100, color = strain_type,
               shape = strain_type) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_est_lwr*100, ymax = boot_est_upr*100),
    alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(alpha = 0.95, position = ggplot2::position_dodge(1)) +
  ggplot2::scale_x_continuous(limits = c(2012.5, 2021.5),
                              breaks = 2013:2021,
                              labels = labelsboth) +
  ggplot2::coord_cartesian(ylim=c(-100, 150)) +
  ggplot2::scale_color_manual(values = ve_comp_colors, name = "Type") +
  ggplot2::labs(x = "Season \n(H1N1 Strain / H3N2 Strain)",
                y = NULL,
                title = "OAHD vs. YASD",
                shape = "Type") +
  ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::annotate("label", x = 2014, y = -70, 
                    label = "Better VE with SD \n vaccine in younger adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::annotate("label", x = 2014, y = 120, 
                    label = "Better VE with HD \n vaccines in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::theme(legend.position = "bottom")

oahd_yasd

ggplot2::ggsave(here::here("results/figures/oahd_yasd.png"), oahd_yasd,
                height = 5, width = 4.5)

oahd_yasd_slides <- 
  ggplot2::ggplot(data = dplyr::filter(ive, 
                                       comparison == "OAHD_YASD",
                                       type == "Absolute",
                                       data == "Main analysis")) +
  ggplot2::aes(x = season, y = boot_est*100, color = strain_type) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_est_lwr*100, ymax = boot_est_upr*100),
    alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(alpha = 0.95, position = ggplot2::position_dodge(1)) +
  ggplot2::scale_x_continuous(limits = c(2012.5, 2021.5),
                              breaks = 2013:2021,
                              labels = labelsboth) +
  ggplot2::coord_cartesian(ylim=c(-100, 150)) +
  ggplot2::scale_color_manual(values = ve_comp_colors, name = "Type") +
  ggplot2::labs(x = "Season \n(H1N1 Strain / H3N2 Strain)",
                y = NULL,
                title = "OAHD vs. YASD") +
  ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::annotate("label", x = 2014, y = -70, 
                    label = "Better VE with SD \n vaccine in younger adults",
                    fill = "#F3F3F370", color = "#BA0C2F", size = 2.5) +
  ggplot2::annotate("label", x = 2014, y = 120, 
                    label = "Better VE with HD \n vaccines in older adults",
                    fill = "#F3F3F370", color = "#BA0C2F", size = 2.5) +
  ggplot2::theme(legend.position = "right")

oahd_yasd_slides

ggplot2::ggsave(here::here("results/figures/oahd_yasd_slides.png"), 
                oahd_yasd_slides,
                height = 5, width = 6.5)

### OASD_YASD-------------------------------------------------------------------
oasd_yasd <- 
  ggplot2::ggplot(data = dplyr::filter(ive, 
                                       comparison == "OASD_YASD",
                                       type == "Absolute",
                                       data == "Main analysis")) +
  ggplot2::aes(x = season, y = boot_est*100, color = strain_type,
               shape = strain_type) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_est_lwr*100, ymax = boot_est_upr*100),
    alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(alpha = 0.95, position = ggplot2::position_dodge(1)) +
  ggplot2::scale_x_continuous(limits = c(2012.5, 2021.5),
                              breaks = 2013:2021,
                              labels = labelsboth) +
  ggplot2::coord_cartesian(ylim=c(-100, 150)) +
  ggplot2::scale_color_manual(values = ve_comp_colors, name = "Type") +
  ggplot2::labs(x = "Season \n(H1N1 Strain / H3N2 Strain)",
                y = NULL,
                title = "OASD vs. YASD",
                shape = 'Type') +
  ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::annotate("label", x = 2014, y = -70, 
                    label = "Better VE with SD \n vaccine in younger adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::annotate("label", x = 2014, y = 120, 
                    label = "Better VE with SD \n vaccine in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::theme(legend.position = "bottom")

oasd_yasd

ggplot2::ggsave(here::here("results/figures/oasd_yasd.png"), oasd_yasd,
                height = 5, width = 4.5)


oasd_yasd_slides <- 
  ggplot2::ggplot(data = dplyr::filter(ive, 
                                       comparison == "OASD_YASD",
                                       type == "Absolute",
                                       data == "Main analysis")) +
  ggplot2::aes(x = season, y = boot_est*100, color = strain_type) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_est_lwr*100, ymax = boot_est_upr*100),
    alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(alpha = 0.95, position = ggplot2::position_dodge(1)) +
  ggplot2::scale_x_continuous(limits = c(2012.5, 2021.5),
                              breaks = 2013:2021,
                              labels = labelsboth) +
  ggplot2::coord_cartesian(ylim=c(-100, 150)) +
  ggplot2::scale_color_manual(values = ve_comp_colors, name = "Type") +
  ggplot2::labs(x = "Season \n(H1N1 Strain / H3N2 Strain)",
                y = NULL,
                title = "OASD vs. YASD") +
  ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::annotate("label", x = 2014, y = -70, 
                    label = "Better VE with SD \n vaccine in younger adults",
                    fill = "#F3F3F370", color = "#BA0C2F", size = 2.5) +
  ggplot2::annotate("label", x = 2014, y = 120, 
                    label = "Better VE with SD \n vaccine in older adults",
                    fill = "#F3F3F370", color = "#BA0C2F", size = 2.5) +
  ggplot2::theme(legend.position = "right")

oasd_yasd_slides

ggplot2::ggsave(here::here("results/figures/oasd_yasd_slides.png"), 
                oasd_yasd_slides,
                height = 5, width = 6.5)
#### Combined-------------------------------------------------------------------

adve <- 
  (oasd_yasd / oahd_yasd / oahd_oasd) /
  patchwork::plot_layout(guides = "collect",
                         axis_titles = "collect",
                         axes = "collect") +
  patchwork::plot_annotation(tag_levels = 'A') &
  ggplot2::theme(legend.position = "bottom",
                 axis.text.x = ggplot2::element_text(size = 7),
                 axis.title.x = ggplot2::element_text(size = 8))

adve

ggplot2::ggsave(here::here("results/figures/adve.png"), plot = adve,
                           height = 8, width = 5)

## With minors------------------------------------------------------------------
### OAHD_YASD-------------------------------------------------------------------
oahd_yasd_kids <-
  ggplot2::ggplot(data = dplyr::filter(ive_kids, 
                                       comparison == "OAHD_YASD",
                                       type == "Absolute")) +
  ggplot2::aes(x = season, y = boot_est*100, color = strain_type,
               shape = data, linetype = data) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_est_lwr*100, ymax = boot_est_upr*100),
    alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(alpha = 0.95, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_errorbar(data = dplyr::filter(ive, type == "Absolute", 
                                              comparison == "OAHD_YASD",
                                              data == "Main analysis"),
                         ggplot2::aes(ymin = boot_est_lwr*100, 
                                      ymax = boot_est_upr*100),
                         alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(data = dplyr::filter(ive, type == "Absolute", 
                                           comparison == "OAHD_YASD",
                                           data == "Main analysis"),
                      alpha = 0.95, position = ggplot2::position_dodge(1)) +
  ggplot2::scale_x_continuous(limits = c(2012.5, 2021.5),
                              breaks = 2013:2021,
                              labels = labelsboth) +
  ggplot2::coord_cartesian(ylim=c(-100, 150)) +
  ggplot2::scale_color_manual(values = ve_comp_colors, name = "Type") +
  ggplot2::labs(x = "Season (H1N1 Strain / H3N2 Strain)",
                y = NULL,
                title = "OAHD vs. YASD, all ages",
                shape = "ADVE Source",
                linetype = "ADVE Source") +
  ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::annotate("label", x = 2014.5, y = -90, 
                    label = "Better VE with SD \n vaccine in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::annotate("label", x = 2014.5, y = 140, 
                    label = "Better VE with HD \n vaccines in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::theme(legend.position = "bottom", legend.box = "vertical")

oahd_yasd_kids


ggplot2::ggsave(filename = here::here("results/figures/oahd_yasd_kids.png"),
                plot = oahd_yasd_kids,
                width = 5, height = 4.5)

### OASD_YASD-------------------------------------------------------------------
oasd_yasd_kids <-
  ggplot2::ggplot(data = dplyr::filter(ive_kids, 
                                       comparison == "OASD_YASD",
                                       type == "Absolute")) +
  ggplot2::aes(x = season, y = boot_est*100, color = strain_type,
               shape = data, linetype = data) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_est_lwr*100, ymax = boot_est_upr*100),
    alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(alpha = 0.95, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_errorbar(data = dplyr::filter(ive, type == "Absolute", 
                                              comparison == "OASD_YASD",
                                              data == "Main analysis"),
                         ggplot2::aes(ymin = boot_est_lwr*100, 
                                      ymax = boot_est_upr*100),
                         alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(data = dplyr::filter(ive, type == "Absolute", 
                                           comparison == "OASD_YASD",
                                           data == "Main analysis"),
                      alpha = 0.95, position = ggplot2::position_dodge(1)) +
  ggplot2::scale_x_continuous(limits = c(2012.5, 2021.5),
                              breaks = 2013:2021,
                              labels = labelsboth) +
  ggplot2::coord_cartesian(ylim=c(-100, 150)) +
  ggplot2::scale_color_manual(values = ve_comp_colors, name = "Type") +
  ggplot2::labs(x = "Season (H1N1 Strain / H3N2 Strain)",
                y = NULL,
                title = "OASD vs. YASD, all ages",
                shape = "ADVE Source",
                linetype = "ADVE Source") +
  ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::annotate("label", x = 2014.5, y = -90, 
                    label = "Better VE with SD \n vaccine in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::annotate("label", x = 2014.5, y = 140, 
                    label = "Better VE with HD \n vaccines in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::theme(legend.position = "bottom",
                 legend.box = "vertical")

oasd_yasd_kids


ggplot2::ggsave(filename = here::here("results/figures/oasd_yasd_kids.png"),
                plot = oasd_yasd_kids,
                width = 5, height = 4.5)

#### Combined-------------------------------------------------------------------

adve_kids <- 
  (oasd_yasd_kids / oahd_yasd_kids ) /
  patchwork::plot_layout(guides = "collect",
                         axis_titles = "collect",
                         axes = "collect") +
  patchwork::plot_annotation(tag_levels = 'A') &
  ggplot2::theme(legend.position = "bottom",
                 axis.text.x = ggplot2::element_text(size = 7),
                 axis.title.x = ggplot2::element_text(size = 8))

adve_kids

ggplot2::ggsave(here::here("results/figures/adve_kids.png"), 
                plot = adve_kids,
                height = 7, width = 5)

## No LoD Values----------------------------------------------------------------
### OAHD_OASD-------------------------------------------------------------------
oahd_oasd_nolod <-
  ggplot2::ggplot(data = dplyr::filter(ive_nolod, 
                                       comparison == "OAHD_OASD",
                                       type == "Absolute")) +
  ggplot2::aes(x = season, y = boot_est*100, color = strain_type,
               shape = data, linetype = data) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_est_lwr*100, ymax = boot_est_upr*100),
    alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_errorbar(data = dplyr::filter(ive, type == "Absolute", 
                                              comparison == "OAHD_OASD",
                                              data == "Main analysis"),
                         ggplot2::aes(ymin = boot_est_lwr*100, 
                                      ymax = boot_est_upr*100),
                         alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(data = dplyr::filter(ive, type == "Absolute", 
                                           comparison == "OAHD_OASD",
                                           data == "Main analysis"),
                      alpha = 0.95, position = ggplot2::position_dodge(1)) +
  ggplot2::scale_x_continuous(limits = c(2012.5, 2021.5),
                              breaks = 2013:2021,
                              labels = labelsboth) +
  ggplot2::coord_cartesian(ylim=c(-100, 150)) +
  ggplot2::scale_color_manual(values = ve_comp_colors, name = "Type") +
  ggplot2::labs(x = "Season (H1N1 Strain / H3N2 Strain)",
                y = NULL,
                title = "OAHD vs. OASD, no LoD values",
                shape = "ADVE Source",
                linetype = "ADVE Source") +
  ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::annotate("label", x = 2014.5, y = -90, 
                    label = "Better VE with SD \n vaccine in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::annotate("label", x = 2014.5, y = 140, 
                    label = "Better VE with HD \n vaccines in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::theme(legend.position = "bottom", legend.box = "vertical")

oahd_oasd_nolod


ggplot2::ggsave(filename = here::here("results/figures/oahd_oasd_nolod.png"),
                plot = oahd_oasd_nolod,
                width = 5, height = 4.5)

### OAHD_YASD-------------------------------------------------------------------
oahd_yasd_nolod <-
  ggplot2::ggplot(data = dplyr::filter(ive_nolod, 
                                       comparison == "OAHD_YASD",
                                       type == "Absolute")) +
  ggplot2::aes(x = season, y = boot_est*100, color = strain_type,
               shape = data, linetype = data) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_est_lwr*100, ymax = boot_est_upr*100),
    alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_errorbar(data = dplyr::filter(ive, type == "Absolute", 
                                              comparison == "OAHD_YASD",
                                              data == "Main analysis"),
                         ggplot2::aes(ymin = boot_est_lwr*100, 
                                      ymax = boot_est_upr*100),
                         alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(data = dplyr::filter(ive, type == "Absolute", 
                                           comparison == "OAHD_YASD",
                                           data == "Main analysis"),
                      alpha = 0.95, position = ggplot2::position_dodge(1)) +
  ggplot2::scale_x_continuous(limits = c(2012.5, 2021.5),
                              breaks = 2013:2021,
                              labels = labelsboth) +
  ggplot2::coord_cartesian(ylim=c(-100, 150)) +
  ggplot2::scale_color_manual(values = ve_comp_colors, name = "Type") +
  ggplot2::labs(x = "Season (H1N1 Strain / H3N2 Strain)",
                y = NULL,
                title = "OAHD vs. YASD, no LoD values",
                shape = "ADVE Source",
                linetype = "ADVE Source") +
  ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::annotate("label", x = 2014.5, y = -90, 
                    label = "Better VE with SD \n vaccine in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::annotate("label", x = 2014.5, y = 140, 
                    label = "Better VE with HD \n vaccines in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::theme(legend.position = "bottom", legend.box = "vertical")

oahd_yasd_nolod


ggplot2::ggsave(filename = here::here("results/figures/oahd_yasd_nolod.png"),
                plot = oahd_yasd_nolod,
                width = 5, height = 4.5)

### OASD_YASD-------------------------------------------------------------------
oasd_yasd_nolod <-
  ggplot2::ggplot(data = dplyr::filter(ive_nolod, 
                                       comparison == "OASD_YASD",
                                       type == "Absolute")) +
  ggplot2::aes(x = season, y = boot_est*100, color = strain_type,
               shape = data, linetype = data) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_est_lwr*100, ymax = boot_est_upr*100),
    alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_errorbar(data = dplyr::filter(ive, type == "Absolute", 
                                              comparison == "OASD_YASD",
                                              data == "Main analysis"),
                         ggplot2::aes(ymin = boot_est_lwr*100, 
                                      ymax = boot_est_upr*100),
                         alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(data = dplyr::filter(ive, type == "Absolute", 
                                           comparison == "OASD_YASD",
                                           data == "Main analysis"),
                      alpha = 0.95, position = ggplot2::position_dodge(1)) +
  ggplot2::scale_x_continuous(limits = c(2012.5, 2021.5),
                              breaks = 2013:2021,
                              labels = labelsboth) +
  ggplot2::coord_cartesian(ylim=c(-100, 150)) +
  ggplot2::scale_color_manual(values = ve_comp_colors, name = "Type") +
  ggplot2::labs(x = "Season (H1N1 Strain / H3N2 Strain)",
                y = NULL,
                title = "OASD vs. YASD, no LoD values",
                shape = "ADVE Source",
                linetype = "ADVE Source") +
  ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::annotate("label", x = 2014.5, y = -90, 
                    label = "Better VE with SD \n vaccine in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::annotate("label", x = 2014.5, y = 140, 
                    label = "Better VE with HD \n vaccines in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::theme(legend.position = "bottom",
                 legend.box = "vertical")

oasd_yasd_nolod


ggplot2::ggsave(filename = here::here("results/figures/oasd_yasd_nolod.png"),
                plot = oasd_yasd_nolod,
                width = 5, height = 4.5)

#### Combined-------------------------------------------------------------------

adve_nolod <- 
  (oasd_yasd_nolod / oahd_yasd_nolod / oahd_oasd_nolod) /
  patchwork::plot_layout(guides = "collect",
                         axis_titles = "collect",
                         axes = "collect") +
  patchwork::plot_annotation(tag_levels = 'A') &
  ggplot2::theme(legend.position = "bottom",
                 axis.text.x = ggplot2::element_text(size = 7),
                 axis.title.x = ggplot2::element_text(size = 8))

adve_nolod

ggplot2::ggsave(here::here("results/figures/adve_nolod.png"), 
                plot = adve_nolod,
                height = 8, width = 5)


## Individual-based VE----------------------------------------------------------
### OAHD_OASD-------------------------------------------------------------------
oahd_oasd_ind <-
  ggplot2::ggplot(data = dplyr::filter(ive_ind, 
                                       comparison == "OAHD_OASD",
                                       type == "Absolute")) +
  ggplot2::aes(x = season, y = boot_est*100, color = strain_type,
               shape = data, linetype = data) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_est_lwr*100, ymax = boot_est_upr*100),
    alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(alpha = 0.95, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_errorbar(data = dplyr::filter(ive, type == "Absolute", 
                                              comparison == "OAHD_OASD",
                                              data == "Main analysis"),
                         ggplot2::aes(ymin = boot_est_lwr*100, 
                                      ymax = boot_est_upr*100),
                         alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(data = dplyr::filter(ive, type == "Absolute", 
                                           comparison == "OAHD_OASD",
                                           data == "Main analysis"),
                      alpha = 0.95, position = ggplot2::position_dodge(1)) +
  ggplot2::scale_x_continuous(limits = c(2012.5, 2021.5),
                              breaks = 2013:2021,
                              labels = labelsboth) +
  ggplot2::coord_cartesian(ylim=c(-100, 150)) +
  ggplot2::scale_color_manual(values = ve_comp_colors, name = "Type") +
  ggplot2::labs(x = "Season (H1N1 Strain / H3N2 Strain)",
                y = NULL,
                title = "OAHD vs. OASD, individual-based VE",
                shape = "ADVE Source",
                linetype = "ADVE Source") +
  ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::annotate("label", x = 2014.5, y = -90, 
                    label = "Better VE with SD \n vaccine in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::annotate("label", x = 2014.5, y = 140, 
                    label = "Better VE with HD \n vaccines in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::theme(legend.position = "bottom", legend.box = "vertical")

oahd_oasd_ind


ggplot2::ggsave(filename = here::here("results/figures/oahd_oasd_ind.png"),
                plot = oahd_oasd_ind,
                width = 5, height = 4.5)

### OAHD_YASD-------------------------------------------------------------------
oahd_yasd_ind <-
  ggplot2::ggplot(data = dplyr::filter(ive_ind, 
                                       comparison == "OAHD_YASD",
                                       type == "Absolute")) +
  ggplot2::aes(x = season, y = boot_est*100, color = strain_type,
               shape = data, linetype = data) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_est_lwr*100, ymax = boot_est_upr*100),
    alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(alpha = 0.95, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_errorbar(data = dplyr::filter(ive, type == "Absolute", 
                                              comparison == "OAHD_YASD",
                                              data == "Main analysis"),
                         ggplot2::aes(ymin = boot_est_lwr*100, 
                                      ymax = boot_est_upr*100),
                         alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(data = dplyr::filter(ive, type == "Absolute", 
                                           comparison == "OAHD_YASD",
                                           data == "Main analysis"),
                      alpha = 0.95, position = ggplot2::position_dodge(1)) +
  ggplot2::scale_x_continuous(limits = c(2012.5, 2021.5),
                              breaks = 2013:2021,
                              labels = labelsboth) +
  ggplot2::coord_cartesian(ylim=c(-100, 150)) +
  ggplot2::scale_color_manual(values = ve_comp_colors, name = "Type") +
  ggplot2::labs(x = "Season (H1N1 Strain / H3N2 Strain)",
                y = NULL,
                title = "OAHD vs. YASD, individual-based VE",
                shape = "ADVE Source",
                linetype = "ADVE Source") +
  ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::annotate("label", x = 2014.5, y = -90, 
                    label = "Better VE with SD \n vaccine in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::annotate("label", x = 2014.5, y = 140, 
                    label = "Better VE with HD \n vaccines in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::theme(legend.position = "bottom", legend.box = "vertical")

oahd_yasd_ind


ggplot2::ggsave(filename = here::here("results/figures/oahd_yasd_ind.png"),
                plot = oahd_yasd_ind,
                width = 5, height = 4.5)

### OASD_YASD-------------------------------------------------------------------
oasd_yasd_ind <-
  ggplot2::ggplot(data = dplyr::filter(ive_ind, 
                                       comparison == "OASD_YASD",
                                       type == "Absolute")) +
  ggplot2::aes(x = season, y = boot_est*100, color = strain_type,
               shape = data, linetype = data) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_est_lwr*100, ymax = boot_est_upr*100),
    alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(alpha = 0.95, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_errorbar(data = dplyr::filter(ive, type == "Absolute", 
                                              comparison == "OASD_YASD",
                                              data == "Main analysis"),
                         ggplot2::aes(ymin = boot_est_lwr*100, ymax = boot_est_upr*100),
                         alpha = 0.75, position = ggplot2::position_dodge(1)) +
  ggplot2::geom_point(data = dplyr::filter(ive, type == "Absolute", 
                                           comparison == "OASD_YASD",
                                           data == "Main analysis"),
                      alpha = 0.95, position = ggplot2::position_dodge(1)) +
  ggplot2::scale_x_continuous(limits = c(2012.5, 2021.5),
                              breaks = 2013:2021,
                              labels = labelsboth) +
  ggplot2::coord_cartesian(ylim=c(-100, 150)) +
  ggplot2::scale_color_manual(values = ve_comp_colors, name = "Type") +
  ggplot2::labs(x = "Season (H1N1 Strain / H3N2 Strain)",
                y = NULL,
                title = "OASD vs. YASD, individual-based VE",
                shape = "ADVE Source",
                linetype = "ADVE Source") +
  ggplot2::geom_hline(ggplot2::aes(yintercept = 0)) +
  ggplot2::geom_vline(xintercept = seq(2012.5, 2021.5, 1),
                      alpha = 0.4) +
  ggplot2::annotate("label", x = 2014.5, y = -90, 
                    label = "Better VE with SD \n vaccine in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::annotate("label", x = 2014.5, y = 140, 
                    label = "Better VE with HD \n vaccines in older adults",
                    fill = "#F3F3F370", color = "black", size = 2.5) +
  ggplot2::theme(legend.position = "bottom",
                 legend.box = "vertical")

oasd_yasd_ind


ggplot2::ggsave(filename = here::here("results/figures/oasd_yasd_ind.png"),
                plot = oasd_yasd_ind,
                width = 5, height = 4.5)

### Combined--------------------------------------------------------------------

adve_ind <- 
  (oasd_yasd_ind / oahd_yasd_ind / oahd_oasd_ind) /
  patchwork::plot_layout(guides = "collect",
                         axis_titles = "collect",
                         axes = "collect") +
  patchwork::plot_annotation(tag_levels = 'A') &
  ggplot2::theme(legend.position = "bottom",
                 axis.text.x = ggplot2::element_text(size = 7),
                 axis.title.x = ggplot2::element_text(size = 8))

adve_ind

ggplot2::ggsave(here::here("results/figures/adve_ind.png"), 
                plot = adve_ind,
                height = 8, width = 5)

## CDC Comparison---------------------------------------------------------------

cdc_seasons <- 
  ve_compare |> 
  dplyr::mutate(age_group = factor(age_group, 
                                   levels = c("18–49", "50–64", "≥65"))) |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = season, y = VE) +
  ggplot2::geom_point(
    ggplot2::aes(color = "CDC"), shape = "square"
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = lwr, ymax = upr, color = "CDC")
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y = boot_ve*100, color = strain_type)
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_ve_lwr*100,
                 ymax = boot_ve_upr*100, color = strain_type),
  ) +
  ggplot2::theme(axis.text.x = 
                   ggplot2::element_text(angle = 30, vjust = 0.55))  +
  ggplot2::scale_color_manual(values = ve_comp_colors, name = "VE Type") +
  ggplot2::facet_wrap(~age_group) +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::labs(x = "Season", y = "VE (%) and 95% CI")


cdc_seasons

ggplot2::ggsave(here::here("results/figures/cdc_compare.png"), cdc_seasons,
                height = 4, width = 7.25)

# CDC comparison plot including information on dominance 
cdc_seasons_circ <- 
  ve_compare |> 
  dplyr::mutate(age_group = factor(age_group, 
                                   levels = c("18–49", "50–64", "≥65"))) |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = season, y = VE) +
  ggplot2::geom_point(
    ggplot2::aes(color = "CDC"), shape = "square"
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = lwr, ymax = upr, color = "CDC")
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y = boot_ve*100, color = strain_type)
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_ve_lwr*100,
                 ymax = boot_ve_upr*100, color = strain_type)
  ) +
  ggplot2::theme(axis.text.x = 
                   ggplot2::element_text(angle = 30, vjust = 0.55))  +
  ggplot2::scale_color_manual(values = ve_comp_colors, name = "VE Type") +
  ggplot2::facet_wrap(~age_group) +
  ggrepel::geom_text_repel(ggplot2::aes(label = paste0(percemt, "%", 
                                                       stringr::str_sub(
                                                         strain_type, 1L, 2L)),
                                        y = -35),
                           nudge_y = ifelse(ve_compare$strain_type == "H1N1",
                                            1, 0)) +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::labs(x = "Season", y = "VE (%) and 95% CI")

cdc_seasons_circ

ggplot2::ggsave(here::here("results/figures/cdc_compare_circ.png"),
                cdc_seasons_circ,
                height = 8, width = 14)

# Just one age group 
cdc_seasons_circ_one <- 
  ve_compare |> 
  dplyr::filter(age_group ==  "18–49") |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = season, y = VE) +
  ggplot2::geom_point(
    ggplot2::aes(color = "CDC"), shape = "square"
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = lwr, ymax = upr, color = "CDC")
  ) +
  ggplot2::geom_point(
    ggplot2::aes(y = boot_ve*100, color = strain_type)
  ) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_ve_lwr*100,
                 ymax = boot_ve_upr*100, color = strain_type)
  ) +
  ggplot2::theme(axis.text.x = 
                   ggplot2::element_text(angle = 30, vjust = 0.55))  +
  ggplot2::scale_color_manual(values = ve_comp_colors, name = "VE Type") +
  ggrepel::geom_text_repel(ggplot2::aes(label = paste0(percemt, "%", 
                                                       stringr::str_sub(
                                                         strain_type, 1L, 2L)),
                                        y = -35),
                           nudge_y = ifelse(ve_compare$strain_type == "H1N1",
                                            1, 0)) +
  ggplot2::theme(legend.position = "bottom") +
  ggplot2::labs(x = "Season", y = "VE (%) and 95% CI",
                title = "Comparison of estimated VE to CDC repots in 18-49 age group") +
  ggplot2::theme(legend.position = "bottom",
                 title = ggplot2::element_text(size = 12),
                 legend.text = ggplot2::element_text(size = 11),
                 axis.text.y = ggplot2::element_text(size = 11),
                 legend.title = ggplot2::element_text(size = 12),
                 axis.title.y = ggplot2::element_text(size = 12),
                 axis.text.x = ggplot2::element_text(size = 11, 
                                                     color = "black")) 

cdc_seasons_circ_one

ggplot2::ggsave(here::here("results/figures/cdc_compare_circ_one.png"),
                cdc_seasons_circ_one,
                height = 8, width = 14)

## Schematic figures------------------------------------------------------------
### Titers - 2018 H1N1----------------------------------------------------------
h1n1_2018_pretiters <- 
  homologous_dat |> 
  dplyr::filter(strain_type=="H1N1",
                age >= 18,
                season == "2018") |> 

  ggplot2::ggplot() +
  ggplot2::aes(x = season, y = log(pretiter), 
               group = interaction(season, group),
               fill = group,
               alpha = group) +
  ggplot2::scale_x_continuous(breaks = 2013:2021,
                              labels = labelsh1n1) +
  ggplot2::geom_boxplot() +
  ggplot2::theme(legend.position = "right",
                 title = ggplot2::element_text(size = 12),
                 legend.text = ggplot2::element_text(size = 11),
                 axis.text.y = ggplot2::element_text(size = 11),
                 legend.title = ggplot2::element_text(size = 12),
                 axis.title.y = ggplot2::element_text(size = 12),
                 axis.text.x = ggplot2::element_text(size = 11, 
                                                     color = "black")) +
  ggplot2::labs(fill = "Group",
                x = NULL,
                y = "Log HAI Titer",
                alpha = "Group") +
  ggplot2::scale_fill_manual(values = group_colors) +
  ggplot2::ggtitle("Pre-Vaccine Log Titers") 

h1n1_2018_pretiters

ggplot2::ggsave(plot = h1n1_2018_pretiters, 
                filename = here::here("results/figures/h1n1_2018_pretiters.png"),
                width = 4, height = 2.25)

h1n1_2018_postiters <- 
  homologous_dat |> 
  dplyr::filter(strain_type=="H1N1",
                age >= 18,
                season == "2018") |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = season, y = log(postiter), 
               group = interaction(season, group),
               fill = group,
               alpha = group) +
  ggplot2::scale_x_continuous(breaks = 2013:2021,
                              labels = labelsh1n1) +
  ggplot2::geom_boxplot() +
  ggplot2::theme(legend.position = "left",
                 title = ggplot2::element_text(size = 12),
                 legend.text = ggplot2::element_text(size = 11),
                 axis.text.y = ggplot2::element_text(size = 11),
                 legend.title = ggplot2::element_text(size = 12),
                 axis.title.y = ggplot2::element_text(size = 12),
                 axis.text.x = ggplot2::element_text(size = 11, 
                                                     color = "black")) +
  ggplot2::labs(fill = "Group",
                x = NULL,
                y = "Log HAI Titer",
                alpha = "Group") +
  ggplot2::scale_fill_manual(values = group_colors) +
  ggplot2::labs(x = NULL) +
  ggplot2::ggtitle("Post-Vaccine Log Titers")


h1n1_2018_postiters

ggplot2::ggsave(plot = h1n1_2018_postiters, 
                filename = 
                  here::here("results/figures/h1n1_2018_postiters.png"),
                width = 4, height = 2.25)

### Protection - H1N1 2018------------------------------------------------------
h1n1_2018_preprotect <- 
  protection |> 
  dplyr::filter(strain_type=="H1N1",
                age >= 18,
                season == "2018") |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = season, y = pre_protect, 
               group = interaction(season, group),
               fill = group,
               alpha = group) +
  ggplot2::scale_x_continuous(breaks = 2013:2021,
                              labels = labelsh1n1) +
  ggplot2::geom_boxplot() +
  ggplot2::theme(legend.position = "right",
                 title = ggplot2::element_text(size = 12),
                 legend.text = ggplot2::element_text(size = 11),
                 axis.text.y = ggplot2::element_text(size = 11),
                 legend.title = ggplot2::element_text(size = 12),
                 axis.title.y = ggplot2::element_text(size = 12),
                 axis.text.x = ggplot2::element_text(size = 11, 
                                                     color = "black")) + 
  ggplot2::labs(fill = "Group",
                x = NULL,
                y = "Probability of Protection",
                alpha = "Group") +
  ggplot2::scale_fill_manual(values = group_colors) +
  ggplot2::ggtitle("Pre-Vaccine Protection")

h1n1_2018_preprotect

ggplot2::ggsave(plot = h1n1_2018_preprotect, 
                filename = 
                  here::here("results/figures/h1n1_2018_preprotect.png"),
                width = 4, height = 2.25)

h1n1_2018_postprotect <- 
  protection |> 
  dplyr::filter(strain_type=="H1N1",
                season == "2018") |> 
  ggplot2::ggplot() +
  ggplot2::aes(x = season, y = post_protect, 
               group = interaction(season, group),
               fill = group,
               alpha = group) +
  ggplot2::scale_x_continuous(breaks = 2013:2021,
                              labels = labelsh1n1) +
  ggplot2::geom_boxplot() +
  ggplot2::theme(legend.position = "left",
                 title = ggplot2::element_text(size = 12),
                 legend.text = ggplot2::element_text(size = 11),
                 axis.text.y = ggplot2::element_text(size = 11),
                 legend.title = ggplot2::element_text(size = 12),
                 axis.title.y = ggplot2::element_text(size = 12),
                 axis.text.x = ggplot2::element_text(size = 11, 
                                                     color = "black")) +
  ggplot2::labs(fill = "Group",
                x = NULL,
                y = "Probability of Protection",
                alpha = "Group") +
  ggplot2::scale_fill_manual(values = group_colors) +
  ggplot2::ggtitle("Post-Vaccine Protection")

h1n1_2018_postprotect

ggplot2::ggsave(plot = h1n1_2018_postprotect, 
                filename = 
                  here::here("results/figures/h1n1_2018_postprotect.png"),
                width = 4, height = 2.25)

### VE - 2018 H1N1--------------------------------------------------------------
h1n1_2018_ve <- 
  ggplot2::ggplot(data = dplyr::filter(ve, strain_type == "H1N1",
                                       data == "Main analysis",
                                       season == "2018")) +
  ggplot2::aes(x = season, 
               group = factor(group, labels =  c("YASD", "OAHD", "OASD"),
                              levels = c("YASD", "OAHD", "OASD")), 
               color = factor(group, labels =  c("YASD", "OAHD", "OASD"),
                              levels = c("YASD", "OAHD", "OASD")),
               fill =  factor(group, labels =  c("YASD", "OAHD", "OASD"),
                              levels = c("YASD", "OAHD", "OASD")),
               shape = factor(group, labels =  c("YASD", "OAHD", "OASD"),
                              levels = c("YASD", "OAHD", "OASD"))) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_ve_lwr*100, ymax = boot_ve_upr*100),
    linewidth = 0.75, alpha = 0.6, 
    position = ggplot2::position_dodge(width = 0.75)) +
  ggplot2::geom_point(size = 2, alpha = 0.85,
                      position = ggplot2::position_dodge(width = 0.75),
                      ggplot2::aes(y = boot_ve*100)) +
  ggplot2::labs(color = "Group",
                shape = "Group",
                fill = "Group",
                x = NULL,
                y = "VE and 95% CI",
                title = "Vaccine Efficacy (VE)") +
  ggplot2::scale_x_continuous(breaks = 2013:2021,
                              labels = labelsh1n1) +
  ggplot2::scale_color_manual(values = group_colors) +
  ggplot2::scale_fill_manual(values = group_colors) +
  ggplot2::scale_shape_manual(values = group_shapes) +
  ggplot2::theme(legend.position = "bottom",
                 title = ggplot2::element_text(size = 12),
                 legend.text = ggplot2::element_text(size = 11),
                 axis.text.y = ggplot2::element_text(size = 11),
                 legend.title = ggplot2::element_text(size = 12),
                 axis.title.y = ggplot2::element_text(size = 12),
                 axis.text.x = ggplot2::element_text(size = 11, 
                                                     color = "black")) +  
  ggplot2::ylim(c(0,100)) 

h1n1_2018_ve

ggplot2::ggsave(plot = h1n1_2018_ve, 
                filename = here::here("results/figures/h1n1_2018_ve.png"),
                width = 3.5, height = 3)

### ADVE - 2018 H1N1------------------------------------------------------------
h1n1_2018_adve <-
  ggplot2::ggplot(data = dplyr::filter(ive, 
                                     season == 2018,
                                     strain_type == "H1N1",
                                     type == "Absolute",
                                     data == "Main analysis")) +
  ggplot2::aes(x = season,
               group = factor(comparison, levels = 
                                c("OASD_YASD", "OAHD_OASD", "OAHD_YASD")),
               shape = factor(comparison, levels =
                                c("OASD_YASD", "OAHD_OASD", "OAHD_YASD"))) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = boot_est_lwr*100, ymax = boot_est_upr*100),
    linewidth = 0.75, alpha = 0.6, 
    position = ggplot2::position_dodge(width = 0.75)) +
  ggplot2::geom_point(size = 2, alpha = 0.85,
                      position = ggplot2::position_dodge(width = 0.75),
                      ggplot2::aes(y = boot_est*100)) +

  ggplot2::labs(shape = NULL,
                fill = NULL,
                x = NULL,
                y = "DVE and 95% CI",
                title = "Difference in Vaccine Efficacy (DVE)") +
  ggplot2::theme(legend.position = "bottom",
                 title = ggplot2::element_text(size = 12),
                 legend.text = ggplot2::element_text(size = 11),
                 axis.text.y = ggplot2::element_text(size = 11),
                 legend.title = ggplot2::element_text(size = 12),
                 axis.title.y = ggplot2::element_text(size = 12),
                 axis.text.x = ggplot2::element_text(size = 11, 
                                                     color = "black")) +
  ggplot2::scale_x_continuous(breaks = 2013:2021,
                              labels = labelsh1n1)
h1n1_2018_adve

ggplot2::ggsave(plot = h1n1_2018_adve, 
                filename = here::here("results/figures/h1n1_2018_adve.png"),
                width = 4.5, height = 3)

