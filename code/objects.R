################################################################################
### TITLE: Objects Script                                                    ###
### AUTHOR: Savannah Hammerton                                               ###
### PURPOSE: This script will store all the various functions to be used in  ###
###          the other coding scripts of this project. Example objects are   ###
###          color and axis vectors for plots, vectors of functions, etc.    ###
###          This script does not need to be run within its own file, but    ###
###          will instead be run at the start of every script using          ###
###          `source(objects.R)`                                             ###
################################################################################

# Package loading 

library(NatParksPalettes)
library(paletteer)


# Color scales------------------------------------------------------------------
# These color scales will allow me to maintain a consistent color scheme in all
# my plots, particularly keeping group colors the same 
NatParksPalettes::natparks.pals("DeathValley")
NatParksPalettes::natparks.pals("Cuyahoga")

timing_colors <- c(PreVaccine = 
                     NatParksPalettes::natparks.pals("DeathValley")[[5]],
                   PostVaccine = 
                     NatParksPalettes::natparks.pals("DeathValley")[[6]])
dose_colors <- c(SD = 
                   NatParksPalettes::natparks.pals("Cuyahoga")[[1]],
                 HD =
                   NatParksPalettes::natparks.pals("Cuyahoga")[[4]],
                 Intradermal = 
                   NatParksPalettes::natparks.pals("Cuyahoga")[[2]])
NatParksPalettes::natparks.pals("Arches", n = 16)

count_colors <- c(NatParksPalettes::natparks.pals("Arches", n = 16)[[1]],
                  NatParksPalettes::natparks.pals("Arches", n = 16)[[2]],
                  NatParksPalettes::natparks.pals("Arches", n = 16)[[3]],
                  NatParksPalettes::natparks.pals("Arches", n = 16)[[4]],
                  NatParksPalettes::natparks.pals("Arches", n = 16)[[5]],
                  NatParksPalettes::natparks.pals("Arches", n = 16)[[6]],
                  NatParksPalettes::natparks.pals("Arches", n = 16)[[7]])
group_colors <- c(OAHD = "#5F7759",
                  OASD = "#0072B2",
                  YASD = "#D55E00"
                    )


group_shapes <- c(
  OAHD = 16,
  OASD = 17,
  YASD = 15
)

comp_shapes <- c(
  OASD_YASD = 16,
  OAHD_YASD = 17,
  OAHD_OASD = 15
)

comp_colors <- c(OASD_YASD = "#5F7759",
                  OAHD_YASD = "#0072B2",
                  OAHD_OASD = "#D55E00"
                  
)


type_colors <- c(
  Absolute = paletteer::paletteer_d("nationalparkcolors::MtRainier")[[1]], 
  Relative = paletteer::paletteer_d("nationalparkcolors::Arches")[[3]]
)

type_shapes = c(
  Absolute = 16,
  Relative = 18
)


# Axis labels 
h1n1_x_labels <- 
  homologous_dat |> 
  dplyr::filter(strain_type == "H1N1") |> 
  dplyr::group_by(season) |> 
  dplyr::summarise(strain = strains_short[1]) |> 
  data.frame()

breaks.h1n1 <- h1n1_x_labels[,1] |> as.vector()
breaks.h1n1.strain <- h1n1_x_labels[,2] |> as.vector()
labelsh1n1 <- paste0(breaks.h1n1, "\n", breaks.h1n1.strain)

# Axis labels 
h3n2_x_labels <- 
  homologous_dat |> 
  dplyr::filter(strain_type == "H3N2") |> 
  dplyr::group_by(season) |> 
  dplyr::summarise(strain = strains_short[1]) |> 
  data.frame()

breaks.h3n2 <- h3n2_x_labels[,1] |> as.vector()
breaks.h3n2.strain <- h3n2_x_labels[,2] |> as.vector()
labelsh3n2 <- paste0(breaks.h3n2, "\n", breaks.h3n2.strain)


labelsboth <- paste0(labelsh1n1, "\n", breaks.h3n2.strain)

# Age groups (based on CDC dataset)
age_groups <- c("9–17", "18–49", "50–64", "≥65")

# Colors for CDC VE comparison plots 
ve_comp_colors <- c(
  CDC = "black",
  H1N1 = "royalblue",
  H3N2 = "red"
)

# Shapes for CDC VE comparison plots indicating which type was matched

