################################################################################
##     "High molecular weight glutenin gene diversity in Aegilops tauschii    ##
##         demonstrates unique origin of superior wheat quality"              ##
##                                                                            ##
##          Figure 3: Cryptic haplotypes within SDS-PAGE alleles              ##
##                                                                            ##
################################################################################
# doi: 

# geographic distribution code derived from Dr. Paula Silva's
# because she is my best friend and her figures are always beautiful
#  - https://github.com/SilvaPaula/scripts/blob/main/script.map.R

# colorblind friendliness checked with davidmathlogic.com/colorblind/

############################
##    Work space setup    ##
############################

library(ggmap)
library(maps)
library(mapdata)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(viridis)
library(ggpubr)
library(ggplot2)
library(dplyr)

##  Create blank world map object
worldmap <- ne_countries(scale = 'medium', type = 'map_units', returnclass = 'sf')

## Read in the population passport information file, provided as Supplemental Data 1
popmap <- read.csv(file="Supplemental-Data-1_popmap.csv", header=T)

#########################################################################
##                  Prepare data for plotting                          ##
#########################################################################
##  Create new column of SDS-PAGE x and y subunits together
popmap$SDSPAGE <- paste(popmap$D1x_SDS_PAGE, "+", popmap$D1y_SDS_PAGE, sep="")

## Ae. tauschii accessions without SDS-PAGE info will have only a "+"
## set those to NA
popmap$SDSPAGE <- recode_factor(popmap$SDSPAGE, "+" = NA_character_)

## Select accessions with the 2+12 SDS-PAGE allele
popmap_cry <- popmap %>% filter(SDSPAGE == "2x+12y")

## Change the factor level name from 2x+12y to 2+12 for cleaner plot labels
popmap_cry$SDSPAGE <- recode_factor(popmap_cry$SDSPAGE, "2x+12y" = "2+12")

#########################################################################
##      Figure 3, Geographic distribution of cryptic haplotpyes        ##
#########################################################################
cyp_geo <- 
  ggplot() + geom_sf(data = worldmap, fill="gray95") +  
  theme_bw() + theme(legend.position="right") +
  coord_sf(xlim = c(44,58), ylim = c(32.75, 43), expand = FALSE) +
  # plot first the accessions that aren't of interest
  geom_point(data=popmap, mapping = aes(x=Longitude, y=Latitude, col=Subclade, 
                                        shape=Lineage),cex=2, stroke=1.5, 
                                        alpha=0.5) +
  # plot accessions of interest on top
  geom_point(data=popmap_cry, mapping = aes(x=Longitude, y=Latitude, 
                                            col=Subclade, shape=Lineage), cex=2,
                                            stroke=1.5, alpha=0.7) +
  scale_color_viridis(alpha = 0.7, begin = 0, end = 1, direction = 1, 
                      discrete = FALSE, option = "D")  +
  scale_shape_manual(values=c(15, 16, 17, 0)) +
  # Plot SDS-PAGE labels at either the top or bottom of map, depending on 
  # coordinates of the corresponding accession
  geom_label_repel(data = subset(popmap_cry, Latitude > 38.5),
                   aes(x=Longitude, y=Latitude, label = SDSPAGE),
                   nudge_y = 44 - subset(popmap_cry, Latitude > 38.5)$Latitude,
                   box.padding = 0.35, point.padding = 0.5, force = 100,
                   segment.color = 'black', alpha=1, direction = "x") +
  geom_label_repel(data = subset(popmap_cry, Latitude < 38.5),
                   aes(x=Longitude, y=Latitude, label = SDSPAGE),
                   nudge_y = 32.75 - subset(popmap_cry, Latitude < 38.5)$Latitude,
                   box.padding = 0.35, point.padding = 0.5,
                   force = 100, segment.color = 'black',
                   alpha=1, direction = "x") 

######################################################################
######                 Write plot to a pdf                    ########
######################################################################

pdf('Figure3.pdf', width = 7, height = 5)
  cyp_geo
dev.off()
