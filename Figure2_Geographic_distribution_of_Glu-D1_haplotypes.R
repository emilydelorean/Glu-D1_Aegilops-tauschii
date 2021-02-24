################################################################################
##     "High molecular weight glutenin gene diversity in Aegilops tauschii    ##
##         demonstrates unique origin of superior wheat quality"              ##
##                                                                            ##
##        Figure 2: Geographic distribution of Glu-D1 haplotypes              ##
##                                                                            ##
################################################################################
# doi: 
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

##  Create blank world map object
worldmap <- ne_countries(scale = 'medium', type = 'map_units', returnclass = 'sf')

## Read in the population passport information file, provided as Supplemental Data 1
popmap <- read.csv(file="Supplemental-Data-1_popmap.csv", header=T)


#########################################################################
##      Figure 2, panel A, Geographic distribution by major clade      ##
#########################################################################

majclade <- 
  ggplot() + geom_sf(data = worldmap) +  theme_bw() + theme(legend.position="right") +
  coord_sf(xlim = c(36,76), ylim = c(30, 42.5), expand = FALSE) +
  geom_point(data=popmap, mapping = aes(x=Longitude, y=Latitude,color=Major_clade, shape = Lineage),cex=2,stroke=1.5, alpha=0.5) +
  scale_color_manual(values = c("#1B91A2","#882255", "#DCAC6E" )) + 
  scale_shape_manual(values=c(15, 16, 17, 0)) +
    annotate(geom="text", x=39, y=39, label="TURKEY", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 38.5, y = 35, label = "SYRIA", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 43, y = 33, label = "IRAQ", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 54, y = 33, label = "IRAN", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 65, y = 33.5, label = "AFGHANISTAN", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 71.5, y = 33.5, label = "PAKISTAN", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 59, y = 39.5, label = "TURKMENISTAN", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 65, y = 40, label = "UZBEKISTAN", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 75, y = 42, label = "KYRGYZSTAN", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 72, y = 38.8, label = "TAJIKISTAN", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 75, y = 39.5, label = "CHINA", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 48, y = 40.1, label = "AZERBAIJAN", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 68.5, y = 42, label = "KAZAKHSTAN", color="grey50", size=3,fontface=2)


######################################################################
##      Figure 2, panel B, Geographic distribution by subclade      ##
######################################################################

subclade <- 
  ggplot() + geom_sf(data = worldmap) +  theme_bw() + theme(legend.position="right") +
  coord_sf(xlim = c(36,76), ylim = c(30, 42.5), expand = FALSE) +
  geom_point(data=popmap, mapping = aes(x=Longitude, y=Latitude,color=Subclade, shape = Lineage),cex=2,stroke=1.5, alpha=0.5) +
  scale_color_viridis(alpha = 0.7, begin = 0, end = 1, direction = 1, discrete = FALSE, option = "D")  +
  scale_shape_manual(values=c(15, 16, 17, 0))+
  scale_colour_steps2() +
    annotate(geom="text", x=39, y=39, label="TURKEY", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 38.5, y = 35, label = "SYRIA", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 43, y = 33, label = "IRAQ", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 54, y = 33, label = "IRAN", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 65, y = 33.5, label = "AFGHANISTAN", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 71.5, y = 33.5, label = "PAKISTAN", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 59, y = 39.5, label = "TURKMENISTAN", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 65, y = 40, label = "UZBEKISTAN", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 75, y = 42, label = "KYRGYZSTAN", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 72, y = 38.8, label = "TAJIKISTAN", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 75, y = 39.5, label = "CHINA", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 48, y = 40.1, label = "AZERBAIJAN", color="grey50", size=3,fontface=2) +
    annotate(geom = "text", x = 68.5, y = 42, label = "KAZAKHSTAN", color="grey50", size=3,fontface=2)


######################################################################
##          Plots panels together and write to a pdf                ##
######################################################################

pdf('Figure2.pdf', width = 11, height = 7)
ggarrange(majclade, subclade, ncol = 1, nrow = 2, labels="AUTO")
dev.off()

