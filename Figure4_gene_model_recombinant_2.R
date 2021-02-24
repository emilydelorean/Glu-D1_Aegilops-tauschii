################################################################################
##     "High molecular weight glutenin gene diversity in Aegilops tauschii    ##
##         demonstrates unique origin of superior wheat quality"              ##
##                                                                            ##
##                     Figure 4: Glu-D1 recombinants                          ##
##                                                                            ##
################################################################################
# doi: 

# geographic distribution code derived from Dr. Paula Silva's
# because she is my best friend and her figures are always beautiful
#  - https://github.com/SilvaPaula/scripts/blob/main/script.map.R

# colorblind friendliness checked with davidmathlogic.com/colorblind/

###################################################
##                Gene model figure              ##
##      for Glu-D1 recombinants (R1 and R2)      ##
###################################################

#  - Glu-D1x length is 2568 bp (419306988 bp to 419309556 bp)
#  - Glu-D1y length is 1980 bp (419364015 bp to 419365995 bp)
#  - Distance between x and y subunit is 54,549 bp (419309556 to 419364015)

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
worldmap <- ne_countries(scale = 'medium', type = 'map_units',returnclass = 'sf')

## Read in the population passport information file, provided as Supplemental Data 1
popmap <- read.csv(file="Supplemental-Data-1_popmap.csv", header=TRUE)

## Read in numeric genotype calls file
## See vcf_to_numeric.R code for generating this file from the vcf
genos <- read.csv(file="Numeric_genotype_calls.csv", header=TRUE)


#########################################################################
##                  Prepare data for plotting                          ##
#########################################################################
popmapR1    <- popmap %>% filter(Recom_1 %in% c('x_subunit', 'y_subunit', 'R1'))
popmapR1_R1 <- popmap %>% filter(Recom_1 %in% c('R1'))
popmapR2    <- popmap %>% filter(Recom_2 %in% c('x_subunit', 'x_subunit', 'y_subunit', 'R2'))
popmapR2_R2 <- popmap %>% filter(recom_2 %in% c('R2'))


######################################
##     Create the Gene models      ###
######################################
#  - Glu-D1x CDS length is 2568 bp (419306988 bp to 419309556 bp)
#  - Glu-D1y CDS length is 1980 bp (419364015 bp to 419365995 bp)
#  - Distance between x and y subunit is 54,549 bp (419309556 to 419364015)
#  - Haplotypes included the 2.5 kb flanking regions around each gene

# The coordinates used for the gene models are (on the x axis) 
# Glu-D1x CDS is 1 to 2569
#  - flanking regions are -2501 to 0 and 2570 to 5070
# Glu-D1y CDS is 12000 to 13980
#  - flanking regions are 9500 to 12000 and 13980 to 16481


buf=0.8     # y axis space between each gene model
bot=2       # y axis bottom of the bottom most gene model     
top=2.5     # y axis top of the bottom most gene model 
mid=top+(bot-top)/2  # y axis midpoint of bottom most gene model
s2=0.05  # x axis width of tick marks 
ht = 0.1 # y axis height of tick marks 


# Create the rectangles for the x and y subunit genes
d1=data.frame(x1=c(1,12000, 1,12000, 1,12000), 
              x2=c(2569,12000+1980, 2569,12000+1980, 2569,12000+1980), 
              y1=c(bot,bot, bot+buf,bot+buf, bot+2*buf,bot+2*buf),
              y2=c(top,top, top+buf,top+buf, top+2*buf,top+2*buf))
d1
# Create the rectangle for x flanking region
d2=data.frame(x1=c(-2501, -2501, -2501), x2=c(2570+2500, 2570+2500, 2570+2500),
              y1=c(bot, bot+buf, bot+2*buf), y2=c(top, top+buf, top+2*buf), )
d2
# Create the thin rectangle to connect the x and y subunit
d3=data.frame(x1=c(2570+2500, 2570+2500, 2570+2500), 
              x2=c(12000-2500, 12000-2500, 12000-2500), 
              y1=c(mid-.01, mid-.01+buf, mid-.01+2*buf), 
              y2=c(mid+0.01, mid+0.01+buf, mid+0.01+2*buf))

# Create the rectangle for y flanking region
d4=data.frame(x1=c(12000-2500, 12000-2500, 12000-2500), 
              x2=c(12000+1980+2501, 12000+1980+2501, 12000+1980+2501), 
              y1=c(top, top+buf, top+2*buf), 
              y2=c(bot, bot+buf, bot+2*buf))


# Hatch marks 
hatch<-data.frame(x=c(7300,7500, 7300,7500, 7300,7500), 
                  xend=c(7400,7600, 7400,7600, 7400,7600),
                  y=c(2.15,2.15, mid+1*buf-0.1, mid+1*buf-0.1, mid+2*buf-0.1, mid+2*buf-0.1), 
                  yend=c(2.35,2.35, mid+1*buf+0.1, mid+1*buf+0.1, mid+2*buf+0.1, mid+2*buf+0.1))

# Tick marks for labels 
# at -2000, -1000, x end, x start, +1000, +2000 
# at -2000, -1000, y end, y start, +1000, +2000 

tick <- data.frame(x=c(-2000, -1000, 0, 2568, 2568+1000, 2568+2000, 
                       12000-2000, 12000-1000, 12000, 12000+1980, 12000+1980+1000, 12000+1980+2000 ),
                   xend=c(-2000-s2, -1000-s2, 0+s2, 2568+s2, 2568+1000+s2, 2568+2000+s2,
                          12000-2000+s2, 12000-1000-s2, 12000+s2, 12000+1980+s2,12000+1980+1000+s2, 12000+1980+2000+s2),
                   label=c("+2 kb", "+1 kb", "end", "start", "-1 kb", "-2 kb",
                           "+2 kb", "+1 kb", "end", "start", "-1 kb", "-2 kb")
                   )

tick <- tick %>% mutate( y = top+2*buf ) %>% mutate(yend = y+ht) %>% mutate(lab_ht = yend+s2)

###############################
##     Mutate positions      ##
## to match the gene models  ##
###############################

dfx <- genos %>% filter(POS >= 419306988-2501 & POS <= 419309556+2501) %>% mutate(POS = POS - 419306988)
dfy <- genos %>% filter(POS >= 419364015-2501 & POS <= 419365995+2501 ) %>%  mutate(POS = POS - 419364015 + 12000)
SNP <- rbind(dfx, dfy)

#keep tidy
rm(dfx, dfy)

# x axis width of the SNP bar
size <- 0.1

############################
##   Variant Positions    ##
##        for R1          ##
############################

TA10081 <- SNP %>% select(TA10081, POS) %>% mutate(POS = POS*TA10081) %>% mutate(POS2 = POS+size) %>% mutate(ymin = bot, ymax = top)             %>% filter(POS != 0)
TA1668  <- SNP %>% select(TA1668, POS)  %>% mutate(POS = POS*TA1668)  %>% mutate(POS2 = POS+size) %>% mutate(ymin = bot+1*buf, ymax = top+1*buf) %>% filter(POS != 0)
TA1678  <- SNP %>% select(TA1678, POS)  %>% mutate(POS = POS*TA1678)  %>% mutate(POS2 = POS+size) %>% mutate(ymin = bot+2*buf, ymax = top+2*buf) %>% filter(POS != 0)


############################
##   Variant Positions    ##
##        for R2          ##
############################

TA2466    <- SNP %>% select(TA2466, POS)   %>% mutate(POS = POS*TA2466)    %>% mutate(POS2 = POS+size) %>% mutate(ymin = bot, ymax = top)             %>% filter(POS != 0)
TA2576    <- SNP %>% select(TA2576, POS)   %>% mutate(POS = POS*TA2576)    %>% mutate(POS2 = POS+size) %>% mutate(ymin = bot+1*buf, ymax = top+1*buf) %>% filter(POS != 0)
BW_01028  <- SNP %>% select(BW_01028, POS) %>% mutate(POS = POS*BW_01028)  %>% mutate(POS2 = POS+size) %>% mutate(ymin = bot+2*buf, ymax = top+2*buf) %>% filter(POS != 0)


######################################
##     Plot Recombinant 1           ##
##        gene models               ##
######################################

R1_mod <- 
  ggplot() + 
  ylim(1.5, 5) + xlim(-5000, 20000) +
  theme_void() +
  geom_rect(data=d1, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color="grey55", fill="gray55", alpha=0.5) +    # make the gene rectangles
  geom_rect(data=d2, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="grey60", fill="gray88", alpha=0.5) +            # make the x flanking line 
  geom_rect(data=d3, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="grey60", fill="gray88", alpha=0.5) +            # make the connecting line between genes
  geom_rect(data=d4, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="grey60", fill="gray88", alpha=0.5) +            # make the y flanking line
  geom_rect(data=TA1678, mapping=aes(xmin=POS, xmax=POS2, ymin=ymin, ymax=ymax), color="purple4", alpha=0.5) +               # make SNP lines
  geom_rect(data=TA10081, mapping=aes(xmin=POS, xmax=POS2, ymin=ymin, ymax=ymax), color="purple4", alpha=0.5) +              # make SNP lines
  geom_rect(data=TA1668, mapping=aes(xmin=POS, xmax=POS2, ymin=ymin, ymax=ymax), color="purple4", alpha=0.5)  +              # make SNP lines
  annotate("segment", x=hatch$x, y=hatch$y, xend=hatch$xend, yend=hatch$yend, col="gray50") +                                 # add hatch marks
  geom_rect(data=tick, mapping=aes(xmin=x, xmax=xend, ymin=y, ymax=yend), color="grey20", alpha=0.5) +                  # make position ticks
  annotate("segment", x=hatch$x, y=hatch$y, xend=hatch$xend, yend=hatch$yend, col="gray50") +                                 # add hatch marks
  annotate(geom="text", x=tick$x, y=tick$lab_ht, label=tick$label, color="grey50", size=2.75,fontface=2) +
  annotate(geom="text", x=c(1300, 13000), y=c(4.5, 4.5), label=c("x subunit", "y subunit"), color="grey20", size=4, fontface=2) +
  annotate(geom="text", x=c(-3800, 7500, 18000), y=c(4.8, 4.8, 4.8), label=c("Lineage", "Haplotypes", "SDS-PAGE"), color="grey20", size=4, fontface=2) +
  annotate(geom="text", x=c(-3800, 7500, 18000, -3800, 7500, 18000, -3800, 7500, 18000), y=c(3.87, 4.05, 3.87, 3.1, 3.3, 3.1, 2.25, 2.5, 2.25), 
           label=c("L1 & L2", "x9a + y9b", "2 + 12" , "L2", "x9a + y5a", "2 + 10.2", "L2", "x5b + y5a", "2 + 10.2"  ), color="grey20", size=4, fontface=1) 



######################################
##     Plot Recombinant 2           ##
##        gene models               ##
######################################

R2_mod <- 
  ggplot() + 
  ylim(1.5, 5) + xlim(-5000, 20000) +
  theme_void() +
  geom_rect(data=d1, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color="grey55", fill="gray55", alpha=0.5) +    # make the gene rectangles
  geom_rect(data=d2, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="grey60", fill="gray88", alpha=0.5) +            # make the x flanking line 
  geom_rect(data=d3, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="grey60", fill="gray88", alpha=0.5) +            # make the connecting line between genes
  geom_rect(data=d4, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="grey60", fill="gray88", alpha=0.5) +            # make the y flanking line
  geom_rect(data=BW_01028, mapping=aes(xmin=POS, xmax=POS2, ymin=ymin, ymax=ymax), color="purple4", alpha=0.5)  +              # make SNP lines
  geom_rect(data=TA2576, mapping=aes(xmin=POS, xmax=POS2, ymin=ymin, ymax=ymax), color="purple4", alpha=0.5) +              # make SNP lines
  geom_rect(data=TA2466, mapping=aes(xmin=POS, xmax=POS2, ymin=ymin, ymax=ymax), color="purple4", alpha=0.5) +               # make SNP lines
  geom_rect(data=tick, mapping=aes(xmin=x, xmax=xend, ymin=y, ymax=yend), color="grey20", alpha=0.5) +                  # make position ticks
  annotate("segment", x=hatch$x, y=hatch$y, xend=hatch$xend, yend=hatch$yend, col="gray50") +                                 # add hatch marks
  annotate(geom="text", x=tick$x, y=tick$lab_ht, label=tick$label, color="grey50", size=2.75,fontface=2) +
  annotate(geom="text", x=c(1300, 13000), y=c(4.5, 4.5), label=c("x subunit", "y subunit"), color="grey20", size=4, fontface=2) +
  annotate(geom="text", x=c(-3800, 7500, 18000), y=c(4.8, 4.8, 4.8), label=c("Lineage", "Haplotypes", "SDS-PAGE"), color="grey20", size=4, fontface=2) +
  annotate(geom="text", x=c(-3800, 7500, 18000, -3800, 7500, 18000, -3800, 7500, 18000), y=c(3.87, 4.05, 3.87, 3.1, 3.3, 3.1, 2.25, 2.5, 2.25), 
           label=c("L3", "x7a + y7a", "" , "L3", "xR2a + yR2a", "", "L2", "x15a + y15a", "2.1 + 12.3"  ), color="grey20", size=4, fontface=1) 
  



##############################################################
##        Geographic distibution plots for R1 and R2        ##
##############################################################
R1_geo <- 
  ggplot() + geom_sf(data = worldmap, fill="gray95") +  theme_bw() + theme(legend.position="none") +
  coord_sf(xlim = c(42,55), ylim = c(36, 42.5), expand = FALSE) +
  geom_point(data=popmap, mapping = aes(x=Longitude, y=Latitude, col=Recom_1, shape=Lineage),cex=1.5,stroke=1.5, alpha=0.7) +
  geom_point(data=popmapR1, mapping = aes(x=Longitude, y=Latitude, col=Recom_1, shape=Lineage),cex=1.5,stroke=1.5, alpha=0.7) +
  geom_point(data=popmapR1_R1, mapping = aes(x=Longitude, y=Latitude, col=Recom_1, shape=Lineage),cex=1.5,stroke=1.5, alpha=0.7) +
  scale_color_manual(values = c("gray65", "lawngreen", "turquoise4", "darkgoldenrod2")) +
  scale_shape_manual(values=c(15, 16, 17, 0)) 

R2_geo <- 
  ggplot() + geom_sf(data = worldmap, fill="gray95") +  theme_bw() + theme(legend.position="none") +
  coord_sf(xlim = c(42,55), ylim = c(36, 42.5), expand = FALSE) +
  geom_point(data=popmap, mapping = aes(x=Longitude, y=Latitude, col=Recom_2,shape=Lineage),cex=2,stroke=1.5, alpha=0.7) +
  geom_point(data=popmapR2, mapping = aes(x=Longitude, y=Latitude, col=Recom_2, shape=Lineage),cex=1.5,stroke=1.5, alpha=0.7) +
  scale_color_manual(values = c("gray65", "lawngreen", "turquoise4", "darkgoldenrod2")) +
  scale_shape_manual(values=c(15, 16, 17, 0)) 


######################################################################
##      Figure 4, plot all panels together and export to pdf        ##
######################################################################

pdf('Figure4.pdf', width = 15, height = 7)
  ggarrange(R1_mod, R1_geo, R2_mod, R2_geo, ncol = 2, nrow = 2, widths = c(1.75,1), labels=c("A", "C", "B", "D"))
dev.off()
