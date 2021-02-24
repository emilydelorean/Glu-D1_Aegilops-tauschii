################################################################################
##     "High molecular weight glutenin gene diversity in Aegilops tauschii    ##
##         demonstrates unique origin of superior wheat quality"              ##
##                                                                            ##
##        Figure 1: Geographic distribution of Glu-D1 haplotypes              ##
##                                                                            ##
################################################################################
# doi: 
# phylogenetic tree code derived from Dr. Paula Silva's with permission
#  - https://github.com/SilvaPaula/scripts/blob/main/script.phylogeny.R
# colorblind friendliness checked with davidmathlogic.com/colorblind/

library(ape)
library(phyclust)
library(tidytree)
library(dplyr)
library(ggplot2)

## Read in the population passport information file, provided as Supplemental Data 1
popmap <- read.csv(file="Supplemental-Data-1_popmap.csv", header=TRUE)

## Read in numeric genotype calls file
## See vcf_to_numeric.R code for generating this file from the vcf
genos <- read.csv(file="Numeric_genotype_calls.csv", header=TRUE, check.names=FALSE)

## Select only the nonredudant lines
popmap <- popmap %>% filter(Dup_Keep == 'y')
popmap_list <- as.vector(unlist(popmap$TA))
genos_temp <- genos %>% select(popmap_list)


###################################
## Hierarchical Cluster Analysis 
###################################
# compute genetic distance
distMat <- dist(t(genos_temp))

# cluster accessions and convert to phylo object for ape
hc2 <- as.phylo(hclust(distMat)) #, method = 'average'

# cluster coloring
edgecols <- cbind('TA'=NA, 1:nrow(hc2$edge), color='#808080') # create data frame
edgecols[,1] <- hc2$tip.label[hc2$edge[,2]] # get labels
edgecols <- as.matrix(merge(edgecols, popmap, by = 'TA', all.x = T)) # get samples info
edgecols <- edgecols[order(as.numeric(edgecols[,2])), ] # get samples in original order

# Painting the edges
#edgecols[,3][edgecols[,6] == 'I'] = "black"
#edgecols[,3][edgecols[,5] == 'II'] = "blue"
#edgecols[,3][edgecols[,5] == 'III'] = "green"

#######################
## Leaf painting ##
#######################
leaves <- as.matrix(merge(hc2$tip.label, edgecols, by = 1, sort = F))

# Painting the tips by lineage
leaves[,3][leaves[,5] == 'L1'] = "#696969"
leaves[,3][leaves[,5] == 'L2'] = "#20B2AA"
leaves[,3][leaves[,5] == 'L3'] = "#A1D209"
leaves[,3][leaves[,5] == 'Wheat'] = "#6B0633"

# Plot phylogenetic tree
plot.phylo(hc2, type = 'p', show.tip.label = T, label.offset = 0.2, cex = 1,
           edge.color = edgecols[,3], edge.width = 2, tip.color = leaves[,3]) 








###################################################
##                Gene model figure              ##
##                                               ##
###################################################

#  - Glu-D1x length is 2568 bp (419306988 bp to 419309556 bp)
#  - Glu-D1y length is 1980 bp (419364015 bp to 419365995 bp)
#  - Distance between x and y subunit is 54,549 bp (419309556 to 419364015)

#   colorblind checked with davidmathlogic.com/colorblind/
#   palettes from color-hex.com


############################
##      Gene models      ###
############################
bot=2       # 
top=2.5
mid=top+(bot-top)/2
ht=0.03

# Create the rectangle for x flanking region
d2=data.frame(x1=c(-2501), 
              x2=c(2570+2500), 
              y1=c(top), 
              y2=c(bot))

# Create the rectangles for the x and y subunit genes
d1=data.frame(x1=c(1,12000), 
              x2=c(2569,12000+1980), 
              y1=c(top+ht,top+ht), 
              y2=c(bot-ht,bot-ht))

# Create the rectangle to connect the x and y subunits
d3=data.frame(x1=c(2570+2500), 
              x2=c(12000-2500), 
              y1=c(mid-.01), 
              y2=c(mid+0.01))

# Create the rectangle for y flanking region
d4=data.frame(x1=c(12000-2500), 
              x2=c(12000+1980+2501), 
              y1=c(top), 
              y2=c(bot))


# Hatch marks 
hatch<-data.frame(x=c(7300,7500), 
                  y=c(2.15,2.15), 
                  xend=c(7400,7600), 
                  yend=c(2.35,2.35))

# Tick marks
# at -2000, -1000, x end, x start, +1000, +2000 
# at -2000, -1000, y end, y start, +1000, +2000 

s2=0.05
ht = 0.1

tick <- data.frame(x=c(-2000, -1000, 0, 2568, 2568+1000, 2568+2000, 
                       12000-2000, 12000-1000, 12000, 12000+1980, 12000+1980+1000, 12000+1980+2000 ),
                   xend=c(-2000-s2, -1000-s2, 0+s2, 2568+s2, 2568+1000+s2, 2568+2000+s2,
                          12000-2000+s2, 12000-1000-s2, 12000+s2, 12000+1980+s2,12000+1980+1000+s2, 12000+1980+2000+s2),
                   y=c(top,top, top, top, top, top,
                       top,top, top, top, top, top), 
                   yend=c(top+ht,top+ht, top+ht, top+ht, top+ht, top+ht,
                          top+ht,top+ht, top+ht, top+ht, top+ht, top+ht),
                   label=c("+2 kb", "+1 kb", "end", "start", "-1 kb", "-2 kb",
                           "+2 kb", "+1 kb", "end", "start", "-1 kb", "-2 kb"))



############################
##   Variant Positions    ##
##                        ##
############################
genos <- read.csv(file="Numeric_genotype_calls.csv", header=TRUE)

dfx <- genos %>% select(POS) %>% filter(POS >= 419306988-2501 & POS <= 419309556+2501) %>% mutate(POS = POS - 419306988) 
dfy <- genos %>% select(POS) %>% filter(POS >= 419364015-2501 & POS <= 419365995+2501 ) %>%  mutate(POS = POS - 419364015 + 12000)
SNP <- rbind(dfx, dfy)

rm(dfx, dfy)

size <- 0.05
SNP <- SNP %>% mutate(POS2 = POS+size) %>% mutate(ymin = bot, ymax = top) %>% filter(POS != 0)

# Create the vectors to connect the SNPs to the haplotype table

d5=data.frame(x1=SNP$POS, x2=seq(-2500, 12000+1980+2501, by=(12000+1980+2501+2500)/nrow(SNP) + 0.001), y1=1.5, y2=2)
d5

########################################################
##          Gene Model x and y subunits               ##
########################################################

ggplot() + 
  ylim(1.4, 2.7) + xlim(-2510, 12000+1980+2501) +
  theme_void() +
  geom_rect(data=d1, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="grey55", fill="gray55", alpha=0.5) +            # make the gene rectangles
  geom_rect(data=d2, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="grey60", fill="gray88", alpha=0.5) +            # make the x flanking line 
  geom_rect(data=d3, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="grey60", fill="gray88", alpha=0.5) +            # make the connecting line between genes
  geom_rect(data=d4, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="grey60", fill="gray88", alpha=0.5) +            # make the y flanking line
  geom_rect(data=SNP, mapping=aes(xmin=POS, xmax=POS2, ymin=ymin, ymax=ymax), color="purple4", alpha=0.5) +                  # make SNP lines
  geom_rect(data=tick, mapping=aes(xmin=x, xmax=xend, ymin=y, ymax=yend), color="grey20", alpha=0.5) +                  # make position ticks
  annotate("segment", x=hatch$x, y=hatch$y, xend=hatch$xend, yend=hatch$yend, col="gray50")   +                               # add hatch marks
  annotate("segment", x=d5$x1, y=d5$y2, xend=d5$x2, yend=d5$y1, color="purple4", alpha=0.1, size=1)  +
  annotate(geom="text", x=c(1300, 13000), y=2.66, label=c("x subunit", "y subunit"), color="grey20", size=4, fontface=2) +
  annotate(geom="text", x=tick$x, y=2.62, label=tick$label, color="grey50", size=2.75,fontface=2) 
  