################################################################################
##     "High molecular weight glutenin gene diversity in Aegilops tauschii    ##
##         demonstrates unique origin of superior wheat quality"              ##
##                                                                            ##
##            converting vcf to numeric date frame for figures                ##
##                                                                            ##
################################################################################
# doi: 
library(vcfR)
library(dplyr)

######################
# Read in data files #
######################

# Population information
popmap <- read.csv(file="Supplemental-Data-1_popmap.csv", header=TRUE)

# Variant calls 
vcf_in <- read.vcfR(file="Supplemental-Data-2_Glu-D1.vcf",  verbose = TRUE)

#######################
# Process vcf into df #
#######################
# Retrieve variant info for each call
SNP_pos<-data.frame(chr='chr1D', POS=as.integer(vcf_in@fix[,2]), QUAL=as.integer(vcf_in@fix[,6]), 
                    REF=as.character(vcf_in@fix[,4]), ALT=as.character(vcf_in@fix[,5]))

# Retrieve numeric genotype calls
genos_temp <- extract.gt(vcf_in, as.numeric = T)

#################################
# Fixing accession designations #
#################################
head(popmap)
# 1) select only the TAs in popmap
popmap_list <- as.vector(unlist(popmap$BW))
genos_temp2 <- as.data.frame(genos_temp)
genos_temp2 <- genos_temp2 %>% select(popmap_list)
colnames(genos_temp2)

# 2) rename accessions from BW to TA number
temp<-as.data.frame(colnames(genos_temp2))
colnames(temp)<-"BW"
temp2<-full_join(temp, popmap, by="BW")
colnames(genos_temp2)<-temp2$TA
colnames(genos_temp2)

# Bind variant info back to genotype calls
genos <- cbind(SNP_pos, genos_temp2)

# Keep tidy
rm(SNP_pos, genos_temp, genos_temp2, temp, temp2, popmap_list, vcf_in)

#################################
# Writing numeric calls to csv  #
#################################
write.csv(genos, file="Numeric_genotype_calls.csv", row.names = FALSE, quote = FALSE)

