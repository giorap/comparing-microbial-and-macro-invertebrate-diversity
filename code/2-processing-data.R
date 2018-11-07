###################################################
################### micro-macro ###################
###################################################
################# DATA PROCESSING #################
###################################################
###############################
##### -- Load packages -- #####
###############################
my_packages <- c("plyr", "dplyr", "magrittr", "XLConnect", "scales", ### data manipulation
                 "SpadeR", "vegan", "SYNCSA", "iNEXT", "MuMIn", ### algorithms
                 "RColorBrewer", "ggplot2", "ggpubr", "GGally" ### visualizations
)
#### Check for packages that are not already installed and isntall them
new_packages <- my_packages[!(my_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
#### Load all packages
lapply(my_packages, require, character.only=T)

#####################################
##### -- Load microbial data -- #####
#####################################
micro_abundance_byLake <- readRDS("data/microbe_abundance_byLake.rds")

#########################################
##### -- Load macro-inverts data -- #####
#########################################
inverts_abundance_byLake <- readRDS("data/inverts_abundance_byLake.rds") 

#####################################################################
##### -- Identify lakes with complete biological information -- #####
#####################################################################
shared_lake_codes <- intersect(unique(micro_abundance_byLake$lake_code), unique(inverts_abundance_byLake$lake_code))

#########################################
##### -- Load environmental data -- #####
#########################################
environment_byLake <- readRDS("data/environment_byLake.rds")
