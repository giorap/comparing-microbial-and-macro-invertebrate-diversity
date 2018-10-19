###################################################
################### micro-macro ###################
###################################################
################# DATA PROCESSING #################
###################################################
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
