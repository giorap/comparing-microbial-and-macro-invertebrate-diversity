#########################################
############## micro-macro ##############
#########################################
############# DATA ANALYSIS ############# 
#########################################
#########################################
##### -- Beta diversity of lakes -- #####
#########################################
#########################################
##### -- Estimate beta diversity -- #####
#########################################
##### -- Generate datasets -- #####
#### Macro-invertebrates
### Filter dataset for shared lakes
inverts_combined_byLake_complete <- subset(inverts_byLake_combined, lake_code %in% shared_lake_codes)
### Generate presence-absence datasets
inverts_combined_byLake_complete_binary <- data.frame(lake_code = inverts_combined_byLake_complete[, 1], decostand(inverts_combined_byLake_complete[, -1], method = "pa"))

#### Microbes
### Filter dataset for shared lakes
micro_byLake_complete <- subset(micro_samples_normalized, lake_code %in% shared_lake_codes)
### Generate presence-absence datasets
micro_byLake_complete_binary <- data.frame(lake_code = micro_byLake_complete[, 1], decostand(micro_byLake_complete[, -1], method = "pa"))

##### Create beta_byLake data frame to include dissimilarity summaries
beta_byLake <- data.frame(lake_code = shared_lake_codes)
#### Add stratified field
beta_byLake$stratified <- alpha_byLake$stratified

##### -- Calculating dissimilarities -- #####
##### Macro-invertebrates
#### Abundance-based
inverts_combined_byLake_horn <- vegdist(inverts_combined_byLake_complete[, -1], method = "horn")
inverts_combined_byLake_horn_df <- inverts_combined_byLake_horn %>% as.matrix() %>% as.data.frame() 
beta_byLake$inverts_combined_horn_median <- apply(inverts_combined_byLake_horn_df, 1, function(x) quantile(x[x > 0], .5, na.rm = TRUE))
beta_byLake$inverts_combined_horn_lower <- apply(inverts_combined_byLake_horn_df, 1, function(x) quantile(x[x > 0], .375, na.rm = TRUE))
beta_byLake$inverts_combined_horn_upper <- apply(inverts_combined_byLake_horn_df, 1, function(x) quantile(x[x > 0], .675, na.rm = TRUE))
#### Binary presence-absence
inverts_combined_byLake_horn_binary <- vegdist(inverts_combined_byLake_complete_binary[, -1], method = "horn")
inverts_combined_byLake_horn_binary_df <- inverts_combined_byLake_horn_binary %>% as.matrix() %>% as.data.frame() 
beta_byLake$inverts_combined_horn_binary_median <- apply(inverts_combined_byLake_horn_binary_df, 1, function(x) quantile(x[x > 0], .5, na.rm = TRUE))
beta_byLake$inverts_combined_horn_binary_lower <- apply(inverts_combined_byLake_horn_binary_df, 1, function(x) quantile(x[x > 0], .375, na.rm = TRUE))
beta_byLake$inverts_combined_horn_binary_upper <- apply(inverts_combined_byLake_horn_binary_df, 1, function(x) quantile(x[x > 0], .675, na.rm = TRUE))

##### Microbes
#### Abundance-based
micro_byLake_horn <- vegdist(micro_byLake_complete[, -1], method = "horn") %>% as.matrix() %>% as.data.frame() 
micro_byLake_horn_df <- micro_byLake_horn %>% as.matrix() %>% as.data.frame() 
beta_byLake$micro_horn_median <- apply(micro_byLake_horn_df, 1, function(x) quantile(x[x > 0], .5, na.rm = TRUE))
beta_byLake$micro_horn_lower <- apply(micro_byLake_horn_df, 1, function(x) quantile(x[x > 0], .375, na.rm = TRUE))
beta_byLake$micro_horn_upper <- apply(micro_byLake_horn_df, 1, function(x) quantile(x[x > 0], .675, na.rm = TRUE))
#### Abundance-based
micro_byLake_horn_binary <- vegdist(micro_byLake_complete_binary[, -1], method = "horn") %>% as.matrix() %>% as.data.frame() 
micro_byLake_horn_binary_df <- micro_byLake_horn_binary %>% as.matrix() %>% as.data.frame() 
beta_byLake$micro_horn_binary_median <- apply(micro_byLake_horn_binary_df, 1, function(x) quantile(x[x > 0], .5, na.rm = TRUE))
beta_byLake$micro_horn_binary_lower <- apply(micro_byLake_horn_binary_df, 1, function(x) quantile(x[x > 0], .375, na.rm = TRUE))
beta_byLake$micro_horn_binary_upper <- apply(micro_byLake_horn_binary_df, 1, function(x) quantile(x[x > 0], .675, na.rm = TRUE))

##### -- Comparing dissimilarities -- #####
##### Mantel test of the correlation between two dissimilarity matrices
#### Check correlation between microbes and macro-invertebrates
## Abundance-based
mantel(micro_byLake_horn, inverts_combined_byLake_horn, permutations = 9999)
## Binary
mantel(micro_byLake_horn_binary, inverts_combined_byLake_horn_binary, permutations = 9999)
##### Run principal coordinates analysis
### Abundance-based
inverts_combined_beta_horn_PCoA <- capscale(inverts_combined_byLake_complete[, -1] ~ 1, distance = "horn")
micro_beta_horn_PCoA <- capscale(micro_byLake_complete[, -1] ~ 1, distance = "horn")
### Binary
inverts_combined_binary_beta_horn_PCoA <- capscale(inverts_combined_byLake_complete_binary[, -1] ~ 1, distance = "horn")
micro_binary_beta_horn_PCoA <- capscale(micro_byLake_complete_binary[, -1] ~ 1, distance = "horn")
#### Check correlation between microbes and macro-invertebrates using procrustes rotation tests
protest(micro_beta_horn_PCoA, inverts_combined_beta_horn_PCoA, permutations = 9999)
protest(micro_binary_beta_horn_PCoA, inverts_combined_binary_beta_horn_PCoA, permutations = 9999)

##### Generate a data frame with pairwise dissimilarities
names(inverts_combined_byLake_horn_df) <- inverts_combined_byLake_complete$lake_code
row.names(inverts_combined_byLake_horn_df) <- inverts_combined_byLake_complete$lake_code
pairwise_horn_dist <- stack(inverts_combined_byLake_horn_df)
pairwise_horn_dist$lake2 <- rep(inverts_combined_byLake_complete$lake_code, times = 12)
names(pairwise_horn_dist) <- c("inverts_combined_horn_dissimilarity", "lake1", "lake2")
pairwise_horn_dist$micro_horn_dissimilarity <- stack(micro_byLake_horn_df)$values
pairwise_horn_dist$inverts_combined_binary_horn_dissimilarity <- stack(inverts_combined_byLake_horn_binary_df)$values
pairwise_horn_dist$micro_binary_horn_dissimilarity <- stack(micro_byLake_horn_binary_df)$values
pairwise_horn_dist$comparison_type <- 1
pairwise_horn_dist$comparison_type <- ifelse(pairwise_horn_dist$lake1 %in% subset(alpha_byLake, stratified == 1)$lake_code & pairwise_horn_dist$lake2 %in% subset(alpha_byLake, stratified == 1)$lake_code, 2, pairwise_horn_dist$comparison_type) 
pairwise_horn_dist$comparison_type <- ifelse(pairwise_horn_dist$lake1 %in% subset(alpha_byLake, stratified == 0)$lake_code & pairwise_horn_dist$lake2 %in% subset(alpha_byLake, stratified == 0)$lake_code, 3, pairwise_horn_dist$comparison_type) 
pairwise_horn_dist$comparison_type <- as.factor(pairwise_horn_dist$comparison_type)
pairwise_horn_dist <- pairwise_horn_dist[c("lake1", "lake2", "comparison_type", "inverts_combined_horn_dissimilarity", "micro_horn_dissimilarity", "inverts_combined_binary_horn_dissimilarity", "micro_binary_horn_dissimilarity")]
pairwise_horn_dist <- pairwise_horn_dist[-unlist(lapply(seq_along(1:12), function(i) (12*(i-1))+(1:i))), ]
############################################################
##### -- Environmental correlates of beta diversity -- #####
############################################################
##### Effect of lake type on pairwise dissimilairty 
summary.lm(lm(micro_dissimilarity ~ comparison_type, data = pairwise_dist))
summary.lm(lm(inverts_combined_dissimilarity ~ comparison_type, data = pairwise_dist))
with(subset(pairwise_horn_dist, comparison_type == 1), c(median(micro_dissimilarity), median(inverts_combined_dissimilarity)))
with(subset(pairwise_horn_dist, comparison_type == 2), c(median(micro_dissimilarity), median(inverts_combined_dissimilarity)))
with(subset(pairwise_horn_dist, comparison_type == 3), c(median(micro_dissimilarity), median(inverts_combined_dissimilarity)))

##### Constrained Redundancy Analysis
##### Macro-invertebrates
##### Binary
#### Combined dataset
#### Join environmental data
inverts_combined_byLake_complete_binary_withEnv <- left_join(inverts_combined_byLake_complete_binary, environment_byLake, by = "lake_code")
### Add row.names
row.names(inverts_combined_byLake_complete_binary_withEnv) <- inverts_combined_byLake_complete_binary_withEnv$lake_code
### RDA MAM model
inverts_combined_beta_rda_horn_binary_intercept <- capscale(inverts_combined_byLake_complete_binary_withEnv[names(inverts_combined_byLake_complete_binary)[-1]] ~ 1, inverts_combined_byLake_complete_binary_withEnv[names(environment_byLake)[-c(1, ncol(environment_byLake))]], distance = "horn")
inverts_combined_beta_rda_horn_binary_full <- inverts_combined_beta_rda_horn_binary_reduced <- capscale(inverts_combined_byLake_complete_binary_withEnv[names(inverts_combined_byLake_complete_binary)[-1]] ~ ., inverts_combined_byLake_complete_binary_withEnv[names(environment_byLake)[-c(1, ncol(environment_byLake))]], distance = "horn")
anova(inverts_combined_beta_rda_horn_binary_full, by="term", permutations=99999)
inverts_combined_beta_rda_horn_binary_step_backward <- ordistep(inverts_combined_beta_rda_horn_binary_full, scope = list(lower = formula(inverts_combined_beta_rda_horn_binary_intercept), upper = formula(inverts_combined_beta_rda_horn_binary_full)), direction = "backward", permutations = 9999)
inverts_combined_beta_rda_horn_binary_reduced <- update(inverts_combined_beta_rda_horn_binary_reduced, formula(paste(". ~ . - ", row.names(inverts_combined_beta_rda_horn_binary_reduced$CCA$biplot)[which.min(anova(inverts_combined_beta_rda_horn_binary_reduced,  by="term", permutations=9999)$F)], sep = "")))
inverts_combined_beta_rda_horn_binary_reduced <- update(inverts_combined_beta_rda_horn_binary_reduced, formula(paste(". ~ . - ", row.names(inverts_combined_beta_rda_horn_binary_reduced$CCA$biplot)[which.min(anova(inverts_combined_beta_rda_horn_binary_reduced,  by="term", permutations=9999)$F)], sep = "")))
inverts_combined_beta_rda_horn_binary_reduced <- update(inverts_combined_beta_rda_horn_binary_reduced, formula(paste(". ~ . - ", row.names(inverts_combined_beta_rda_horn_binary_reduced$CCA$biplot)[which.min(anova(inverts_combined_beta_rda_horn_binary_reduced,  by="term", permutations=9999)$F)], sep = "")))
inverts_combined_beta_rda_horn_binary_reduced <- update(inverts_combined_beta_rda_horn_binary_reduced, formula(paste(". ~ . - ", row.names(inverts_combined_beta_rda_horn_binary_reduced$CCA$biplot)[which.min(anova(inverts_combined_beta_rda_horn_binary_reduced,  by="term", permutations=9999)$F)], sep = "")))
anova(inverts_combined_beta_rda_horn_binary_reduced, by="term", permutations=99999)

##### Abundance-based
#### Combined dataset
#### Join environmental data
inverts_combined_byLake_complete_withEnv <- left_join(inverts_combined_byLake_complete, environment_byLake, by = "lake_code")
### Add row.names
row.names(inverts_combined_byLake_complete_withEnv) <- inverts_combined_byLake_complete_withEnv$lake_code
### RDA MAM model
inverts_combined_beta_rda_horn_intercept <- capscale(inverts_combined_byLake_complete_withEnv[names(inverts_combined_byLake_complete)[-1]] ~ 1, inverts_combined_byLake_complete_withEnv[names(environment_byLake)[-c(1, ncol(environment_byLake))]], distance = "horn")
inverts_combined_beta_rda_horn_full <- inverts_combined_beta_rda_horn_reduced <- capscale(inverts_combined_byLake_complete_withEnv[names(inverts_combined_byLake_complete)[-1]] ~ ., inverts_combined_byLake_complete_withEnv[names(environment_byLake)[-c(1, ncol(environment_byLake))]], distance = "horn")
anova(inverts_combined_beta_rda_horn_full, by="term", permutations=99999)
inverts_combined_beta_rda_horn_step_backward <- ordistep(inverts_combined_beta_rda_horn_full, scope = list(lower = formula(inverts_combined_beta_rda_horn_intercept), upper = formula(inverts_combined_beta_rda_horn_full)), direction = "backward", permutations = 9999)
inverts_combined_beta_rda_horn_reduced <- update(inverts_combined_beta_rda_horn_reduced, formula(paste(". ~ . - ", row.names(inverts_combined_beta_rda_horn_reduced$CCA$biplot)[which.min(anova(inverts_combined_beta_rda_horn_reduced,  by="term", permutations=9999)$F)], sep = "")))
inverts_combined_beta_rda_horn_reduced <- update(inverts_combined_beta_rda_horn_reduced, formula(paste(". ~ . - ", row.names(inverts_combined_beta_rda_horn_reduced$CCA$biplot)[which.min(anova(inverts_combined_beta_rda_horn_reduced,  by="term", permutations=9999)$F)], sep = "")))
inverts_combined_beta_rda_horn_reduced <- update(inverts_combined_beta_rda_horn_reduced, formula(paste(". ~ . - ", row.names(inverts_combined_beta_rda_horn_reduced$CCA$biplot)[which.min(anova(inverts_combined_beta_rda_horn_reduced,  by="term", permutations=9999)$F)], sep = "")))
inverts_combined_beta_rda_horn_reduced <- update(inverts_combined_beta_rda_horn_reduced, formula(paste(". ~ . - ", row.names(inverts_combined_beta_rda_horn_reduced$CCA$biplot)[which.min(anova(inverts_combined_beta_rda_horn_reduced,  by="term", permutations=9999)$F)], sep = "")))
anova(inverts_combined_beta_rda_horn_reduced, by="term", permutations=99999)

##### Microbes
##### Binary
#### Join environmental data
micro_byLake_complete_binary_withEnv <- left_join(micro_byLake_complete_binary, environment_byLake, by = "lake_code")
### Add row.names
row.names(micro_byLake_complete_binary_withEnv) <- micro_byLake_complete_binary_withEnv$lake_code
### RDA MAM model
micro_beta_rda_horn_binary_intercept <- capscale(micro_byLake_complete_binary_withEnv[names(micro_byLake_complete_binary)[-1]] ~ 1, micro_byLake_complete_binary_withEnv[names(environment_byLake)[-c(1, ncol(environment_byLake))]], distance = "horn")
micro_beta_rda_horn_binary_full <- micro_beta_rda_horn_binary_reduced <- capscale(micro_byLake_complete_binary_withEnv[names(micro_byLake_complete_binary)[-1]] ~ ., micro_byLake_complete_binary_withEnv[names(environment_byLake)[-c(1, ncol(environment_byLake))]], distance = "horn")
anova(micro_beta_rda_horn_binary_full, by="term", permutations=9999)
micro_beta_rda_horn_binary_backward <- ordistep(micro_beta_rda_horn_binary_full, scope = list(lower = formula(micro_beta_rda_horn_binary_intercept), upper = formula(micro_beta_rda_horn_binary_full)), direction = "backward", permutations = 9999)
micro_beta_rda_horn_binary_reduced <- update(micro_beta_rda_horn_binary_reduced, formula(paste(". ~ . - ", row.names(micro_beta_rda_horn_binary_reduced$CCA$biplot)[which.min(anova(micro_beta_rda_horn_binary_reduced,  by="term", permutations=9999)$F)], sep = "")))
micro_beta_rda_horn_binary_reduced <- update(micro_beta_rda_horn_binary_reduced, formula(paste(". ~ . - ", row.names(micro_beta_rda_horn_binary_reduced$CCA$biplot)[which.min(anova(micro_beta_rda_horn_binary_reduced,  by="term", permutations=9999)$F)], sep = "")))
micro_beta_rda_horn_binary_reduced <- update(micro_beta_rda_horn_binary_reduced, formula(paste(". ~ . - ", row.names(micro_beta_rda_horn_binary_reduced$CCA$biplot)[which.min(anova(micro_beta_rda_horn_binary_reduced,  by="term", permutations=9999)$F)], sep = "")))
micro_beta_rda_horn_binary_reduced <- update(micro_beta_rda_horn_binary_reduced, formula(paste(". ~ . - ", row.names(micro_beta_rda_horn_binary_reduced$CCA$biplot)[which.min(anova(micro_beta_rda_horn_binary_reduced,  by="term", permutations=9999)$F)], sep = "")))
anova(micro_beta_rda_horn_binary_reduced, by="term", permutations=99999)

##### Abundance-based 
#### Join environmental data
micro_byLake_complete_withEnv <- left_join(micro_byLake_complete, environment_byLake, by = "lake_code")
### Add row.names
row.names(micro_byLake_complete_withEnv) <- micro_byLake_complete_withEnv$lake_code
### RDA MAM model
micro_beta_rda_horn_intercept <- capscale(micro_byLake_complete_withEnv[names(micro_byLake_complete)[-1]] ~ 1, micro_byLake_complete_withEnv[names(environment_byLake)[-c(1, ncol(environment_byLake))]], distance = "horn")
micro_beta_rda_horn_full <- micro_beta_rda_horn_reduced <- capscale(micro_byLake_complete_withEnv[names(micro_byLake_complete)[-1]] ~ ., micro_byLake_complete_withEnv[names(environment_byLake)[-c(1, ncol(environment_byLake))]], distance = "horn")
anova(micro_beta_rda_horn_full, by="term", permutations=9999)
micro_beta_rda_horn_backward <- ordistep(micro_beta_rda_horn_full, scope = list(lower = formula(micro_beta_rda_horn_intercept), upper = formula(micro_beta_rda_horn_full)), direction = "backward", permutations = 9999)
micro_beta_rda_horn_intercept <- capscale(micro_byLake_complete_withEnv[names(micro_byLake_complete)[-1]] ~ 1, micro_byLake_complete_withEnv[names(environment_byLake)[-c(1, ncol(environment_byLake))]], distance = "horn")
micro_beta_rda_horn_full <- micro_beta_rda_horn_reduced <- capscale(micro_byLake_complete_withEnv[names(micro_byLake_complete)[-1]] ~ ., micro_byLake_complete_withEnv[names(environment_byLake)[-c(1, ncol(environment_byLake))]], distance = "horn")
anova(micro_beta_rda_horn_full, by="term", permutations=9999)
micro_beta_rda_horn_backward <- ordistep(micro_beta_rda_horn_full, scope = list(lower = formula(micro_beta_rda_horn_intercept), upper = formula(micro_beta_rda_horn_full)), direction = "backward", permutations = 9999)
micro_beta_rda_horn_reduced <- update(micro_beta_rda_horn_reduced, formula(paste(". ~ . - ", row.names(micro_beta_rda_horn_reduced$CCA$biplot)[which.min(anova(micro_beta_rda_horn_reduced,  by="term", permutations=9999)$F)], sep = "")))
micro_beta_rda_horn_reduced <- update(micro_beta_rda_horn_reduced, formula(paste(". ~ . - ", row.names(micro_beta_rda_horn_reduced$CCA$biplot)[which.min(anova(micro_beta_rda_horn_reduced,  by="term", permutations=9999)$F)], sep = "")))
micro_beta_rda_horn_reduced <- update(micro_beta_rda_horn_reduced, formula(paste(". ~ . - ", row.names(micro_beta_rda_horn_reduced$CCA$biplot)[which.min(anova(micro_beta_rda_horn_reduced,  by="term", permutations=9999)$F)], sep = "")))
micro_beta_rda_horn_reduced <- update(micro_beta_rda_horn_reduced, formula(paste(". ~ . - ", row.names(micro_beta_rda_horn_reduced$CCA$biplot)[which.min(anova(micro_beta_rda_horn_reduced,  by="term", permutations=9999)$F)], sep = "")))
anova(micro_beta_rda_horn_reduced, by="term", permutations=99999)
