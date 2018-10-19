##########################################
############## micro-macro ###############
##########################################
############# DATA ANALYSIS ############## 
##########################################
##########################################
##### -- Alpha diversity of lakes -- #####
##########################################
##########################################
##### -- Estimate alpha diversity -- #####
##########################################
##### Macro-invertebrates
inverts_abundance_byLake_alpha <- calc_div(inverts_abundance_byLake, label = "inverts_combined")
##### Microbes
microbiota_alpha_normalized_byLake <- calc_div(micro_abundance_byLake, label = "microbiota_normalized")
##### Merge micro and macro
#### Start with shared lake codes
alpha_byLake <- data.frame(lake_code = shared_lake_codes)
#### Merge combined data
alpha_byLake <- join_all(list(alpha_byLake, microbiota_alpha_normalized_byLake, inverts_abundance_byLake_alpha), by = "lake_code")

#############################################################
##### -- Environmental correlates of alpha diversity -- #####
#############################################################
##### Merge alpha diversity with environmental data
alpha_byLake <- left_join(alpha_byLake, environment_byLake, by = "lake_code")
##### Calculate linear correlations
micro_alpha_ace_corr <- alpha_byLake[c("microbiota_normalized_ace_estimate", "oxygen_median", "oxygen_sd", "conductivity_median", "CHL_median", "temperature_median", "PAR_sd", "surface_area_m2", "distance_to_ocean_min_m")] %>% cor()
micro_alpha_ace_corr <- round(micro_alpha_ace_corr[1, -1], 3)
micro_alpha_shannon_corr <- alpha_byLake[c("microbiota_normalized_shannon_diversity_estimate", "oxygen_median", "oxygen_sd", "conductivity_median", "CHL_median", "temperature_median", "PAR_sd", "surface_area_m2", "distance_to_ocean_min_m")] %>% cor()
micro_alpha_shannon_corr <- round(micro_alpha_shannon_corr[1, -1], 3)
micro_alpha_simpson_corr <- alpha_byLake[c("microbiota_normalized_simpson_evenness", "oxygen_median", "oxygen_sd", "conductivity_median", "CHL_median", "temperature_median", "PAR_sd", "surface_area_m2", "distance_to_ocean_min_m")] %>% cor()
micro_alpha_simpson_corr <- round(micro_alpha_simpson_corr[1, -1], 3)
micro_alpha_dominance_corr <- alpha_byLake[c("microbiota_normalized_dominance", "oxygen_median", "oxygen_sd", "conductivity_median", "CHL_median", "temperature_median", "PAR_sd", "surface_area_m2", "distance_to_ocean_min_m")] %>% cor()
micro_alpha_dominance_corr <- round(micro_alpha_dominance_corr[1, -1], 3)
inverts_alpha_ace_corr <- alpha_byLake[c("inverts_combined_ace_estimate", "oxygen_median", "oxygen_sd", "conductivity_median", "CHL_median", "temperature_median", "PAR_sd", "surface_area_m2", "distance_to_ocean_min_m")] %>% cor()
inverts_alpha_ace_corr <- round(inverts_alpha_ace_corr[1, -1], 3)
inverts_alpha_shannon_corr <- alpha_byLake[c("inverts_combined_shannon_diversity_estimate", "oxygen_median", "oxygen_sd", "conductivity_median", "CHL_median", "temperature_median", "PAR_sd", "surface_area_m2", "distance_to_ocean_min_m")] %>% cor()
inverts_alpha_shannon_corr <- round(inverts_alpha_shannon_corr[1, -1], 3)
inverts_alpha_simpson_corr <- alpha_byLake[c("inverts_combined_simpson_evenness", "oxygen_median", "oxygen_sd", "conductivity_median", "CHL_median", "temperature_median", "PAR_sd", "surface_area_m2", "distance_to_ocean_min_m")] %>% cor()
inverts_alpha_simpson_corr <- round(inverts_alpha_simpson_corr[1, -1], 3)
inverts_alpha_dominance_corr <- alpha_byLake[c("inverts_combined_dominance", "oxygen_median", "oxygen_sd", "conductivity_median", "CHL_median", "temperature_median", "PAR_sd", "surface_area_m2", "distance_to_ocean_min_m")] %>% cor()
inverts_alpha_dominance_corr <- round(inverts_alpha_dominance_corr[1, -1], 3)

write.csv(data.frame(taxon = rep(c("Microbes", "Macro-invertebrates"), each = 4),
                     measure = rep(c("ACE", "Shannon diversity", "Evenness", "Dominance"), times = 2),
                     rbind(micro_alpha_ace_corr, micro_alpha_shannon_corr, micro_alpha_simpson_corr, micro_alpha_dominance_corr, 
                           inverts_alpha_ace_corr, inverts_alpha_shannon_corr, inverts_alpha_simpson_corr, inverts_alpha_dominance_corr)
),
"output/alpha_correlates_table.csv", row.names = FALSE
)



