###################################################
################### micro-macro ###################
###################################################
##################### FIGURES #####################
###################################################
##########################################
##### -- Alpha diversity of lakes -- #####
##########################################
##### -- Figure 2: Alpha diversity of macro-invertebrates vs microbes -- #####
#### ACE diversity
alpha_fig_ACE <- alpha_div_fig(y_var = "microbiota_normalized_ace_estimate", x_var = "inverts_combined_ace_estimate", 
                           fit_line = FALSE, 
                           error_bars = TRUE, y_upper = "microbiota_normalized_ace_upper", y_lower = "microbiota_normalized_ace_lower",
                           x_upper = "inverts_combined_ace_upper", x_lower = "inverts_combined_ace_lower",
                           x_min = 0, x_max = 200, y_min = 40000, y_max = 132000, 
                           text_v = c(rep(-1, 3), 1.6, rep(-1, 2), 1.6, 1.3, rep(-1, 3), 1.6), text_h = c(rep(1.4, 10), -0.5, 1.4), 
                           y_lab = "Richness of microbes", x_lab = "Richness of macro-invertebrates")
#### Simpson diversity
alpha_fig_Evenness <- alpha_div_fig(y_var = "microbiota_normalized_simpson_evenness", x_var = "inverts_combined_simpson_evenness", 
                                fit_line = FALSE, 
                                error_bars = FALSE, 
                                x_min = 0, x_max = 0.10, y_min = 0.0002, y_max = 0.0009,
                                text_v = c(-1, 1.6, -0.2, -0.5, 1.6, -1, -1, 1.6, -1, -1, -0.2, -1), text_h = c(rep(1.4, 10), 1.8, 1.4), 
                                y_lab = "Evenness of microbes", x_lab = "Evenness of macro-invertebrates")


#### Shannon diversity
alpha_fig_Obs <- alpha_div_fig(y_var = "microbiota_normalized_shannon_diversity_estimate", x_var = "inverts_combined_shannon_diversity_estimate", 
                                   fit_line = FALSE, 
                                   error_bars = TRUE, y_upper = "microbiota_normalized_shannon_diversity_upper", y_lower = "microbiota_normalized_shannon_diversity_lower",
                                   x_upper = "inverts_combined_shannon_diversity_upper", x_lower = "inverts_combined_shannon_diversity_lower",
                                   x_min = 0, x_max = 65, y_min = 0, y_max = 3500, 
                                   text_v = c(-1, 1.6, -1, -1, 1.6, -1, -1, 1.6, -1, -1, -0.2, -1), text_h = c(rep(1.4, 10), 1.8, 1.4), 
                                   y_lab = "OTU richness microbes", x_lab = "Species richness of macro-invertebrates")
#### Shannon diversity
alpha_fig_Shannon <- alpha_div_fig(y_var = "microbiota_normalized_shannon_diversity_estimate", x_var = "inverts_combined_shannon_diversity_estimate", 
                                   fit_line = FALSE, 
                                   error_bars = TRUE, y_upper = "microbiota_normalized_shannon_diversity_upper", y_lower = "microbiota_normalized_shannon_diversity_lower",
                                   x_upper = "inverts_combined_shannon_diversity_upper", x_lower = "inverts_combined_shannon_diversity_lower",
                                   x_min = 0, x_max = 65, y_min = 0, y_max = 3500, 
                                   text_v = c(-1, 1.6, -1, -1, 1.6, -1, -1, 1.6, -1, -1, -0.2, -1), text_h = c(rep(1.4, 10), 1.8, 1.4), 
                                   y_lab = "Shannon diversity of microbes", x_lab = "Shannon diversity of macro-invertebrates")
#### Pielou
alpha_fig_Pielou <- alpha_div_fig(y_var = "microbiota_normalized_pielou", x_var = "inverts_combined_pielou", 
                                   fit_line = FALSE, 
                                   x_min = 0.4, x_max = 1, y_min = 0.4, y_max = 1, 
                                   text_v = c(-1, 1.6, -1, -1, 1.6, -1, -1, 1.6, -1, -1, -0.2, -1), text_h = c(rep(1.4, 10), 1.8, 1.4), 
                                   y_lab = "Pielou's evenness of microbes", x_lab = "Pielou's evenness of macro-invertebrates")
#### Dominance diversity
alpha_fig_Dominance <- alpha_div_fig(y_var = "microbiota_normalized_dominance", x_var = "inverts_combined_dominance", 
                                 fit_line = FALSE, 
                                 error_bars = FALSE,
                                 x_min = 0, x_max = 0.6, y_min = 0.025, y_max = 0.09, 
                                 text_v = c(-1, 1.6, -1, -1, 1.6, -1, -1, 1.6, -1, -1, -0.2, -1), text_h = c(rep(1.4, 10), 1.8, 1.4), 
                                 y_lab = "Dominance of microbes", x_lab = "Dominance of macro-invertebrates")

tiff("figures/fig2-alpha_diversity.tiff", width = 12.0, height = 6.0, units = 'in', res = 300)
ggarrange(alpha_fig_ACE, alpha_fig_Evenness, labels = c("A", "B"), font.label = list(size = 18, face = "bold"))
dev.off()

bmp("figures/fig2-alpha_diversity.bmp", width = 12.0, height = 6.0, units = 'in', res = 220)
ggarrange(alpha_fig_ACE, alpha_fig_Evenness, labels = c("A", "B"), font.label = list(size = 18, face = "bold"))
dev.off()

png("figures/fig2-alpha_diversity.png", width = 12.0, height = 6.0, units = 'in', res = 220)
ggarrange(alpha_fig_ACE, alpha_fig_Evenness, labels = c("A", "B"), font.label = list(size = 18, face = "bold"))
dev.off()

jpeg("figures/fig2-alpha_diversity.jpg", width = 12.0, height = 6.0, units = 'in', res = 300)
ggarrange(alpha_fig_ACE, alpha_fig_Evenness, labels = c("A", "B"), font.label = list(size = 18, face = "bold"))
dev.off()

tiff("figures/figS2-alpha_diversity.tif", width = 12.0, height = 12.0, units = 'in', res = 600)
ggarrange(alpha_fig_Obs, alpha_fig_Shannon, alpha_fig_Pielou, alpha_fig_Dominance, labels = c("A", "B", "C", "D"))
dev.off()

#### Correlations
with(alpha_byLake, cor.test(microbiota_normalized_observed_richness, inverts_combined_observed_richness, method = "spearman"))
with(alpha_byLake, cor.test(microbiota_normalized_ace_estimate, inverts_combined_ace_estimate, method = "spearman"))
with(alpha_byLake, cor.test(microbiota_normalized_chao1_estimate, inverts_combined_chao1_estimate, method = "spearman"))
with(alpha_byLake, cor.test(microbiota_normalized_shannon_diversity_estimate, inverts_combined_shannon_diversity_estimate, method = "spearman"))
with(alpha_byLake, cor.test(microbiota_normalized_simpson_diversity_estimate, inverts_combined_simpson_diversity_estimate, method = "spearman"))
with(alpha_byLake, cor.test(microbiota_normalized_pielou, inverts_combined_pielou, method = "spearman"))
with(alpha_byLake, cor.test(microbiota_normalized_simpson_evenness, inverts_combined_simpson_evenness, method = "spearman"))
with(alpha_byLake, cor.test(microbiota_normalized_dominance, inverts_combined_dominance, method = "spearman"))

#############################################################
##### -- Environmental correlates of alpha diversity -- #####
#############################################################
##### -- Figure 2: Environmental correlates of alpha diversity -- #####
top_col1 <- ggarrange(div_byEnv_factor(dat = alpha_byLake, response = "microbiota_normalized_ace_estimate", predictor = "stratified", y_lab = "Richness", x_lab = "", annot = c(1, 60000), plot.margin=margin(l = 14, r = 5, b = 4)),
                      div_byEnv_factor(dat = alpha_byLake, response = "microbiota_normalized_simpson_evenness", predictor = "stratified", y_lab = "Evenness", x_lab = "", annot = c(2, 0.0003), y_min = 0, y_max = 0.0009, plot.margin=margin(l = 14, r = 5, b = 4)),
                      nrow = 2, ncol = 1, align = "v"
)
top_col2 <- ggarrange(div_byEnv(dat = alpha_byLake, response = "microbiota_normalized_ace_estimate", predictor = "oxygen_median", model_type = "lm", y_lab = "", x_lab = "", annot = c(3, 60000), axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin=margin(l = 0, r = 5, b = 4)),
                      div_byEnv(dat = alpha_byLake, response = "microbiota_normalized_simpson_evenness", predictor = "oxygen_median", y_lab = "", x_lab = "", annot = c(3, 0.0001), y_min = 0, y_max = 0.0009, axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin=margin(l = 0, r = 5, b = 4)),
                      nrow = 2, ncol = 1, align = "v"
)
top_col3 <- ggarrange(div_byEnv(dat = alpha_byLake, response = "microbiota_normalized_ace_estimate", predictor = "oxygen_sd", y_lab = "", x_lab = "", annot = c(0, 60000), axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin=margin(l = 0, r = 5, b = 4)),
                      div_byEnv(dat = alpha_byLake, response = "microbiota_normalized_simpson_evenness", predictor = "oxygen_sd", y_lab = "", x_lab = "", annot = c(-0.5, 0.0001), y_min = 0, y_max = 0.0009, axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin=margin(l = 0, r = 5, b = 4)),
                      nrow = 2, ncol = 1, align = "v"
)                    
top_col4 <- ggarrange(div_byEnv(dat = alpha_byLake, response = "microbiota_normalized_ace_estimate", predictor = "distance_to_ocean_min_m", y_lab = "", x_lab = "", annot = c(200, 60000), axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin=margin(l = 0, r = 5, b = 4)),
                      div_byEnv(dat = alpha_byLake, response = "microbiota_normalized_simpson_evenness", predictor = "distance_to_ocean_min_m", y_lab = "", x_lab = "", annot = c(200, 0.0001), y_min = 0, y_max = 0.0009, axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin=margin(l = 0, r = 5, b = 4)),
                      nrow = 2, ncol = 1, align = "v"
)
top_rows <- ggarrange(top_col1, top_col2, top_col3, top_col4, nrow = 1, ncol = 4, widths = c(1.25, 1, 1, 1))

bottom_col1 <- ggarrange(div_byEnv_factor(dat = alpha_byLake, response = "inverts_combined_ace_estimate", predictor = "stratified", y_lab = "Richness", x_lab = "", annot = c(1, 40), plot.margin=margin(l = 14, r = 5, b = 4)),
                         div_byEnv_factor(dat = alpha_byLake, response = "inverts_combined_simpson_evenness", predictor = "stratified", y_lab = "Evenness", x_lab = "Lake type", annot = c(2, 0.05), plot.margin=margin(l = 14, r = 5, b = 4)),
                         nrow = 2, ncol = 1, align = "v"
)
bottom_col2 <- ggarrange(div_byEnv(dat = alpha_byLake, response = "inverts_combined_ace_estimate", predictor = "oxygen_median", y_lab = "", x_lab = "", annot = c(3.5, 40), axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin=margin(l = 0, r = 5, b = 4)),
                         div_byEnv(dat = alpha_byLake, response = "inverts_combined_simpson_evenness", predictor = "oxygen_median", y_lab = "", x_lab = "Oxygen median", annot = c(3, 0.05), axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin=margin(l = 0, r = 5, b = 4)),
                         nrow = 2, ncol = 1, align = "v"
)
bottom_col3 <- ggarrange(div_byEnv(dat = alpha_byLake, response = "inverts_combined_ace_estimate", predictor = "oxygen_sd", y_lab = "", x_lab = "", annot = c(-0.5, 40), axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin=margin(l = 0, r = 0, b = 4)),
                         div_byEnv(dat = alpha_byLake, response = "inverts_combined_simpson_evenness", predictor = "oxygen_sd", y_lab = "", x_lab = "Oxygen variation", annot = c(1, 0.05), axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin=margin(l = 0, r = 5, b = 4)),
                         nrow = 2, ncol = 1, align = "v"
)                    
bottom_col4 <- ggarrange(div_byEnv(dat = alpha_byLake, response = "inverts_combined_ace_estimate", predictor = "distance_to_ocean_min_m", y_lab = "", x_lab = "", annot = c(200, 90), axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin=margin(l = 0, r = 0, b = 4)),
                         div_byEnv(dat = alpha_byLake, response = "inverts_combined_simpson_evenness", predictor = "distance_to_ocean_min_m", y_lab = "", x_lab = "Distance to ocean", annot = c(200, 0.05), axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin=margin(l = 0, r = 5, b = 4)),
                         nrow = 2, ncol = 1, align = "v"
)

bottom_rows <- ggarrange(bottom_col1, bottom_col2, bottom_col3, bottom_col4, nrow = 1, ncol = 4, widths = c(1.25, 1, 1, 1))

tiff("figures/fig3-alpha_correlates.tiff", width = 12.0, height = 12.0, units = 'in', res = 600)
ggarrange(annotate_figure(top_rows, left = text_grob("      Microbes", rot = 90, size = 19)),
          annotate_figure(bottom_rows, left = text_grob("      Macro-invertebrates", rot = 90, size = 19)),
          nrow = 2, ncol = 1
)
dev.off()

##### Correlations
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

############################################
##### -- Beta-diversity among lakes -- #####
############################################
#### Check correlation between macrobiota and microbiota
### Mantel test uses dissimilarity matrices
binary_diss_cor <- mantel(micro_byLake_horn_binary, inverts_combined_byLake_horn_binary, permutations = 9999, method = "spearman")
abun_diss_cor <- mantel(micro_byLake_horn, inverts_combined_byLake_horn, permutations = 9999, method = "spearman")

##### -- Figure 4: Beta diversity of macro-invertebrates vs microbes -- #####
##### Binary dissimilarity
#### Plot correlation between pairwise distances
fig4_pairwise_diss_bin <- ggplot(pairwise_horn_dist, aes(y = micro_binary_horn_dissimilarity, x = inverts_combined_binary_horn_dissimilarity, group = comparison_type, color = comparison_type)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(grey(.7), "green", "deepskyblue3")) +
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.25, color = grey(0.35)) +
  geom_smooth(aes(y = micro_binary_horn_dissimilarity, x = inverts_combined_binary_horn_dissimilarity), pairwise_horn_dist, method=lm, se = FALSE, color="black", inherit.aes = FALSE) +  # Add linear regression line
  theme_classic() +
  theme(legend.position = "none", text = element_text(size=17)) + 
  ylab("Microbial dissimilarity") +
  xlab("") +
  annotate("text", 0.80, 0.10, label = paste("r = ", round(binary_diss_cor$statistic, 2), "*", sep = ""), parse = FALSE, size = 6)


fig4_beta_diversity_bin <- beta_div_fig(beta_byLake, y_var = "micro_horn_binary_median", x_var = "inverts_combined_horn_binary_median",
                                          error_bars = TRUE, y_upper = "micro_horn_binary_upper", y_lower = "micro_horn_binary_lower",
                                          x_upper = "inverts_combined_horn_binary_upper", x_lower = "inverts_combined_horn_binary_lower",
                                          y_lab = "",
                                          x_lab = "",
                                          x_min = 0.70,
                                          y_min = 0.70,
                                          median_lines = FALSE,
                                          text_v = c(-0.7, -0.7, -0.7, 1.7, 1.7, -0.7, 1.7, -0.7, -0.7, -0.7, 1.7, 1.7), 
                                          text_h = c(1.5, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, 1.5, 1.5, -0.3, -0.3, -0.3)
)

##### Abundance-based dissimilarity
#### Plot correlation between pairwise distances
fig4_pairwise_diss_abun <- ggplot(pairwise_horn_dist, aes(y = micro_horn_dissimilarity, x = inverts_combined_horn_dissimilarity, group = comparison_type, color = comparison_type)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(grey(.7), "green", "deepskyblue3")) +
  xlim(0, 1) +
  ylim(0, 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.25, color = grey(0.35)) +
  geom_smooth(aes(y = micro_horn_dissimilarity, x = inverts_combined_horn_dissimilarity), pairwise_horn_dist, method=lm, se = FALSE, color="black", inherit.aes = FALSE) +  # Add linear regression line
  theme_classic() +
  theme(legend.position = "none", text = element_text(size=17)) + 
  ylab("Microbial dissimilarity") +
  xlab("Macro-invertebrate dissimilarity") +
  annotate("text", 0.80, 0.10, label = paste("r = ", round(abun_diss_cor$statistic, 2), "*", sep = ""), parse = FALSE, size = 6)

#### Morisita-Horn medians +- 50%
fig4_beta_diversity_abun <- beta_div_fig(beta_byLake, y_var = "micro_horn_median", x_var = "inverts_combined_horn_median",
                                               error_bars = TRUE, y_upper = "micro_horn_upper", y_lower = "micro_horn_lower",
                                               x_upper = "inverts_combined_horn_upper", x_lower = "inverts_combined_horn_lower",
                                               y_lab = "",
                                               x_lab = "Macro-invertebrate dissimilarity",
                                               x_min = 0.80,
                                               y_min = 0.50,
                                               median_lines = FALSE,
                                               text_v = c(2, -0.7, -0.3, 1.7, 2.2, 0, 0, 1.7, -0.7, -0.7, 1.7, 0), 
                                               text_h = c(-0.3, 1.5, -0.3, 1.5, 0.3, -0.3, -0.6, -0.3, 1.5, -0.3, -0.3, -0.3)
)

#### Combine panels
tiff("figures/fig4-beta_diversity.tiff", width = 12.0, height = 12.0, units = 'in', res = 600)
ggarrange(fig4_pairwise_diss_bin, fig4_beta_diversity_bin,
          fig4_pairwise_diss_abun, fig4_beta_diversity_abun,
          labels = c("A", "B", "C", "D"), font.label = list(size = 18, face = "bold"))
dev.off()

##### -- Figure 5: RDA and Procrustes analysis -- #####
tiff("figures/fig5-RDA_analysis.tiff", width = 14.0, height = 9.0, units = 'in', res = 600)
par(mfrow = c(2, 3), mar = c(2, 2, 1, 2))
#### Panel A:
plot(micro_beta_rda_horn_binary_reduced, type = "n", xlim = c(-2, 2), ylim = c(-2, 2))
text(micro_beta_rda_horn_binary_reduced, labels = c("", "", "", ""), display = "bp", head.arrow = 0.08, col = grey(.3), cex = 1.4)
text(c(0, 1.8, 1.85, -1), c(-1.5, -1, 0.3, -1.6), c("Temperature \nmedian", "Salinity", "Oxygen \nmedian", "Surface \narea"), col = grey(.3), cex = 1.4, lwd = 2)
lake_col <- c("green", "deepskyblue3", rep("green", 3), rep("deepskyblue3", 3), rep("green", 3), "deepskyblue3")
points(micro_beta_rda_horn_binary_reduced, pch = 16, cex = 1.6, col = lake_col)
text(scores(micro_beta_rda_horn_binary_reduced, display = "sites") - c(0.3, -0.3, rep(0.3, 4), -0.3, rep(0.3, 2), 0.3, -0.3, -0.3, rep(0, 12)), 
     labels = dimnames(scores(micro_beta_rda_horn_binary_reduced, display = "sites"))[[1]], cex = 1.4)
text(-2, 1.9, "A", cex = 2.7)
#### Panel B:
plot(inverts_combined_beta_rda_horn_binary_reduced, type = "n", xlim = c(-2, 2), ylim = c(2, -2))
text(inverts_combined_beta_rda_horn_binary_reduced, labels = c("", "", "", ""), display = "bp", head.arrow = 0.08, col = grey(.3), cex = 1.4)
text(c(0.2, 1.8, 1.8, -1.75), c(1.6, 0.65, -0.8, 0.8), c("Temperature \nmedian", "Salinity", "Oxygen \nmedian", "Oxygen \nvariation"), col = grey(.3), cex = 1.4)
lake_col <- c("green", "deepskyblue3", rep("green", 3), rep("deepskyblue3", 3), rep("green", 3), "deepskyblue3")
points(inverts_combined_beta_rda_horn_binary_reduced, pch = 16, cex = 1.6, col = lake_col)
text(scores(inverts_combined_beta_rda_horn_binary_reduced, display = "sites") - c(0.3, -0.3, 0.3, -0.3, rep(-0.3, 2), 0.3, -0.3, rep(0.3, 3), 0.3, rep(0, 9), -0.2, rep(0, 2)), 
     labels = dimnames(scores(inverts_combined_beta_rda_horn_binary_reduced, display = "sites"))[[1]], cex = 1.4)
text(-2, -1.9, "B", cex = 2.7)
#### Panel C:
pro_test <- protest(micro_beta_rda_horn_binary_reduced, inverts_combined_beta_rda_horn_binary_reduced, permutations = 9999)
plot(pro_test, main = "", xlab = "CAP1", ylab = "CAP2")
text(pro_test$X, labels = shared_lake_codes, adj = 1.2, col = lake_col, cex = 1.4)
text(-0.41, 0.52, "C", cex = 2.2)
text(0.25, 0.3, paste("Rotation r = ", round(pro_test$scale, 3), sep = ""), cex = 1.6)
#### Panel D:
plot(micro_beta_rda_horn_reduced, type = "n", xlim = c(-2, 2), ylim = c(-2, 2))
text(micro_beta_rda_horn_reduced, labels = c("", "", "", ""), display = "bp", head.arrow = 0.08, col = grey(.3), cex = 1.4)
text(c(0, 1.7, 1.8, -1.1), c(-1.4, -1.2, -0.05, -1.65), c("Temperature \nmedian", "Salinity", "Oxygen \nmedian", "Surface \narea"), col = grey(.3), cex = 1.4)
lake_col <- c("green", "deepskyblue3", rep("green", 3), rep("deepskyblue3", 3), rep("green", 3), "deepskyblue3")
points(micro_beta_rda_horn_reduced, pch = 16, cex = 1.6, col = lake_col)
text(scores(micro_beta_rda_horn_reduced, display = "sites") - c(0.3, -0.3, 0.3, 0.3, 0.25, rep(-0.3, 3), rep(0.3, 3), -0.3, rep(0, 4), 0.05, rep(0, 7)),
     labels = dimnames(scores(micro_beta_rda_horn_reduced, display = "sites"))[[1]], cex = 1.4)
text(-2, 1.9, "D", cex = 2.7)
#### Panel E:
plot(inverts_combined_beta_rda_horn_reduced, type = "n", xlim = c(-2, 2), ylim = c(-2, 2))
text(inverts_combined_beta_rda_horn_reduced, labels = c("", "", "", ""), display = "bp", head.arrow = 0.08, col = grey(.3), cex = 1.4)
text(c(0.95, 1.88, -1.7, 0.4), c(-1.5, -0.5, 0.7, 0.6), c("Salinity", "Oxygen \nmedian", "Oxygen \nvariation", "Radiation \nvariation"), col = grey(.3), cex = 1.4)
lake_col <- c("green", "deepskyblue3", rep("green", 3), rep("deepskyblue3", 3), rep("green", 3), "deepskyblue3")
points(inverts_combined_beta_rda_horn_reduced, pch = 16, cex = 1.6, col = lake_col)
text(scores(inverts_combined_beta_rda_horn_reduced, display = "sites") - c(0.3, -0.25, 0.3, 0.13, 0.3, -0.3, 0.25, -0.3, -0.3, 0.3, -0.25, -0.20, 0, 0.10, 0, -0.23, rep(0, 2), -0.08, 0, -0.15, 0, 0.10, 0.10),
     labels = dimnames(scores(inverts_combined_beta_rda_horn_reduced, display = "sites"))[[1]], cex = 1.4)
text(-2, 1.9, "E", cex = 2.7)
#### Panel F:
pro_test2 <- protest(micro_beta_rda_horn_reduced, inverts_combined_beta_rda_horn_reduced, permutations = 9999)
plot(pro_test2, main = "", xlab = "CAP1", ylab = "CAP2")
text(pro_test2$X, labels = shared_lake_codes, adj = 1.2, col = lake_col, cex = 1.4)
text(-0.31, 0.42, "F", cex = 2.2)
text(0.2, 0.3, paste("Rotation r = ", round(pro_test2$scale, 3), sep = ""), cex = 1.6)
dev.off()

#### Correlation in Procrustes test
pro_test
pro_test2

#### Test of difference in median dissimilarity between micro- and macro-organisms
kruskal.test(pairwise_horn_dist$micro_horn_dissimilarity, pairwise_horn_dist$inverts_combined_horn_dissimilarity)
