#########################
###### micro-macro ######
#########################
###### FUNCTIONS ########
#########################
##### -- Data manipulation -- #####
lake_composition <- function(abun_dat, lake, min_abun = 1){
  subset(abun_dat, lake_code == lake)[which(subset(abun_dat, lake_code == lake) >= min_abun)][, -1] %>% 
    sort(decreasing = TRUE)
}

##### -- Alpha diversity of lakes -- #####
#### Calculate Simpson's evenness index from a vector of abundances
simpson_evenness <- function(x, zeroes=TRUE) {
  if (!zeroes) {
    x[x > 0]
  }
  # Species richness (number of species)
  S <- length(x)
  # Relative abundances
  p <- x/sum(x)
  # Simpson index
  lambda <- sum(p^2)
  # Simpson evenness (Simpson diversity per richness)
  (1/lambda)/S
}
#### Calculate Pielou's evenness index from vector of abundances
pielou <- function(x) {
  # Remove zeroes
  x <- x[x > 0]
  # Species richness (number of species)
  S <- length(x)
  # Relative abundances
  p <- x/sum(x)
  # Shannon index
  H <- (-sum(p * log(p)))
  # Simpson evenness
  H/log(S)
}

#### Wrapper for Diversity function in library(SpadeR)
calc_div <- function(abun_data, label = ""){
  abun_data_div <- abun_data[, -1]
  div_list <- vector("list", nrow(abun_data_div))
  for (i in 1:nrow(abun_data_div)){
    print(i)
    div <- Diversity(t(abun_data_div[i, ]), datatype = "abundance", q = c(0, 1))
    div_list[[i]] <- summarize_Diversity_output(div)
    div_list[[i]]$dominance <- max(abun_data_div[i, ]/sum(abun_data_div[i, ]))
    div_list[[i]]$simpson_evenness <- simpson_evenness(abun_data_div[i, ])
    div_list[[i]]$pielou <- pielou(abun_data_div[i, ])
  }
  abun_data_div <- do.call("rbind", div_list)
  abun_data_div <- data.frame(lake_code = abun_data[, 1], abun_data_div)
  names(abun_data_div)[-1] <- paste(label, names(abun_data_div)[-1], sep = "_")
  return(abun_data_div)
}
#### Generate data.frame from output of SpadeR:Diversity 
summarize_Diversity_output <- function(Diversity_output){
  out <- data.frame(observed_richness = as.numeric(as.character(Diversity_output$Basic_data$Value[2])),
             estimated_sample_coverage = as.numeric(as.character(Diversity_output$Basic_data$Value[3])),
             chao1_estimate = as.numeric(as.character(Diversity_output$Species_richness))[1],
             chao1_lower = as.numeric(as.character(Diversity_output$Species_richness))[11],
             chao1_upper = as.numeric(as.character(Diversity_output$Species_richness))[16],
             ichao1_estimate = as.numeric(as.character(Diversity_output$Species_richness))[1+2],
             ichao1_lower = as.numeric(as.character(Diversity_output$Species_richness))[11+2],
             ichao1_upper = as.numeric(as.character(Diversity_output$Species_richness))[16+2],
             ace_estimate = as.numeric(as.character(Diversity_output$Species_richness))[1+3],
             ace_lower = as.numeric(as.character(Diversity_output$Species_richness))[11+3],
             ace_upper = as.numeric(as.character(Diversity_output$Species_richness))[16+3],
             ace1_estimate = as.numeric(as.character(Diversity_output$Species_richness))[1+4],
             ace1_lower = as.numeric(as.character(Diversity_output$Species_richness))[11+4],
             ace1_upper = as.numeric(as.character(Diversity_output$Species_richness))[16+4],
             shannon_diversity_estimate = as.numeric(as.character(Diversity_output$Shannon_diversity))[4],
             shannon_diversity_lower = as.numeric(as.character(Diversity_output$Shannon_diversity))[12],
             shannon_diversity_upper = as.numeric(as.character(Diversity_output$Shannon_diversity))[16],
             simpson_diversity_estimate = as.numeric(as.character(Diversity_output$Simpson_diversity))[2],
             simpson_diversity_lower = as.numeric(as.character(Diversity_output$Simpson_diversity))[6],
             simpson_diversity_upper = as.numeric(as.character(Diversity_output$Simpson_diversity))[8],
             gini_simpson_estimate = 1 - as.numeric(as.character(Diversity_output$Simpson_index))[2],
             gini_simpson_lower = 1 - as.numeric(as.character(Diversity_output$Simpson_index))[6],
             gini_simpson_upper = 1 - as.numeric(as.character(Diversity_output$Simpson_index))[8]
  )
  return(out)
}

#### Wrapper function for diversity measure functions in iNEXT
calc_Chao <- function(abun_data, FUN = ChaoRichness, type = c("abundance", "incidence_freq"), prefix = "chao1", label = "", ...){
  type <- match.arg(type)
  abun_data_div <- abun_data[, -1]
  if (type == "incidence_freq") abun_data_div[abun_data_div > 0] <- 1
  abun_data_div <- apply(abun_data_div, 1, function(x){
    if (sum(x) > 0) FUN(x, datatype = type, ...)
    else data.frame("Observed" = 0, "Estimator" = 0, "Est_s.e." = 0, "95% Lower" = 0, "95% Upper" = 0)
  })
  abun_data_div <- lapply(abun_data_div, function(x){
    names(x) <- paste(prefix, c("Observed", "Estimator", "Est_s.e.", "95%_Lower", "95%_Upper"), sep = "_")
    x })
  abun_data_div <- do.call("rbind", abun_data_div)
  names(abun_data_div) <- paste(label, names(abun_data_div), sep = "_")
  abun_data_div <- data.frame(lake_code = abun_data[, 1], abun_data_div)
  return(abun_data_div)
}

##### Species accumulation curves
sacc_curves <- function(abun_data = macrobiota_.7, lakes = NULL, col_to_rm = 1:3){
  if (is.null(lakes)) lakes <- unique(abun_data$lake_code)
  specaccum_list <- lapply(lakes, function(x){
    out <- subset(abun_data, lake_code == x)[, -col_to_rm] %>% specaccum()
    data.frame(lake_code = x, site = out$sites, richness = out$richness, sd = out$sd)
  })
  specaccum_df <- do.call("rbind", specaccum_list)
  specaccum_df$lake_code <- as.factor(specaccum_df$lake_code)
  return(specaccum_df)
}

#########################
##### -- FIGURES -- #####
#########################
##########################################
##### -- Alpha diversity of lakes -- #####
##########################################
##### alpha_div_fig()
alpha_div_fig <- function(dat = alpha_byLake, y_var, x_var, y_min = NULL, y_max = NULL, x_min = NULL, x_max = NULL, fit_line = FALSE,
                          error_bars = FALSE, y_upper = NULL, y_lower = NULL, x_upper = NULL, x_lower = NULL,
                          y_lab = "Microbial diversity", x_lab = "Macrobial diversity", text_v = 1.7, text_h = -0.3){
  if (is.null(y_min)) y_min <- min(dat[, y_var], na.rm = TRUE)
  if (is.null(y_max)) y_max <- max(dat[, y_var], na.rm = TRUE)
  if (is.null(x_min)) x_min <- min(dat[, x_var], na.rm = TRUE)
  if (is.null(x_max)) x_max <- max(dat[, x_var], na.rm = TRUE)
  p <- ggplot(dat, aes_string(y = y_var, x = x_var, group = "stratified", color = "stratified")) +
    geom_point(size = 3) +
    scale_y_continuous(labels=function(n){format(n, scientific = FALSE)}, limits = c(y_min, y_max)) +
    scale_color_manual(values = c("deepskyblue3", "green")) +
    geom_text(aes(label=lake_code), size = 4, vjust = text_v, hjust = text_h, col = "black") +
    xlim(x_min, x_max) +
    theme_classic() +
    theme(legend.position = "none", text = element_text(size=18), 
          axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
    ylab(y_lab) +
    xlab(x_lab)
  if (isTRUE(fit_line)) p <- p + geom_smooth(method=lm, se=FALSE)
  if (isTRUE(error_bars)) {
    p <- p + geom_errorbar(aes_string(ymax = y_upper, ymin = y_lower), size = 0.2, width = 0) +
      geom_errorbarh(aes_string(xmax = x_upper, xmin = x_lower), size = 0.2, height = 0)
  }
  p
  return(p)
}
############################################
##### -- Beta-diversity among lakes -- #####
############################################
##### beta_div_fig()
beta_div_fig <- function(dat, y_var, x_var, y_min = 0, y_max = 1.01, x_min = 0, x_max = 1.01, point_labels = TRUE, median_lines = TRUE,
                         error_bars = FALSE, y_upper = NULL, y_lower = NULL, x_upper = NULL, x_lower = NULL,
                         y_lab = "Beta diversity", x_lab = "Beta diversity", text_v = 1.7, text_h = -0.3){
  p <- ggplot(dat, aes_string(y = y_var, x = x_var, color = "stratified")) +
    geom_point(size = 10* (abs(dat[, y_var] - dat[, x_var]))) +
    scale_color_manual(values = c("deepskyblue3", "green")) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 0.25, color = grey(0.35)) +
    xlim(x_min, x_max) +
    ylim(y_min, y_max) +
    theme_classic() +
    theme(legend.position = "none", text = element_text(size=17)) +
    ylab(y_lab) +
    xlab(x_lab)
    if (isTRUE(error_bars)) {
      p <- p + geom_errorbar(aes_string(ymax = y_upper, ymin = y_lower), size = 0.2, width = 0) +
        geom_errorbarh(aes_string(xmax = x_upper, xmin = x_lower), size = 0.2, height = 0)
    }
    if (isTRUE(point_labels)) {
      p <- p + geom_text(aes(label=lake_code), color = "black", size = 3, vjust = text_v, hjust = text_h)
    }
  if (isTRUE(median_lines)) {
    p <- p + geom_hline(yintercept = median(dat[, y_var]), color = grey(0.75), size = 0.3, linetype = "solid") +
      geom_vline(xintercept = median(dat[, x_var]), color = grey(0.75), size = 0.3, linetype = "solid") 
  } 
  p
  return(p)
}

##############################################################
##### -- Alpha-diversity as a function of environment -- #####
##############################################################
##### For continuous variables
div_byEnv <- function(dat, response, level = .5, type = "spearman", predictor = NULL, model_type = c("lm", "gam"), y_lab = "Microbial diversity", x_lab = "", x_min = NULL, x_max = NULL, y_min = NULL, y_max = NULL, annot = c(Inf, Inf), text_v = 1.7, text_h = -0.3, ...){
  model_type <- match.arg(model_type)
  high_cor <- NULL
  if (is.null(predictor)){
    high_cor <- which(abs(cor(dat[c(response, names(dat)[(which(names(dat) == "stratified")+1):ncol(dat)])], method = type)[1, -1]) >= level)
    if (length(high_cor) > 0){
      for (i in 1:length(high_cor)){
        p <- ggplot(dat,
                    aes_string(y = response, x = names(high_cor)[i])) +
          geom_smooth(method=lm, se = FALSE, color="black") +  # Add linear regression lines +
          geom_point(size = 3) +
          theme_classic() +
          theme(legend.position = "none") +    
          ylab(y_lab) +
          xlab(names(high_cor)[i])
        print(p)
      }
    }
  } else {
    if (is.null(y_min)) y_min <- min(dat[, response], na.rm = TRUE)
    if (is.null(y_max)) y_max <- max(dat[, response], na.rm = TRUE)
    if (is.null(x_min)) x_min <- min(dat[, predictor], na.rm = TRUE)
    if (is.null(x_max)) x_max <- max(dat[, predictor], na.rm = TRUE)
    p <- ggplot(dat,
                aes_string(y = response, x = predictor, color = "stratified")) +
      geom_point(size = 3) +
      scale_color_manual(values = c("deepskyblue3", "green")) +
      #geom_text(aes(label=lake_code), size = 4, vjust = text_v, hjust = text_h, col = "black") +
      #scale_y_continuous(labels = comma) +
      xlim(x_min, x_max) +
      ylim(y_min, y_max) +
      theme_classic() +
      theme(legend.position = "none", text = element_text(size=15), ...) +
      ylab(y_lab) +
      xlab(x_lab)
    if (model_type == "gam"){
      mod <- gam(as.formula(paste(response, " ~ s(", predictor, ", k = 3)", sep = "")), data = dat)
      p <- p + geom_smooth(aes_string(y = response, x = predictor), method="gam", formula = y ~ s(x, k = 3), se = TRUE, inherit.aes = FALSE, color = "black") + 
        annotate("text", annot[1], annot[2], label = paste("R^2 == ", round(summary(mod)$r.sq, 2)), parse = TRUE, size = 5)
    } else { 
      #mod <- lm(as.formula(paste(response, " ~ ", predictor, sep = "")), data = dat)
      mod_cor <- cor.test(dat[, response], dat[, predictor])
      if (mod_cor$p.value < 0.05) { annot_label <- paste("r = ", round(mod_cor$estimate, 2), "*", sep = "") }
      else {annot_label <- paste("r = ", round(mod_cor$estimate, 2))}
      p <- p + geom_smooth(aes_string(y = response, x = predictor), method="lm", formula = y ~ x, se = FALSE, inherit.aes = FALSE, color = "black") + 
        annotate("text", annot[1], annot[2], label = annot_label, parse = FALSE, size = 5)
    }
  }
  if (length(high_cor) > 0) return(high_cor)
  else return(p)
}

##### For factor variables
div_byEnv_factor <- function(dat, response, predictor, y_lab = "Microbial diversity", x_lab = "", annot = c(Inf, Inf), y_min = NULL, y_max = NULL, ...){
  mod <- lm(as.formula(paste(response, " ~ ", predictor, sep = "")), data = dat)
  if (is.null(y_min)) y_min <- min(dat[, response], na.rm = TRUE)
  if (is.null(y_max)) y_max <- max(dat[, response], na.rm = TRUE)
  p <- ggplot(dat,
              aes_string(y = response, x = predictor, color = "stratified")) +
    geom_boxplot() +
    scale_color_manual(values = c("deepskyblue3", "green")) +
    scale_x_discrete(labels = c("mixed", "stratified")) +
    theme_classic() +
    ylim(y_min, y_max) +
    theme(legend.position = "none", text = element_text(size=15), ...) +
    ylab(y_lab) +
    xlab(x_lab) +
    annotate("text", annot[1], annot[2], label = paste("R^2 == ", round(summary(mod)$r.sq, 2)), parse = TRUE, size = 5)
  if (summary(mod)$coefficients[2, 4] < 0.05) p + annotate("text", annot[1]+0.1, annot[2], label = "*", size = 5) 
  return(p)
}

