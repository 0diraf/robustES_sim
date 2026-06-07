library(boot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(knitr)

n <- 10000        # Sample size
m <- floor(n^0.7) # Subsample size
R <- 10000        # Bootstrap Replicates
set.seed(1)

#### Helpers ####

#### Estimators and Data Generation Helpers ####

calc_cohens_d <- function(g1, g2) {
  n1 <- length(g1)
  n2 <- length(g2)
  
  mean_diff <- mean(g2) - mean(g1)
  
  sp <- sqrt(
    ((n1 - 1) * var(g1) + (n2 - 1) * var(g2)) /
      (n1 + n2 - 2)
  )
  
  return(mean_diff / sp)
}

calc_delta_mad <- function(g1, g2) {
  diff_med <- abs(median(g1) - median(g2))
  n1 <- length(g1)
  n2 <- length(g2)
  
  pooled <- ((n1 - 1) * mad(g1) + (n2 - 1) * mad(g2))/(n1+n2-2)
  return(diff_med/pooled)
}

# calc_hedges_g <- function(g1, g2) {
#   d <- calc_cohens_d(g1, g2)
#   
#   n1 <- length(g1)
#   n2 <- length(g2)
#   df <- n1 + n2 - 2
#   J <- 1 - (3 / (4 * df - 1))
#   g <- d * J
#   
#   return(g)
# }

calc_robust_cohens_d <- function(g1, g2, trim = 0.2) {
  n1 <- length(g1); n2 <- length(g2)
  
  m1 <- mean(g1, trim = trim)
  m2 <- mean(g2, trim = trim)
  
  wvar <- function(x, t) {
    lo <- quantile(x, t, names=FALSE)
    hi <- quantile(x, 1-t, names=FALSE)
    x[x<lo] <- lo; x[x>hi] <- hi
    var(x)
  }
  v1 <- wvar(g1, trim)
  v2 <- wvar(g2, trim)
  
  spw <- sqrt(((n1 - 1)*v1 + (n2 - 1)*v2) / (n1 + n2 - 2))
  
  return((m2 - m1) / (spw / 0.642)) # 0.642 is the scale correction for 20% trimmed mean
}

rpareto <- function(n, alpha, xm, double = FALSE) {
  if (n == 0) return(numeric(0))
  
  r <- runif(n, 0, 1)
  mag <- (xm * (r)^(-1 / (alpha)))
  
  if (double) {
    signs <- sample(c(-1,1), n, replace=TRUE)
    return(signs * mag)
  } else {
    return(mag)
  }
}

normal_contaminated <- function(n, contamination_rate, dist_type = "pareto", alpha = 1.5,
                                mean = 0, sd = 1, bulk_sd = 1, xmin = 1, double = FALSE) {

  n_contaminated <- rbinom(1, n, contamination_rate)
  n_normal <- n - n_contaminated

  # bulk_sd controls the normal bulk; sd controls only the contaminant's spread
  # (the lognormal sdlog). Decoupling these keeps the treatment bulk at N(mean, bulk_sd)
  # while the lognormal sweep varies sigma_log, instead of the two moving together.
  data_normal <- rnorm(n_normal, mean = mean, sd = bulk_sd)

  if (dist_type == "pareto") {
    data_heavy <- rpareto(n_contaminated, xm = xmin, alpha = alpha, double = double)
  } else if (dist_type == "lognormal") {
    data_heavy <- rlnorm(n_contaminated, meanlog = 0, sdlog = sd)
  }
  
  combined_data <- sample(c(data_normal, data_heavy))
  return(combined_data)
}


#### Contamination Simulation Helpers ####

contamination_sweep <- function(n, alpha, n_iter = 2000, cmean = 0, csd = 1,
                                tmean = 0.5, tsd = 1, double=FALSE, treat_only = FALSE, dist_type='pareto') {
  
  fracs <- seq(0.01, 0.10, by = 0.01)
  results <- data.frame()
  
  for (f in fracs) {
    d_vals   <- numeric(n_iter)
    mad_vals <- numeric(n_iter)
    wd_vals  <- numeric(n_iter)
    # g_vals <- numeric(n_iter)  # Hedges' g removed
    
    for(i in 1:n_iter) {
      if (treat_only) {
        control <- rnorm(n, mean = cmean, sd = csd)
        treat <- normal_contaminated(n, contamination_rate = f, alpha = alpha, double = double,
                                     mean=tmean, sd=tsd, dist_type=dist_type)
      } else {
        control <- normal_contaminated(n, contamination_rate=f, mean=cmean, sd=csd, alpha=alpha, double=double, dist_type=dist_type)
        treat <- normal_contaminated(n, contamination_rate=f, mean=tmean, sd=tsd, alpha=alpha, double=double, dist_type=dist_type)
      }
      
      d_vals[i]   <- calc_cohens_d(control, treat)
      mad_vals[i] <- calc_delta_mad(control, treat)
      # g_vals[i] <- calc_hedges_g(control, treat)  # Hedges' g removed
      wd_vals[i]  <- calc_robust_cohens_d(control, treat)
    }
    
    results <- rbind(results, data.frame(
      Fraction      = f,
      Cohen_d       = mean(d_vals, na.rm=TRUE),
      Cohen_d_sd    = sd(d_vals, na.rm=TRUE),
      Delta_MAD     = mean(mad_vals, na.rm=TRUE),
      Delta_MAD_sd  = sd(mad_vals, na.rm=TRUE),
      # Hedges_g    = mean(g_vals, na.rm=TRUE),  #Hedges' g removed
      Winsorized_d     = mean(wd_vals, na.rm=TRUE),
      Winsorized_d_sd  = sd(wd_vals, na.rm=TRUE)
    ))
  }
  return(results)
}

run_alpha_sweep <- function(n=2000, n_iter=1000, cmean=0, csd=1,
                            tmean=0.5, tsd=1, double=FALSE, dist_type='pareto') {
  alphas <- c(0.5, 1.5, 2.5)
  sds <- c(1, 1.8, 3)
  all_results <- data.frame()
  
  if (dist_type=='pareto') {
    print(paste("Starting Alpha Sweep"))
    
    for (a in alphas) {
      cat(paste("\nRunning Alpha =", a, "... "))
      
      res <- contamination_sweep(n=n, alpha=a, n_iter=n_iter, double=double, treat_only=TRUE, cmean=cmean, tmean=tmean, csd = csd, tsd = tsd, dist_type=dist_type)
      
      res$Alpha <- as.factor(a)
      all_results <- rbind(all_results, res)
    }
  } else {
    print(paste("Starting Sweep"))
    
    for (s in sds) {
      cat(paste("\nRunning SD =", s, "... "))
      
      res <- contamination_sweep(n=n, alpha=0.5, n_iter=n_iter, double=double, treat_only=TRUE, 
                                 cmean=cmean, tmean=tmean, csd = csd, tsd = s, dist_type=dist_type)
      
      res$SD <- as.factor(s)
      all_results <- rbind(all_results, res)
    }
  }
  return(all_results)
}

#### Summary Table Helper ####

# Computes offset (mean estimator - true_delta) and SD at selected contamination levels.

bias_table <- function(sweep_result, true_delta, group_var, f_vals = c(0.05, 0.10)) {
  sweep_result %>%
    filter(round(Fraction, 10) %in% round(f_vals, 10)) %>%
    mutate(
      Bias_Cohen_d    = round(Cohen_d - true_delta, 4),
      Bias_Delta_MAD  = round(Delta_MAD - true_delta, 4),
      Bias_Winsorized = round(Winsorized_d - true_delta, 4),
      SD_Cohen_d      = round(Cohen_d_sd, 4),
      SD_Delta_MAD    = round(Delta_MAD_sd, 4),
      SD_Winsorized   = round(Winsorized_d_sd, 4)
    ) %>%
    select(all_of(group_var), Fraction,
           Bias_Cohen_d, SD_Cohen_d,
           Bias_Delta_MAD, SD_Delta_MAD,
           Bias_Winsorized, SD_Winsorized)
}


##### Exps ####


##### CONTAMINATION EFFECT SIMULATIONS ####


##### NORMAL + PARETO ####
#### Small ####
sweep_results_poss <- run_alpha_sweep(n=n, n_iter=10000, cmean=0, tmean=0.2, double=F)

sweep_long_poss <- sweep_results_poss %>%
  pivot_longer(cols = c("Cohen_d", "Delta_MAD", "Winsorized_d"), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         "Cohen_d" = "Cohen's d",
                         "Delta_MAD" = "Blaine's d",
                         "Winsorized_d" = "AKP robust d"))

p1_poss <- ggplot(sweep_long_poss, aes(x = Fraction, y = Value, color = Metric)) +
  
  geom_line(aes(linetype = Metric), linewidth = 1.2, alpha = 0.8) + 
  
  geom_point(aes(shape = Metric), size = 3) +
  
  geom_hline(yintercept = 0.2, linetype = "dotted", color = "gray20", linewidth = 1) +
  
  facet_wrap(~Alpha, labeller = label_both) + 
  
  scale_color_manual(values = c("Cohen's d" = "black", 
                                "Blaine's d" = "#1f78b4",  
                                "AKP robust d" = "#d95f02")) + 
  
  scale_shape_manual(values = c("Cohen's d" = 16,
                                "Blaine's d" = 15,
                                "AKP robust d" = 18)) +
  
  scale_linetype_manual(values = c("Cohen's d" = "solid", 
                                   "Blaine's d" = "solid", 
                                   "AKP robust d" = "dotdash")) +
  
  labs(y = "Estimated Effect Size",
       x = "Contamination Fraction") +
  
  theme_classic(base_size = 16, base_family = "sans") + 
  theme(
    legend.position = "bottom",          
    legend.title = element_blank(),           
    
    panel.grid.major = element_blank(),       
    panel.grid.minor = element_blank(),       
    
    strip.background = element_blank(),       
    strip.text = element_text(size = 12, face = "bold"),
    
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 11, color = "black"), 
    #axis.line = element_line(color = "black"),             
    
    axis.line = element_blank(),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5), 
    panel.spacing = unit(1.5, "lines") 
  )


print(p1_poss)

# kable(
#   bias_table(sweep_results_poss, true_delta = 0.2, group_var = "Alpha", f_vals = c(0.01, 0.05, 0.10)) %>%
#     rename("Bias (Cohen's d)" = Bias_Cohen_d, "SD (Cohen's d)" = SD_Cohen_d,
#            "Bias (Blaine's d)" = Bias_Delta_MAD, "SD (Blaine's d)" = SD_Delta_MAD,
#            "Bias (AKP robust d)" = Bias_Winsorized, "SD (AKP robust d)" = SD_Winsorized),
#   format = "simple",
#   caption = "Bias and SD of estimators under Pareto contamination (small effect, delta = 0.2)"
# )
# 
# kable(bias_table(sweep_results_poss, true_delta = 0.2, group_var = "Alpha", f_vals = c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11)) %>%
#         rename("Bias (Cohen's d)" = Bias_Cohen_d, "SD (Cohen's d)" = SD_Cohen_d,
#                "Bias (Blaine's d)" = Bias_Delta_MAD, "SD (Blaine's d)" = SD_Delta_MAD,
#                "Bias (AKP robust d)" = Bias_Winsorized, "SD (AKP robust d)" = SD_Winsorized),
#       format = "simple",
#       caption = "Bias and SD of estimators under Pareto contamination (small effect, delta = 0.2)")

#### Medium ####

sweep_results_posm <- run_alpha_sweep(n=n, n_iter=10000, cmean=0, tmean=0.5, double=F)

sweep_long_posm <- sweep_results_posm  %>%
  pivot_longer(cols = c("Cohen_d", "Delta_MAD", "Winsorized_d"), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         "Cohen_d" = "Cohen's d",
                         "Delta_MAD" = "Blaine's d",
                         "Winsorized_d" = "AKP robust d"))

p1_posm <- ggplot(sweep_long_posm, aes(x = Fraction, y = Value, color = Metric)) +
  
  geom_line(aes(linetype = Metric), linewidth = 1.2, alpha = 0.8) + 
  
  geom_point(aes(shape = Metric), size = 3) +
  
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "gray20", linewidth = 1) +
  
  facet_wrap(~Alpha, labeller = label_both) + 
  
  scale_color_manual(values = c("Cohen's d" = "black", 
                                "Blaine's d" = "#1f78b4",  
                                "AKP robust d" = "#d95f02")) + 
  
  scale_shape_manual(values = c("Cohen's d" = 16,
                                "Blaine's d" = 15,
                                "AKP robust d" = 18)) +
  
  scale_linetype_manual(values = c("Cohen's d" = "solid", 
                                   "Blaine's d" = "solid", 
                                   "AKP robust d" = "dotdash")) +
  
  labs(y = "Estimated Effect Size",
       x = "Contamination Fraction") +
  
  theme_classic(base_size = 16, base_family = "sans") + 
  theme(
    legend.position = "bottom",          
    legend.title = element_blank(),           
    legend.text = element_text(size = 14),    
    
    panel.grid.major = element_blank(),       
    panel.grid.minor = element_blank(),       
    
    strip.background = element_blank(),       
    strip.text = element_text(size = 16, face = "bold"), 
    
    axis.title = element_text(size = 16),                
    axis.text = element_text(size = 14, color = "black"),
    
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.spacing = unit(1.5, "lines") 
  )

print(p1_posm)

#### Large ####

sweep_results_posl <- run_alpha_sweep(n=n, n_iter=10000, cmean=0, tmean=0.8, double=F)

sweep_long_posl <- sweep_results_posl %>%
  pivot_longer(cols = c("Cohen_d", "Delta_MAD", "Winsorized_d"), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         "Cohen_d" = "Cohen's d",
                         "Delta_MAD" = "Blaine's d",
                         "Winsorized_d" = "AKP robust d"))

p1_posl <- ggplot(sweep_long_posl, aes(x = Fraction, y = Value, color = Metric)) +
  
  geom_line(aes(linetype = Metric), linewidth = 1.2, alpha = 0.8) + 
  
  geom_point(aes(shape = Metric), size = 3) +
  
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray20", linewidth = 1) +
  
  facet_wrap(~Alpha, labeller = label_both) + 
  
  scale_color_manual(values = c("Cohen's d" = "black", 
                                "Blaine's d" = "#1f78b4",  
                                "AKP robust d" = "#d95f02")) + 
  
  scale_shape_manual(values = c("Cohen's d" = 16,
                                "Blaine's d" = 15,
                                "AKP robust d" = 18)) +
  
  scale_linetype_manual(values = c("Cohen's d" = "solid", 
                                   "Blaine's d" = "solid", 
                                   "AKP robust d" = "dotdash")) +
  
  labs(y = "Estimated Effect Size",
       x = "Contamination Fraction") +
  
  theme_classic(base_size = 16, base_family = "sans") + 
  theme(
    legend.position = "bottom",          
    legend.title = element_blank(),           
    legend.text = element_text(size = 14),    
    
    panel.grid.major = element_blank(),       
    panel.grid.minor = element_blank(),       
    
    strip.background = element_blank(),       
    strip.text = element_text(size = 16, face = "bold"), 
    
    axis.title = element_text(size = 16),                
    axis.text = element_text(size = 14, color = "black"),
    
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.spacing = unit(1.5, "lines") 
  )

print(p1_posl)

pareto_table <- bind_rows(
  bias_table(sweep_results_poss, true_delta = 0.2, group_var = "Alpha", f_vals = c(0.01, 0.05, 0.10)) %>%
    mutate(Effect = "Small (δ = 0.2)"),
  bias_table(sweep_results_posm, true_delta = 0.5, group_var = "Alpha", f_vals = c(0.01, 0.05, 0.10)) %>%
    mutate(Effect = "Medium (δ = 0.5)"),
  bias_table(sweep_results_posl, true_delta = 0.8, group_var = "Alpha", f_vals = c(0.01, 0.05, 0.10)) %>%
    mutate(Effect = "Large (δ = 0.8)")
) %>%
  select(Effect, Alpha, Fraction, everything()) %>%
  rename(
    "Offset (Cohen's d)"    = Bias_Cohen_d,    "SD (Cohen's d)"    = SD_Cohen_d,
    "Offset (Blaine's d)"   = Bias_Delta_MAD,  "SD (Blaine's d)"   = SD_Delta_MAD,
    "Offset (AKP robust d)" = Bias_Winsorized,  "SD (AKP robust d)" = SD_Winsorized
  )

kable(pareto_table, format = "simple", digits = 4,
      caption = "Offset and SD of estimators under Pareto contamination at f = 0.01, 0.05, 0.10 across small, medium, and large effects")



#### NORMAL + LOGNORMAL ####
#### Small ####

lsweep_results_poss <- run_alpha_sweep(n=n, n_iter=10000, cmean=0, 
                                       tmean=0.2, double=F, dist_type = 'lognormal')

lsweep_long_poss <- lsweep_results_poss %>%  
  pivot_longer(cols = c("Cohen_d", "Delta_MAD", "Winsorized_d"), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         "Cohen_d" = "Cohen's d",
                         "Delta_MAD" = "Blaine's d",
                         "Winsorized_d" = "AKP robust d"))

lp1_poss <- ggplot(lsweep_long_poss, aes(x = Fraction, y = Value, color = Metric)) +
  
  geom_line(aes(linetype = Metric), linewidth = 1.2, alpha = 0.8) + 
  
  geom_point(aes(shape = Metric), size = 3) +
  
  geom_hline(yintercept = 0.2, linetype = "dotted", color = "gray20", linewidth = 1) +
  
  facet_wrap(~SD, labeller = label_both) + 
  
  scale_color_manual(values = c("Cohen's d" = "black", 
                                "Blaine's d" = "#1f78b4",  
                                "AKP robust d" = "#d95f02")) + 
  
  scale_shape_manual(values = c("Cohen's d" = 16,
                                "Blaine's d" = 15,
                                "AKP robust d" = 18)) +
  
  scale_linetype_manual(values = c("Cohen's d" = "solid", 
                                   "Blaine's d" = "solid", 
                                   "AKP robust d" = "dotdash")) +
  
  labs(y = "Estimated Effect Size",
       x = "Contamination Fraction") +
  
  theme_classic(base_size = 16, base_family = "sans") + 
  theme(
    legend.position = "bottom",          
    legend.title = element_blank(),           
    legend.text = element_text(size = 14),    
    
    panel.grid.major = element_blank(),       
    panel.grid.minor = element_blank(),       
    
    strip.background = element_blank(),       
    strip.text = element_text(size = 16, face = "bold"), 
    
    axis.title = element_text(size = 16),                
    axis.text = element_text(size = 14, color = "black"),
    
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.spacing = unit(1.5, "lines") 
  )
print(lp1_poss)

#### Medium ####

lsweep_results_posm <- run_alpha_sweep(n=n, n_iter=10000, cmean=0, tmean=0.5, double=F, dist_type='lognormal')

lsweep_long_posm <- lsweep_results_posm %>%
  pivot_longer(cols = c("Cohen_d", "Delta_MAD", "Winsorized_d"), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         "Cohen_d" = "Cohen's d",
                         "Delta_MAD" = "Blaine's d",
                         "Winsorized_d" = "AKP robust d"))

lsweep_long_posm$SD <- factor(lsweep_long_posm$SD, levels = c(3, 1.8, 1))

lp1_posm <- ggplot(lsweep_long_posm, aes(x = Fraction, y = Value, color = Metric)) +
  
  geom_line(aes(linetype = Metric), linewidth = 1.2, alpha = 0.8) + 
  
  geom_point(aes(shape = Metric), size = 3) +
  
  geom_hline(yintercept = 0.5, linetype = "dotted", color = "gray20", linewidth = 1) +
  
  facet_wrap(~SD, labeller = label_both) + 
  
  scale_color_manual(values = c("Cohen's d" = "black", 
                                "Blaine's d" = "#1f78b4",  
                                "AKP robust d" = "#d95f02")) + 
  
  scale_shape_manual(values = c("Cohen's d" = 16,
                                "Blaine's d" = 15,
                                "AKP robust d" = 18)) +
  
  scale_linetype_manual(values = c("Cohen's d" = "solid", 
                                   "Blaine's d" = "solid", 
                                   "AKP robust d" = "dotdash")) +
  
  labs(y = "Estimated Effect Size",
       x = "Contamination Fraction") +
  
  theme_classic(base_size = 16, base_family = "sans") + 
  theme(
    legend.position = "bottom",          
    legend.title = element_blank(),           
    legend.text = element_text(size = 14),    
    
    panel.grid.major = element_blank(),       
    panel.grid.minor = element_blank(),       
    
    strip.background = element_blank(),       
    strip.text = element_text(size = 16, face = "bold"), 
    
    axis.title = element_text(size = 16),                
    axis.text = element_text(size = 14, color = "black"),
    
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.spacing = unit(1.5, "lines") 
  )
print(lp1_posm)

#### Large ####

lsweep_results_posl <- run_alpha_sweep(n=n, n_iter=10000, cmean=0, tmean=0.8, double=F, dist_type='lognormal')

lsweep_long_posl <- lsweep_results_posl %>%
  pivot_longer(cols = c("Cohen_d", "Delta_MAD", "Winsorized_d"), names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(Metric,
                         "Cohen_d" = "Cohen's d",
                         "Delta_MAD" = "Blaine's d",
                         "Winsorized_d" = "AKP robust d"))

lp1_posl <- ggplot(lsweep_long_posl, aes(x = Fraction, y = Value, color = Metric)) +
  
  geom_line(aes(linetype = Metric), linewidth = 1.2, alpha = 0.8) + 
  
  geom_point(aes(shape = Metric), size = 3) +
  
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray20", linewidth = 1) +
  
  facet_wrap(~SD, labeller = label_both) + 
  
  scale_color_manual(values = c("Cohen's d" = "black", 
                                "Blaine's d" = "#1f78b4",  
                                "AKP robust d" = "#d95f02")) + 
  
  scale_shape_manual(values = c("Cohen's d" = 16,
                                "Blaine's d" = 15,
                                "AKP robust d" = 18)) +
  
  scale_linetype_manual(values = c("Cohen's d" = "solid", 
                                   "Blaine's d" = "solid", 
                                   "AKP robust d" = "dotdash")) +
  
  labs(y = "Estimated Effect Size",
       x = "Contamination Fraction") +
  
  theme_classic(base_size = 16, base_family = "sans") + 
  theme(
    legend.position = "bottom",          
    legend.title = element_blank(),           
    legend.text = element_text(size = 14),    
    
    panel.grid.major = element_blank(),       
    panel.grid.minor = element_blank(),       
    
    strip.background = element_blank(),       
    strip.text = element_text(size = 16, face = "bold"), 
    
    axis.title = element_text(size = 16),                
    axis.text = element_text(size = 14, color = "black"),
    
    axis.line = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.spacing = unit(1.5, "lines") 
  )

print(lp1_posl)

lognormal_table <- bind_rows(
  bias_table(lsweep_results_poss, true_delta = 0.2, group_var = "SD", f_vals = c(0.01, 0.05, 0.10)) %>%
    mutate(Effect = "Small (δ = 0.2)"),
  bias_table(lsweep_results_posm, true_delta = 0.5, group_var = "SD", f_vals = c(0.01, 0.05, 0.10)) %>%
    mutate(Effect = "Medium (δ = 0.5)"),
  bias_table(lsweep_results_posl, true_delta = 0.8, group_var = "SD", f_vals = c(0.01, 0.05, 0.10)) %>%
    mutate(Effect = "Large (δ = 0.8)")
) %>%
  select(Effect, SD, Fraction, everything()) %>%
  rename(
    "Offset (Cohen's d)"    = Bias_Cohen_d,    "SD (Cohen's d)"    = SD_Cohen_d,
    "Offset (Blaine's d)"   = Bias_Delta_MAD,  "SD (Blaine's d)"   = SD_Delta_MAD,
    "Offset (AKP robust d)" = Bias_Winsorized,  "SD (AKP robust d)" = SD_Winsorized
  )

kable(lognormal_table, format = "simple", digits = 4,
      caption = "Bias and SD of estimators under lognormal contamination at f = 0.01, 0.05, 0.10 across small, medium, and large effects")
