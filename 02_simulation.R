### Libraries ##################################################################
library(dplyr)
library(tidyr)
library(MASS)
library(ggplot2)
library(patchwork)

okabe_ito <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

# ------------------------------------------------------------------------------
# 1. Social/Family Well-Being (SWB) 0-28, 7 questions
# ------------------------------------------------------------------------------
set.seed(2323)

# Timepoints and SWB specific data
timepoints <- c(0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10)
swb_max <- 28

# Coefficients (Race and ethnicity)
race_coefs <- c(
  "American Indian or Alaska Native" = -0.65,
  "Asian" = -0.7,
  "Black or African American" = -0.6,
  "White" = 0
)

eth_coefs <- c(
  "Hispanic or Latino" = -1.0,
  "Not Hispanic or Latino" = 0
)

# Subset only required variables
req_cols_swb <- c(
  # Required cases.~ variables
  "cases.case_id",
  # Required demographic.~ variables
  "demographic.age_at_index", "demographic.ethnicity", "demographic.race"
)

df_swb1 <- df_final_1 %>%
  dplyr::select(all_of(req_cols_swb))

# Weighted average for missing race
race_counts <- df_swb1 %>%
  filter(demographic.race != "Missing") %>%
  count(demographic.race)

# join counts to coef values
race_counts <- race_counts %>%
  mutate(coef = race_coefs)

# Calculate average effect for race
avg_race_effect <- sum(race_counts$n * race_counts$coef) / sum(race_counts$n)

# Weighted average for missing ethnicity
eth_counts <- df_swb1 %>%
  filter(demographic.ethnicity != "Missing") %>%
  count(demographic.ethnicity)

# join counts to coef values
eth_counts <- eth_counts %>%
  mutate(coef = eth_coefs)

# Calculate average effect for ethnicity
avg_eth_effect <- sum(eth_counts$n * eth_counts$coef) / sum(eth_counts$n)

# Map the values
df_swb2 <- df_swb1 %>%
  mutate(
    race_effect = ifelse(demographic.race == "Missing",
                         avg_race_effect,
                         race_coefs[as.character(demographic.race)]),
    eth_effect = ifelse(demographic.ethnicity == "Missing",
                        avg_eth_effect,
                        eth_coefs[as.character(demographic.ethnicity)])
  )

# Randomise trajectory #########################################################
df_swb3 <- df_swb2 %>%
  mutate(
    trajectory = ifelse(runif(n()) < 0.70, "Linear", "Quadratic")
  )

### Random effects definition ##################################################

# For both linear and quadratic
sd_intercept <- 2 # standard deviation of random intercept
sd_slope <- 0.5 # standard deviation of the average change

# For linear (2x2)
cor_int_slope <- 0.2
Sigma_lin <- matrix(c(sd_intercept^2,
                      cor_int_slope * sd_intercept * sd_slope,
                      cor_int_slope * sd_intercept * sd_slope,
                      sd_slope^2),
                    nrow = 2, byrow = TRUE)

# For quadratic (3x3)
sd_quad <- 0.08
cor_int_quad <- 0.1
cor_slope_quad <- 0.1
Sigma_quad <- matrix(c(
  sd_intercept^2,
  cor_int_slope * sd_intercept * sd_slope,
  cor_int_quad  * sd_intercept * sd_quad,
  
  cor_int_slope * sd_intercept * sd_slope,
  sd_slope^2,
  cor_slope_quad * sd_slope * sd_quad,
  
  cor_int_quad  * sd_intercept * sd_quad,
  cor_slope_quad * sd_slope * sd_quad,
  sd_quad^2
), nrow = 3, byrow = TRUE)

### Random effects assignment ##################################################

N <- nrow(df_swb3)

df_swb3$re_int <- NA
df_swb3$re_slope <- NA
df_swb3$re_quad <- 0  # default for linear subjects

# For linear subjects
idx_lin <- df_swb3$trajectory == "Linear"

re_lin <- mvrnorm(sum(idx_lin), mu = c(0,0), Sigma = Sigma_lin)
df_swb3$re_int[idx_lin] <- re_lin[,1]
df_swb3$re_slope[idx_lin] <- re_lin[,2]

# For quadratic subjects
idx_quad <- df_swb3$trajectory == "Quadratic"

re_quad <- mvrnorm(sum(idx_quad), mu = c(0,0,0), Sigma = Sigma_quad)
df_swb3$re_int[idx_quad] <- re_quad[,1]
df_swb3$re_slope[idx_quad] <- re_quad[,2]
df_swb3$re_quad[idx_quad] <- re_quad[,3]

### Fixed effects ##############################################################

pop_intercept <- 21
age_coef <- -0.02
slope_intercept <- -0.05
age_on_slope <- 0.0005
fixed_quad <- -0.01
resid_sd <- 1.8

# Compute scores and expand to long format #####################################
n_time <- length(timepoints)

req_cols_swb2 <- c(req_cols_swb, "time", "swb_total")

df_swb_long <- df_swb3 %>%
  slice(rep(1:n(), each = n_time)) %>%
  mutate(time = rep(timepoints, times = N)) %>%
  arrange(cases.case_id, time) %>%
  mutate(
    fixed_baseline = pop_intercept + age_coef * demographic.age_at_index +
      race_effect + eth_effect,
    fixed_slope = slope_intercept +
      age_on_slope *
      (demographic.age_at_index - mean(demographic.age_at_index, na.rm = TRUE)),
    latent = fixed_baseline + re_int + (fixed_slope + re_slope) * time +
      (fixed_quad + re_quad) * time^2,
    observed_cont = latent + rnorm(n(), mean = 0, sd = resid_sd),
    observed_raw = round(observed_cont),
    observed_clamped = pmin(pmax(observed_raw, 0), swb_max),
    swb_total = as.integer(observed_clamped)
  ) %>%
  dplyr::select(all_of(req_cols_swb2))

summary_by_time_swb <- df_swb_long %>%
  group_by(time) %>%
  summarise(n = n(), mean_swb = mean(swb_total), sd_swb = sd(swb_total))

baseline <- df_swb_long %>% filter(time == 0)

group_means_race <- baseline %>%
  group_by(demographic.race) %>%
  summarise(n = n(), mean_swb = mean(swb_total), sd_swb = sd(swb_total))
group_means_eth <- baseline %>%
  group_by(demographic.ethnicity) %>%
  summarise(n = n(), mean_swb = mean(swb_total), sd_swb = sd(swb_total))

#print(summary_by_time_swb)
#print(group_means_race)
#print(group_means_eth)

# Summary: mean and sd of swb_total by time and race ###########################
df_swb_long <- df_swb_long %>%
  mutate(
    time_num = as.numeric(as.character(time)) # for plotting
  )

summary_by_stage_time <- df_swb_long %>%
  group_by(time, demographic.race) %>%
  summarise(
    n = n(),
    mean_swb = mean(swb_total, na.rm = TRUE),
    sd_swb = sd(swb_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(time, demographic.race)

# Wide table showing mean (sd)
summary_wide <- summary_by_stage_time %>%
  mutate(mean_sd = sprintf("%.2f (%.2f)", mean_swb, sd_swb)) %>%
  dplyr::select(time, demographic.race, mean_sd) %>%
  pivot_wider(names_from = time, values_from = mean_sd, names_prefix = "time_")


# Plot mean trajectories by stage (with 95% CI)
plot_df <- df_swb_long %>%
  group_by(time_num, demographic.race) %>%
  summarise(
    mean_swb = mean(swb_total, na.rm = TRUE),
    se_swb = sd(swb_total, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

p_swb_race <- ggplot(plot_df, aes(x = time_num, y = mean_swb,
                         color = demographic.race,
                         group = demographic.race)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_swb - 1.96 * se_swb,
                    ymax = mean_swb + 1.96 * se_swb),
                width = 0.08) +
  scale_x_continuous(breaks = sort(unique(plot_df$time_num))) +
  labs(x = "Time (years)", y = "Mean SWB total", color = "Race",
       title = "Mean SWB by Race across Timepoints") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#print(p_swb_race)

# ------------------------------------------------------------------------------
# 2. Emotional Well-Being (EWB) 0-24, 6 questions
# ------------------------------------------------------------------------------

# Timepoints and EWB specific data
timepoints <- c(0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10)
ewb_max <- 24

# Coefficients (Race, ethnicity, stage, metastasis and prior malingancy status)
race_coefs <- c(
  "American Indian or Alaska Native" = -0.4,
  "Asian" = -0.1,
  "Black or African American" = -0.3,
  "White" = 0.0
)

eth_coefs <- c(
  "Hispanic or Latino" = -0.3,
  "Not Hispanic or Latino" = 0.0
)

stage_coefs <- c(
  "I" = 0.8,
  "II" = -2,
  "III" = -3.5,
  "IV" = -6
)

prior_malig_coefs <- c(
  "N" = 0.0,
  "Y" = -1.5
)

# Subset only required variables
req_cols_ewb <- c(
  # Required cases.~ variables
  "cases.case_id",
  # Required demographic.~ variables
  "demographic.age_at_index", "demographic.ethnicity", "demographic.race",
  # Required diagnoses.~ variables
  "diagnoses.ajcc_pathologic_stage", "diagnoses.prior_malignancy"
)

df_ewb1 <- df_final_1 %>%
  dplyr::select(all_of(req_cols_ewb))

# Get only stage value
df_ewb1 <- df_ewb1 %>%
  mutate(
    stage_raw = as.character(diagnoses.ajcc_pathologic_stage),
    stage_clean = trimws(gsub("(?i)^Stage\\s+", "", stage_raw, perl = TRUE)),
    stage_group = case_when(
      is.na(stage_clean) ~ "X",
      grepl("^IV", stage_clean) ~ "IV",
      grepl("^III", stage_clean) ~ "III",
      grepl("^II", stage_clean) ~ "II",
      grepl("^I", stage_clean) ~ "I",
      TRUE ~ "X"
    )
  )

# Weighted average for missing race
race_counts <- df_ewb1 %>%
  filter(demographic.race != "Missing") %>%
  count(demographic.race)

# join counts to coef values
race_counts <- race_counts %>%
  mutate(coef = race_coefs)

# Calculate average effect for race
avg_race_effect <- sum(race_counts$n * race_counts$coef) / sum(race_counts$n)

# Weighted average for missing ethnicity
eth_counts <- df_ewb1 %>%
  filter(demographic.ethnicity != "Missing") %>%
  count(demographic.ethnicity)

# join counts to coef values
eth_counts <- eth_counts %>%
  mutate(coef = eth_coefs)

# Calculate average effect for ethnicity
avg_eth_effect <- sum(eth_counts$n * eth_counts$coef) / sum(eth_counts$n)

# Weighted average for missing stage
stage_counts <- df_ewb1 %>%
  filter(stage_group != "X") %>%
  count(stage_group)

# join counts to coef values
stage_counts <- stage_counts %>%
  mutate(coef = stage_coefs)

# Calculate average effect for stage
avg_stage_effect <-
  sum(stage_counts$n * stage_counts$coef) / sum(stage_counts$n)

# Weighted average for prior malignancy
prior_malig_counts <- df_ewb1 %>%
  filter(!is.na(diagnoses.prior_malignancy)) %>%
  count(diagnoses.prior_malignancy)

# join counts to coef values
prior_malig_counts <- prior_malig_counts %>%
  mutate(coef = prior_malig_coefs)

# Calculate average effect for prior malignancy
avg_prior_malig_effect <-
  sum(prior_malig_counts$n * prior_malig_counts$coef) /
  sum(prior_malig_counts$n)

# Map the values
df_ewb2 <- df_ewb1 %>%
  mutate(
    race_effect = ifelse(demographic.race == "Missing",
                         avg_race_effect,
                         race_coefs[as.character(demographic.race)]),
    eth_effect = ifelse(demographic.ethnicity == "Missing",
                        avg_eth_effect,
                        eth_coefs[as.character(demographic.ethnicity)]),
    stage_effect = ifelse(stage_group == "X",
                        avg_stage_effect,
                        stage_coefs[as.character(stage_group)]),
    prior_effect = ifelse(is.na(diagnoses.prior_malignancy),
                          avg_prior_malig_effect,
                          prior_malig_coefs[as.character(
                            diagnoses.prior_malignancy)])
  )

# Randomise trajectory #########################################################
df_ewb3 <- df_ewb2 %>%
  mutate(
    trajectory = ifelse(runif(n()) < 0.70, "Linear", "Quadratic")
  )

### Random effects definition ##################################################

# For both linear and quadratic
sd_intercept <- 4 # standard deviation of random intercept
sd_slope <- 0.1 # standard deviation of the average change

# For linear (2x2)
cor_int_slope <- -0.3
Sigma_lin <- matrix(c(sd_intercept^2,
                      cor_int_slope * sd_intercept * sd_slope,
                      cor_int_slope * sd_intercept * sd_slope,
                      sd_slope^2),
                    nrow = 2, byrow = TRUE)

# For quadratic (3x3)
sd_quad <- 0.03
cor_int_quad <- 0.25
cor_slope_quad <- 0.2
Sigma_quad <- matrix(c(
  sd_intercept^2,
  cor_int_slope * sd_intercept * sd_slope,
  cor_int_quad  * sd_intercept * sd_quad,
  
  cor_int_slope * sd_intercept * sd_slope,
  sd_slope^2,
  cor_slope_quad * sd_slope * sd_quad,
  
  cor_int_quad  * sd_intercept * sd_quad,
  cor_slope_quad * sd_slope * sd_quad,
  sd_quad^2
), nrow = 3, byrow = TRUE)

### Random effects assignment ##################################################

N <- nrow(df_ewb3)

df_ewb3$re_int <- NA
df_ewb3$re_slope <- NA
df_ewb3$re_quad <- 0  # default for linear subjects

# For linear subjects
idx_lin <- df_ewb3$trajectory == "Linear"

re_lin <- mvrnorm(sum(idx_lin), mu = c(0,0), Sigma = Sigma_lin)
df_ewb3$re_int[idx_lin] <- re_lin[,1]
df_ewb3$re_slope[idx_lin] <- re_lin[,2]

# For quadratic subjects
idx_quad <- df_ewb3$trajectory == "Quadratic"

re_quad <- mvrnorm(sum(idx_quad), mu = c(0,0,0), Sigma = Sigma_quad)
df_ewb3$re_int[idx_quad] <- re_quad[,1]
df_ewb3$re_slope[idx_quad] <- re_quad[,2]
df_ewb3$re_quad[idx_quad] <- re_quad[,3]

### Fixed effects ##############################################################

pop_intercept <- 17
age_coef <- -0.01
slope_intercept <- -0.05
age_on_slope <- 0.0005
fixed_quad <- -0.01
resid_sd <- 1.4

# Compute scores and expand to long format #####################################
n_time <- length(timepoints)

req_cols_ewb2 <- c(req_cols_ewb, "time", "ewb_total", "stage_group")

df_ewb_long <- df_ewb3 %>%
  slice(rep(1:n(), each = n_time)) %>%
  mutate(time = rep(timepoints, times = N)) %>%
  arrange(cases.case_id, time) %>%
  mutate(
    fixed_baseline = pop_intercept +
      age_coef * demographic.age_at_index +
      race_effect + eth_effect + stage_effect + prior_effect,
    fixed_slope = slope_intercept +
      age_on_slope * (demographic.age_at_index - mean(demographic.age_at_index,
                                                      na.rm = TRUE)),
    latent = fixed_baseline + re_int + (fixed_slope + re_slope) * time +
      (fixed_quad + re_quad) * time^2,
    observed_cont = latent + rnorm(n(), mean = 0, sd = resid_sd),
    observed_raw = round(observed_cont),
    observed_clamped = pmin(pmax(observed_raw, 0), ewb_max),
    ewb_total = as.integer(observed_clamped)
  ) %>%
  dplyr::select(all_of(req_cols_ewb2))

summary_by_time_ewb <- df_ewb_long %>%
  group_by(time) %>%
  summarise(n = n(), mean_ewb = mean(ewb_total), sd_ewb = sd(ewb_total))

baseline <- df_ewb_long %>% filter(time == 0)

group_means_race <- baseline %>%
  group_by(demographic.race) %>%
  summarise(n = n(), mean_ewb = mean(ewb_total), sd_ewb = sd(ewb_total))
group_means_eth <- baseline %>%
  group_by(demographic.ethnicity) %>%
  summarise(n = n(), mean_ewb = mean(ewb_total), sd_ewb = sd(ewb_total))
group_means_stage <- baseline %>%
  group_by(stage_group) %>%
  summarise(n = n(), mean_ewb = mean(ewb_total), sd_ewb = sd(ewb_total))
group_means_prior_malig <- baseline %>%
  group_by(diagnoses.prior_malignancy) %>%
  summarise(n = n(), mean_ewb = mean(ewb_total), sd_ewb = sd(ewb_total))

#print(summary_by_time_ewb)
#print(group_means_race)
#print(group_means_eth)
#print(group_means_stage)
#print(group_means_prior_malig)

# Summary: mean and sd of ewb_total by time and stage ##########################
df_ewb_long <- df_ewb_long %>%
  mutate(
    time_num = as.numeric(as.character(time)) # for plotting
  )

summary_by_stage_time <- df_ewb_long %>%
  group_by(time, stage_group) %>%
  summarise(
    n = n(),
    mean_ewb = mean(ewb_total, na.rm = TRUE),
    sd_ewb = sd(ewb_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(time, stage_group)

# Wide table showing mean (sd)
summary_wide <- summary_by_stage_time %>%
  mutate(mean_sd = sprintf("%.2f (%.2f)", mean_ewb, sd_ewb)) %>%
  dplyr::select(time, stage_group, mean_sd) %>%
  pivot_wider(names_from = time, values_from = mean_sd, names_prefix = "time_")


# Plot mean trajectories by stage (with 95% CI)
plot_df <- df_ewb_long %>%
  group_by(time_num, stage_group) %>%
  summarise(
    mean_ewb = mean(ewb_total, na.rm = TRUE),
    se_ewb = sd(ewb_total, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

p_ewb_stage <- ggplot(plot_df, aes(x = time_num, y = mean_ewb,
                         color = stage_group,
                         group = stage_group)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_ewb - 1.96 * se_ewb,
                    ymax = mean_ewb + 1.96 * se_ewb),
                width = 0.08) +
  scale_x_continuous(breaks = sort(unique(plot_df$time_num))) +
  labs(x = "Time (years)", y = "Mean SWB total", color = "Stage",
       title = "Mean EWB by Stage across Timepoints") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#print(p_ewb_stage)

# ------------------------------------------------------------------------------
# 3. Functional Well-Being (FWB) 0-28, 7 questions
# ------------------------------------------------------------------------------

# Timepoints and FWB specific data
timepoints <- c(0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10)
fwb_max <- 28

# Subject information coefficients
# (Race, ethnicity, stage, metastasis, prior malingancy, prior treatment,
#  menopausal status)
race_coefs <- c(
  "American Indian or Alaska Native" = -0.4,
  "Asian" = -0.1,
  "Black or African American" = -0.3,
  "White" = 0.0
)

eth_coefs <- c(
  "Hispanic or Latino" = -0.3,
  "Not Hispanic or Latino" = 0.0
)

stage_coefs <- c(
  "I" = 0.0,
  "II" = -1.5,
  "III" = -3.0,
  "IV" = -5.0
)

prior_malig_coefs <- c(
  "N" = 0.0,
  "Y" = -1.5
)

prior_treat_coefs <- c(
  "N" = 0.0,
  "Y" = -1.2
)

meno_coefs <- c(
  "Perimenopausal" = -0.8,
  "Postmenopausal" = -0.7,
  "Premenopausal" = 0.0
)

# Treatment information coefficients
# (lhrh only for Premonopausal subjects)
treat_coefs <- list(
  rad   = -0.8,
  chemo = -2.0,
  immuno = -1.2,
  serm  = -0.6,
  lhrh  = -0.8,
  ai    = -1.0,
  serd  = -1.0
)

# Subset only required variables
req_cols_fwb <- c(
  # Required cases.~ variables
  "cases.case_id",
  # Required demographic.~ variables
  "demographic.age_at_index", "demographic.ethnicity", "demographic.race",
  # Required diagnoses.~ variables
  "diagnoses.ajcc_pathologic_stage", "diagnoses.prior_malignancy",
  "diagnoses.prior_treatment",
  # Menopausal status
  "meno.status",
  # Derived treatment flags
  "rad.fl", "serm.fl", "lhrh_a.fl", "ai.fl", "serd.fl", "immuno.fl", "chemo.fl"
)

df_fwb1 <- df_final_1 %>%
  dplyr::select(all_of(req_cols_fwb))

# Get only stage value
df_fwb1 <- df_fwb1 %>%
  mutate(
    stage_raw = as.character(diagnoses.ajcc_pathologic_stage),
    stage_clean = trimws(gsub("(?i)^Stage\\s+", "", stage_raw, perl = TRUE)),
    stage_group = case_when(
      is.na(stage_clean) ~ "X",
      grepl("^IV", stage_clean) ~ "IV",
      grepl("^III", stage_clean) ~ "III",
      grepl("^II", stage_clean) ~ "II",
      grepl("^I", stage_clean) ~ "I",
      TRUE ~ "X"
    )
  )

# Weighted average for missing race
race_counts <- df_fwb1 %>%
  filter(demographic.race != "Missing") %>%
  count(demographic.race)

# join counts to coef values
race_counts <- race_counts %>%
  mutate(coef = race_coefs)

# Calculate average effect for race
avg_race_effect <- sum(race_counts$n * race_counts$coef) / sum(race_counts$n)

# Weighted average for missing ethnicity
eth_counts <- df_fwb1 %>%
  filter(demographic.ethnicity != "Missing") %>%
  count(demographic.ethnicity)

# join counts to coef values
eth_counts <- eth_counts %>%
  mutate(coef = eth_coefs)

# Calculate average effect for ethnicity
avg_eth_effect <- sum(eth_counts$n * eth_counts$coef) / sum(eth_counts$n)

# Weighted average for missing stage
stage_counts <- df_fwb1 %>%
  filter(stage_group != "X") %>%
  count(stage_group)

# join counts to coef values
stage_counts <- stage_counts %>%
  mutate(coef = stage_coefs)

# Calculate average effect for stage
avg_stage_effect <-
  sum(stage_counts$n * stage_counts$coef) / sum(stage_counts$n)

# Weighted average for prior malignancy
prior_malig_counts <- df_fwb1 %>%
  filter(!is.na(diagnoses.prior_malignancy)) %>%
  count(diagnoses.prior_malignancy)

# join counts to coef values
prior_malig_counts <- prior_malig_counts %>%
  mutate(coef = prior_malig_coefs)

# Calculate average effect for prior malignancy
avg_prior_malig_effect <-
  sum(prior_malig_counts$n * prior_malig_counts$coef) /
  sum(prior_malig_counts$n)

# Weighted average for prior treatment
prior_treat_counts <- df_fwb1 %>%
  filter(!is.na(diagnoses.prior_treatment)) %>%
  count(diagnoses.prior_treatment)

# join counts to coef values
prior_treat_counts <- prior_treat_counts %>%
  mutate(coef = prior_treat_coefs)

# Calculate average effect for prior treatment
avg_prior_treat_effect <-
  sum(prior_treat_counts$n * prior_treat_counts$coef) /
  sum(prior_treat_counts$n)

# Weighted average for menopausal status
meno_counts <- df_fwb1 %>%
  filter(meno.status != "Unknown") %>%
  count(meno.status)

# join counts to coef values
meno_counts <- meno_counts %>%
  mutate(coef = meno_coefs)

# Calculate average effect for menopausal status
avg_meno_effect <-
  sum(meno_counts$n * meno_counts$coef) /
  sum(meno_counts$n)

# Map the values
df_fwb2 <- df_fwb1 %>%
  mutate(
    # With different values
    race_effect = ifelse(demographic.race == "Missing",
                         avg_race_effect,
                         race_coefs[as.character(demographic.race)]),
    eth_effect = ifelse(demographic.ethnicity == "Missing",
                        avg_eth_effect,
                        eth_coefs[as.character(demographic.ethnicity)]),
    stage_effect = ifelse(stage_group == "X",
                          avg_stage_effect,
                          stage_coefs[as.character(stage_group)]),
    prior_malig_effect = ifelse(is.na(diagnoses.prior_malignancy),
                                avg_prior_malig_effect,
                                prior_malig_coefs[as.character(
                                  diagnoses.prior_malignancy)]),
    prior_treat_effect = ifelse(is.na(diagnoses.prior_treatment),
                                avg_prior_treat_effect,
                                prior_treat_coefs[as.character(
                                  diagnoses.prior_treatment)]),
    meno_effect = ifelse(meno.status == "Unknown",
                         avg_meno_effect,
                         meno_coefs[as.character(meno.status)]),
    # On/Off effect
    rad_effect  = ifelse(rad.fl  == "Y", treat_coefs$rad, 0),
    chemo_effect = ifelse(chemo.fl == "Y", treat_coefs$chemo, 0),
    immuno_effect = ifelse(immuno.fl == "Y", treat_coefs$immuno, 0),
    serm_effect = ifelse(serm.fl == "Y", treat_coefs$serm, 0),
    ai_effect   = ifelse(ai.fl   == "Y", treat_coefs$ai, 0),
    serd_effect = ifelse(serd.fl == "Y", treat_coefs$serd, 0),
    lhrh_effect = ifelse(lhrh_a.fl == "Y" & meno.status == "Premenopausal",
                         treat_coefs$lhrh, 0),
    total_treat_effect = rad_effect + chemo_effect + immuno_effect +
      serm_effect + lhrh_effect + ai_effect + serd_effect
  )

# Randomise trajectory #########################################################
df_fwb3 <- df_fwb2 %>%
  mutate(
    trajectory = ifelse(runif(n()) < 0.70, "Linear", "Quadratic")
  )

### Random effects definition ##################################################

# For both linear and quadratic
sd_intercept <- 4 # standard deviation of random intercept
sd_slope <- 0.4 # standard deviation of the average change

# For linear (2x2)
cor_int_slope <- -0.4
Sigma_lin <- matrix(c(sd_intercept^2,
                      cor_int_slope * sd_intercept * sd_slope,
                      cor_int_slope * sd_intercept * sd_slope,
                      sd_slope^2),
                    nrow = 2, byrow = TRUE)

# For quadratic (3x3)
sd_quad <- 0.03
cor_int_quad <- 0.25
cor_slope_quad <- 0.25
Sigma_quad <- matrix(c(
  sd_intercept^2,
  cor_int_slope * sd_intercept * sd_slope,
  cor_int_quad  * sd_intercept * sd_quad,
  
  cor_int_slope * sd_intercept * sd_slope,
  sd_slope^2,
  cor_slope_quad * sd_slope * sd_quad,
  
  cor_int_quad  * sd_intercept * sd_quad,
  cor_slope_quad * sd_slope * sd_quad,
  sd_quad^2
), nrow = 3, byrow = TRUE)

### Random effects assignment ##################################################

N <- nrow(df_ewb3)

df_fwb3$re_int <- NA
df_fwb3$re_slope <- NA
df_fwb3$re_quad <- 0  # default for linear subjects

# For linear subjects
idx_lin <- df_fwb3$trajectory == "Linear"

re_lin <- mvrnorm(sum(idx_lin), mu = c(0,0), Sigma = Sigma_lin)
df_fwb3$re_int[idx_lin] <- re_lin[,1]
df_fwb3$re_slope[idx_lin] <- re_lin[,2]

# For quadratic subjects
idx_quad <- df_fwb3$trajectory == "Quadratic"

re_quad <- mvrnorm(sum(idx_quad), mu = c(0,0,0), Sigma = Sigma_quad)
df_fwb3$re_int[idx_quad] <- re_quad[,1]
df_fwb3$re_slope[idx_quad] <- re_quad[,2]
df_fwb3$re_quad[idx_quad] <- re_quad[,3]

### Fixed effects ##############################################################

pop_intercept <- 20
age_coef <- -0.01
slope_intercept <- -0.05
age_on_slope <- 0.0005
fixed_quad <- -0.012
resid_sd <- 1.8

# Compute scores and expand to long format #####################################
n_time <- length(timepoints)

req_cols_fwb2 <- c(req_cols_fwb, "time", "fwb_total", "stage_group")

df_fwb_long <- df_fwb3 %>%
  slice(rep(1:n(), each = n_time)) %>%
  mutate(time = rep(timepoints, times = N)) %>%
  arrange(cases.case_id, time) %>%
  mutate(
    fixed_baseline = pop_intercept +
      age_coef * demographic.age_at_index +
      race_effect + eth_effect + stage_effect +
      prior_malig_effect + prior_treat_effect + meno_effect +
      total_treat_effect,
    
    fixed_slope = slope_intercept +
      age_on_slope * (demographic.age_at_index - mean(demographic.age_at_index,
                                                      na.rm = TRUE)),
    
    latent = fixed_baseline + re_int + (fixed_slope + re_slope) * time +
      (fixed_quad + re_quad) * time^2,
    
    observed_cont = latent + rnorm(n(), mean = 0, sd = resid_sd),
    
    observed_raw = round(observed_cont),
    observed_clamped = pmin(pmax(observed_raw, 0), fwb_max),
    fwb_total = as.integer(observed_clamped)
  ) %>%
  dplyr::select(all_of(req_cols_fwb2))

summary_by_time_fwb <- df_fwb_long %>%
  group_by(time) %>%
  summarise(n = n(), mean_fwb = mean(fwb_total), sd_fwb = sd(fwb_total))

baseline <- df_fwb_long %>% filter(time == 0)

group_means_race <- baseline %>%
  group_by(demographic.race) %>%
  summarise(n = n(), mean_fwb = mean(fwb_total), sd_fwb = sd(fwb_total))
group_means_eth <- baseline %>%
  group_by(demographic.ethnicity) %>%
  summarise(n = n(), mean_fwb = mean(fwb_total), sd_fwb = sd(fwb_total))
group_means_stage <- baseline %>%
  group_by(stage_group) %>%
  summarise(n = n(), mean_fwb = mean(fwb_total), sd_fwb = sd(fwb_total))
group_means_prior_malig <- baseline %>%
  group_by(diagnoses.prior_malignancy) %>%
  summarise(n = n(), mean_fwb = mean(fwb_total), sd_fwb = sd(fwb_total))
group_means_prior_treat <- baseline %>%
  group_by(diagnoses.prior_treatment) %>%
  summarise(n = n(), mean_fwb = mean(fwb_total), sd_fwb = sd(fwb_total))
group_means_meno <- baseline %>%
  group_by(meno.status) %>%
  summarise(n = n(), mean_fwb = mean(fwb_total), sd_fwb = sd(fwb_total))
group_means_rad <- baseline %>%
  group_by(rad.fl) %>%
  summarise(n = n(), mean_fwb = mean(fwb_total), sd_fwb = sd(fwb_total))
group_means_serm <- baseline %>%
  group_by(serm.fl) %>%
  summarise(n = n(), mean_fwb = mean(fwb_total), sd_fwb = sd(fwb_total))
group_means_lhrh <- baseline %>%
  group_by(lhrh_a.fl) %>%
  summarise(n = n(), mean_fwb = mean(fwb_total), sd_fwb = sd(fwb_total))
group_means_ai <- baseline %>%
  group_by(ai.fl) %>%
  summarise(n = n(), mean_fwb = mean(fwb_total), sd_fwb = sd(fwb_total))
group_means_serd <- baseline %>%
  group_by(serd.fl) %>%
  summarise(n = n(), mean_fwb = mean(fwb_total), sd_fwb = sd(fwb_total))
group_means_immuno <- baseline %>%
  group_by(immuno.fl) %>%
  summarise(n = n(), mean_fwb = mean(fwb_total), sd_fwb = sd(fwb_total))
group_means_chemo <- baseline %>%
  group_by(chemo.fl) %>%
  summarise(n = n(), mean_fwb = mean(fwb_total), sd_fwb = sd(fwb_total))

#print(summary_by_time_fwb)
#print(group_means_race)
#print(group_means_eth)
#print(group_means_stage)
#print(group_means_prior_malig)
#print(group_means_prior_treat)
#print(group_means_meno)
#print(group_means_rad)
#print(group_means_serm)
#print(group_means_lhrh)
#print(group_means_ai)
#print(group_means_serd)
#print(group_means_immuno)
#print(group_means_chemo)

# Summary: mean and sd of fwb_total by time and stage ##########################
df_fwb_long <- df_fwb_long %>%
  mutate(
    time_num = as.numeric(as.character(time)) # for plotting
  )

summary_by_stage_time <- df_fwb_long %>%
  group_by(time, stage_group) %>%
  summarise(
    n = n(),
    mean_fwb = mean(fwb_total, na.rm = TRUE),
    sd_fwb = sd(fwb_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(time, stage_group)

# Wide table showing mean (sd)
summary_wide <- summary_by_stage_time %>%
  mutate(mean_sd = sprintf("%.2f (%.2f)", mean_fwb, sd_fwb)) %>%
  dplyr::select(time, stage_group, mean_sd) %>%
  pivot_wider(names_from = time, values_from = mean_sd, names_prefix = "time_")


# Plot mean trajectories by stage (with 95% CI)
plot_df <- df_fwb_long %>%
  group_by(time_num, stage_group) %>%
  summarise(
    mean_fwb = mean(fwb_total, na.rm = TRUE),
    se_fwb = sd(fwb_total, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

p_fwb_stage <- ggplot(plot_df, aes(x = time_num, y = mean_fwb,
                         color = stage_group,
                         group = stage_group)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_fwb - 1.96 * se_fwb,
                    ymax = mean_fwb + 1.96 * se_fwb),
                width = 0.08) +
  scale_x_continuous(breaks = sort(unique(plot_df$time_num))) +
  labs(x = "Time (years)", y = "Mean FWB total", color = "Stage",
       title = "Mean FWB by Stage across Timepoints") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#print(p_fwb_stage)

# ------------------------------------------------------------------------------
# 4. Breast Cancer Subscale (BCS) 0-40, 10 questions
# ------------------------------------------------------------------------------

# Timepoints and BCS specific data
timepoints <- c(0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10)
bcs_max <- 40

# Subject information coefficients
# (Race, ethnicity, stage, menopausal status, margin status, molecular subtype)
race_coefs <- c(
  "American Indian or Alaska Native" = -0.4,
  "Asian" = -0.1,
  "Black or African American" = -0.3,
  "White" = 0.0
)

eth_coefs <- c(
  "Hispanic or Latino" = -0.3,
  "Not Hispanic or Latino" = 0.0
)

stage_coefs <- c(
  "I" = 0.0,
  "II" = -1.0,
  "III" = -2.5,
  "IV" = -4.0
)

meno_coefs <- c(
  "Perimenopausal" = -0.8,
  "Postmenopausal" = -0.7,
  "Premenopausal" = 0.0
)

margin_coefs <- c(
  "Involved" = -1.0,
  "Uninvolved" = 0.0
)

subtype_coefs <- c(
  "HER2-enriched" = -1.0,
  "Luminal A-like" = 0.0,
  "Luminal B-like (HER2-neg)" = -0.6,
  "Luminal B-like (HER2-pos)" = -0.8,
  "Triple-negative" = -2.0
)

# Treatment information coefficients
# (lhrh only for Premonopausal subjects)
treat_coefs <- list(
  rad   = -1.5,
  chemo = -3.5,
  immuno = -1.5,
  serm  = -1.0,
  lhrh  = -1.2,
  ai    = -1.6,
  serd  = -1.4
)


# Subset only required variables
req_cols_bcs <- c(
  # Required cases.~ variables
  "cases.case_id",
  # Required demographic.~ variables
  "demographic.age_at_index", "demographic.ethnicity", "demographic.race",
  # Required diagnoses.~ variables
  "diagnoses.ajcc_pathologic_stage",
  # Require first_line.~ variables
  "first_line.margin_status",
  # Menopausal status, molecular subtype, margin status
  "meno.status", "molecular_subtype",
  # Derived treatment flags
  "rad.fl", "serm.fl", "lhrh_a.fl", "ai.fl", "serd.fl", "immuno.fl", "chemo.fl"
)

df_bcs1 <- df_final_1 %>%
  dplyr::select(all_of(req_cols_bcs))

# Get only stage value
df_bcs1 <- df_bcs1 %>%
  mutate(
    stage_raw = as.character(diagnoses.ajcc_pathologic_stage),
    stage_clean = trimws(gsub("(?i)^Stage\\s+", "", stage_raw, perl = TRUE)),
    stage_group = case_when(
      is.na(stage_clean) ~ "X",
      grepl("^IV", stage_clean) ~ "IV",
      grepl("^III", stage_clean) ~ "III",
      grepl("^II", stage_clean) ~ "II",
      grepl("^I", stage_clean) ~ "I",
      TRUE ~ "X"
    )
  )

# Weighted average for missing race
race_counts <- df_bcs1 %>%
  filter(demographic.race != "Missing") %>%
  count(demographic.race)

# join counts to coef values
race_counts <- race_counts %>%
  mutate(coef = race_coefs)

# Calculate average effect for race
avg_race_effect <- sum(race_counts$n * race_counts$coef) / sum(race_counts$n)

# Weighted average for missing ethnicity
eth_counts <- df_bcs1 %>%
  filter(demographic.ethnicity != "Missing") %>%
  count(demographic.ethnicity)

# join counts to coef values
eth_counts <- eth_counts %>%
  mutate(coef = eth_coefs)

# Calculate average effect for ethnicity
avg_eth_effect <- sum(eth_counts$n * eth_counts$coef) / sum(eth_counts$n)

# Weighted average for missing stage
stage_counts <- df_bcs1 %>%
  filter(stage_group != "X") %>%
  count(stage_group)

# join counts to coef values
stage_counts <- stage_counts %>%
  mutate(coef = stage_coefs)

# Calculate average effect for stage
avg_stage_effect <-
  sum(stage_counts$n * stage_counts$coef) / sum(stage_counts$n)

# Weighted average for menopausal status
meno_counts <- df_bcs1 %>%
  filter(meno.status != "Unknown") %>%
  count(meno.status)

# join counts to coef values
meno_counts <- meno_counts %>%
  mutate(coef = meno_coefs)

# Calculate average effect for menopausal status
avg_meno_effect <-
  sum(meno_counts$n * meno_counts$coef) /
  sum(meno_counts$n)

# Weighted average for margin status
margin_counts <- df_bcs1 %>%
  filter(first_line.margin_status != "Indeterminate" &
           first_line.margin_status != "Missing") %>%
  count(first_line.margin_status)

# join counts to coef values
margin_counts <- margin_counts %>%
  mutate(coef = margin_coefs)

# Calculate average effect for margin status
avg_margin_effect <-
  sum(margin_counts$n * margin_counts$coef) /
  sum(margin_counts$n)

# Weighted average for molecular subtype
subtype_counts <- df_bcs1 %>%
  filter(molecular_subtype != "Unknown") %>%
  count(molecular_subtype)

# join counts to coef values
subtype_counts <- subtype_counts %>%
  mutate(coef = subtype_coefs)

# Calculate average effect for molecular subtype
avg_subtype_effect <-
  sum(subtype_counts$n * subtype_counts$coef) /
  sum(subtype_counts$n)

# Map the values
df_bcs2 <- df_bcs1 %>%
  mutate(
    # With different values
    race_effect = ifelse(demographic.race == "Missing",
                         avg_race_effect,
                         race_coefs[as.character(demographic.race)]),
    eth_effect = ifelse(demographic.ethnicity == "Missing",
                        avg_eth_effect,
                        eth_coefs[as.character(demographic.ethnicity)]),
    stage_effect = ifelse(stage_group == "X",
                          avg_stage_effect,
                          stage_coefs[as.character(stage_group)]),
    meno_effect = ifelse(meno.status == "Unknown",
                         avg_meno_effect,
                         meno_coefs[as.character(meno.status)]),
    margin_effect = ifelse(first_line.margin_status == "Indeterminate" |
                             first_line.margin_status == "Missing",
                           avg_margin_effect,
                           margin_coefs[as.character(
                             first_line.margin_status)]),
    subtype_effect = ifelse(molecular_subtype == "Unknown",
                            avg_subtype_effect,
                            subtype_coefs[as.character(molecular_subtype)]),
    # On/Off effect
    rad_effect  = ifelse(rad.fl  == "Y", treat_coefs$rad, 0),
    chemo_effect = ifelse(chemo.fl == "Y", treat_coefs$chemo, 0),
    immuno_effect = ifelse(immuno.fl == "Y", treat_coefs$immuno, 0),
    serm_effect = ifelse(serm.fl == "Y", treat_coefs$serm, 0),
    ai_effect   = ifelse(ai.fl   == "Y", treat_coefs$ai, 0),
    serd_effect = ifelse(serd.fl == "Y", treat_coefs$serd, 0),
    lhrh_effect = ifelse(lhrh_a.fl == "Y" & meno.status == "Premenopausal",
                         treat_coefs$lhrh, 0),
    total_treat_effect = rad_effect + chemo_effect + immuno_effect +
      serm_effect + lhrh_effect + ai_effect + serd_effect
  )

# Randomise trajectory #########################################################
df_bcs3 <- df_bcs2 %>%
  mutate(
    trajectory = ifelse(runif(n()) < 0.70, "Linear", "Quadratic")
  )

### Random effects definition ##################################################

# For both linear and quadratic
sd_intercept <- 7 # standard deviation of random intercept
sd_slope <- 1 # standard deviation of the average change

# For linear (2x2)
cor_int_slope <- -0.4
Sigma_lin <- matrix(c(sd_intercept^2,
                      cor_int_slope * sd_intercept * sd_slope,
                      cor_int_slope * sd_intercept * sd_slope,
                      sd_slope^2),
                    nrow = 2, byrow = TRUE)

# For quadratic (3x3)
sd_quad <- 0.03
cor_int_quad <- 0.25
cor_slope_quad <- 0.25
Sigma_quad <- matrix(c(
  sd_intercept^2,
  cor_int_slope * sd_intercept * sd_slope,
  cor_int_quad  * sd_intercept * sd_quad,
  
  cor_int_slope * sd_intercept * sd_slope,
  sd_slope^2,
  cor_slope_quad * sd_slope * sd_quad,
  
  cor_int_quad  * sd_intercept * sd_quad,
  cor_slope_quad * sd_slope * sd_quad,
  sd_quad^2
), nrow = 3, byrow = TRUE)

### Random effects assignment ##################################################

N <- nrow(df_bcs3)

df_bcs3$re_int <- NA
df_bcs3$re_slope <- NA
df_bcs3$re_quad <- 0  # default for linear subjects

# For linear subjects
idx_lin <- df_bcs3$trajectory == "Linear"

re_lin <- mvrnorm(sum(idx_lin), mu = c(0,0), Sigma = Sigma_lin)
df_bcs3$re_int[idx_lin] <- re_lin[,1]
df_bcs3$re_slope[idx_lin] <- re_lin[,2]

# For quadratic subjects
idx_quad <- df_bcs3$trajectory == "Quadratic"

re_quad <- mvrnorm(sum(idx_quad), mu = c(0,0,0), Sigma = Sigma_quad)
df_bcs3$re_int[idx_quad] <- re_quad[,1]
df_bcs3$re_slope[idx_quad] <- re_quad[,2]
df_bcs3$re_quad[idx_quad] <- re_quad[,3]

### Fixed effects ##############################################################

pop_intercept <- 28
age_coef <- -0.01
slope_intercept <- -0.05
age_on_slope <- 0.0005
fixed_quad <- -0.01
resid_sd <- 2.4

# Compute scores and expand to long format #####################################
n_time <- length(timepoints)

req_cols_bcs2 <- c(req_cols_bcs, "time", "bcs_total", "stage_group")

df_bcs_long <- df_bcs3 %>%
  slice(rep(1:n(), each = n_time)) %>%
  mutate(time = rep(timepoints, times = N)) %>%
  arrange(cases.case_id, time) %>%
  mutate(
    fixed_baseline = pop_intercept +
      age_coef * demographic.age_at_index +
      race_effect + eth_effect + stage_effect +
      meno_effect + margin_effect + subtype_effect +
      total_treat_effect,
    
    fixed_slope = slope_intercept +
      age_on_slope * (demographic.age_at_index - mean(demographic.age_at_index,
                                                      na.rm = TRUE)),
    
    latent = fixed_baseline + re_int + (fixed_slope + re_slope) * time +
      (fixed_quad + re_quad) * time^2,
    
    observed_cont = latent + rnorm(n(), mean = 0, sd = resid_sd),
    
    observed_raw = round(observed_cont),
    observed_clamped = pmin(pmax(observed_raw, 0), bcs_max),
    bcs_total = as.integer(observed_clamped)
  ) %>%
  dplyr::select(all_of(req_cols_bcs2))

summary_by_time_bcs <- df_bcs_long %>%
  group_by(time) %>%
  summarise(n = n(), mean_bcs = mean(bcs_total), sd_bcs = sd(bcs_total))

baseline <- df_bcs_long %>% filter(time == 0)

group_means_race <- baseline %>%
  group_by(demographic.race) %>%
  summarise(n = n(), mean_bcs = mean(bcs_total), sd_bcs = sd(bcs_total))
group_means_eth <- baseline %>%
  group_by(demographic.ethnicity) %>%
  summarise(n = n(), mean_bcs = mean(bcs_total), sd_bcs = sd(bcs_total))
group_means_stage <- baseline %>%
  group_by(stage_group) %>%
  summarise(n = n(), mean_bcs = mean(bcs_total), sd_bcs = sd(bcs_total))
group_means_meno <- baseline %>%
  group_by(meno.status) %>%
  summarise(n = n(), mean_bcs = mean(bcs_total), sd_bcs = sd(bcs_total))
group_means_margin <- baseline %>%
  group_by(first_line.margin_status) %>%
  summarise(n = n(), mean_bcs = mean(bcs_total), sd_bcs = sd(bcs_total))
group_means_subtype <- baseline %>%
  group_by(molecular_subtype) %>%
  summarise(n = n(), mean_bcs = mean(bcs_total), sd_bcs = sd(bcs_total))

group_means_rad <- baseline %>%
  group_by(rad.fl) %>%
  summarise(n = n(), mean_bcs = mean(bcs_total), sd_bcs = sd(bcs_total))
group_means_serm <- baseline %>%
  group_by(serm.fl) %>%
  summarise(n = n(), mean_bcs = mean(bcs_total), sd_bcs = sd(bcs_total))
group_means_lhrh <- baseline %>%
  group_by(lhrh_a.fl) %>%
  summarise(n = n(), mean_bcs = mean(bcs_total), sd_bcs = sd(bcs_total))
group_means_ai <- baseline %>%
  group_by(ai.fl) %>%
  summarise(n = n(), mean_bcs = mean(bcs_total), sd_bcs = sd(bcs_total))
group_means_serd <- baseline %>%
  group_by(serd.fl) %>%
  summarise(n = n(), mean_bcs = mean(bcs_total), sd_bcs = sd(bcs_total))
group_means_immuno <- baseline %>%
  group_by(immuno.fl) %>%
  summarise(n = n(), mean_bcs = mean(bcs_total), sd_bcs = sd(bcs_total))
group_means_chemo <- baseline %>%
  group_by(chemo.fl) %>%
  summarise(n = n(), mean_bcs = mean(bcs_total), sd_bcs = sd(bcs_total))

#print(summary_by_time_bcs)
#print(group_means_race)
#print(group_means_eth)
#print(group_means_stage)
#print(group_means_meno)
#print(group_means_margin)
#print(group_means_subtype)

#print(group_means_rad)
#print(group_means_serm)
#print(group_means_lhrh)
#print(group_means_ai)
#print(group_means_serd)
#print(group_means_immuno)
#print(group_means_chemo)

# Summary: mean and sd of bcs_total by time and stage ##########################
df_bcs_long <- df_bcs_long %>%
  mutate(
    time_num = as.numeric(as.character(time)) # for plotting
  )

summary_by_stage_time <- df_bcs_long %>%
  group_by(time, stage_group) %>%
  summarise(
    n = n(),
    mean_bcs = mean(bcs_total, na.rm = TRUE),
    sd_bcs = sd(bcs_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(time, stage_group)

# Wide table showing mean (sd)
summary_wide <- summary_by_stage_time %>%
  mutate(mean_sd = sprintf("%.2f (%.2f)", mean_bcs, sd_bcs)) %>%
  dplyr::select(time, stage_group, mean_sd) %>%
  pivot_wider(names_from = time, values_from = mean_sd, names_prefix = "time_")


# Plot mean trajectories by stage (with 95% CI)
plot_df <- df_bcs_long %>%
  group_by(time_num, stage_group) %>%
  summarise(
    mean_bcs = mean(bcs_total, na.rm = TRUE),
    se_bcs = sd(bcs_total, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

p_bcs_stage <- ggplot(plot_df, aes(x = time_num, y = mean_bcs,
                         color = stage_group, group = stage_group)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_bcs - 1.96 * se_bcs,
                    ymax = mean_bcs + 1.96 * se_bcs),
                width = 0.08) +
  scale_x_continuous(breaks = sort(unique(plot_df$time_num))) +
  scale_color_manual(values = okabe_ito) +
  labs(x = "Time (years)", y = "Mean BCS total (Max = 40)", color = "Stage",
       title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1))

#print(p_bcs_stage)

# ------------------------------------------------------------------------------
# 5. Physical Well-Being (PWB) 0-28, 7 questions
# ------------------------------------------------------------------------------

# Timepoints and PWB specific data
timepoints <- c(0, 0.25, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10)
pwb_max <- 28

# Subject information coefficients
# (Race, ethnicity, stage, menopausal status, margin status, molecular subtype,
#  prior malignancy, prior treatment)
race_coefs <- c(
  "American Indian or Alaska Native" = -0.6,
  "Asian" = -0.2,
  "Black or African American" = -0.8,
  "White" = 0.0
)

eth_coefs <- c(
  "Hispanic or Latino" = -0.4,
  "Not Hispanic or Latino" = 0.0
)

stage_coefs <- c(
  "I" = 0.0,
  "II" = -1.2,
  "III" = -2.8,
  "IV" = -4.5
)

meno_coefs <- c(
  "Perimenopausal" = -0.8,
  "Postmenopausal" = -0.7,
  "Premenopausal" = 0.0
)

margin_coefs <- c(
  "Involved" = -1.0,
  "Uninvolved" = 0.0
)

subtype_coefs <- c(
  "HER2-enriched" = -1.0,
  "Luminal A-like" = 0.0,
  "Luminal B-like (HER2-neg)" = -0.7,
  "Luminal B-like (HER2-pos)" = -0.9,
  "Triple-negative" = -2.2
)

# Treatment information coefficients
# (lhrh only for Premonopausal subjects)
treat_coefs <- list(
  rad   = -1.5,
  chemo = -4.0,
  immuno = -1.5,
  serm  = -1.0,
  lhrh  = -1.2,
  ai    = -1.6,
  serd  = -1.4
)

prior_malig_coefs <- c(
  "N" = 0.0,
  "Y" = -1.5
)

prior_treat_coefs <- c(
  "N" = 0.0,
  "Y" = -1.2
)

# Subset only required variables
req_cols_pwb <- c(
  # Required cases.~ variables
  "cases.case_id",
  # Required demographic.~ variables
  "demographic.age_at_index", "demographic.ethnicity", "demographic.race",
  # Required diagnoses.~ variables
  "diagnoses.ajcc_pathologic_stage", "diagnoses.prior_malignancy",
  "diagnoses.prior_treatment",
  # Require first_line.~ variables
  "first_line.margin_status",
  # Menopausal status, molecular subtype, margin status
  "meno.status", "molecular_subtype",
  # Derived treatment flags
  "rad.fl", "serm.fl", "lhrh_a.fl", "ai.fl", "serd.fl", "immuno.fl", "chemo.fl"
)

df_pwb1 <- df_final_1 %>%
  dplyr::select(all_of(req_cols_pwb))

# Get only stage value
df_pwb1 <- df_pwb1 %>%
  mutate(
    stage_raw = as.character(diagnoses.ajcc_pathologic_stage),
    stage_clean = trimws(gsub("(?i)^Stage\\s+", "", stage_raw, perl = TRUE)),
    stage_group = case_when(
      is.na(stage_clean) ~ "X",
      grepl("^IV", stage_clean) ~ "IV",
      grepl("^III", stage_clean) ~ "III",
      grepl("^II", stage_clean) ~ "II",
      grepl("^I", stage_clean) ~ "I",
      TRUE ~ "X"
    )
  )

# Weighted average for missing race
race_counts <- df_pwb1 %>%
  filter(demographic.race != "Missing") %>%
  count(demographic.race)

# join counts to coef values
race_counts <- race_counts %>%
  mutate(coef = race_coefs)

# Calculate average effect for race
avg_race_effect <- sum(race_counts$n * race_counts$coef) / sum(race_counts$n)

# Weighted average for missing ethnicity
eth_counts <- df_pwb1 %>%
  filter(demographic.ethnicity != "Missing") %>%
  count(demographic.ethnicity)

# join counts to coef values
eth_counts <- eth_counts %>%
  mutate(coef = eth_coefs)

# Calculate average effect for ethnicity
avg_eth_effect <- sum(eth_counts$n * eth_counts$coef) / sum(eth_counts$n)

# Weighted average for missing stage
stage_counts <- df_pwb1 %>%
  filter(stage_group != "X") %>%
  count(stage_group)

# join counts to coef values
stage_counts <- stage_counts %>%
  mutate(coef = stage_coefs)

# Calculate average effect for stage
avg_stage_effect <-
  sum(stage_counts$n * stage_counts$coef) / sum(stage_counts$n)

# Weighted average for menopausal status
meno_counts <- df_pwb1 %>%
  filter(meno.status != "Unknown") %>%
  count(meno.status)

# join counts to coef values
meno_counts <- meno_counts %>%
  mutate(coef = meno_coefs)

# Calculate average effect for menopausal status
avg_meno_effect <-
  sum(meno_counts$n * meno_counts$coef) /
  sum(meno_counts$n)

# Weighted average for margin status
margin_counts <- df_pwb1 %>%
  filter(first_line.margin_status != "Indeterminate" &
           first_line.margin_status != "Missing") %>%
  count(first_line.margin_status)

# join counts to coef values
margin_counts <- margin_counts %>%
  mutate(coef = margin_coefs)

# Calculate average effect for margin status
avg_margin_effect <-
  sum(margin_counts$n * margin_counts$coef) /
  sum(margin_counts$n)

# Weighted average for molecular subtype
subtype_counts <- df_pwb1 %>%
  filter(molecular_subtype != "Unknown") %>%
  count(molecular_subtype)

# join counts to coef values
subtype_counts <- subtype_counts %>%
  mutate(coef = subtype_coefs)

# Calculate average effect for molecular subtype
avg_subtype_effect <-
  sum(subtype_counts$n * subtype_counts$coef) /
  sum(subtype_counts$n)

# Weighted average for prior malignancy
prior_malig_counts <- df_pwb1 %>%
  filter(!is.na(diagnoses.prior_malignancy)) %>%
  count(diagnoses.prior_malignancy)

# join counts to coef values
prior_malig_counts <- prior_malig_counts %>%
  mutate(coef = prior_malig_coefs)

# Calculate average effect for prior malignancy
avg_prior_malig_effect <-
  sum(prior_malig_counts$n * prior_malig_counts$coef) /
  sum(prior_malig_counts$n)

# Weighted average for prior treatment
prior_treat_counts <- df_pwb1 %>%
  filter(!is.na(diagnoses.prior_treatment)) %>%
  count(diagnoses.prior_treatment)

# join counts to coef values
prior_treat_counts <- prior_treat_counts %>%
  mutate(coef = prior_treat_coefs)

# Calculate average effect for prior treatment
avg_prior_treat_effect <-
  sum(prior_treat_counts$n * prior_treat_counts$coef) /
  sum(prior_treat_counts$n)

# Map the values
df_pwb2 <- df_pwb1 %>%
  mutate(
    # With different values
    race_effect = ifelse(demographic.race == "Missing",
                         avg_race_effect,
                         race_coefs[as.character(demographic.race)]),
    eth_effect = ifelse(demographic.ethnicity == "Missing",
                        avg_eth_effect,
                        eth_coefs[as.character(demographic.ethnicity)]),
    stage_effect = ifelse(stage_group == "X",
                          avg_stage_effect,
                          stage_coefs[as.character(stage_group)]),
    meno_effect = ifelse(meno.status == "Unknown",
                         avg_meno_effect,
                         meno_coefs[as.character(meno.status)]),
    margin_effect = ifelse(first_line.margin_status == "Indeterminate" |
                             first_line.margin_status == "Missing",
                           avg_margin_effect,
                           margin_coefs[as.character(
                             first_line.margin_status)]),
    subtype_effect = ifelse(molecular_subtype == "Unknown",
                            avg_subtype_effect,
                            subtype_coefs[as.character(molecular_subtype)]),
    prior_malig_effect = ifelse(is.na(diagnoses.prior_malignancy),
                                avg_prior_malig_effect,
                                prior_malig_coefs[as.character(
                                  diagnoses.prior_malignancy)]),
    prior_treat_effect = ifelse(is.na(diagnoses.prior_treatment),
                                avg_prior_treat_effect,
                                prior_treat_coefs[as.character(
                                  diagnoses.prior_treatment)]),    
    # On/Off effect
    rad_effect  = ifelse(rad.fl  == "Y", treat_coefs$rad, 0),
    chemo_effect = ifelse(chemo.fl == "Y", treat_coefs$chemo, 0),
    immuno_effect = ifelse(immuno.fl == "Y", treat_coefs$immuno, 0),
    serm_effect = ifelse(serm.fl == "Y", treat_coefs$serm, 0),
    ai_effect   = ifelse(ai.fl   == "Y", treat_coefs$ai, 0),
    serd_effect = ifelse(serd.fl == "Y", treat_coefs$serd, 0),
    lhrh_effect = ifelse(lhrh_a.fl == "Y" & meno.status == "Premenopausal",
                         treat_coefs$lhrh, 0),
    total_treat_effect = rad_effect + chemo_effect + immuno_effect +
      serm_effect + lhrh_effect + ai_effect + serd_effect
  )

# Randomise trajectory #########################################################
df_pwb3 <- df_pwb2 %>%
  mutate(
    trajectory = ifelse(runif(n()) < 0.70, "Linear", "Quadratic")
  )

### Random effects definition ##################################################

# For both linear and quadratic
sd_intercept <- 5 # standard deviation of random intercept
sd_slope <- 0.15 # standard deviation of the average change

# For linear (2x2)
cor_int_slope <- -0.4
Sigma_lin <- matrix(c(sd_intercept^2,
                      cor_int_slope * sd_intercept * sd_slope,
                      cor_int_slope * sd_intercept * sd_slope,
                      sd_slope^2),
                    nrow = 2, byrow = TRUE)

# For quadratic (3x3)
sd_quad <- 0.03
cor_int_quad <- 0.25
cor_slope_quad <- 0.25
Sigma_quad <- matrix(c(
  sd_intercept^2,
  cor_int_slope * sd_intercept * sd_slope,
  cor_int_quad  * sd_intercept * sd_quad,
  
  cor_int_slope * sd_intercept * sd_slope,
  sd_slope^2,
  cor_slope_quad * sd_slope * sd_quad,
  
  cor_int_quad  * sd_intercept * sd_quad,
  cor_slope_quad * sd_slope * sd_quad,
  sd_quad^2
), nrow = 3, byrow = TRUE)

### Random effects assignment ##################################################

N <- nrow(df_bcs3)

df_pwb3$re_int <- NA
df_pwb3$re_slope <- NA
df_pwb3$re_quad <- 0  # default for linear subjects

# For linear subjects
idx_lin <- df_pwb3$trajectory == "Linear"

re_lin <- mvrnorm(sum(idx_lin), mu = c(0,0), Sigma = Sigma_lin)
df_pwb3$re_int[idx_lin] <- re_lin[,1]
df_pwb3$re_slope[idx_lin] <- re_lin[,2]

# For quadratic subjects
idx_quad <- df_pwb3$trajectory == "Quadratic"

re_quad <- mvrnorm(sum(idx_quad), mu = c(0,0,0), Sigma = Sigma_quad)
df_pwb3$re_int[idx_quad] <- re_quad[,1]
df_pwb3$re_slope[idx_quad] <- re_quad[,2]
df_pwb3$re_quad[idx_quad] <- re_quad[,3]

### Fixed effects ##############################################################

pop_intercept <- 18
age_coef <- -0.01
slope_intercept <- 0.06 # average increase
age_on_slope <- 0.0005
fixed_quad <- -0.008
resid_sd <- 2.2

# Compute scores and expand to long format #####################################
n_time <- length(timepoints)

req_cols_pwb2 <- c(req_cols_pwb, "time", "pwb_total", "stage_group")

df_pwb_long <- df_pwb3 %>%
  slice(rep(1:n(), each = n_time)) %>%
  mutate(time = rep(timepoints, times = N)) %>%
  arrange(cases.case_id, time) %>%
  mutate(
    fixed_baseline = pop_intercept +
      age_coef * demographic.age_at_index +
      race_effect + eth_effect + stage_effect +
      meno_effect + margin_effect + subtype_effect +
      prior_malig_effect + prior_treat_effect + total_treat_effect,
    
    fixed_slope = slope_intercept +
      age_on_slope * (demographic.age_at_index - mean(demographic.age_at_index,
                                                      na.rm = TRUE)),
    
    latent = fixed_baseline + re_int + (fixed_slope + re_slope) * time +
      (fixed_quad + re_quad) * time^2,
    
    observed_cont = latent + rnorm(n(), mean = 0, sd = resid_sd),
    
    observed_raw = round(observed_cont),
    observed_clamped = pmin(pmax(observed_raw, 0), pwb_max),
    pwb_total = as.integer(observed_clamped)
  ) %>%
  dplyr::select(all_of(req_cols_pwb2))

summary_by_time_pwb <- df_pwb_long %>%
  group_by(time) %>%
  summarise(n = n(), mean_pwb = mean(pwb_total), sd_pwb = sd(pwb_total))

baseline <- df_pwb_long %>% filter(time == 0)

group_means_race <- baseline %>%
  group_by(demographic.race) %>%
  summarise(n = n(), mean_pwb = mean(pwb_total), sd_pwb = sd(pwb_total))
group_means_eth <- baseline %>%
  group_by(demographic.ethnicity) %>%
  summarise(n = n(), mean_pwb = mean(pwb_total), sd_pwb = sd(pwb_total))
group_means_stage <- baseline %>%
  group_by(stage_group) %>%
  summarise(n = n(), mean_pwb = mean(pwb_total), sd_pwb = sd(pwb_total))
group_means_meno <- baseline %>%
  group_by(meno.status) %>%
  summarise(n = n(), mean_pwb = mean(pwb_total), sd_pwb = sd(pwb_total))
group_means_margin <- baseline %>%
  group_by(first_line.margin_status) %>%
  summarise(n = n(), mean_pwb = mean(pwb_total), sd_pwb = sd(pwb_total))
group_means_subtype <- baseline %>%
  group_by(molecular_subtype) %>%
  summarise(n = n(), mean_pwb = mean(pwb_total), sd_pwb = sd(pwb_total))
group_means_prior_malig <- baseline %>%
  group_by(diagnoses.prior_malignancy) %>%
  summarise(n = n(), mean_fwb = mean(pwb_total), sd_fwb = sd(pwb_total))
group_means_prior_treat <- baseline %>%
  group_by(diagnoses.prior_treatment) %>%
  summarise(n = n(), mean_fwb = mean(pwb_total), sd_fwb = sd(pwb_total))

group_means_rad <- baseline %>%
  group_by(rad.fl) %>%
  summarise(n = n(), mean_pwb = mean(pwb_total), sd_pwb = sd(pwb_total))
group_means_serm <- baseline %>%
  group_by(serm.fl) %>%
  summarise(n = n(), mean_pwb = mean(pwb_total), sd_pwb = sd(pwb_total))
group_means_lhrh <- baseline %>%
  group_by(lhrh_a.fl) %>%
  summarise(n = n(), mean_pwb = mean(pwb_total), sd_pwb = sd(pwb_total))
group_means_ai <- baseline %>%
  group_by(ai.fl) %>%
  summarise(n = n(), mean_pwb = mean(pwb_total), sd_pwb = sd(pwb_total))
group_means_serd <- baseline %>%
  group_by(serd.fl) %>%
  summarise(n = n(), mean_pwb = mean(pwb_total), sd_pwb = sd(pwb_total))
group_means_immuno <- baseline %>%
  group_by(immuno.fl) %>%
  summarise(n = n(), mean_pwb = mean(pwb_total), sd_pwb = sd(pwb_total))
group_means_chemo <- baseline %>%
  group_by(chemo.fl) %>%
  summarise(n = n(), mean_pwb = mean(pwb_total), sd_pwb = sd(pwb_total))

#print(summary_by_time_pwb)
#print(group_means_race)
#print(group_means_eth)
#print(group_means_stage)
#print(group_means_meno)
#print(group_means_margin)
#print(group_means_subtype)
#print(group_means_prior_malig)
#print(group_means_prior_treat)

#print(group_means_rad)
#print(group_means_serm)
#print(group_means_lhrh)
#print(group_means_ai)
#print(group_means_serd)
#print(group_means_immuno)
#print(group_means_chemo)

# Summary: mean and sd of pwb_total by time and stage ##########################
df_pwb_long <- df_pwb_long %>%
  mutate(
    time_num = as.numeric(as.character(time)) # for plotting
  )

summary_by_stage_time <- df_pwb_long %>%
  group_by(time, stage_group) %>%
  summarise(
    n = n(),
    mean_pwb = mean(pwb_total, na.rm = TRUE),
    sd_pwb = sd(pwb_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(time, stage_group)

# Wide table showing mean (sd)
summary_wide <- summary_by_stage_time %>%
  mutate(mean_sd = sprintf("%.2f (%.2f)", mean_pwb, sd_pwb)) %>%
  dplyr::select(time, stage_group, mean_sd) %>%
  pivot_wider(names_from = time, values_from = mean_sd, names_prefix = "time_")


# Plot mean trajectories by stage (with 95% CI)
plot_df <- df_pwb_long %>%
  group_by(time_num, stage_group) %>%
  summarise(
    mean_pwb = mean(pwb_total, na.rm = TRUE),
    se_pwb = sd(pwb_total, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop")

p_pwb_stage <- ggplot(plot_df, aes(x = time_num, y = mean_pwb,
                                   color = stage_group, group = stage_group)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_pwb - 1.96 * se_pwb,
                    ymax = mean_pwb + 1.96 * se_pwb),
                width = 0.08) +
  scale_x_continuous(breaks = sort(unique(plot_df$time_num))) +
  scale_color_manual(values = okabe_ito) +
  labs(x = "Time (years)", y = "Mean PWB total (Max = 28)", color = "Stage",
       title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1))

#print(p_pwb_stage)

# ------------------------------------------------------------------------------
# 6. Join all and calculate total
# ------------------------------------------------------------------------------

df_tot <- df_pwb_long %>%
  left_join(df_swb_long %>% dplyr::select(cases.case_id, time, swb_total),
            by = c("cases.case_id","time")) %>%
  left_join(df_ewb_long %>% dplyr::select(cases.case_id, time, ewb_total),
            by = c("cases.case_id","time")) %>%
  left_join(df_fwb_long %>% dplyr::select(cases.case_id, time, fwb_total),
            by = c("cases.case_id","time")) %>%
  left_join(df_bcs_long %>% dplyr::select(cases.case_id, time, bcs_total),
            by = c("cases.case_id","time")) %>%
  mutate(fact_b_total = swb_total + ewb_total + fwb_total + bcs_total
         + pwb_total)

summary_by_time_fact_b <- df_tot %>%
  group_by(time) %>%
  summarise(n = n(),
            mean_fact_b_total = mean(fact_b_total),
            sd_fact_b_total = sd(fact_b_total))

#print(summary_by_time_fact_b)

# Summary: mean and sd of fact_b_total by time and stage #######################
df_fact_b_long <- df_tot %>%
  mutate(
    time_num = as.numeric(as.character(time)) # for plotting
  )

summary_by_stage_time <- df_fact_b_long %>%
  group_by(time, stage_group) %>%
  summarise(
    n = n(),
    mean_FACT_B = mean(fact_b_total, na.rm = TRUE),
    sd_FACT_B = sd(fact_b_total, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(time, stage_group)

# Wide table showing mean (sd)
summary_wide <- summary_by_stage_time %>%
  mutate(mean_sd = sprintf("%.2f (%.2f)", mean_FACT_B, sd_FACT_B)) %>%
  dplyr::select(time, stage_group, mean_sd) %>%
  pivot_wider(names_from = time, values_from = mean_sd, names_prefix = "time_")

# Plot mean trajectories by stage (with 95% CI)
plot_df <- df_fact_b_long %>%
  group_by(time_num, stage_group) %>%
  summarise(
    mean_FACT_B = mean(fact_b_total, na.rm = TRUE),
    se_FACT_B = sd(fact_b_total, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop")

p_fact_b_stage <- ggplot(plot_df, aes(x = time_num, y = mean_FACT_B,
                                   color = stage_group, group = stage_group)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_FACT_B - 1.96 * se_FACT_B,
                    ymax = mean_FACT_B + 1.96 * se_FACT_B),
                width = 0.08) +
  scale_x_continuous(breaks = sort(unique(plot_df$time_num))) +
  scale_color_manual(values = okabe_ito) +
  labs(x = "Time (years)", y = "Mean FACT-B total (Max = 148)", color = "Stage",
       title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#print(p_fact_b_stage)

### Visualisation for BCS, PWB, and FACT-B total by stage and time #############
p1 <- p_bcs_stage + theme(legend.position = "none")
p2 <- p_pwb_stage + theme(legend.position = "none")
p3 <- p_fact_b_stage + theme(legend.position = "right")

(p1 | p2) / (p3)


### Mean % and SD % visualisation for all scores (all data) ####################
summary_by_time_all <- summary_by_time_swb %>%
  left_join(summary_by_time_ewb, by = c("time", "n")) %>%
  left_join(summary_by_time_fwb, by = c("time", "n")) %>%
  left_join(summary_by_time_bcs, by = c("time", "n")) %>%
  left_join(summary_by_time_pwb, by = c("time", "n")) %>%
  left_join(summary_by_time_fact_b, by = c("time", "n")) %>%
  mutate(
    max_SWB = 28,
    max_EWB = 24,
    max_FWB = 28,
    max_BCS = 40,
    max_PWB = 28,
    max_FACT_B = 148,
    pct_SWB = (mean_swb/max_SWB)*100,
    pct_EWB = (mean_ewb/max_EWB)*100,
    pct_FWB = (mean_fwb/max_FWB)*100,
    pct_BCS = (mean_bcs/max_BCS)*100,
    pct_PWB = (mean_pwb/max_PWB)*100,
    pct_FACT_B = (mean_fact_b_total/max_FACT_B)*100,
    pct_sd_SWB = (sd_swb/max_SWB)*100,
    pct_sd_EWB = (sd_ewb/max_EWB)*100,
    pct_sd_FWB = (sd_fwb/max_FWB)*100,
    pct_sd_BCS = (sd_bcs/max_BCS)*100,
    pct_sd_PWB = (sd_pwb/max_PWB)*100,
    pct_sd_FACT_B = (sd_fact_b_total/max_FACT_B)*100
  )

plot_df <- summary_by_time_all %>%
  dplyr::select(
    time,
    starts_with("pct_"),
  ) %>%
  pivot_longer(
    cols = -time,
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    type = ifelse(grepl("^pct_sd_", variable), "sd", "mean"),
    domain = gsub("^pct_sd_|^pct_", "", variable)
  ) %>%
  dplyr::select(time, domain, type, value) %>%
  pivot_wider(
    names_from = type,
    values_from = value
  )

plot_df$domain <- factor(
  plot_df$domain,
  levels = c("PWB", "SWB", "EWB", "FWB", "BCS", "FACT_B")
)

p_all_scores <-
ggplot(plot_df, aes(x = time, y = mean, color = domain)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    width = 0.08,
    alpha = 0.7
  ) +
  scale_y_continuous(
    limits = c(15, 95),
    breaks = seq(15, 95, 10)
  ) +
  scale_x_continuous(breaks = sort(unique(plot_df$time))) +
  scale_color_manual(values = okabe_ito) +
  labs(
    title = "",
    x = "Time (years)",
    y = "% of maximum score",
    color = "Domain"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 12)
  )

#p_all_scores

# ------------------------------------------------------------------------------
# 7. Discard measurements after day of death or last follow up
# ------------------------------------------------------------------------------

df_keep <- df_final_1 %>%
  dplyr::select(cases.case_id, demographic.days_to_death,
                diagnoses.days_to_last_follow_up) %>%
  mutate(years_to_death = demographic.days_to_death / 365.25,
         years_to_lost_to_follow_up =
           diagnoses.days_to_last_follow_up / 365.25,
         year_last = pmin(years_to_death, years_to_lost_to_follow_up,
                         na.rm = TRUE))

df_keep2 <- df_keep %>%
  slice(rep(1:n(), each = n_time)) %>%
  mutate(time = rep(timepoints, times = N)) %>%
  arrange(cases.case_id, time) %>%
  mutate(structural_miss = time > year_last) %>%
  dplyr::select(cases.case_id, time, structural_miss)

df_tot2 <- df_tot %>%
  left_join(df_keep2, by = c("cases.case_id","time")) %>%
  mutate(
    across(
      c(swb_total, ewb_total, fwb_total, bcs_total, pwb_total, fact_b_total),
      ~ ifelse(structural_miss, NA, .x)))

df_qol <- df_tot2

### Mean % and SD % visualisation for all scores (NOT structurally missing) ####
summary_by_time_all_str_miss <- df_qol %>%
  group_by(time) %>%
  summarise(
    n = n(),
    mean_PWB = mean(pwb_total, na.rm = TRUE),
    sd_PWB = sd(pwb_total, na.rm = TRUE),
    se_PWB = sd_PWB / sqrt(n),
    mean_SWB = mean(swb_total, na.rm = TRUE),
    sd_SWB = sd(swb_total, na.rm = TRUE),
    se_SWB = sd_SWB / sqrt(n),
    mean_EWB = mean(ewb_total, na.rm = TRUE),
    sd_EWB = sd(ewb_total, na.rm = TRUE),
    se_EWB = sd_EWB / sqrt(n),
    mean_FWB = mean(fwb_total, na.rm = TRUE),
    sd_FWB = sd(fwb_total, na.rm = TRUE),
    se_FWB = sd_FWB / sqrt(n),
    mean_BCS = mean(bcs_total, na.rm = TRUE),
    sd_BCS = sd(bcs_total, na.rm = TRUE),
    se_BCS = sd_BCS / sqrt(n),
    mean_FACT_B = mean(fact_b_total, na.rm = TRUE),
    sd_FACT_B = sd(fact_b_total, na.rm = TRUE),
    se_FACT_B = sd_FACT_B / sqrt(n),
  ) %>%
  mutate(
    max_SWB = 28,
    max_EWB = 24,
    max_FWB = 28,
    max_BCS = 40,
    max_PWB = 28,
    max_FACT_B = 148,
    pct_SWB = (mean_SWB/max_SWB)*100,
    pct_EWB = (mean_EWB/max_EWB)*100,
    pct_FWB = (mean_FWB/max_FWB)*100,
    pct_BCS = (mean_BCS/max_BCS)*100,
    pct_PWB = (mean_PWB/max_PWB)*100,
    pct_FACT_B = (mean_FACT_B/max_FACT_B)*100,
    pct_sd_SWB = (sd_SWB/max_SWB)*100,
    pct_sd_EWB = (sd_EWB/max_EWB)*100,
    pct_sd_FWB = (sd_FWB/max_FWB)*100,
    pct_sd_BCS = (sd_BCS/max_BCS)*100,
    pct_sd_PWB = (sd_PWB/max_PWB)*100,
    pct_sd_FACT_B = (sd_FACT_B/max_FACT_B)*100    
  )

plot_df <- summary_by_time_all_str_miss %>%
  dplyr::select(
    time,
    starts_with("pct_"),
  ) %>%
  pivot_longer(
    cols = -time,
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    type = ifelse(grepl("^pct_sd_", variable), "sd", "mean"),
    domain = gsub("^pct_sd_|^pct_", "", variable)
  ) %>%
  dplyr::select(time, domain, type, value) %>%
  pivot_wider(
    names_from = type,
    values_from = value
  )

plot_df$domain <- factor(
  plot_df$domain,
  levels = c("PWB", "SWB", "EWB", "FWB", "BCS", "FACT_B")
)

p_all_scores_str_miss <-
  ggplot(plot_df, aes(x = time, y = mean, color = domain)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    width = 0.08,
    alpha = 0.7
  ) +
  scale_y_continuous(
    limits = c(15, 95),
    breaks = seq(15, 95, 10)
  ) +
  scale_x_continuous(breaks = sort(unique(plot_df$time))) +
  scale_color_manual(values = okabe_ito) +
  labs(
    title = "",
    x = "Time (years)",
    y = "% of maximum score",
    color = "Domain"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 12)
  )

#p_all_scores_str_miss

# Visualise the proportion on measurements that are structurally missing
df_struct_summary <- df_qol %>%
  group_by(time) %>%
  summarise(
    n_total = n(),
    n_struct_missing = sum(structural_miss),
    n_available = sum(!structural_miss),
    prop_struct_missing = n_struct_missing / n_total,
    prop_available = n_available / n_total
  )

df_struct_long <- df_struct_summary %>%
  dplyr::select(time, prop_available, prop_struct_missing) %>%
  tidyr::pivot_longer(
    cols = c(prop_available, prop_struct_missing),
    names_to = "type",
    values_to = "proportion"
  )

p_miss <- ggplot(df_struct_long, aes(x = factor(time),
                                     y = proportion, fill = type)) +
  geom_col(position = "stack") +
  geom_text(
    data = df_struct_summary,
    aes(
      x = factor(time),
      y = 1.02,
      label = paste0("n = ", n_available)
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(
    values = c(
      prop_available = "green",
      prop_struct_missing = "red"
    ),
    labels = c(
      prop_available = "Available",
      prop_struct_missing = "Structurally Missing"
    )
  ) +
  labs(
    x = "Timepoint",
    y = "Proportion",
    fill = "Status",
    title = "Structured Missingness Over Time",
    subtitle = "Percentage of available vs structurally missing QoL measurements"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#print(p_miss)
