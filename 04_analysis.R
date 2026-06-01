### Libraries ##################################################################
library(survival)
library(survminer)
library(dplyr)
library(ggplot2)
library(broom)
library(nlme)
library(JMbayes2)
library(bayesplot)

# ------------------------------------------------------------------------------
# 1. PFS analysis
# ------------------------------------------------------------------------------

### Prepare data for Survival analysis
df_pfs <- df_final_1 %>%
  mutate(
    pfs_event = ifelse(!is.na(pfs_event.d), 1, 0),
    pfs_time = ifelse(
      pfs_event == 1,
      pfs_event.d,
      diagnoses.days_to_last_follow_up
    ),
    pfs_time_years = as.numeric(pfs_time) / 365.25,
    # Get only main stage
    stage_raw = as.character(diagnoses.ajcc_pathologic_stage),
    stage_clean = trimws(gsub("(?i)^Stage\\s+", "", stage_raw, perl = TRUE)),
    stage_group = case_when(
      is.na(stage_clean) ~ "X",
      grepl("^IV", stage_clean) ~ "IV",
      grepl("^III", stage_clean) ~ "III",
      grepl("^II", stage_clean) ~ "II",
      grepl("^I", stage_clean) ~ "I",
      TRUE ~ "X"
    ),
    # Group ages
    age_group = case_when(
      demographic.age_at_index < 50 ~ "<50",
      demographic.age_at_index < 65 ~ "50–64",
      TRUE ~ "65+"
    ),
    # "Missing" and "Indeterminate" coded as "Missing"
    margin_status = case_when(
      first_line.margin_status == "Uninvolved" ~ "Uninvolved",
      first_line.margin_status == "Involved" ~ "Involved",
      TRUE ~ "Missing"
    ),
    # Add "American Indian or Alaska Native" to "Missing" due to n<3
    demographic.race = case_when(
      demographic.race == "American Indian or Alaska Native" ~ "Missing",
      TRUE ~ demographic.race
    )
  ) %>%
  dplyr::select(-stage_raw, -stage_clean)

### var_dict specification #####################################################

### OVERALL                     Variable name                   Number
# overall                       1                               0

### CORE PROGNOSTIC FACTORS     Variable name                   Number
# stage_group                   stage_group                     1
# age_group                     age_group                       2
# subtype                       molecular_subtype               3
# margin_status                 margin_status                   4
# prior_malignancy              diagnoses.prior_malignancy      5

### TREATMENT EXPOSURES         Variable name                    Number
# chemo.fl                      chemo.fl                        6
# rad.fl                        rad.fl                          7
# serm.fl                       serm.fl                         8
# ai.fl                         ai.fl                           9
# serd.fl                       serd.fl                         10
# lhrh_a.fl                     lhrh_a.fl                       11
# immuno.fl                     immuno.fl                       12

### ADDITIONAL FACTORS          Variable name                   Number
# ethnicity                     demographic.ethnicity           13
# race                          demographic.race                14
# meno.status                   meno.status                     15

var_dict <- list(
  
  # 0
  overall = list(
    var = "1",
    label = "Overall PFS",
    palette_n = 1
  ),
  
  # 1
  age_group = list(
    var = "age_group",
    label = "Age group",
    palette_n = 3
  ),
  
  # 2
  stage_group = list(
    var = "stage_group",
    label = "AJCC Stage",
    palette_n = 5
  ),
  
  # 3
  subtype = list(
    var = "molecular_subtype",
    label = "Molecular subtype",
    palette_n = 6
  ),
  
  # 4
  margin_status = list(
    var = "margin_status",
    label = "Margin status",
    palette_n = 3
  ),
  
  # 5
  prior_malignancy = list(
    var = "diagnoses.prior_malignancy",
    label = "Prior malignancy",
    palette_n = 2
  ),

  # 6
  chemo.fl = list(
    var = "chemo.fl",
    label = "Chemotherapy",
    palette_n = 2
  ),
  
  # 7
  rad.fl = list(
    var = "rad.fl",
    label = "Radiation therapy",
    palette_n = 2
  ), 
  
  # 8
  serm.fl = list(
    var = "serm.fl",
    label = "Selective estrogen receptor modulators",
    palette_n = 2
  ), 
  
  # 9
  ai.fl = list(
    var = "ai.fl",
    label = "Aromatase inhibitors",
    palette_n = 2
  ), 
  
  # 10
  serd.fl = list(
    var = "serd.fl",
    label = "Estrogen receptor downregulator",
    palette_n = 2
  ), 
  
  # 11
  lhrh_a.fl = list(
    var = "lhrh_a.fl",
    label = "LHRH (GnRH) agonists",
    palette_n = 2
  ), 
  
  # 12
  immuno.fl = list(
    var = "immuno.fl",
    label = "Immunotherapy",
    palette_n = 2
  ),
  
  # 13
  ethnicity = list(
    var = "demographic.ethnicity",
    label = "Ethnicity",
    palette_n = 3
  ),
  
  # 14
  race = list(
    var = "demographic.race",
    label = "Race",
    palette_n = 4
  ),
  
  # 15
  meno.status = list(
    var = "meno.status",
    label = "Menopausal status",
    palette_n = 4
  )
)

### General check for selected variables and KP survival curve #################

plot_pfs_km <- function(v) {
  
  if (v$var == "1") {
    fml <- Surv(pfs_time_years, pfs_event) ~ 1
  } else {
    fml <- as.formula(
      paste0("Surv(pfs_time_years, pfs_event) ~ ", v$var)
    )
  }
  
  fit <- eval(bquote(survfit(.(fml), data = df_pfs)))
  
  okabe_ito <- c(
    "#000000", "#E69F00", "#56B4E9", "#009E73",
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
  )
  
  pval_flag <- v$var != "1"
  
  ggsurvplot(
    fit,
    data = df_pfs,
    conf.int = FALSE,
    censor.shape = "|",
    censor.size = 1.5,
    palette = okabe_ito[1:v$palette_n],
    ggtheme = theme_minimal(base_size = 12),
    xlab = "Time (years)",
    ylab = "PFS probability",
    ylim = c(0, 1),
    surv.scale = "percent",
    risk.table = TRUE,
    risk.table.height = 0.30,
    risk.table.y.text.col = TRUE,
    risk.table.y.text = FALSE,
    pval = pval_flag,
    pval.coord = c(0, 0.1)
  )
}

### OVERALL                     Variable name                   Number
# overall                       1                               0
plot_pfs_km(var_dict$overall)

### CORE PROGNOSTIC FACTORS     Variable name                   Number
# stage_group                   stage_group                     1
# age_group                     age_group                       2
# subtype                       molecular_subtype               3
# margin_status                 margin_status                   4
# prior_malignancy              diagnoses.prior_malignancy      5
plot_pfs_km(var_dict$stage_group)
plot_pfs_km(var_dict$age_group)
plot_pfs_km(var_dict$subtype)
plot_pfs_km(var_dict$margin_status)
plot_pfs_km(var_dict$prior_malignancy)

### TREATMENT EXPOSURES         Variable name                    Number
# chemo.fl                      chemo.fl                        6
# rad.fl                        rad.fl                          7
# serm.fl                       serm.fl                         8
# ai.fl                         ai.fl                           9
# serd.fl                       serd.fl                         10
# lhrh_a.fl                     lhrh_a.fl                       11
# immuno.fl                     immuno.fl                       12
plot_pfs_km(var_dict$chemo.fl)
plot_pfs_km(var_dict$rad.fl)
plot_pfs_km(var_dict$serm.fl)
plot_pfs_km(var_dict$ai.fl)
plot_pfs_km(var_dict$serd.fl)
plot_pfs_km(var_dict$lhrh_a.fl)
plot_pfs_km(var_dict$immuno.fl)

### ADDITIONAL FACTORS          Variable name                   Number
# ethnicity                     demographic.ethnicity           13
# race                          demographic.race                14
# meno.status                   meno.status                     15
plot_pfs_km(var_dict$ethnicity)
plot_pfs_km(var_dict$race)
plot_pfs_km(var_dict$meno.status)

### Log-rank test for selected variables ####################################### 

# Log-rank test (PFS ~ variable)
# H_0: PFS curves are identical across groups
# H_1: At least one stage group differs

logrank_pfs <- function(v) {
  fml <- as.formula(paste0("Surv(pfs_time_years, pfs_event) ~ ", v$var))
  survdiff(fml, data = df_pfs)
}

### CORE PROGNOSTIC FACTORS     Variable name                   Number
# stage_group                   stage_group                     1
# age_group                     age_group                       2
# subtype                       molecular_subtype               3
# margin_status                 margin_status                   4
# prior_malignancy              diagnoses.prior_malignancy      5
logrank_pfs(var_dict$stage_group)       # p<0.01    diff
logrank_pfs(var_dict$age_group)         # p<0.01    diff
logrank_pfs(var_dict$subtype)           # p=0.03    diff
logrank_pfs(var_dict$margin_status)     # p<0.01    diff
logrank_pfs(var_dict$prior_malignancy)  # p=0.2     NO diff

### TREATMENT EXPOSURES         Variable name                    Number
# chemo.fl                      chemo.fl                        6
# rad.fl                        rad.fl                          7
# serm.fl                       serm.fl                         8
# ai.fl                         ai.fl                           9
# serd.fl                       serd.fl                         10
# lhrh_a.fl                     lhrh_a.fl                       11
# immuno.fl                     immuno.fl                       12
logrank_pfs(var_dict$chemo.fl)          # p<0.01    diff
logrank_pfs(var_dict$rad.fl)            # p<0.01    diff
logrank_pfs(var_dict$serm.fl)           # p=0.06    diff
logrank_pfs(var_dict$ai.fl)             # p=0.02    diff
logrank_pfs(var_dict$serd.fl)           # p<0.01    diff
logrank_pfs(var_dict$lhrh_a.fl)         # p=0.3     NO diff
logrank_pfs(var_dict$immuno.fl)         # p=0.4     NO diff

### ADDITIONAL FACTORS          Variable name                   Number
# ethnicity                     demographic.ethnicity           13
# race                          demographic.race                14
# meno.status                   meno.status                     15
logrank_pfs(var_dict$ethnicity)         # p=0.01    diff
logrank_pfs(var_dict$race)              # p=0.7     NO diff
logrank_pfs(var_dict$meno.status)       # p<0.01    diff

### Univariate Cox proportional hazards models #################################

cox_uni <- function(v, data = df_pfs) {
  fml <- as.formula(paste0("Surv(pfs_time_years, pfs_event) ~ ", v$var))
  fit <- survival::coxph(fml, data = data)
  summary(fit)
}

### CORE PROGNOSTIC FACTORS     Variable name                   Number
# stage_group                   stage_group                     1
# age_group                     age_group                       2
# subtype                       molecular_subtype               3
# margin_status                 margin_status                   4
# prior_malignancy              diagnoses.prior_malignancy      5

# Relevel the stage_group so that I is a reference
df_pfs$stage_group <- factor(
  df_pfs$stage_group,
  levels = c("I", "II", "III", "IV", "X")
)
cox_uni(var_dict$stage_group)
# Reference group - Stage I
# Stage II v ref: HR = 1.50, 95% CI = 0.89, 2.51, p=0.126
# Stage III v ref: HR = 2.97, 95% CI = 1.73, 5.10, p<0.01
# Stage IV v ref: HR = 10.33, 95% CI = 5.13, 20.8, p<0.01
# Stage X v ref: HR = 6.11, 95% CI = 2.99, 12.47, p<0.01
# Missing group closer to advanced stage

# Check without missing stage (X)
df_temp <- df_pfs[!df_pfs$stage_group %in% c("Missing", "Unknown", "X"), ]
# Relevel without "X"
df_temp$stage_group <- relevel(
  factor(df_temp$stage_group),
  ref = "I"
)
cox_uni(var_dict$stage_group, data=df_temp)
# Reference group - Stage I
# Stage II v ref: HR = 1.50, 95% CI = 0.90, 2.52, p=0.12
# Stage III v ref: HR = 2.98, 95% CI = 1.74, 5.11, p<0.01
# Stage IV v ref: HR = 10.29, 95% CI = 5.11, 20.71, p<0.01

cox_uni(var_dict$age_group)
# Reference group - <50
# 50-65 v ref: HR = 0.93, 95% CI = 0.64, 1.36, p=0.71
# 65+ v ref: HR = 1.71, 95% CI = 1.17, 2.50, p<0.01

# Relevel the molecular subtype so that Luminal A-Like is a reference
df_pfs$molecular_subtype <- relevel(
  factor(df_pfs$molecular_subtype),
  ref = "Luminal A-like"
)

cox_uni(var_dict$subtype)
# Reference group: Luminal A-like
# HER2-enriched v ref: HR = 2.01, 95% CI = 0.96, 4.22, p=0.06
# Luminal B–like (HER2–neg) v ref: HR = 1.24, 95% CI = 0.65, 2.37, p=0.51
# Luminal B–like (HER2–pos) v ref: HR = 1.59, 95% CI = 0.98, 2.58, p=0.06
# Triple-negative v ref: HR = 1.41, 95% CI = 0.88, 2.25 , p=0.15
# Unknown v ref: HR = 1.89, 95% CI = 1.27, 2.82, p<0.01
# Unknown group worst outcome

# Check without missing subtype (Unknown)
df_temp <- df_pfs[!df_pfs$molecular_subtype %in% c("Missing", "Unknown", "X"), ]
# Relevel without "Unknown"
df_temp$molecular_subtype <- relevel(
  factor(df_temp$molecular_subtype),
  ref = "Luminal A-like"
)

cox_uni(var_dict$subtype, data=df_temp)
# Reference group: Luminal A-like
# HER2-enriched v ref: HR = 2.04, 95% CI = 0.97, 4.28, p=0.06
# Luminal B–like (HER2–neg) v ref: HR = 1.23, 95% CI = 0.65, 2.35, p=0.53
# Luminal B–like (HER2–pos) v ref: HR = 1.59, 95% CI = 0.98, 2.58, p=0.06
# Triple-negative v ref: HR = 1.42, 95% CI = 0.89, 2.27 , p=0.15

cox_uni(var_dict$margin_status)
# Reference group - Involved
# Missing v ref: HR = 1.58, 95% CI = 0.94, 2.67, p=0.09
# Uninvolved v ref: HR = 0.52, 95% CI = 0.33, 0.84, p<0.01

# Check without missing margin status (Missing)
df_temp <- df_pfs[!df_pfs$margin_status %in% c("Missing", "Unknown", "X"), ]

cox_uni(var_dict$margin_status, data=df_temp)
# Reference group - Involved
# Uninvolved v ref: HR = 0.52, 95% CI = 0.32, 0.83, p<0.01

cox_uni(var_dict$prior_malignancy)
# Reference group - N
# Y v ref: HR = 1.62, 95% CI = 0.71, 3.66, p=0.25

### ADDITIONAL FACTORS          Variable name                   Number
# ethnicity                     demographic.ethnicity           13
# race                          demographic.race                14
# meno.status                   meno.status                     15

cox_uni(var_dict$ethnicity)
# Reference group - Hispanic or Latino
# Missing v ref: HR = 1.82, 95% CI = 0.40, 8.23, p=0.44
# Not Hispanic or Latino v ref: HR = 3.75, 95% CI = 0.93, 15.13, p=0.06

# Relevel the race so that White is a reference
df_pfs$demographic.race <- relevel(
  factor(df_pfs$demographic.race),
  ref = "White"
)

cox_uni(var_dict$race)
# Reference group - White
# Asian v ref: HR = 0.85, 95% CI = 0.31, 2.30, p=0.75
# Black or African American v ref: HR = 1.26, 95% CI = 0.86, 1.84, p=0.24
# Missing v ref: HR = 0.95, 95% CI = 0.44, 2.05, p=0.89

# Relevel the menopausal status so that Premenopausal is a reference
df_pfs$meno.status <- relevel(
  factor(df_pfs$meno.status),
  ref = "Premenopausal"
)

cox_uni(var_dict$meno.status)
# Reference group - Premenopausal
# Perimenopausal v ref: HR = 0.23, 95% CI = 0.03, 1.68, p=0.15
# Postmenopausal v ref: HR = 1.50, 95% CI = 0.99, 2.28, p=0.05
# Unknown v ref: HR = 2.38, 95% CI = 1.47, 3.85, p<0.01

# Check without missing menopausal status (Unknown)
df_temp <- df_pfs[!df_pfs$meno.status %in% c("Missing", "Unknown", "X"), ]
# Relevel without "Unknown"
df_temp$meno.status <- relevel(
  factor(df_temp$meno.status),
  ref = "Premenopausal"
)

cox_uni(var_dict$meno.status, data=df_temp)
# Reference group - Premenopausal
# Perimenopausal v ref: HR = 0.23, 95% CI = 0.03, 1.68, p=0.15
# Postmenopausal v ref: HR = 1.48, 95% CI = 0.98, 2.25, p=0.06

### Multivariate Cox proportional hazards models ###############################

df_pfs <- droplevels(df_pfs)

# Initial model
fml_mv <- Surv(pfs_time_years, pfs_event) ~ 
  stage_group +
  age_group +
  molecular_subtype +
  margin_status +
  diagnoses.prior_malignancy +
  demographic.race +
  meno.status

cox_mv <- survival::coxph(fml_mv, data = df_pfs)
summary(cox_mv)

# Reduced model (-race)
fml_mv2 <- Surv(pfs_time_years, pfs_event) ~ 
  stage_group +
  age_group +
  molecular_subtype +
  margin_status +
  diagnoses.prior_malignancy +
  meno.status

cox_mv2 <- survival::coxph(fml_mv2, data = df_pfs)
summary(cox_mv2)

# Reduced model (-prior malignancy status)
fml_mv3 <- Surv(pfs_time_years, pfs_event) ~ 
  stage_group +
  age_group +
  molecular_subtype +
  margin_status +
  meno.status

cox_mv3 <- survival::coxph(fml_mv3, data = df_pfs)
summary(cox_mv3)

# Reduced model (-prior malignancy status)
fml_mv4 <- Surv(pfs_time_years, pfs_event) ~ 
  stage_group +
  age_group +
  molecular_subtype +
  meno.status

cox_mv4 <- survival::coxph(fml_mv4, data = df_pfs)
summary(cox_mv4)

# Reduced model (-menopausal status)
fml_mv5 <- Surv(pfs_time_years, pfs_event) ~ 
  stage_group +
  age_group +
  molecular_subtype

cox_mv5 <- survival::coxph(fml_mv5, data = df_pfs)
summary(cox_mv5)

# Reduced model (-age)
fml_mv6 <- Surv(pfs_time_years, pfs_event) ~ 
  stage_group +
  molecular_subtype

cox_mv6 <- survival::coxph(fml_mv6, data = df_pfs)
summary(cox_mv6)

# Time dependence check
cox.zph(cox_mv6)
plot(cox.zph(cox_mv6))

# Final model
fml_final <- Surv(pfs_time_years, pfs_event) ~ 
  stage_group +
  molecular_subtype
cox_final <- survival::coxph(fml_final, data = df_pfs)


tbl_final <- tidy(cox_final, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    term = gsub("stage_group", "Stage ", term),
    term = gsub("molecular_subtype", "Subtype: ", term)
  )

tbl_final <- tbl_final %>%
  mutate(
    term = factor(
      term,
      levels = c(
        "Stage II", "Stage III", "Stage IV", "Stage X",
        "Subtype: HER2-enriched",
        "Subtype: Luminal B-like (HER2-neg)",
        "Subtype: Luminal B-like (HER2-pos)",
        "Subtype: Triple-negative",
        "Subtype: Unknown"
      )
    )
  )

# Forest plot
tbl_final$label <- sprintf("p = %.3f", tbl_final$p.value)

ggplot(tbl_final, aes(x = term, y = estimate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.15) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_text(aes(label = label), hjust = 0, vjust = -1, size = 3.5) +
  coord_flip(clip = "off") +
  labs(
    x = "",
    y = "Hazard Ratio (95% CI)",
    title = "Multivariable Cox Model: Stage and Molecular Subtype"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(hjust = 0),
    plot.margin = margin(5.5, 40, 5.5, 5.5)
  )

# Model summary table
summary_tbl <- tidy(cox_final, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    term = case_when(
      grepl("stage_group", term) ~ gsub("stage_group", "Stage ", term),
      grepl("molecular_subtype", term) ~ gsub("molecular_subtype", "Subtype: ", term),
      TRUE ~ term
    ),
    HR = round(estimate, 2),
    CI = sprintf("%.2f–%.2f", conf.low, conf.high),
    p = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
  ) %>%
  dplyr::select(term, HR, CI, p)

summary_tbl

# ------------------------------------------------------------------------------
# 2. Analysis with QoL summary statistic
# ------------------------------------------------------------------------------

### Summary statistic - mean over the first year of FU #########################
# Values to be tested - PWB, BCS, FACT-B total

run_mean_models <- function(df_in) {
  
  # Summary derivation  
  df_qol_summary_mean <- df_in %>%
    filter(time_num <= 1 & !is.na(fact_b_total)) %>%
    group_by(cases.case_id) %>%
    summarise(
      PWB_mean_0_1         = mean(pwb_total, na.rm = TRUE),
      BCS_mean_0_1         = mean(bcs_total, na.rm = TRUE),
      FACTB_total_mean_0_1 = mean(fact_b_total, na.rm = TRUE),
      n_obs_0_1            = n(),
      .groups = "drop"
    )
  
  # Merge with PFS data
  df_surv_qol_mean <- df_pfs %>%
    left_join(df_qol_summary_mean, by = "cases.case_id") %>%
    filter(!is.na(FACTB_total_mean_0_1)) %>%
    dplyr::select(
      cases.case_id, pfs_time_years, pfs_event, stage_group,
      molecular_subtype, PWB_mean_0_1, BCS_mean_0_1,
      FACTB_total_mean_0_1, n_obs_0_1
    )
  
  # Formulae and names lists 
  formulas_mean <- list(
    # With stage and molecular subtype
    Surv(pfs_time_years, pfs_event) ~
      stage_group + molecular_subtype + PWB_mean_0_1 + BCS_mean_0_1,
    Surv(pfs_time_years, pfs_event) ~
      stage_group + molecular_subtype + PWB_mean_0_1,
    Surv(pfs_time_years, pfs_event) ~
      stage_group + molecular_subtype + BCS_mean_0_1,
    Surv(pfs_time_years, pfs_event) ~
      stage_group + molecular_subtype + FACTB_total_mean_0_1,
    # With stage
    Surv(pfs_time_years, pfs_event) ~
      stage_group + PWB_mean_0_1 + BCS_mean_0_1,
    Surv(pfs_time_years, pfs_event) ~
      stage_group + PWB_mean_0_1,
    Surv(pfs_time_years, pfs_event) ~
      stage_group + BCS_mean_0_1,
    Surv(pfs_time_years, pfs_event) ~
      stage_group + FACTB_total_mean_0_1,
    # Only QoL
    Surv(pfs_time_years, pfs_event) ~
      PWB_mean_0_1 + BCS_mean_0_1,
    Surv(pfs_time_years, pfs_event) ~
      PWB_mean_0_1,
    Surv(pfs_time_years, pfs_event) ~
      BCS_mean_0_1,
    Surv(pfs_time_years, pfs_event) ~
      FACTB_total_mean_0_1
  )
  
  model_names <- c(
    # With stage and molecular subtype
    "Model1: PWB + BCS",
    "Model2: PWB",
    "Model3: BCS",
    "Model4: FACT-B total",
    # With stage
    "Model5: PWB + BCS (no subtype)",
    "Model6: PWB (no subtype)",
    "Model7: BCS (no subtype)",
    "Model8: FACT-B total (no subtype)",
    # Only QoL
    "Model9: PWB + BCS (only)",
    "Model10: PWB (only)",
    "Model11: BCS (only)",
    "Model12: FACT-B total (only)"
  )
  
  # Models fit
  models <- lapply(formulas_mean, function(f) {
    fit <- survival::coxph(f, data = df_surv_qol_mean)
    list(model = fit, summary = summary(fit))
  })
  names(models) <- model_names
  
  return(models)
}

res_mean <- run_mean_models(df_qol)
res_mean_mcar <- run_mean_models(df_mcar)
res_mean_mar  <- run_mean_models(df_mar)
res_mean_mnar <- run_mean_models(df_mnar)

### Summary statistic - slope over the first year of FU ########################
# Values to be tested - PWB, BCS, FACT-B total

run_slope_models <- function(df_in) {
  
  # Summary derivation
  df_qol_summary_slope <- df_in %>%
    filter(time_num <= 1 & !is.na(fact_b_total)) %>%
    group_by(cases.case_id) %>%
    summarise(
      PWB_slope_0_1 = {
        d <- pick(time_num, pwb_total) 
        if (nrow(d) < 2) {
          NA_real_
        } else {
          coef(lm(pwb_total ~ time_num, data = d))[2]
        }
      },
      BCS_slope_0_1 = {
        d <- pick(time_num, bcs_total) 
        if (nrow(d) < 2) {
          NA_real_
        } else {
          coef(lm(bcs_total ~ time_num, data = d))[2]
        }
      },
      FACTB_total_slope_0_1 = {
        d <- pick(time_num, fact_b_total) 
        if (nrow(d) < 2) {
          NA_real_
        } else {
          coef(lm(fact_b_total ~ time_num, data = d))[2]
        }
      },
      n_obs_0_1           = n(),
      .groups = "drop"
    )
  
  # PFS data merge
  df_surv_qol_slope <- df_pfs %>%
    left_join(df_qol_summary_slope, by ="cases.case_id") %>%
    filter(!is.na(FACTB_total_slope_0_1)) %>%
    dplyr::select(cases.case_id, pfs_time_years, pfs_event, stage_group,
                  molecular_subtype, PWB_slope_0_1, BCS_slope_0_1,
                  FACTB_total_slope_0_1, n_obs_0_1)
  
  # Formulae and names lists 
  formulas_slope <- list(
    # With stage and molecular subtype
    Surv(pfs_time_years, pfs_event) ~
      stage_group + molecular_subtype + PWB_slope_0_1 + BCS_slope_0_1,
    Surv(pfs_time_years, pfs_event) ~
      stage_group + molecular_subtype + PWB_slope_0_1,
    Surv(pfs_time_years, pfs_event) ~
      stage_group + molecular_subtype + BCS_slope_0_1,
    Surv(pfs_time_years, pfs_event) ~
      stage_group + molecular_subtype + FACTB_total_slope_0_1,
    # With stage
    Surv(pfs_time_years, pfs_event) ~
      stage_group + PWB_slope_0_1 + BCS_slope_0_1,
    Surv(pfs_time_years, pfs_event) ~
      stage_group + PWB_slope_0_1,
    Surv(pfs_time_years, pfs_event) ~
      stage_group + BCS_slope_0_1,
    Surv(pfs_time_years, pfs_event) ~
      stage_group + FACTB_total_slope_0_1,
    # Only QoL
    Surv(pfs_time_years, pfs_event) ~
      PWB_slope_0_1 + BCS_slope_0_1,
    Surv(pfs_time_years, pfs_event) ~
      PWB_slope_0_1,
    Surv(pfs_time_years, pfs_event) ~
      BCS_slope_0_1,
    Surv(pfs_time_years, pfs_event) ~
      FACTB_total_slope_0_1
  )
  
  model_names <- c(
    # With stage and molecular subtype
    "Model1: PWB + BCS",
    "Model2: PWB",
    "Model3: BCS",
    "Model4: FACT-B total",
    # With stage
    "Model5: PWB + BCS (no subtype)",
    "Model6: PWB (no subtype)",
    "Model7: BCS (no subtype)",
    "Model8: FACT-B total (no subtype)",
    # Only QoL
    "Model9: PWB + BCS (only)",
    "Model10: PWB (only)",
    "Model11: BCS (only)",
    "Model12: FACT-B total (only)"
  )
  
  # Models fit
  models <- lapply(formulas_slope, function(f) {
    fit <- survival::coxph(f, data = df_surv_qol_slope)
    list(model = fit, summary = summary(fit))
  })
  names(models) <- model_names
  
  return(models)
}

res_slope <- run_slope_models(df_qol)
res_slope_mcar <- run_slope_models(df_mcar)
res_slope_mar  <- run_slope_models(df_mar)
res_slope_mnar <- run_slope_models(df_mnar)

### Summary statistic - AUC over the first year of FU ##########################
# Values to be tested - PWB, BCS, FACT-B total

run_auc_models <- function(df_in) {
  
  # Summary derivation
  df_qol_summary_auc <- df_in %>%
    filter(time_num <= 1 & !is.na(fact_b_total)) %>%
    group_by(cases.case_id) %>%
    summarise(
      PWB_auc_0_1 = {
        d <- pick(time_num, pwb_total)
        if (nrow(d) < 2) NA_real_
        else
          sum((d$pwb_total[-1] + d$pwb_total[-nrow(d)]) * diff(d$time_num) / 2)
      },
      BCS_auc_0_1 = {
        d <- pick(time_num, bcs_total)
        if (nrow(d) < 2) NA_real_
        else
          sum((d$bcs_total[-1] + d$bcs_total[-nrow(d)]) * diff(d$time_num) / 2)
      },
      FACTB_total_auc_0_1 = {
        d <- pick(time_num, fact_b_total)
        if (nrow(d) < 2) NA_real_
        else
          sum((d$fact_b_total[-1] +
                 d$fact_b_total[-nrow(d)]) * diff(d$time_num) / 2)
      },
      n_obs_0_1 = n(),
      .groups = "drop"
    )
  
  # PFS data merge
  df_surv_qol_auc <- df_pfs %>%
    left_join(df_qol_summary_auc, by ="cases.case_id") %>%
    filter(!is.na(FACTB_total_auc_0_1)) %>%
    dplyr::select(cases.case_id, pfs_time_years, pfs_event, stage_group,
                  molecular_subtype, PWB_auc_0_1, BCS_auc_0_1,
                  FACTB_total_auc_0_1, n_obs_0_1)
  
  # Formulae and names lists 
  formulas_auc <- list(
    # With stage and molecular subtype
    Surv(pfs_time_years, pfs_event) ~
      stage_group + molecular_subtype + PWB_auc_0_1 + BCS_auc_0_1,
    Surv(pfs_time_years, pfs_event) ~
      stage_group + molecular_subtype + PWB_auc_0_1,
    Surv(pfs_time_years, pfs_event) ~
      stage_group + molecular_subtype + BCS_auc_0_1,
    Surv(pfs_time_years, pfs_event) ~
      stage_group + molecular_subtype + FACTB_total_auc_0_1,
    # With stage
    Surv(pfs_time_years, pfs_event) ~
      stage_group + PWB_auc_0_1 + BCS_auc_0_1,
    Surv(pfs_time_years, pfs_event) ~
      stage_group + PWB_auc_0_1,
    Surv(pfs_time_years, pfs_event) ~
      stage_group + BCS_auc_0_1,
    Surv(pfs_time_years, pfs_event) ~
      stage_group + FACTB_total_auc_0_1,
    # Only QoL
    Surv(pfs_time_years, pfs_event) ~
      PWB_auc_0_1 + BCS_auc_0_1,
    Surv(pfs_time_years, pfs_event) ~
      PWB_auc_0_1,
    Surv(pfs_time_years, pfs_event) ~
      BCS_auc_0_1,
    Surv(pfs_time_years, pfs_event) ~
      FACTB_total_auc_0_1
  )
  
  model_names <- c(
    # With stage and molecular subtype
    "Model1: PWB + BCS",
    "Model2: PWB",
    "Model3: BCS",
    "Model4: FACT-B total",
    # With stage
    "Model5: PWB + BCS (no subtype)",
    "Model6: PWB (no subtype)",
    "Model7: BCS (no subtype)",
    "Model8: FACT-B total (no subtype)",
    # Only QoL
    "Model9: PWB + BCS (only)",
    "Model10: PWB (only)",
    "Model11: BCS (only)",
    "Model12: FACT-B total (only)"
  )

  # Models fit
  models <- lapply(formulas_auc, function(f) {
    fit <- survival::coxph(f, data = df_surv_qol_auc)
    list(model = fit, summary = summary(fit))
  })
  names(models) <- model_names
  
  return(models)
}

res_auc <- run_auc_models(df_qol)
res_auc_mcar <- run_auc_models(df_mcar)
res_auc_mar  <- run_auc_models(df_mar)
res_auc_mnar <- run_auc_models(df_mnar)

### Time dependence diagnostics - Summary statistics models ####################
run_ph_diagnostics <- function(res_in) {
  
  model_names <- c(
    # With stage and molecular subtype
    "Model1: PWB + BCS",
    "Model2: PWB",
    "Model3: BCS",
    "Model4: FACT-B total",
    # With stage
    "Model5: PWB + BCS (no subtype)",
    "Model6: PWB (no subtype)",
    "Model7: BCS (no subtype)",
    "Model8: FACT-B total (no subtype)",
    # Only QoL
    "Model9: PWB + BCS (only)",
    "Model10: PWB (only)",
    "Model11: BCS (only)",
    "Model12: FACT-B total (only)"
  ) 
  
  for (m in model_names) {
    cat("Proportional Hazards Check for:", m)
    
    fit <- res_in[[m]]$model
    ph_test <- cox.zph(fit)
    print(ph_test)
    plot(ph_test, main = paste("PH Diagnostics:", m))
  }
}

# Mean
# run_ph_diagnostics(res_mean)
# run_ph_diagnostics(res_mean_mcar)
# run_ph_diagnostics(res_mean_mar)
# run_ph_diagnostics(res_mean_mnar)

# Slope
# run_ph_diagnostics(res_slope)
# run_ph_diagnostics(res_slope_mcar)
# run_ph_diagnostics(res_slope_mar)
# run_ph_diagnostics(res_slope_mnar)

# AUC
# run_ph_diagnostics(res_auc)
# run_ph_diagnostics(res_auc_mcar)
# run_ph_diagnostics(res_auc_mar)
# run_ph_diagnostics(res_auc_mnar)


# ------------------------------------------------------------------------------
# 3. Analysis with time-varying model
# ------------------------------------------------------------------------------

# Prepare time-varying data
df_pfs_1 <- df_pfs %>%
  filter(pfs_time_years > 0)

tv_base <- tmerge(
  data1 = df_pfs_1,
  data2 = df_pfs_1,
  id = cases.case_id,
  tstart = 0,
  tstop = pfs_time_years,
  event = event(pfs_time_years, pfs_event)
)

run_tv_models <- function(df_in) {
  
  add_lag <- function(qol_name) {
    tmerge(
      data1 = tv_base,
      data2 = df_in,
      id = cases.case_id,
      var_lag = tdc(time_num, df_in[[qol_name]])
    )
  }
  
  tv_pwb  <- add_lag("pwb_total") %>%
    rename(PWB_lag = var_lag)
 
  tv_bcs  <- add_lag("bcs_total") %>%
    rename(BCS_lag = var_lag) %>%
    dplyr::select(cases.case_id, tstart, tstop, BCS_lag)
  
  tv_factb  <- add_lag("fact_b_total") %>%
    rename(FACTB_lag = var_lag) %>%
    dplyr::select(cases.case_id, tstart, tstop, FACTB_lag)
  
  tv_all <- tv_pwb %>%
    left_join(tv_bcs, by = c("cases.case_id","tstart","tstop")) %>%
    left_join(tv_factb, by = c("cases.case_id","tstart","tstop"))
  
  # Formulae and names lists 
  formulas_tv <- list(
    # With stage and molecular subtype
    Surv(tstart, tstop, event) ~
      stage_group + molecular_subtype + PWB_lag + BCS_lag,
    Surv(tstart, tstop, event) ~
      stage_group + molecular_subtype + PWB_lag,
    Surv(tstart, tstop, event) ~
      stage_group + molecular_subtype + BCS_lag,
    Surv(tstart, tstop, event) ~
      stage_group + molecular_subtype + FACTB_lag,
    # With stage
    Surv(tstart, tstop, event) ~
      stage_group + PWB_lag + BCS_lag,
    Surv(tstart, tstop, event) ~
      stage_group + PWB_lag,
    Surv(tstart, tstop, event) ~
      stage_group + BCS_lag,
    Surv(tstart, tstop, event) ~
      stage_group + FACTB_lag,
    # Only QoL
    Surv(tstart, tstop, event) ~
      PWB_lag + BCS_lag,
    Surv(tstart, tstop, event) ~
      PWB_lag,
    Surv(tstart, tstop, event) ~
      BCS_lag,
    Surv(tstart, tstop, event) ~
      FACTB_lag
  )
  
  model_names <- c(
    # With stage and molecular subtype
    "Model1: PWB + BCS",
    "Model2: PWB",
    "Model3: BCS",
    "Model4: FACT-B total",
    # With stage
    "Model5: PWB + BCS (no subtype)",
    "Model6: PWB (no subtype)",
    "Model7: BCS (no subtype)",
    "Model8: FACT-B total (no subtype)",
    # Only QoL
    "Model9: PWB + BCS (only)",
    "Model10: PWB (only)",
    "Model11: BCS (only)",
    "Model12: FACT-B total (only)"
  )
  
  # Models fit
  models <- lapply(formulas_tv, function(f) {
    fit <- survival::coxph(f, data = tv_all)
    list(model = fit, summary = summary(fit))
  })
  names(models) <- model_names
  
  return(models)
}

res_tv <- run_tv_models(df_qol)
res_tv_mcar <- run_tv_models(df_mcar)
res_tv_mar <- run_tv_models(df_mar)
res_tv_mnar <- run_tv_models(df_mnar)


# time dependence check
run_ph_diagnostics <- function(res_in) {
  
  model_names <- c(
    "Model1: PWB + BCS",
    "Model2: PWB",
    "Model3: BCS",
    "Model4: PWB (no subtype)",
    "Model5: BCS (no subtype)",
    "Model6: FACT-B total",
    "Model7: FACT-B total (no subtype)"
  )  
  
  for (m in model_names) {
    cat("Proportional Hazards Check for:", m, "\n")
    
    fit <- res_in[[m]]$model
    ph_test <- cox.zph(fit)
    print(ph_test)
    plot(ph_test, main = paste("PH Diagnostics:", m))
  }
}

# run_ph_diagnostics(res_tv)
# run_ph_diagnostics(res_tv_mcar)
# run_ph_diagnostics(res_tv_mar)
# run_ph_diagnostics(res_tv_mnar)

# ------------------------------------------------------------------------------
# 4. Basic joint modelling
# ------------------------------------------------------------------------------

run_jm_model <- function(data_in, qol_var, rmv_subtype = "N") {
  
  # Longitudinal model: QoL over time
  df_qol_lme <- data_in %>%
    filter(!is.na(.data[[qol_var]]))
  
  ids_to_keep <- unique(df_qol_lme$cases.case_id)
  
  lme_fit <- lme(as.formula(paste0(qol_var, "~ time")),
                 data = df_qol_lme,
                 random = ~ time | cases.case_id,
                 control = lmeControl(opt = "optim"))
  
  # PFS model
  df_pfs_sync <- df_pfs %>%
    filter(cases.case_id %in% ids_to_keep)
  
  if (rmv_subtype == "Y") {
    fml_final <- Surv(pfs_time_years, pfs_event) ~
      stage_group
  } else {
    fml_final <- Surv(pfs_time_years, pfs_event) ~
      stage_group + molecular_subtype
  }
  
  cox_fit <- survival::coxph(fml_final, data = df_pfs_sync)
  
  # Joint model
  set.seed(2323)
  
  jm_fit <- jm(
    Surv_object = cox_fit,
    Mixed_object = lme_fit,
    time_var = "time",
    functional_forms = ~ value(qol_var),
    n_chains = 3,
    n_burnin = 500,
    n_iter = 3500,
    n_thin = 1
  )
  
  return(jm_fit)
  
}

jm_fit_bcs <- run_jm_model(data_in = df_qol, qol_var = "bcs_total")
jm_fit_pwb <- run_jm_model(data_in = df_qol, qol_var = "pwb_total")
jm_fit_fact_b <- run_jm_model(data_in = df_qol, qol_var = "fact_b_total")

jm_fit_bcs_no_sub <- run_jm_model(data_in = df_qol, qol_var = "bcs_total",
                                  rmv_sub = "Y")
jm_fit_bcs_mcar <- run_jm_model(data_in = df_mcar, qol_var = "bcs_total")
jm_fit_bcs_mar <- run_jm_model(data_in = df_mar, qol_var = "bcs_total")
jm_fit_bcs_mnar <- run_jm_model(data_in = df_mnar, qol_var = "bcs_total")

# Diagnostics
diagnose_jm <- function(jm_fit, thin_by = 20) {
  
  thin_data <- function(block) {
    keep_idx <- seq(1, nrow(block[[1]]), by = thin_by)
    lapply(block, function(x) x[keep_idx, , drop = FALSE])
  }
  
  # association parameter
  print(mcmc_trace(thin_data(jm_fit$mcmc$alphas)))
  # survival covariates
  print(mcmc_trace(thin_data(jm_fit$mcmc$gammas)))
  # longitudinal fixed effects
  print(mcmc_trace(thin_data(jm_fit$mcmc$betas1)))
  # residual SD of the longitudinal model
  print(mcmc_trace(thin_data(jm_fit$mcmc$sigmas)))
  # random-effects covariance params
  print(mcmc_trace(thin_data(jm_fit$mcmc$D)))
  # baseline hazard spline coefficients
  print(mcmc_trace(thin_data(jm_fit$mcmc$bs_gammas)))
  
  print(summary(jm_fit))
}

diagnose_jm(jm_fit_bcs)
diagnose_jm(jm_fit_pwb)
diagnose_jm(jm_fit_fact_b)
