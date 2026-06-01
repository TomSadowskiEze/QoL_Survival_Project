### Libraries ##################################################################
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

### Read in data (clinical) ####################################################
# Clinical data set
read_df_init <- read_tsv("tcga-brca-clinical/clinical.tsv")

# Columns to be kept
req_cols <- c(
  # Required cases.~ variables
  "cases.case_id", "cases.lost_to_followup",
  # Required demographic.~ variables
  "demographic.age_at_index", "demographic.days_to_death",
  "demographic.ethnicity", "demographic.race", "demographic.vital_status",
  # Required diagnoses.~ variables
  "diagnoses.ajcc_pathologic_m", "diagnoses.ajcc_pathologic_n",
  "diagnoses.ajcc_pathologic_stage", "diagnoses.ajcc_pathologic_t",
  "diagnoses.days_to_last_follow_up", "diagnoses.metastasis_at_diagnosis",
  "diagnoses.primary_diagnosis", "diagnoses.prior_malignancy",
  "diagnoses.prior_treatment", "diagnoses.synchronous_malignancy",
  "diagnoses.year_of_diagnosis",
  # Required treatments.~ variables 
  "treatments.margin_status", "treatments.therapeutic_agents",
  "treatments.treatment_type",
  # For selection only
  "select.treatment_intent_type"
  )

# Subset (female subjects, Ductal and Lobular Neoplasms, primary tumor,
# first-line and adjuvant therapies)
df_init <- read_df_init %>%
  filter(grepl("Ductal and Lobular Neoplasms", cases.disease_type),) %>%
  filter(diagnoses.classification_of_tumor == "primary" &
           demographic.gender == "female" &
           (treatments.treatment_intent_type == "First-Line Therapy" |
           treatments.treatment_intent_type == "Adjuvant")) %>%
  rename(select.treatment_intent_type = treatments.treatment_intent_type) %>%
  dplyr::select(all_of(req_cols))

### Get cases data separately and check for inconsistencies ####################
# Only cases.~ columns
df_cases <- df_init %>%
  dplyr::select(starts_with("cases"))

# Subset inconsistent cases
inconsistent_cases_cases <- df_cases %>%
  group_by(cases.case_id) %>%
  summarise(across(everything(), ~ n_distinct(.x)), .groups = "drop") %>%
  filter(if_any(-cases.case_id, ~ .x > 1))

# Print outcome
if (nrow(inconsistent_cases_cases) == 0) {
  print("Cases - no inconsistencies")
} else {
  print("Cases - inconsistencies detected")
}

### Get demographics data separately and check for inconsistencies #############
# Only demographic.~ columns
df_demo <- df_init %>%
  dplyr::select(cases.case_id, starts_with("demographic"))

# Subset inconsistent cases
inconsistent_cases_demo <- df_demo %>%
  group_by(cases.case_id) %>%
  summarise(across(everything(), ~ n_distinct(.x)), .groups = "drop") %>%
  filter(if_any(-cases.case_id, ~ .x > 1))

# Print outcome
if (nrow(inconsistent_cases_demo) == 0) {
  print("Demographics - no inconsistencies")
} else {
  print("Demographics - inconsistencies detected")
}

### Get diagnoses data separately and check for inconsistencies ################
# Only diagnoses.~ columns
df_diag <- df_init %>%
  dplyr::select(cases.case_id, starts_with("diagnoses"))

# Subset inconsistent cases
inconsistent_cases_diag <- df_diag %>%
  group_by(cases.case_id) %>%
  summarise(across(everything(), ~ n_distinct(.x)), .groups = "drop") %>%
  filter(if_any(-cases.case_id, ~ .x > 1))

# Print outcome
if (nrow(inconsistent_cases_diag) == 0) {
  print("Diagnoses - no inconsistencies")
} else {
  print("Diagnoses - inconsistencies detected")
}

### Data checks ################################################################
# Unique cases (1050 but only 1040 with first-line)
length(unique(df_init$cases.case_id))

# Count by treatment intent (first-line: 1040, adjuvant: 2262)
table(df_init$select.treatment_intent_type)

# Count overall
count(df_init, treatments.treatment_type) %>%
  arrange(desc(n)) -> trt_counts

# Count unique
df_init %>%
  group_by(treatments.treatment_type) %>%
  summarise(unique_subjects = n_distinct(cases.case_id)) %>%
  arrange(desc(unique_subjects)) -> trt_counts_unique

### Get first-line therapy data separately #####################################
# Get first-line only along with all cases.~, demographic.~, and diagnoses.~
df_fl <- df_init %>%
  filter(grepl("First-Line Therapy", select.treatment_intent_type)) %>%
  dplyr::select(starts_with("cases"), starts_with("demographic"),
         starts_with("diagnoses"), "treatments.margin_status") %>%
  rename_with(~ sub("^treatments", "first_line", .x),
              .cols = starts_with("treatments"))

# Get unique case_id to subset other treatments
unique_ids <- unique(df_fl$cases.case_id)

### Subset initial data to include only subjects with first line ###############
df_init_with_fl <- df_init %>%
  filter(cases.case_id %in% unique_ids)

### Edit the TNM notation system ###############################################
df_fl <- df_fl %>%
  mutate(
    t_suffix =
      str_sub(diagnoses.ajcc_pathologic_t, 3),
    diagnoses.ajcc_pathologic_t =
      str_sub(diagnoses.ajcc_pathologic_t, 1, 2),
    diagnoses.ajcc_pathologic_t_supp =
      if_else(is.na(t_suffix) | t_suffix == "", "", t_suffix),
    diagnoses.ajcc_pathologic_t_supp =
      str_trim(diagnoses.ajcc_pathologic_t_supp)
  ) %>%
  dplyr::select(-t_suffix) %>%
  mutate(
    n_suffix =
      str_sub(diagnoses.ajcc_pathologic_n, 3),
    diagnoses.ajcc_pathologic_n =
      str_sub(diagnoses.ajcc_pathologic_n, 1, 2),
    diagnoses.ajcc_pathologic_n_supp =
      if_else(is.na(n_suffix) | n_suffix == "", "", n_suffix),
    diagnoses.ajcc_pathologic_n_supp =
      str_trim(diagnoses.ajcc_pathologic_n_supp)
  ) %>%
  dplyr::select(-n_suffix) %>%
  mutate(
    diagnoses.ajcc_pathologic_m_supp =
      if_else(diagnoses.ajcc_pathologic_m == "cM0 (i+)",
              "c~(i+)", ""),
    diagnoses.ajcc_pathologic_m  =
      if_else(diagnoses.ajcc_pathologic_m == "cM0 (i+)",
              "M0", diagnoses.ajcc_pathologic_m)
  )

### Get all Radiation Therapy data separately and merge with main data #########
# List all considered radiation therapy
list_rad <- c("Radiation Therapy, NOS", "Radiation, External Beam",
              "Radiation, Stereotactic/Gamma Knife/SRS",
              "Brachytherapy, High Dose", "Brachytherapy, NOS",
              "Radiation, Radioisotope")

# Subset and create flags for radiation therapy
df_rad <- df_init_with_fl %>%
  filter(treatments.treatment_type %in% list_rad) %>%
  dplyr::select("cases.case_id") %>%
  distinct(cases.case_id, .keep_all = TRUE) %>%
  mutate(rad.fl = "Y")

# Merge with main data
df_work1 <- df_fl %>%
  left_join(df_rad, by ="cases.case_id") %>%
  mutate(rad.fl = replace_na(rad.fl, "N"))

### Get Hormone Therapy data separately ########################################
# Subset all cases with hormone therapy
df_horm <- df_init_with_fl %>%
  filter(grepl("Adjuvant", select.treatment_intent_type)) %>%
  filter(grepl("Hormone Therapy", treatments.treatment_type)) %>%
  dplyr::select("cases.case_id", starts_with("treatments")) %>%
  dplyr::select(-treatments.margin_status) %>%
  rename_with(~ sub("^treatments", "horm", .x),
              .cols = starts_with("treatments"))

### Get Chemotherapy data separately ###########################################
# Subset all cases with chemotherapy
pre_df_chemo <- df_init_with_fl %>%
  filter(grepl("Adjuvant", select.treatment_intent_type)) %>%
  filter(grepl("Chemotherapy", treatments.treatment_type)) %>%
  dplyr::select("cases.case_id", starts_with("treatments")) %>%
  dplyr::select(-treatments.margin_status) %>% 
  rename_with(~ sub("^treatments", "chemo", .x),
              .cols = starts_with("treatments"))

# Tamoxifen, Exemestaneand, and Anastrozole reclassified as hormone therapy
list_chg <- c("Tamoxifen", "Exemestane", "Anastrozole")

df_horm2 <- df_init_with_fl %>%
  filter(grepl("Chemotherapy", treatments.treatment_type)) %>%
  filter(treatments.therapeutic_agents %in% list_chg) %>%
  dplyr::select("cases.case_id", starts_with("treatments")) %>%
  dplyr::select(-treatments.margin_status) %>%
  rename_with(~ sub("^treatments", "horm", .x),
              .cols = starts_with("treatments"))

# Same cases removed from the chemotherapy data
df_chemo <- pre_df_chemo %>%
  filter(!chemo.therapeutic_agents %in% list_chg)

### Get Pharmaceutical Therapy, NOS data separately ############################
# Subset
df_phar <- df_init_with_fl %>%
  filter(grepl("Adjuvant", select.treatment_intent_type)) %>%
  filter(grepl("Pharmaceutical Therapy, NOS", treatments.treatment_type)) %>%
  dplyr::select("cases.case_id", starts_with("treatments")) %>%
  rename_with(~ sub("^treatments", "phar", .x),
              .cols = starts_with("treatments"))

# Triptorelin and Anastrozole case to be reclassified as hormone therapy
list_chg <- c("Triptorelin", "Anastrozole")

df_horm3 <- df_init_with_fl %>%
  filter(grepl("Pharmaceutical Therapy, NOS", treatments.treatment_type)) %>%
  filter(treatments.therapeutic_agents %in% list_chg) %>%
  dplyr::select("cases.case_id", starts_with("treatments")) %>%
  dplyr::select(-treatments.margin_status) %>%
  rename_with(~ sub("^treatments", "horm", .x),
              .cols = starts_with("treatments"))

### Get Targeted Molecular Therapy data separately #############################

df_targ <- df_init_with_fl %>%
  filter(grepl("Adjuvant", select.treatment_intent_type)) %>%
  filter(grepl("Targeted Molecular Therapy", treatments.treatment_type)) %>%
  dplyr::select("cases.case_id", starts_with("treatments")) %>%
  rename_with(~ sub("^treatments", "targ", .x),
              .cols = starts_with("treatments"))

# Bevacizumab, Trastuzumab, and Lapatinib cases grouped with chemotherapy
# (Effectively all observations)

df_chemo2 <- df_init_with_fl %>%
  filter(grepl("Targeted Molecular Therapy", treatments.treatment_type)) %>%
  dplyr::select("cases.case_id", starts_with("treatments")) %>%
  dplyr::select(-treatments.margin_status) %>%
  rename_with(~ sub("^treatments", "chemo", .x),
              .cols = starts_with("treatments"))

### Get Immunotherapy (Including Vaccines) data separately #####################

df_immuno_data_subset <- df_init_with_fl %>%
  filter(grepl("Adjuvant", select.treatment_intent_type)) %>%
  filter(grepl("Immunotherapy (Including Vaccines)",
               treatments.treatment_type, fixed=TRUE)) %>%
  dplyr::select("cases.case_id", starts_with("treatments")) %>%
  rename_with(~ sub("^treatments", "immuno_data_subset", .x),
              .cols = starts_with("treatments"))

# Trastuzumab cases grouped with chemotherapy (effectively all observations)
df_chemo3 <- df_init_with_fl %>%
  filter(grepl("Immunotherapy (Including Vaccines)",
               treatments.treatment_type, fixed=TRUE)) %>%
  dplyr::select("cases.case_id", starts_with("treatments")) %>%
  dplyr::select(-treatments.margin_status) %>%
  rename_with(~ sub("^treatments", "chemo", .x),
              .cols = starts_with("treatments"))

### Merge all hormone therapy data and classify by mechanism ###################

bind_rows(df_horm, df_horm2, df_horm3) %>%
  dplyr::select(-horm.treatment_type) -> df_all_horm

# List - Selective estrogen receptor modulators (SERMs) (group serm)
serm_lst <- c("Tamoxifen", "Tamoxifen Citrate", "Toremifene Citrate")

# List - LHRH agonists (group lhrh_a)
lhrh_a_lst <- c("Leuprolide", "Leuprolide Acetate", "Triptorelin",
                "Goserelin", "Goserelin Acetate")

# List - Aromatase inhibitors (AIs) (group ai)
ai_lst <- c("Anastrozole", "Letrozole", "Exemestane")

# List - Selective estrogen receptor degraders (SERDs) (group serd)
serd_lst <- c("Fulvestrant")

# Complete list of all
all_lst <- c(serm_lst, lhrh_a_lst, ai_lst, serd_lst)

lst1 <- c("df_horm_serm", "df_horm_lhrh_a", "df_horm_ai", "df_horm_serd",
          "df_horm_any")
lst2 <- list(serm_lst, lhrh_a_lst, ai_lst, serd_lst, all_lst)
lst3 <- c("serm.fl", "lhrh_a.fl", "ai.fl", "serd.fl", "any_horm.fl")

for (i in seq_along(lst1)) {
  temp_result <- df_all_horm %>%
    filter(horm.therapeutic_agents %in% lst2[[i]]) %>%
    dplyr::select("cases.case_id") %>%
    distinct(cases.case_id, .keep_all = TRUE) %>%
    mutate(!!lst3[i] := "Y")
  
  assign(lst1[i], temp_result)
}

# Merge with main data
df_work2 <- df_work1

for (i in seq_along(lst1)) {
  
  df_to_join <- get(lst1[i])
  col_name <- lst3[i]
  
  df_work2 <- df_work2 %>%
    left_join(df_to_join, by ="cases.case_id") %>%
    mutate(!!sym(col_name) := replace_na(.data[[col_name]], "N"))
}

### Merge all chemotherapy data and classify by type ###########################

bind_rows(df_chemo, df_chemo2, df_chemo3) %>%
  dplyr::select(-chemo.treatment_type) %>%
  mutate(chemo.therapeutic_agents = if_else(
    chemo.therapeutic_agents == "'--",
    "Unknown", chemo.therapeutic_agents)
  ) -> df_all_chemo

# List - Trastuzumab and Bevacizumab separately due to different effect on QoL
immuno_lst <- c("Trastuzumab", "Bevacizumab")

immuno_temp_df <- df_all_chemo %>%
  filter(chemo.therapeutic_agents %in% immuno_lst) %>%
  dplyr::select("cases.case_id") %>%
  distinct(cases.case_id, .keep_all = TRUE) %>%
  mutate(immuno.fl = "Y")

chemo_temp_df <- df_all_chemo %>%
  filter(!chemo.therapeutic_agents %in% immuno_lst) %>%
  dplyr::select("cases.case_id") %>%
  distinct(cases.case_id, .keep_all = TRUE) %>%
  mutate(chemo.fl = "Y")

df_work3 <- df_work2 %>%
  left_join(immuno_temp_df, by ="cases.case_id") %>%
  mutate(immuno.fl = replace_na(immuno.fl, "N")) %>%
  left_join(chemo_temp_df, by ="cases.case_id") %>%
  mutate(chemo.fl = replace_na(chemo.fl, "N"))

### Read in data (clinical) ####################################################
# Clinical data set
read_df_init_fu <- read_tsv("tcga-brca-clinical/follow_up.tsv")

# Subset initial data to include only subjects with first line
df_init_fu <- read_df_init_fu %>%
  filter(cases.case_id %in% unique_ids)

### Follow-ups data ############################################################
# Columns to be kept - follow ups
req_cols_fu <- c(
  # Required cases.~ variables
  "cases.case_id",
  # Required follow_ups.~ variables
  "follow_ups.days_to_follow_up", "follow_ups.days_to_progression",
  "follow_ups.days_to_recurrence", "follow_ups.disease_response",
  "follow_ups.timepoint_category"
)

# Subset follow up data
df_init_fu_data <- df_init_fu %>%
  dplyr::select(all_of(req_cols_fu)) %>%
  filter(!if_all(req_cols_fu[-1], ~ .x == "'--"))

# Subset progression subjects data, keep first progression date
df_prog <- df_init_fu_data %>%
  filter(follow_ups.timepoint_category == "Post Initial Treatment") %>%
  filter(follow_ups.days_to_progression != "'--") %>%
  arrange(cases.case_id, follow_ups.days_to_progression) %>%
  distinct(cases.case_id, .keep_all = TRUE) %>%
  mutate(prog.fl = "Y") %>%
  mutate(prog.d = follow_ups.days_to_progression) %>%
  dplyr::select(cases.case_id, prog.fl, prog.d)

# Subset recurrence subjects data, keep first recurrence date
df_rec <- df_init_fu_data %>%
  filter(follow_ups.timepoint_category == "Post Initial Treatment") %>%
  filter(follow_ups.days_to_recurrence != "'--") %>%
  arrange(cases.case_id, follow_ups.days_to_recurrence) %>%
  distinct(cases.case_id, .keep_all = TRUE) %>%
  mutate(rec.fl = "Y") %>%
  mutate(rec.d = follow_ups.days_to_recurrence) %>%
  dplyr::select(cases.case_id, rec.fl, rec.d)

# Join all follow-up data and add pfs_event
df_fu_final <- df_prog %>%
  full_join(df_rec, by = "cases.case_id") %>%
  mutate(pfs_event.fl = "Y") %>%
  mutate(pfs_event.d = pmin(prog.d, rec.d, na.rm = TRUE))

### Molecular data #############################################################
# Columns to be kept - molecular
req_cols_mol <- c(
  # Required cases.~ variables
  "cases.case_id",
  # Required molecular_tests.~ variables
  "molecular_tests.gene_symbol", "molecular_tests.test_result"
)

# Results that are informative
req_test_results <- c("Positive", "Negative")

# Subset molecular data
df_init_mol <- read_df_init_fu %>%
  dplyr::select(all_of(req_cols_mol)) %>%
  filter(!if_all(req_cols_mol[-1], ~ .x == "'--")) %>%
  filter(molecular_tests.gene_symbol != "Not Applicable") %>%
  filter(molecular_tests.test_result %in% req_test_results)

df_mol <- df_init_mol %>%
  group_by(cases.case_id, molecular_tests.gene_symbol) %>%
  summarise(
    status = case_when(
      any(molecular_tests.test_result == "Positive") ~ "Positive",
      all(molecular_tests.test_result == "Negative") ~ "Negative",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = molecular_tests.gene_symbol,
    values_from = status,
    names_glue = "{tolower(molecular_tests.gene_symbol)}.status"
  ) %>%
  mutate(
    molecular_subtype = case_when(
      esr1.status == "Positive" & pgr.status == "Positive" &
        erbb2.status == "Negative" ~ "Luminal A-like",
      esr1.status == "Positive" & pgr.status == "Negative" &
        erbb2.status == "Negative" ~ "Luminal B-like (HER2-neg)",
      esr1.status == "Positive" & erbb2.status == "Positive" ~
        "Luminal B-like (HER2-pos)",
      esr1.status == "Negative" & pgr.status == "Negative" &
        erbb2.status == "Positive" ~ "HER2-enriched",
      esr1.status == "Negative" & pgr.status == "Negative" &
        erbb2.status == "Negative" ~ "Triple-negative",
      TRUE ~ NA_character_
    )
  )

### Menopausal data ############################################################
# Columns to be kept - attr
req_cols_meno <- c(
  # Required cases.~ variables
  "cases.case_id",
  # Required other_clinical_attributes.~ variables 
  "other_clinical_attributes.menopause_status"
)

# Subset menopausal data
df_meno <- read_df_init_fu %>%
  dplyr::select(all_of(req_cols_meno)) %>%
  filter(!if_all(req_cols_meno[-1], ~ .x == "'--")) %>%
  rename(meno.status=other_clinical_attributes.menopause_status)

### Merge all follow-up data with the original data set ########################

df_work4 <- df_work3 %>%
  left_join(df_fu_final, by ="cases.case_id") %>%
  mutate(prog.fl = replace_na(prog.fl, "N")) %>%
  mutate(rec.fl = replace_na(rec.fl, "N")) %>%
  mutate(pfs_event.fl = replace_na(pfs_event.fl, "N")) %>%
  left_join(df_mol, by ="cases.case_id") %>%
  left_join(df_meno, by = "cases.case_id")

### Additional data clean up ###################################################

# Uniform flags and categories
df_work5 <- df_work4 %>%
  mutate(demographic.age_at_index = as.integer(demographic.age_at_index)) %>%
  mutate(cases.lost_to_followup = if_else(
    demographic.vital_status == "Dead",
    "No", cases.lost_to_followup)) %>%
  mutate(cases.lost_to_followup = if_else(
    cases.lost_to_followup == "'--",
    NA, cases.lost_to_followup)) %>%  
  mutate(demographic.days_to_death = if_else(
    demographic.vital_status == "Dead" & demographic.days_to_death == "'--",
    diagnoses.days_to_last_follow_up, demographic.days_to_death)) %>%
  mutate(demographic.ethnicity = if_else(
    demographic.ethnicity %in% c("not reported","Unknown"),
    "Missing", demographic.ethnicity)) %>%
  mutate(demographic.race = if_else(
    demographic.race %in% c("not reported"),
    "Missing", demographic.race)) %>%
  mutate(diagnoses.ajcc_pathologic_stage = if_else(
    diagnoses.ajcc_pathologic_stage %in% c("'--"),
    NA, diagnoses.ajcc_pathologic_stage)) %>%
  mutate(diagnoses.metastasis_at_diagnosis = case_when(
    diagnoses.metastasis_at_diagnosis == "Metastasis, NOS" ~ "Y",
    diagnoses.metastasis_at_diagnosis == "No Metastasis" ~ "N",
    TRUE ~ NA_character_)) %>%
  mutate(diagnoses.prior_malignancy = case_when(
    diagnoses.prior_malignancy == "yes" ~ "Y",
    diagnoses.prior_malignancy == "no" ~ "N",
    TRUE ~ NA_character_)) %>%
  mutate(diagnoses.prior_treatment = case_when(
    diagnoses.prior_treatment == "Yes" ~ "Y",
    diagnoses.prior_treatment == "No" ~ "N",
    TRUE ~ NA_character_)) %>%
  mutate(diagnoses.synchronous_malignancy = case_when(
    diagnoses.synchronous_malignancy == "Yes" ~ "Y",
    diagnoses.synchronous_malignancy == "No" ~ "N",
    TRUE ~ NA_character_)) %>%
  mutate(diagnoses.year_of_diagnosis = if_else(
    diagnoses.year_of_diagnosis %in% c("'--"),
    NA, diagnoses.year_of_diagnosis)) %>%
  mutate(first_line.margin_status = if_else(
    first_line.margin_status %in% c("Not Reported", "Unknown"),
    "Missing", first_line.margin_status)) %>%
  mutate(erbb2.status = if_else(
    is.na(erbb2.status), "Unknown", erbb2.status)) %>%
  mutate(esr1.status = if_else(
    is.na(esr1.status), "Unknown", esr1.status)) %>%
  mutate(pgr.status = if_else(
    is.na(pgr.status), "Unknown", pgr.status)) %>%
  mutate(molecular_subtype = if_else(
    is.na(molecular_subtype), "Unknown", molecular_subtype)) %>%
  mutate(meno.status = if_else(
    is.na(meno.status), "Unknown", meno.status)) %>%
  mutate(demographic.days_to_death =
      na_if(as.character(demographic.days_to_death), "'--") %>%
        as.integer(),
    diagnoses.days_to_last_follow_up =
      na_if(as.character(diagnoses.days_to_last_follow_up), "'--") %>%
        as.integer()) %>%
  # Propcase for ethnicity and race
  mutate(demographic.ethnicity = case_when(
    demographic.ethnicity == "hispanic or latino" ~ "Hispanic or Latino",
    demographic.ethnicity == "not hispanic or latino" ~
      "Not Hispanic or Latino",
    TRUE ~ "Missing"
  )) %>%
  mutate(demographic.race = case_when(
    demographic.race == "american indian or alaska native" ~
      "American Indian or Alaska Native",
    demographic.race == "asian" ~ "Asian",
    demographic.race == "black or african american" ~
      "Black or African American",
    demographic.race == "white" ~ "White",
    TRUE ~ "Missing"
  ))
  
# Censoring at 10 years
df_work6 <- df_work5 %>%
  # Days to last follow up - if greater replace with 3653
  mutate(diagnoses.days_to_last_follow_up = if_else(
    diagnoses.days_to_last_follow_up > 3653, 3653,
    diagnoses.days_to_last_follow_up)) %>%
  # Days to death - if greater replace with NA
  mutate(demographic.days_to_death = if_else(
    demographic.days_to_death > 3653, NA,
    demographic.days_to_death)) %>%
  # If days to death removed - amend vital status
  mutate(demographic.vital_status = if_else(
    is.na(demographic.days_to_death), "Alive", "Dead")) %>%
  # Adjust the pfs_event flag and day using the days to death
  mutate(pfs_event.fl = if_else(
    demographic.vital_status == "Dead", "Y", pfs_event.fl)) %>%
  mutate(pfs_event.d = if_else(
    pfs_event.fl == "Y", pmin(prog.d, rec.d, demographic.days_to_death,
                              na.rm = TRUE), NA))

### Output final data set ready for the simulation and analysis ################
df_final_1 <- df_work6

### Summary statistics #########################################################
df_final_1_sum <- df_final_1 %>%
  mutate(age_group = case_when(
    demographic.age_at_index < 50 ~ "<50",
    demographic.age_at_index < 65 ~ "50-64",
    TRUE ~ "65+"),
    cases.lost_to_followup = factor(
      cases.lost_to_followup,
      levels = c("Yes", "No", NA)),
    demographic.ethnicity = factor(
      demographic.ethnicity,
      levels = c("Not Hispanic or Latino", "Hispanic or Latino", "Missing")),
    demographic.race = factor(
      demographic.race,
      levels = c("White", "Black or African American", "Asian",
                 "American Indian or Alaska Native", "Missing")),
    meno.status = factor(
      meno.status,
      levels = c("Premenopausal", "Perimenopausal",
                 "Postmenopausal", "Unknown")),
    stage_raw = as.character(diagnoses.ajcc_pathologic_stage),
    stage_clean = trimws(gsub("(?i)^Stage\\s+", "", stage_raw, perl = TRUE)),
    stage_group = case_when(
      is.na(stage_clean) ~ "Unknown",
      grepl("^IV", stage_clean) ~ "Stage IV",
      grepl("^III", stage_clean) ~ "Stage III",
      grepl("^II", stage_clean) ~ "Stage II",
      grepl("^I", stage_clean) ~ "Stage I",
      TRUE ~ "Unknown"
    ),
    stage_group = factor(
      stage_group,
      levels = c("Stage I", "Stage II", "Stage III", "Stage IV", "Unknown")),
    diagnoses.metastasis_at_diagnosis = case_when(
      is.na(diagnoses.metastasis_at_diagnosis) ~ "Unknown",
      diagnoses.metastasis_at_diagnosis == "N" ~ "No",
      diagnoses.metastasis_at_diagnosis == "Y" ~ "Yes"
    ),
    diagnoses.metastasis_at_diagnosis = factor(
      diagnoses.metastasis_at_diagnosis,
      levels = c("No", "Yes", "Unknown")),
    diagnoses.primary_diagnosis = case_when(
      # IDC (Invasive Ductal Carcinoma)
      diagnoses.primary_diagnosis %in% c(
        "Infiltrating duct carcinoma, NOS",
        "Intraductal papillary adenocarcinoma with invasion",
        "Paget disease and infiltrating duct carcinoma of breast"
      ) ~ "IDC (Invasive Ductal Carcinoma)",
      # ILC (Invasive Lobular Carcinoma)
      diagnoses.primary_diagnosis %in% c(
        "Lobular carcinoma, NOS"
      ) ~ "ILC (Invasive Lobular Carcinoma)",
      # Mixed ductal/lobular
      diagnoses.primary_diagnosis %in% c(
        "Infiltrating duct and lobular carcinoma",
        "Infiltrating duct mixed with other types of carcinoma",
        "Infiltrating lobular mixed with other types of carcinoma"
      ) ~ "Mixed ductal/lobular",
      # Other
      TRUE ~ "Other"
    ),
    diagnoses.primary_diagnosis = factor(
      diagnoses.primary_diagnosis,
      levels = c(
        "IDC (Invasive Ductal Carcinoma)",
        "ILC (Invasive Lobular Carcinoma)",
        "Mixed ductal/lobular",
        "Other"
      )),
    diagnoses.prior_malignancy = case_when(
      is.na(diagnoses.prior_malignancy) ~ "Unknown",
      diagnoses.prior_malignancy == "N" ~ "No",
      diagnoses.prior_malignancy == "Y" ~ "Yes"
    ),
    diagnoses.prior_malignancy = factor(
      diagnoses.prior_malignancy,
      levels = c("No", "Yes", "Unknown")),
    diagnoses.synchronous_malignancy = case_when(
      is.na(diagnoses.synchronous_malignancy) ~ "Unknown",
      diagnoses.synchronous_malignancy == "N" ~ "No",
      diagnoses.synchronous_malignancy == "Y" ~ "Yes"
    ),
    diagnoses.synchronous_malignancy = factor(
      diagnoses.synchronous_malignancy,
      levels = c("No", "Yes", "Unknown")),
    molecular_subtype = factor(
      molecular_subtype,
      levels = c("Luminal A-like", "Luminal B-like (HER2-neg)",
                 "Luminal B-like (HER2-pos)", "HER2-enriched",
                 "Triple-negative", "Unknown")),
    first_line.margin_status = case_when(
      first_line.margin_status == "Uninvolved" ~ "Uninvolved",
      first_line.margin_status == "Involved" ~ "Involved",
      TRUE ~ "Unknown"
    ),
    first_line.margin_status = factor(
      first_line.margin_status,
      levels = c("Uninvolved", "Involved", "Unknown")),
    rad.fl = case_when(
      rad.fl == "Y" ~ "Yes",
      rad.fl == "N" ~ "No"
    ),
    rad.fl = factor(
      rad.fl,
      levels = c("Yes", "No")),
    chemo.fl = case_when(
      chemo.fl == "Y" ~ "Yes",
      chemo.fl == "N" ~ "No"
    ),
    chemo.fl = factor(
      chemo.fl,
      levels = c("Yes", "No")),
    any_horm.fl = case_when(
      any_horm.fl == "Y" ~ "Yes",
      any_horm.fl == "N" ~ "No"
    ),
    any_horm.fl = factor(
      any_horm.fl,
      levels = c("Yes", "No")),
    immuno.fl = case_when(
      immuno.fl == "Y" ~ "Yes",
      immuno.fl == "N" ~ "No"
    ),
    immuno.fl = factor(
      immuno.fl,
      levels = c("Yes", "No")),
    diagnoses.prior_treatment = case_when(
      is.na(diagnoses.prior_treatment) ~ "Unknown",
      diagnoses.prior_treatment == "Y" ~ "Yes",
      diagnoses.prior_treatment == "N" ~ "No",
      TRUE ~ "Unknown"
    ),
    diagnoses.prior_treatment = factor(
      rad.fl,
      levels = c("Yes", "No"))
    )

tabulate_var <- function(df, var) {
  
  # Overall
  x <- df[[var]]
  tb_overall <- table(x, useNA = "ifany")
  
  overall_total <- nrow(df)
  pct_overall <- (as.vector(tb_overall) / overall_total) * 100
  
  overall_df <- data.frame(
    group = "Overall",
    level = names(tb_overall),
    n = as.vector(tb_overall),
    pct = pct_overall
  )
  
  # By age group
  g <- df[["age_group"]]
  tb_list <- tapply(x, g, table, useNA = "ifany")
  
  by_df <- do.call(rbind, lapply(names(tb_list), function(grp) {
    tb <- tb_list[[grp]]
    
    group_total <- sum(df$age_group == grp, na.rm = TRUE)
    pct <- (as.vector(tb) / group_total) * 100
    
    data.frame(
      group = grp,
      level = names(tb),
      n = as.vector(tb),
      pct = pct
    )
  }))
  
  rownames(by_df) <- NULL
  
  # Merge long format
  long_df <- rbind(overall_df, by_df)
  
  # Put values in "n (%)" format
  long_df$value <- ifelse(
    long_df$n == 0 & long_df$pct == 0,
    "0",
    paste0(long_df$n, " (", sprintf("%.2f", long_df$pct), ")")
  )
  
  # Pivot wider
  wide_df <- long_df %>%
    dplyr::select(group, level, value) %>%
    tidyr::pivot_wider(
      names_from  = group,
      values_from = value
    )
  
  # Add the group count (keep non‑syntactic names!)
  n_row <- data.frame(
    level   = "N",
    `<50`   = as.character(sum(df$age_group == "<50",   na.rm = TRUE)),
    `50-64` = as.character(sum(df$age_group == "50-64", na.rm = TRUE)),
    `65+`   = as.character(sum(df$age_group == "65+",   na.rm = TRUE)),
    Overall = as.character(nrow(df)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  info_row <- data.frame(
    level = "info",
    `<50`   = "n (%)",
    `50-64` = "n (%)",
    `65+`   = "n (%)",
    Overall = "n (%)",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Bind n-row on top
  final_df <- dplyr::bind_rows(n_row, info_row, wide_df)
  
  # Enforce final column order
  desired_cols <- c("level", "<50", "50-64", "65+", "Overall")
  final_df <- final_df %>%
    dplyr::select(dplyr::any_of(desired_cols))

  final_df

}

### Table 1 - Baseline demographics ############################################
eth_sum <- tabulate_var(df = df_final_1_sum,
                        var = "demographic.ethnicity")
race_sum <- tabulate_var(df = df_final_1_sum,
                        var = "demographic.race")
meno_sum <- tabulate_var(df = df_final_1_sum,
                         var = "meno.status")
eth_sum
race_sum
meno_sum

summarise_age <- function(df) {
  
  age_stats <- function(x) {
    x <- as.numeric(x)
    data.frame(
      n = length(x),
      mean = round(mean(x, na.rm = TRUE), 1),
      sd = round(sd(x, na.rm = TRUE), 1),
      median = round(median(x, na.rm = TRUE), 1),
      IQR = paste0(
        round(quantile(x, 0.25, na.rm = TRUE), 1), "-",
        round(quantile(x, 0.75, na.rm = TRUE), 1)
      ),
      range = paste0(
        round(min(x, na.rm = TRUE), 1), "-",
        round(max(x, na.rm = TRUE), 1)
      ),
      check.names = FALSE
    )
  }
  
  groups <- unique(df$age_group)
  
  out_list <- lapply(groups, function(g) {
    d_sub <- df[df$age_group == g, ]
    stats <- age_stats(d_sub$demographic.age_at_index)
    stats$AgeGroup <- g
    stats
  })
  
  out <- do.call(rbind, out_list)
  overall <- age_stats(df$demographic.age_at_index)
  overall$AgeGroup <- "Overall"
  out <- rbind(out, overall)
  
  out <- out[, c("AgeGroup", "n", "mean", "sd", "median", "IQR", "range")]
  rownames(out) <- NULL
  
  return(out)
}

age_summary <- summarise_age(df_final_1_sum)
age_summary

### Table 2 - Disease characteristics ##########################################
stage_sum <- tabulate_var(df = df_final_1_sum,
                          var = "stage_group")
meta_sum <- tabulate_var(df = df_final_1_sum,
                         var = "diagnoses.metastasis_at_diagnosis")
primary_sum <- tabulate_var(df = df_final_1_sum,
                            var = "diagnoses.primary_diagnosis")
prior_malig_sum <- tabulate_var(df = df_final_1_sum,
                                var = "diagnoses.prior_malignancy")
synch_malig_sum <- tabulate_var(df = df_final_1_sum,
                                var = "diagnoses.synchronous_malignancy")
HER_sum <- tabulate_var(df = df_final_1_sum,
                                var = "erbb2.status")
ER_sum <- tabulate_var(df = df_final_1_sum,
                                var = "esr1.status")
PR_sum <- tabulate_var(df = df_final_1_sum,
                                var = "pgr.status")
molecular_sum <- tabulate_var(df = df_final_1_sum,
                              var = "molecular_subtype")

primary_sum
stage_sum
ER_sum
PR_sum
HER_sum
molecular_sum
meta_sum
prior_malig_sum
synch_malig_sum

### Table 3 - Treatment summary ################################################
margin_sum <- tabulate_var(df = df_final_1_sum,
                           var = "first_line.margin_status")
rad_sum <- tabulate_var(df = df_final_1_sum,
                        var = "rad.fl")
chemo_sum <- tabulate_var(df = df_final_1_sum,
                          var = "chemo.fl")
horm_sum <- tabulate_var(df = df_final_1_sum,
                         var = "any_horm.fl")
immuno_sum <- tabulate_var(df = df_final_1_sum,
                           var = "immuno.fl")
prior_trt_sum <- tabulate_var(df = df_final_1_sum,
                              var = "diagnoses.prior_treatment")

prior_trt_sum
margin_sum
rad_sum
chemo_sum
horm_sum
immuno_sum


### Table 4 - Summarise outcomes ###############################################
summarise_outcomes <- function(df, mode = "O") {
  
  summarise_event <- function(flag, time) {
    time <- as.numeric(time)
    
    n_event <- sum(flag == "Y", na.rm = TRUE)
    pct_event <- mean(flag == "Y", na.rm = TRUE) * 100
    
    if (n_event > 0) {
      times <- time[flag == "Y"]
      med <- median(times, na.rm = TRUE)
      q <- quantile(times, probs = c(0.25, 0.75), na.rm = TRUE)
      time_str <- sprintf("%.0f (%.0f-%.0f)", med, q[1], q[2])
    } else {
      time_str <- NA
    }
    
    list(n = n_event, pct = pct_event, time = time_str)
  }
  
  # Build table
  build_table <- function(d) {
    prog  <- summarise_event(d$prog.fl, d$prog.d)
    rec   <- summarise_event(d$rec.fl,  d$rec.d)
    death_flag <- ifelse(d$demographic.vital_status == "Dead", "Y", "N")
    death <- summarise_event(death_flag, d$demographic.days_to_death)
    
    data.frame(
      Outcome = c("Progression", "Recurrence", "Death"),
      `n (%)` = c(
        sprintf("%d (%.1f%%)", prog$n, prog$pct),
        sprintf("%d (%.1f%%)", rec$n, rec$pct),
        sprintf("%d (%.1f%%)", death$n, death$pct)
      ),
      `Median time (IQR)` = c(prog$time, rec$time, death$time),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  
  # MODE: Overall
  if (mode == "O") {
    return(build_table(df))
  }
  
  # MODE: By Age Group
  if (mode == "A") {
    
    age_groups <- c("<50","50-64","65+")
    
    out_list <- lapply(age_groups, function(g) {
      d_sub <- df[df$age_group == g, ]
      tab <- build_table(d_sub)
      tab$AgeGroup <- g
      tab
    })
    
    out <- do.call(rbind, out_list)
    out <- out[, c("AgeGroup", "Outcome", "n (%)", "Median time (IQR)")]
    rownames(out) <- NULL
    return(out)
  }
  
  # MODE: By Molecular Subtype
  if (mode == "MS") {
    
    mol_sub_groups <- c(
      "Luminal A-like",
      "Luminal B-like (HER2-neg)",
      "Luminal B-like (HER2-pos)",
      "HER2-enriched",
      "Triple-negative",
      "Unknown"
    )
    
    out_list <- lapply(mol_sub_groups, function(g) {
      d_sub <- df[df$molecular_subtype == g, ]
      tab <- build_table(d_sub)
      tab$MolSubGroup <- g
      tab
    })
    
    out <- do.call(rbind, out_list)
    out <- out[, c("MolSubGroup", "Outcome", "n (%)", "Median time (IQR)")]
    rownames(out) <- NULL
    return(out)
  }
  
  # MODE: By Stage
  if (mode == "S") {
    
    stage_groups <- c(
      "Stage I",
      "Stage II",
      "Stage III",
      "Stage IV",
      "Unknown"
    )
    
    
    out_list <- lapply(stage_groups, function(g) {
      d_sub <- df[df$stage_group == g, ]
      tab <- build_table(d_sub)
      tab$StageGroup <- g
      tab
    })
    
    out <- do.call(rbind, out_list)
    out <- out[, c("StageGroup", "Outcome", "n (%)", "Median time (IQR)")]
    rownames(out) <- NULL
    return(out)
  }
  
  stop("mode must be 'O' (overall) or 'A' (by age)
       or 'MS' (by molecular subtype) or 'S' (by stage)")
}

outcome_overall <- summarise_outcomes(df_final_1_sum, mode = "O")
outcome_age <- summarise_outcomes(df_final_1_sum, mode = "A")
outcome_molsub <- summarise_outcomes(df_final_1_sum, mode = "MS")
outcome_stage <- summarise_outcomes(df_final_1_sum, mode = "S")

# Visualisation
okabe_ito <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

clean_outcome_table <- function(df, group_col, group_type_label) {
  
  if (group_col == "AgeGroup") {
    group_levels <- c("<50", "50-64", "65+")
  } else {
    group_levels <- df[[group_col]] %>% unique()
  }
  
  out <- df %>%
    separate(`n (%)`, into = c("n", "pct"), sep = " ", remove = FALSE) %>%
    mutate(
      pct = as.numeric(gsub("[()%]", "", pct)),
      Outcome = factor(Outcome,
                       levels = c("Progression", "Recurrence", "Death")),
      GroupType = group_type_label,
      Group = factor(.data[[group_col]], levels = group_levels)
    ) %>%
    select(GroupType, Group, Outcome, pct)
  
  return(out)
  
}

age_clean    <- clean_outcome_table(outcome_age,
                                    "AgeGroup", "Age group")
molsub_clean <- clean_outcome_table(outcome_molsub,
                                    "MolSubGroup", "Molecular subtype")
stage_clean  <- clean_outcome_table(outcome_stage,
                                    "StageGroup",  "Stage")
combined_rates <- bind_rows(age_clean, molsub_clean, stage_clean) %>%
  mutate(
    Group = case_when(
      GroupType == "Molecular subtype" &
        Group == "Unknown"~ "Unknown (Subtype)",
      GroupType == "Stage" &
        Group == "Unknown" ~ "Unknown (Stage)",
      TRUE ~ Group
    )
  )


combined_rates <- combined_rates %>%
  group_by(GroupType) %>%
  mutate(Group = factor(Group, levels = unique(Group))) %>%
  ungroup()

ggplot(combined_rates, aes(x = Group, y = pct, fill = Outcome)) +
  geom_col(position = position_dodge(width = 0.8)) +
  facet_wrap(~ GroupType, scales = "free_x") +
  scale_fill_manual(values = okabe_ito[1:3]) +
  labs(
    x = "",
    y = "Event rate (%)",
    fill = "Outcome"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    legend.position = c(0, -0.42),
    legend.justification = c(0, 0)
  )

