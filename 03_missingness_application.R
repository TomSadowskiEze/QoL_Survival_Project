### Libraries ##################################################################
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)

# ------------------------------------------------------------------------------
# 1. Missing completely at random (MCAR)
# ------------------------------------------------------------------------------
set.seed(2323)

df_mcar <- df_qol %>%
  mutate(mcar_miss = ifelse(!structural_miss,
                            rbinom(n(), 1, 0.10) == 1, FALSE)) %>%
  mutate(
    across(c(swb_total, ewb_total, fwb_total, bcs_total, pwb_total,
             fact_b_total),
      ~ ifelse(mcar_miss, NA, .x)))

# ------------------------------------------------------------------------------
# 2. Missing at random (MAR)
# ------------------------------------------------------------------------------
df_mar <- df_qol %>%
  group_by(cases.case_id) %>%
  mutate(
    lag_factb = lag(fact_b_total),
    p_miss = ifelse(is.na(lag_factb), 0, plogis(-1.2 + -0.01 * lag_factb)),
    mar_miss = ifelse(!structural_miss, rbinom(n(), 1, p_miss) == 1, FALSE)) %>%
  ungroup() %>%
  mutate(across(
    c(swb_total, ewb_total, fwb_total, bcs_total, pwb_total, fact_b_total),
    ~ ifelse(mar_miss, NA, .x)))

# ------------------------------------------------------------------------------
# 3. Missing not at random (MNAR)
# ------------------------------------------------------------------------------
df_mnar <- df_qol %>%
  group_by(cases.case_id) %>%
  mutate(
    p_miss = ifelse(is.na(fact_b_total),
                    0, plogis(-1.2 + -0.01 * fact_b_total)),
    mnar_miss = ifelse(!structural_miss, rbinom(n(), 1, p_miss) == 1, FALSE)) %>%
  ungroup() %>%
  mutate(across(
    c(swb_total, ewb_total, fwb_total, bcs_total, pwb_total, fact_b_total),
    ~ ifelse(mnar_miss, NA, .x)))

# ------------------------------------------------------------------------------
# 4. Plot missing
# ------------------------------------------------------------------------------

plot_missingness <- function(df, miss_var, title_text) {
  
  miss_sym <- rlang::sym(miss_var)
  df_summary <- df %>%
    group_by(time) %>%
    summarise(
      Structured = sum(structural_miss),
      Simulated = sum( (!!miss_sym) & !structural_miss ),
      Available = n() - Structured - Simulated,
      .groups = "drop"
    )
  
  table_df <- as.data.frame(t(df_summary))
  colnames(table_df) <- table_df[1, ]
  table_df <- table_df[-1, , drop = FALSE]
  table_df <- cbind(Timepoint = rownames(table_df), table_df)
  rownames(table_df) <- NULL
  
  tg <- tableGrob(
    table_df,
    rows = NULL,
    theme = ttheme_minimal(
      core = list(
        fg_params = list(col = "black"),
        bg_params = list(fill = NA),
        padding = unit(c(4, 4), "mm")
      ),
      colhead = list(
        fg_params = list(col = "black"),
        bg_params = list(fill = NA)
      )
    )
  )

  tg$grobs[tg$layout$name == "core-bg"] <- lapply(
    tg$grobs[tg$layout$name == "core-bg"],
    function(x) {
      x$gp <- gpar(col = "black", fill = NA, lwd = 0.5)
      x
    }
  )
  
  tg$grobs[tg$layout$name == "colhead-bg"] <- lapply(
    tg$grobs[tg$layout$name == "colhead-bg"],
    function(x) {
      x$gp <- gpar(col = "black", fill = NA, lwd = 0.8)
      x
    }
  )
  
  df_long <- df %>%
    group_by(time) %>%
    summarise(
      prop_struct = sum(structural_miss) / n(),
      prop_sim = sum( (!!miss_sym) & !structural_miss ) / n(),
      prop_available = 1 - prop_struct - prop_sim,
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(
      cols = c(prop_available, prop_struct, prop_sim),
      names_to = "type",
      values_to = "proportion"
    )
  
  fill_colors <- c(
    "#000000", "#E69F00", "#56B4E9", "#009E73",
    "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
  )
  
  
  fill_labels <- c(
    prop_available = "Available",
    prop_struct = "Structurally Missing",
    prop_sim = title_text
  )
  
  p <- ggplot(df_long, aes(x = factor(time), y = proportion, fill = type)) +
    geom_col(position = "stack") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = fill_colors, labels = fill_labels) +
    labs(
      x = "Timepoint",
      y = "Proportion",
      fill = "Status",
      title = title_text,
      subtitle = "Available vs Structured vs Simulated Missingness"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(hjust = 1)
    )
  
  left_table <- arrangeGrob(
    tg,
    top = NULL,
    padding = unit(0, "mm")
  )
  
  layout <- rbind(
    c(1, 1),
    c(2, 2)
  )
  
  grid.arrange(
    p,
    left_table,
    layout_matrix = layout,
    widths = c(1, 0.0001),
    heights = c(0.8, 0.2)
  )
}

# MCAR
plot_missingness(
  df = df_mcar,
  miss_var = "mcar_miss",
  title_text = "MCAR"
)

# MAR
plot_missingness(
  df = df_mar,
  miss_var = "mar_miss",
  title_text = "MAR"
)

# MNAR
plot_missingness(
  df = df_mnar,
  miss_var = "mnar_miss",
  title_text = "MNAR"
)

### Summary of given domain or total ###########################################

mean_by_time_all <- function(data) {
  
  vars <- c(
    "pwb_total", "swb_total", "ewb_total",
    "fwb_total", "bcs_total", "fact_b_total"
  )
  
  data %>%
    group_by(time) %>%
    summarise(
      across(
        all_of(vars),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd   = ~ sd(.x, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
}

qol_summary <- mean_by_time_all(df_qol)
mcar_summary <- mean_by_time_all(df_mcar)
mar_summary <- mean_by_time_all(df_mar)
mnar_summary <- mean_by_time_all(df_mnar)

gather_var_all <- function(var_name) {
  
  mean_col <- paste0(var_name, "_mean")
  sd_col   <- paste0(var_name, "_sd")
  
  df1 <- qol_summary  %>% dplyr::select(time,
                                 mean = all_of(mean_col),
                                 sd = all_of(sd_col)) %>%
    mutate(dataset = "QOL")
  df2 <- mcar_summary %>% dplyr::select(time,
                                 mean = all_of(mean_col),
                                 sd = all_of(sd_col)) %>%
    mutate(dataset = "MCAR")
  df3 <- mar_summary  %>% dplyr::select(time,
                                 mean = all_of(mean_col),
                                 sd = all_of(sd_col)) %>%
    mutate(dataset = "MAR")
  df4 <- mnar_summary %>% dplyr::select(time,
                                 mean = all_of(mean_col),
                                 sd = all_of(sd_col)) %>%
    mutate(dataset = "MNAR")
  
  # Combine
  bind_rows(df1, df2, df3, df4) %>%
    arrange(time, dataset)
}

okabe_ito <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

pwb_plot <- gather_var_all("pwb_total")

p_pwb_miss <- ggplot(pwb_plot, aes(x = time, y = mean, colour = dataset)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = okabe_ito) +
  scale_y_continuous(limits = c(10, 13)) +
  labs(
    title = "",
    x = "Time (years)",
    y = "Mean PWB Score"
  ) +
  theme_minimal()

bcs_plot <- gather_var_all("bcs_total")

p_bcs_miss <- ggplot(bcs_plot, aes(x = time, y = mean, colour = dataset)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = okabe_ito) +
  scale_y_continuous(limits = c(18, 23)) +
  labs(
    title = "",
    x = "Time (years)",
    y = "Mean BCS Score"
  ) +
  theme_minimal()

fact_b_plot <- gather_var_all("fact_b_total")

p_fact_b_miss <- ggplot(fact_b_plot, aes(x = time, y = mean, colour = dataset)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_color_manual(values = okabe_ito) +
  scale_y_continuous(limits = c(75, 85), breaks = seq(75, 85, by = 1)) +
  labs(
    title = "",
    x = "Time (years)",
    y = "Mean FACT-B Score"
  ) +
  theme_minimal()

p1 <- p_bcs_miss + theme(legend.position = "none")
p2 <- p_pwb_miss + theme(legend.position = "none")
p3 <- p_fact_b_miss + theme(legend.position = "right")

(p1 | p2) / (p3)
