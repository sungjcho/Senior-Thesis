# Bootstrap for GWR

library(sp)
library(spgwr)
library(dplyr)
library(purrr)
library(furrr)
library(tidyr)
library(readr)
library(knitr)
library(tibble)
library(tinytex)

# Data
data_path <- "fertility_rate_data - final2.csv"
fertility_rate_data <- read.csv(data_path)

coordinates(fertility_rate_data) <- ~ X_coord + Y_coord

gwr_formula <- TFR ~ childcare_centers + population_density + net_migration_rate +
  local_population_growth + foreign_marriage_rate + female_workers +
  foreigner_households

# Select optimal bandwidth for GWR
gwr_bw <- gwr.sel(gwr_formula, data = fertility_rate_data)

# Running the GWR
gwr_original <- gwr(
  gwr_formula,
  data = fertility_rate_data,
  bandwidth = gwr_bw,
  hatmatrix = TRUE,
  se.fit = TRUE
)

vars <- all.vars(gwr_formula)[-1]

# ---------------------------------------------
# BOOTSTRAP ITERATIONS
# ---------------------------------------------

plan(multisession)

# Output directory
dir.create("boot_results_5000", showWarnings = FALSE)

# Number of iterations
n_boot <- 5000

future_walk(1:n_boot, function(i) {
  set.seed(i)
  boot_indices <- sample(1:nrow(fertility_rate_data), replace = TRUE)
  boot_data <- fertility_rate_data[boot_indices, ]
  df <- as.data.frame(boot_data)
  coord <- coordinates(boot_data)
  spdf <- SpatialPointsDataFrame(coord, data = df, match.ID = FALSE)
  
  boot_gwr <- gwr(
    gwr_formula,
    data = spdf,
    bandwidth = gwr_bw,
    hatmatrix = TRUE,
    se.fit = TRUE
  )
  
  coeffs <- as.data.frame(boot_gwr$SDF@data[, c("(Intercept)", vars)])
  coeffs$location <- 1:nrow(coeffs)
  
  write_rds(coeffs, sprintf("boot_results_5000/boot_%04d.rds", i))
}, .progress = TRUE)

# ---------------------------------------------
# SUMMARY STATISTICS
# ---------------------------------------------

# Combine bootstrap results
bootstrap_long <- list.files("boot_results_5000", full.names = TRUE) %>%
  set_names() %>%
  imap_dfr(function(path, filename) {
    df <- read_rds(path)
    df$iteration <- as.integer(gsub("[^0-9]", "", filename))
    df
  }) %>%
  pivot_longer(
    cols = -c(iteration, location),
    names_to = "variable",
    values_to = "coefficient"
  )

# Extract original GWR coefficients
original_long <- gwr_original$SDF@data %>%
  select(all_of(vars)) %>%
  mutate(location = 1:n()) %>%
  pivot_longer(cols = all_of(vars), names_to = "variable", values_to = "original_coef")

# Join region names
location_name_lookup <- fertility_rate_data@data %>%
  mutate(location = 1:n()) %>%
  select(location, region_name)

# Compute summary stats
bootstrap_stats <- bootstrap_long %>%
  group_by(location, variable) %>%
  summarise(
    boot_mean = mean(coefficient, na.rm = TRUE),
    ci_lower = quantile(coefficient, 0.025, na.rm = TRUE),
    ci_upper = quantile(coefficient, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(original_long, by = c("location", "variable")) %>%
  left_join(
    bootstrap_long %>%
      left_join(original_long, by = c("location", "variable")) %>%
      group_by(location, variable) %>%
      summarise(
        pct_opposite_sign = mean(sign(coefficient) != sign(original_coef)) * 100,
        .groups = "drop"
      ),
    by = c("location", "variable")
  ) %>%
  left_join(location_name_lookup, by = "location") %>%
  mutate(across(where(is.numeric), ~ round(., 4)))

# ---------------------------------------------
# STABILITY SUMMARIES (ACROSS LOCATIONS)
# ---------------------------------------------

# Count of CIs excluding 0 per variable
ci_summary <- bootstrap_stats %>%
  mutate(ci_excludes_0 = ifelse(ci_lower > 0 | ci_upper < 0, 1, 0)) %>%
  group_by(variable) %>%
  summarise(`CI Excludes 0` = sum(ci_excludes_0), .groups = "drop")

# Binned opposite sign % counts
sign_bin_summary <- bootstrap_stats %>%
  mutate(sign_bin = case_when(
    pct_opposite_sign < 5 ~ "<5%",
    pct_opposite_sign >= 5 & pct_opposite_sign < 25 ~ "[5%–25%)",
    pct_opposite_sign >= 25 & pct_opposite_sign < 50 ~ "[25%–50%)",
    pct_opposite_sign >= 50 ~ "≥50%"
  )) %>%
  group_by(variable, sign_bin) %>%
  summarise(n_locations = n(), .groups = "drop") %>%
  pivot_wider(names_from = sign_bin, values_from = n_locations, values_fill = 0)

# Final summary table of stats
final_summary <- left_join(ci_summary, sign_bin_summary, by = "variable")

print(final_summary)

write.csv(final_summary, "bootstrap_summary_table.csv", row.names = FALSE)
