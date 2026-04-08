# =============================================================================
# DLNM Analysis of Continuous Glucose / Ketone Sensor Data
# Framework: mgcv tensor-product DLNM (Economou-style)
# Outcome: raw ketone levels 
# =============================================================================

library(mgcv)
library(dplyr)
library(tidyr)
library(lubridate)
library(haven)
library(ggplot2)
library(readr)

# -----------------------------------------------------------------------------
# User-defined parameters
# -----------------------------------------------------------------------------

id_col        <- "SUBJECTID"
eff_col       <- "GLUCOSE"
out_col       <- "KETONES"
day_col       <- "DAYS_WORN"
hour_col      <- "HOURS_WORN"
value_col     <- "VALUE"
analyte_col   <- "ANALYTE"
condition_col <- "CONDITION"

L          <- 72
epoch_mins <- 5

group_var  <- condition_col

# -----------------------------------------------------------------------------
# Load / generate data
# -----------------------------------------------------------------------------

num_subjects <- 25
days_per_subject <- 10
readings_per_hour <- 12
total_hours <- 24 * days_per_subject

set.seed(42)

generate_data <- function(subject_id, condition) {
  n_rows <- total_hours * readings_per_hour
  hours_worn <- seq(0, by = 0.0833, length.out = n_rows)
  days_worn <- floor(hours_worn / 24) + 1
  
  # Glucose random walk
  base_g <- if(condition == "Diabetic") 8.5 else 5.5
  glucose <- cumsum(rnorm(n_rows, 0, 0.3)) + base_g
  glucose <- pmax(1.5, pmin(glucose, 20))
  
  # Ketones (continuous)
  ketones <- numeric(n_rows)
  for(i in 1:n_rows) {
    val <- runif(1, 0, 0.1)
    spike_chance <- if(glucose[i] < 3.0) 0.15 else 0.02
    if(runif(1) < spike_chance) {
      val <- rgamma(1, shape = 1.2, scale = 0.8)
      val <- pmin(val + 0.3, 6.0)
    }
    ketones[i] <- val
  }
  
  df_glucose <- data.frame(
    SUBJECTID = subject_id, ANALYTE = "GLUCOSE", DAYS_WORN = days_worn,
    HOURS_WORN = round(hours_worn,4), VALUE = glucose, CONDITION = condition
  )
  df_ketones <- data.frame(
    SUBJECTID = subject_id, ANALYTE = "KETONES", DAYS_WORN = days_worn,
    HOURS_WORN = round(hours_worn,4), VALUE = ketones, CONDITION = condition
  )
  
  bind_rows(df_glucose, df_ketones) %>% arrange(HOURS_WORN, ANALYTE)
}

subject_list <- paste0("S", sprintf("%03d", 1:num_subjects))
conditions <- sample(c("Diabetic", "NonDiabetic"), num_subjects, replace = TRUE)

df_raw <- mapply(generate_data, subject_list, conditions, SIMPLIFY = FALSE) %>% 
  bind_rows()

# -----------------------------------------------------------------------------
# Reshape to wide format
# -----------------------------------------------------------------------------

df_wide <- df_raw %>%
  filter(.data[[analyte_col]] %in% c(eff_col, out_col)) %>%
  pivot_wider(
    id_cols = c(id_col, condition_col, day_col, hour_col),
    names_from = all_of(analyte_col),
    values_from = all_of(value_col)
  ) %>%
  rename(
    subjectid = all_of(id_col),
    condition = all_of(condition_col),
    day       = all_of(day_col),
    hours     = all_of(hour_col),
    glucose   = all_of(eff_col),
    ketones   = all_of(out_col)
  ) %>%
  arrange(subjectid, hours)

# -----------------------------------------------------------------------------
# Construct epoch index + time-of-day
# -----------------------------------------------------------------------------

df_wide <- df_wide %>%
  group_by(subjectid) %>%
  mutate(
    hours_rounded = round(hours * 60 / epoch_mins) * epoch_mins / 60,
    epoch = dense_rank(hours_rounded),
    tod_min = (hours_rounded %% 24) * 60
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------
# Optional: add tiny offset to ketones for Gamma regression
# -----------------------------------------------------------------------------
# df_wide <- df_wide %>% mutate(Y = ketones + 1e-6)

# -----------------------------------------------------------------------------
# Create lagged glucose variables
# -----------------------------------------------------------------------------

df_lagged <- df_wide %>%
  group_by(subjectid) %>%
  arrange(epoch, .by_group = TRUE) %>%
  {
    dat <- .
    for(l in 0:L) {
      dat <- dat %>%
        mutate(!!paste0("glucose_lag", l) := dplyr::lag(glucose, n = l))
    }
    dat
  } %>%
  ungroup() %>%
  drop_na(paste0("glucose_lag", 0:L)) %>%
  mutate(
    subjectid_f = factor(subjectid),
    condition_f = factor(condition),
    day_f       = as.numeric(day)
  )

# -----------------------------------------------------------------------------
# Build exposure + lag matrices
# -----------------------------------------------------------------------------

GLUIndex <- grep("^glucose_lag", colnames(df_lagged))
df_lagged$GLU_mat <- I(as.matrix(df_lagged[, GLUIndex]))
df_lagged$LAG_mat <- I(matrix(
  rep(0:L, each = nrow(df_lagged)), nrow = nrow(df_lagged),
  ncol = L+1
))

# -----------------------------------------------------------------------------
# Fit pooled DLNM (tw regression)
# -----------------------------------------------------------------------------

df_lagged <- df_lagged %>%
  mutate(Y = ketones + 0 )

fit_pooled <- gam(
  Y ~ te(GLU_mat, LAG_mat, k = c(4, 4), bs = c("tp", "tp")) +
    s(tod_min, bs = "cc", k = 10) +
    s(day_f, bs = "tp", k = 5) +
    s(subjectid_f, bs = "re"),
  family = tw(link = "log"),   
  data   = df_lagged,
  method = "REML",
  knots  = list(tod_min = c(0, 1440))
)

#  family = Gamma(link="log"),

summary(fit_pooled)

# -----------------------------------------------------------------------------
# Stratified models (commented out for now)
# -----------------------------------------------------------------------------
# groups <- unique(df_lagged[[group_var]])
# fit_by_group <- list()
# for(grp in groups){
#   sub_data <- df_lagged %>% filter(.data[[group_var]] == grp)
#   sub_data$GLU_mat <- I(as.matrix(sub_data[, GLUIndex]))
#   sub_data$LAG_mat <- I(matrix(rep(0:L, each = nrow(sub_data)), nrow = nrow(sub_data), ncol = L+1))
#   fit_by_group[[grp]] <- gam(
#     Y ~ te(GLU_mat, LAG_mat, k = c(6,6), bs = c("tp","tp")) +
#       s(tod_min, bs = "cc", k = 10) +
#       s(day_f, bs = "tp", k = 5) +
#       s(subjectid_f, bs = "re"),
#     family = Gamma(link="log"),
#     data = sub_data,
#     method = "REML"
#   )
#   cat("Fitted:", grp, "\n")
# }

# -----------------------------------------------------------------------------
# Prediction grid
# -----------------------------------------------------------------------------

glucose_grid <- seq(3, 15, length.out = 100)
lag_grid <- 0:L
dlnm_grid <- expand.grid(glucose = glucose_grid, lag = lag_grid)

example_row <- df_lagged[1, ]
n <- nrow(dlnm_grid)
newdata_grid <- data.frame(
  subjectid_f = example_row$subjectid_f,
  day_f       = example_row$day_f,
  tod_min     = example_row$tod_min
)
newdata_grid <- newdata_grid[rep(1, n), ]

GLU_new <- matrix(dlnm_grid$glucose, nrow = n, ncol = L+1)
LAG_new <- matrix(rep(0:L, each = n), nrow = n, ncol = L+1)
newdata_grid$GLU_mat <- I(GLU_new)
newdata_grid$LAG_mat <- I(LAG_new)

# -----------------------------------------------------------------------------
# Predict
# -----------------------------------------------------------------------------

Xmat <- predict(fit_pooled, type = "lpmatrix", newdata = newdata_grid)
eta <- drop(Xmat %*% coef(fit_pooled))
dlnm_grid$RR <- exp(eta)  # multiplicative effect on expected ketones

# -----------------------------------------------------------------------------
# Visualisation
# -----------------------------------------------------------------------------

ggplot(dlnm_grid, aes(x = glucose, y = lag, fill = log(RR))) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(
    title = "DLNM Surface (Tweddie regression, raw ketones)",
    x = "Glucose (mmol/L)",
    y = paste0("Lag (", epoch_mins, " min epochs)"),
    fill = "log effect"
  ) +
  theme_minimal()

# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------

gam.check(fit_pooled)
concurvity(fit_pooled, full = TRUE)