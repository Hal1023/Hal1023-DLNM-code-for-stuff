# =============================================================================
# GLUCOSE -> KETONE LAGGED RELATIONSHIP: DLNM-GAM MODELLING PIPELINE
# =============================================================================
# Premise: two time-stamped numeric streams per sensor-subject pair.
#   - ketones    : response, zero-heavy, right-skewed, non-negative
#   - glucose    : predictor, continuous ~4-15 mmol/L
#   - uid        : unique sensor-subject identifier (subject_id + sensor)
#   - time       : POSIXct, nominally every 5 minutes
#
# Goal: estimate the lagged effect of glucose on ketones in a way that is
#   (a) distributional honest about the many zeros in ketones
#   (b) flexible about the shape of the glucose effect at each lag
#   (c) accounts for within-subject clustering
#   (d) computationally tractable at 100s of subjects x 10 days x 5-min
#
# Model family: GAM with DLNM cross-basis for glucose lags
# Engine: mgcv (primary), with dlnm for cross-basis construction
# =============================================================================


# -----------------------------------------------------------------------------
# 0. PACKAGES
# -----------------------------------------------------------------------------

pkgs <- c(
  "tidyverse",    # data wrangling
  "lubridate",    # timestamp handling
  "mgcv",         # GAM engine (gam, bam, gamm)
  "dlnm",         # distributed lag non-linear models (cross-basis)
  "glmmTMB",      # hurdle/ZI GAM alternative with mixed effects
  "gratia",       # ggplot2-based GAM diagnostics (draw, appraise)
  "patchwork",    # plot composition
  "data.table",   # fast rolling operations for lag matrix construction
  "zoo"           # rollmean, na.locf for irregular series
)

new_pkgs <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)
invisible(lapply(pkgs, library, character.only = TRUE))


# -----------------------------------------------------------------------------
# 1. DATA PREPARATION
# -----------------------------------------------------------------------------
# Input: df with columns — uid, glucose, ketones, time (POSIXct)
# We construct the lagged glucose matrix here.
# Lags are in units of 5-minute intervals (the nominal sensor period).
# We use up to L = 36 lags = 3 hours of glucose history.
#
# Key decision: we use ROWS as lag units rather than interpolating to a
# fixed grid. For sensors that are truly ~5 min apart this is correct.
# If your timestamps are irregular (dropouts, restarts), interpolate first
# to a regular 5-min grid per uid before running this script.

# --- Parameters you should tune ---
MAX_LAG_INTERVALS  <- 36    # 36 x 5min = 3 hours lag window
N_KNOTS_EXPOSURE   <- 4     # knots over glucose values (exposure dimension)
N_KNOTS_LAG        <- 4     # knots over lag dimension
# Basis type for each dimension (see Section 3 for discussion)
BASIS_EXPOSURE     <- "cr"  # cubic regression splines for glucose values
BASIS_LAG          <- "ps"  # P-splines for lag (better for large n, ordered lags)

# -------------------------------------------------------
# 1.1  Build the lagged exposure matrix (Q matrix in DLNM notation)
# -------------------------------------------------------
# For each observation at time t, Q[t, l] = glucose at time t - l intervals.
# This is the core input to the cross-basis.

build_lag_matrix <- function(df, max_lag = MAX_LAG_INTERVALS) {
  # Work series by series so lags don't bleed across subjects/sensors
  df <- df |>
    arrange(uid, time) |>
    group_by(uid) |>
    mutate(obs_index = row_number()) |>
    ungroup()
  
  uid_list  <- unique(df$uid)
  lag_cols  <- paste0("gl_lag", 0:max_lag)  # lag0 = contemporaneous glucose
  
  result_list <- lapply(uid_list, function(u) {
    sub <- df |> filter(uid == u) |> arrange(time)
    n   <- nrow(sub)
    
    # Build Q: n x (max_lag+1) matrix
    Q <- matrix(NA_real_, nrow = n, ncol = max_lag + 1)
    for (l in 0:max_lag) {
      # lag l means we look back l rows
      Q[, l + 1] <- c(rep(NA_real_, l), sub$glucose[1:(n - l)])
    }
    colnames(Q) <- lag_cols
    
    cbind(sub, as.data.frame(Q))
  })
  
  bind_rows(result_list)
}

df_lagged <- build_lag_matrix(df)

# Drop rows where any lag is missing (the first max_lag rows of each series)
df_model <- df_lagged |>
  filter(complete.cases(select(., starts_with("gl_lag")))) |>
  # Scale glucose for numerical stability in cross-basis
  mutate(glucose_sc = scale(glucose)[, 1])

message(sprintf(
  "Modelling dataset: %d observations | %d series | %d subjects",
  nrow(df_model),
  n_distinct(df_model$uid),
  n_distinct(df_model$subject_id)
))


# -----------------------------------------------------------------------------
# 2. CROSS-BASIS CONSTRUCTION (DLNM)
# -----------------------------------------------------------------------------
# The cross-basis B is a matrix where each column is a basis function over
# the joint (exposure value x lag) space. Its dimension is:
#   n_obs  x  (N_KNOTS_EXPOSURE * N_KNOTS_LAG)
#
# The GAM then estimates coefficients on these basis functions.
# The fitted surface f(exposure, lag) is the bilinear combination.
#
# Exposure dimension (glucose values):
#   We use natural cubic splines ("cr") with knots placed at quantiles of the
#   observed glucose distribution. This allows a nonlinear dose-response at
#   each lag value.
#
# Lag dimension:
#   We use P-splines ("ps") — these penalise differences between adjacent
#   spline coefficients, producing smooth lag curves without end-knot
#   boundary artefacts. Better than "cr" here because:
#     1. Lags are inherently ordered (you care about smoothness along the lag axis)
#     2. P-splines work well at large n (penalty is on coefficient differences,
#        not a dense penalty matrix)
#     3. Natural cubic splines at the lag dimension can produce artefacts
#        near lag=0 and lag=max_lag

# Extract the glucose lag matrix as a plain matrix (required by dlnm::crossbasis)
Q_matrix <- df_model |>
  select(starts_with("gl_lag")) |>
  as.matrix()

# Knot placement: quantiles of the observed glucose distribution
# (avoids placing knots in data-sparse regions)
gluc_knots <- quantile(df_model$glucose, probs = seq(0.1, 0.9,
                        length.out = N_KNOTS_EXPOSURE), na.rm = TRUE)

# Build the cross-basis
cb_glucose <- crossbasis(
  Q_matrix,
  lag    = c(0, MAX_LAG_INTERVALS),
  argvar = list(fun   = BASIS_EXPOSURE,
                knots = gluc_knots),
  arglag = list(fun   = BASIS_LAG,
                df    = N_KNOTS_LAG)   # P-splines: df controls flexibility
)

message(sprintf(
  "Cross-basis constructed: %d columns (= %d exposure knots x %d lag df)",
  ncol(cb_glucose), N_KNOTS_EXPOSURE, N_KNOTS_LAG
))


# =============================================================================
# 3. MODEL FITTING: THREE CANDIDATE SPECIFICATIONS
# =============================================================================
# We fit three models, increasing in complexity, and compare.
#
# M1: bam() + Tweedie, subject random intercept — your current approach improved
# M2: bam() + Tweedie, subject random intercept + AR(1) residual correlation
# M3: Two-part hurdle via glmmTMB (separate binary and gamma GAMs)
#
# WHY bam() OVER gam() AT THIS SCALE:
#   bam() is mgcv's "big data" GAM. It uses:
#   - Discretisation of smooth covariates (reduces effective n for penalty calc)
#   - Sparse matrix methods for the random effects
#   - Parallel computation via nthreads
#   For 100k+ rows, gam() with REML becomes very slow; bam() with fREML
#   (fast REML) is the right call.

# -------------------------------------------------------
# M1: Tweedie GAM with subject random intercept
# -------------------------------------------------------
# Formula components:
#   cb_glucose           : the DLNM cross-basis (fixed smooth)
#   s(uid, bs="re")      : subject-level random intercept (bs="re" in mgcv)
#                          equivalent to a Gaussian random effect b_i ~ N(0, sigma^2)
#   family=tw()          : Tweedie with estimated power p in (1,2)
#   method="fREML"       : fast REML smoothness selection (required for bam)
#   discrete=TRUE        : discretise covariates for speed (key at large n)
#   nthreads             : set to number of available cores

n_cores <- max(1L, parallel::detectCores() - 1L)

m1_tweedie <- bam(
  ketones ~ cb_glucose + s(uid, bs = "re"),
  family  = tw(),          # Tweedie: estimates power p internally
  method  = "fREML",
  discrete = TRUE,
  nthreads = n_cores,
  data    = df_model
)

summary(m1_tweedie)
cat("\nTweedie power p:", m1_tweedie$family$getTheta(TRUE), "\n")
# p close to 1 -> Poisson-like; p close to 2 -> Gamma-like
# p around 1.5 is typical for count-like continuous data with many zeros


# -------------------------------------------------------
# M2: Tweedie GAM + AR(1) within-series correlation
# -------------------------------------------------------
# The observations within each uid are temporally autocorrelated.
# Ignoring this makes the effective sample size look larger than it is,
# which inflates the apparent precision of all smooth estimates.
#
# bam() supports AR(1) via rho argument.
# We estimate rho from the M1 residuals, then refit.
#
# AR(1) structure: epsilon_t = rho * epsilon_{t-1} + nu_t
# This is applied WITHIN each uid (series), not across series.
# The AR.start argument marks the start of each new series so
# the AR structure doesn't bleed across uid boundaries.

# Estimate rho from M1 residuals
resid_m1  <- residuals(m1_tweedie, type = "deviance")
rho_hat   <- acf(resid_m1, lag.max = 1, plot = FALSE)$acf[2]
message(sprintf("Estimated AR(1) rho from M1 residuals: %.4f", rho_hat))

# Mark series starts for the AR structure
df_model <- df_model |>
  group_by(uid) |>
  mutate(is_series_start = row_number() == 1) |>
  ungroup()

m2_tweedie_ar <- bam(
  ketones  ~ cb_glucose + s(uid, bs = "re"),
  family   = tw(),
  method   = "fREML",
  discrete = TRUE,
  nthreads = n_cores,
  rho      = rho_hat,
  AR.start = df_model$is_series_start,
  data     = df_model
)

summary(m2_tweedie_ar)

# Compare M1 and M2 — is the AR correction material?
cat("\nM1 AIC:", AIC(m1_tweedie))
cat("\nM2 AIC:", AIC(m2_tweedie_ar), "\n")
cat("If M2 AIC is substantially lower (>4 units), the AR structure matters.\n")


# -------------------------------------------------------
# M3: Two-part hurdle (explicit separation of zero process)
# -------------------------------------------------------
# The Tweedie lumps the zero-generating process and the continuous process
# into a single parameter p. A hurdle model explicitly says:
#   Part 1: Is this observation zero?  -> logistic smooth GAM
#   Part 2: Given nonzero, how large?  -> gamma GAM
#
# This is more transparent and lets you ask separately:
#   "Does glucose predict WHETHER ketones are elevated?"
#   "Does glucose predict HOW HIGH ketones are when elevated?"
# These may have different lag structures, which is substantively interesting.
#
# We fit each part as a separate bam(), sharing the same cross-basis.

# Part 1: binary (zero vs nonzero)
df_model <- df_model |> mutate(keto_pos = as.integer(ketones > 0))

m3a_binary <- bam(
  keto_pos ~ cb_glucose + s(uid, bs = "re"),
  family   = binomial(link = "logit"),
  method   = "fREML",
  discrete = TRUE,
  nthreads = n_cores,
  data     = df_model
)

# Part 2: gamma for nonzero observations only
df_pos <- df_model |> filter(ketones > 0)

# Rebuild cross-basis for the subset (row subset of Q_matrix)
Q_pos      <- Q_matrix[df_model$ketones > 0, ]
cb_pos     <- crossbasis(
  Q_pos,
  lag    = c(0, MAX_LAG_INTERVALS),
  argvar = list(fun = BASIS_EXPOSURE, knots = gluc_knots),
  arglag = list(fun = BASIS_LAG, df = N_KNOTS_LAG)
)

m3b_gamma <- bam(
  ketones ~ cb_pos + s(uid, bs = "re"),
  family   = Gamma(link = "log"),
  method   = "fREML",
  discrete = TRUE,
  nthreads = n_cores,
  data     = df_pos
)

summary(m3a_binary)
summary(m3b_gamma)


# =============================================================================
# 4. SMOOTHNESS / BASIS SENSITIVITY CHECK
# =============================================================================
# With large data, the smoothness penalty drives the effective degrees of
# freedom (edf). We check whether our knot choices are adequate:
#   - If edf ~= k' (basis dimension), the smooth is under-specified
#     (wants more flexibility than we've allowed)
#   - If edf << k', the smooth is over-specified but penalised correctly
#
# k.check() is the formal test in mgcv.

cat("\n--- Basis dimension check for M1 ---\n")
k.check(m1_tweedie)
# Look at k' and edf columns. If p-value is small and edf ~ k',
# increase N_KNOTS_EXPOSURE or N_KNOTS_LAG and refit.


# =============================================================================
# 5. MODEL DIAGNOSTICS
# =============================================================================

# gratia::appraise() gives a ggplot2 version of gam.check() plots:
#   - QQ plot of deviance residuals
#   - Residuals vs linear predictor
#   - Histogram of residuals
#   - Response vs fitted values

appraise(m2_tweedie_ar, point_col = "#378ADD", point_alpha = 0.15,
         ci_col = "#0C447C") +
  plot_annotation(title = "M2 Tweedie-AR(1) diagnostics")

appraise(m3a_binary, point_col = "#534AB7", point_alpha = 0.15) +
  plot_annotation(title = "M3a Binary (zero vs nonzero) diagnostics")

appraise(m3b_gamma, point_col = "#D85A30", point_alpha = 0.15) +
  plot_annotation(title = "M3b Gamma (nonzero ketones) diagnostics")


# =============================================================================
# 6. VISUALISING THE LAGGED RELATIONSHIP (THE KEY RESULT)
# =============================================================================
# crosspred() from dlnm computes the predicted surface f(exposure, lag)
# over a grid of glucose values and lag values, with confidence intervals.
#
# This is what you're ultimately after: a quantified, smooth description
# of how glucose at various lags predicts ketones.

predict_dlnm_surface <- function(model, crossbasis_obj,
                                 gluc_range = seq(3, 15, by = 0.5),
                                 coef_name  = NULL) {
  # If model has random effects, we need to supply the coefficients
  # corresponding just to the cross-basis (not the RE terms)
  pred <- crosspred(
    basis      = crossbasis_obj,
    model      = model,
    at         = gluc_range,
    cen        = median(df_model$glucose, na.rm = TRUE),  # centre = reference
    cumul      = FALSE   # FALSE = lag-specific; TRUE = cumulative effect
  )
  pred
}

pred_m1 <- predict_dlnm_surface(m1_tweedie, cb_glucose)
pred_m2 <- predict_dlnm_surface(m2_tweedie_ar, cb_glucose)


# -------------------------------------------------------
# 6.1  3D surface plot: f(glucose, lag)
# -------------------------------------------------------
# The surface shows the estimated effect on ketones (on the response scale)
# as a function of both glucose level and lag in 5-min intervals.

plot_dlnm_surface_3d <- function(pred, title = "Glucose-ketone lag surface") {
  # Use base R persp for the 3D surface (no ggplot equivalent for smooth 3D)
  par(mar = c(2, 2, 3, 1))
  plot(pred,
       xlab = "Glucose (mmol/L)",
       zlab = "Predicted ketones (centred)",
       ylab = "Lag (5-min intervals)",
       phi  = 25, theta = -50,
       col  = colorRampPalette(c("#E6F1FB", "#378ADD", "#042C53"))(100),
       main = title)
}

plot_dlnm_surface_3d(pred_m1, "M1: Tweedie GAM — glucose-ketone lag surface")
plot_dlnm_surface_3d(pred_m2, "M2: Tweedie-AR(1) GAM — glucose-ketone lag surface")


# -------------------------------------------------------
# 6.2  Slice plots: effect at specific lag values
# -------------------------------------------------------
# More readable than the 3D surface: show the glucose dose-response curve
# at specific lag times (e.g. lag 0, lag 6, lag 12, lag 24 = immediate,
# 30min, 1hr, 2hr).

plot_lag_slices <- function(pred, lags = c(0, 6, 12, 24, 36),
                             title = "Effect at specific lags") {
  lag_names  <- paste0("Lag ", lags * 5, " min")
  gluc_vals  <- as.numeric(rownames(pred$matfit))
  
  slice_df <- lapply(seq_along(lags), function(i) {
    l <- lags[i]
    tibble(
      glucose  = gluc_vals,
      estimate = pred$matfit[, l + 1],
      lower    = pred$matlow[,  l + 1],
      upper    = pred$mathigh[, l + 1],
      lag      = lag_names[i]
    )
  }) |> bind_rows() |>
    mutate(lag = factor(lag, levels = lag_names))
  
  ggplot(slice_df, aes(x = glucose, y = estimate,
                        colour = lag, fill = lag)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    scale_colour_viridis_d(option = "plasma", end = 0.85) +
    scale_fill_viridis_d(option = "plasma", end = 0.85) +
    labs(title   = title,
         x       = "Glucose (mmol/L)",
         y       = "Effect on ketones (relative to reference)",
         colour  = "Lag",
         fill    = "Lag") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
}

plot_lag_slices(pred_m1, title = "M1: Glucose effect on ketones at each lag")
plot_lag_slices(pred_m2, title = "M2 (AR corrected): Glucose effect at each lag")


# -------------------------------------------------------
# 6.3  Lag dimension slice: effect over time for specific glucose levels
# -------------------------------------------------------
# This answers: "For a glucose value of X, how does the effect on ketones
# unfold over the next 3 hours?"

plot_exposure_slices <- function(pred, glucose_levels = c(4, 6, 8, 10, 12),
                                  max_lag = MAX_LAG_INTERVALS,
                                  title = "Lag profile at fixed glucose levels") {
  lag_seq    <- 0:max_lag
  lag_times  <- lag_seq * 5  # in minutes
  gluc_vals  <- as.numeric(rownames(pred$matfit))
  
  exp_df <- lapply(glucose_levels, function(g) {
    # Find closest glucose row in the prediction grid
    idx <- which.min(abs(gluc_vals - g))
    tibble(
      lag_min  = lag_times,
      estimate = pred$matfit[idx, ],
      lower    = pred$matlow[idx, ],
      upper    = pred$mathigh[idx, ],
      glucose  = paste0(g, " mmol/L")
    )
  }) |> bind_rows() |>
    mutate(glucose = factor(glucose,
                            levels = paste0(sort(glucose_levels), " mmol/L")))
  
  ggplot(exp_df, aes(x = lag_min, y = estimate,
                      colour = glucose, fill = glucose)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    scale_x_continuous(breaks = seq(0, max_lag * 5, by = 30)) +
    scale_colour_viridis_d(option = "viridis", end = 0.85) +
    scale_fill_viridis_d(option = "viridis", end = 0.85) +
    labs(title  = title,
         x      = "Lag (minutes)",
         y      = "Effect on ketones (relative to reference)",
         colour = "Glucose level",
         fill   = "Glucose level") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right")
}

plot_exposure_slices(pred_m1, title = "M1: How glucose level affects ketones over time")
plot_exposure_slices(pred_m2, title = "M2 (AR corrected): Lag profile by glucose level")


# -------------------------------------------------------
# 6.4  Cumulative effect (overall impact integrated over the full lag window)
# -------------------------------------------------------
# crosspred with cumul=TRUE gives the effect of glucose *accumulated*
# over all lags up to L. This collapses the surface to a single curve:
# "overall effect of glucose level X on ketones over the lag window"

pred_cumul <- crosspred(
  basis = cb_glucose,
  model = m2_tweedie_ar,
  at    = seq(3, 15, by = 0.25),
  cen   = median(df_model$glucose, na.rm = TRUE),
  cumul = TRUE
)

cumul_df <- tibble(
  glucose  = as.numeric(rownames(pred_cumul$allfit)),
  estimate = pred_cumul$allfit[, ncol(pred_cumul$allfit)],
  lower    = pred_cumul$alllow[, ncol(pred_cumul$alllow)],
  upper    = pred_cumul$allhigh[, ncol(pred_cumul$allhigh)]
)

p_cumul <- ggplot(cumul_df, aes(x = glucose, y = estimate)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "#534AB7", alpha = 0.2) +
  geom_line(colour = "#534AB7", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title   = "Cumulative effect of glucose on ketones (M2, integrated over all lags)",
       subtitle = "Relative to median glucose; positive = higher ketones",
       x       = "Glucose (mmol/L)",
       y       = "Cumulative effect on ketones") +
  theme_minimal(base_size = 12)

print(p_cumul)


# =============================================================================
# 7. COMPARING THE THREE MODEL FAMILIES
# =============================================================================

compare_models <- function(m1, m2, m3a, m3b) {
  # AIC is comparable within a family but NOT between Tweedie and the
  # two-part model (different likelihoods). We report them separately.
  
  cat("\n", strrep("=", 60), "\n")
  cat("MODEL COMPARISON SUMMARY\n")
  cat(strrep("=", 60), "\n\n")
  
  cat("--- Tweedie models (directly comparable via AIC) ---\n")
  tweedie_comp <- tibble(
    Model    = c("M1: Tweedie, no AR",
                 "M2: Tweedie + AR(1)"),
    EDF      = c(sum(m1$edf), sum(m2$edf)),
    AIC      = c(AIC(m1), AIC(m2)),
    Dev_Expl = c(summary(m1)$dev.expl,
                 summary(m2)$dev.expl) * 100
  ) |> mutate(delta_AIC = AIC - min(AIC))
  print(tweedie_comp)
  
  cat("\n--- Two-part hurdle (separate likelihoods) ---\n")
  hurdle_comp <- tibble(
    Component = c("M3a: Binary (logit)", "M3b: Gamma (nonzero)"),
    EDF       = c(sum(m3a$edf), sum(m3b$edf)),
    AIC       = c(AIC(m3a), AIC(m3b)),
    Dev_Expl  = c(summary(m3a)$dev.expl,
                  summary(m3b)$dev.expl) * 100
  )
  print(hurdle_comp)
  
  cat("\n--- Tweedie p (power parameter) ---\n")
  cat("M1 p:", m1$family$getTheta(TRUE), "\n")
  cat("M2 p:", m2$family$getTheta(TRUE), "\n")
  cat("p in (1,2); closer to 1 = Poisson-like, 2 = Gamma-like\n")
  cat("If p < 1.1 or > 1.9, the Tweedie may not be the right family.\n")
  
  cat("\n--- AR(1) rho estimate ---\n")
  cat(sprintf("Residual autocorrelation (rho) from M1: %.4f\n", rho_hat))
  cat("If rho > 0.3, the AR correction (M2) is substantively important.\n")
}

compare_models(m1_tweedie, m2_tweedie_ar, m3a_binary, m3b_gamma)


# =============================================================================
# 8. CONCURVITY CHECK
# =============================================================================
# In a GAM, concurvity is the analogue of collinearity — if two smooth
# terms are highly related, coefficient estimates become unstable.
# Here the cross-basis columns are correlated by construction (lagged values
# of the same series). concurvity() quantifies this.
# Values > 0.8 in the "estimate" column suggest a problem.

cat("\n--- Concurvity check (M2) ---\n")
conc <- concurvity(m2_tweedie_ar, full = TRUE)
print(round(conc, 3))
# If high: consider reducing N_KNOTS_LAG or constraining the lag basis


# =============================================================================
# 9. PREDICTION FOR A NEW SERIES (INFERENCE USE CASE)
# =============================================================================
# Given a new glucose series, predict the expected ketone trajectory.
# This is what you'd do for a new experimental subject.

predict_new_series <- function(model, crossbasis_obj, new_glucose_series,
                                max_lag = MAX_LAG_INTERVALS) {
  n   <- length(new_glucose_series)
  
  # Build Q matrix for the new series
  Q_new <- matrix(NA_real_, nrow = n, ncol = max_lag + 1)
  for (l in 0:max_lag) {
    Q_new[, l + 1] <- c(rep(NA_real_, l), new_glucose_series[1:(n - l)])
  }
  colnames(Q_new) <- paste0("gl_lag", 0:max_lag)
  
  # Cross-basis for new data
  cb_new <- crossbasis(
    Q_new,
    lag    = c(0, max_lag),
    argvar = list(fun = BASIS_EXPOSURE, knots = gluc_knots),
    arglag = list(fun = BASIS_LAG, df = N_KNOTS_LAG)
  )
  
  # Need a representative uid for the RE term (population average = exclude RE)
  # Using exclude = "s(uid)" gives population-level prediction
  df_new <- as.data.frame(Q_new) |>
    mutate(uid = levels(df_model$uid)[1])  # placeholder
  
  pred <- predict(
    model,
    newdata = df_new,
    exclude = "s(uid)",       # population-level (average over subjects)
    type    = "response",
    se.fit  = TRUE
  )
  
  tibble(
    t        = seq_len(n),
    glucose  = new_glucose_series,
    keto_pred = pred$fit,
    keto_se  = pred$se.fit,
    keto_lwr = pred$fit - 1.96 * pred$se.fit,
    keto_upr = pred$fit + 1.96 * pred$se.fit
  ) |> filter(t > max_lag)
}

# Example usage: create a synthetic 6-hour glucose trace
set.seed(1)
new_glucose <- c(
  seq(5, 5, length.out = 24),   # 2hr baseline
  seq(5, 10, length.out = 6),   # meal spike rising
  seq(10, 6, length.out = 18),  # return to baseline
  seq(6, 4.5, length.out = 24)  # post-meal dip
)

pred_new <- predict_new_series(m2_tweedie_ar, cb_glucose, new_glucose)

p_pred <- ggplot(pred_new, aes(x = t * 5)) +  # convert to minutes
  geom_ribbon(aes(ymin = keto_lwr, ymax = keto_upr),
              fill = "#D85A30", alpha = 0.2) +
  geom_line(aes(y = keto_pred), colour = "#D85A30", linewidth = 0.8) +
  geom_line(aes(y = glucose / 10), colour = "#378ADD",
            linewidth = 0.6, linetype = "dashed") +  # rescale for dual axis
  scale_y_continuous(
    name     = "Predicted ketones (mmol/L)",
    sec.axis = sec_axis(~ . * 10, name = "Glucose (mmol/L) [dashed]")
  ) +
  labs(title = "Population-level prediction for a new glucose trace",
       x     = "Time (minutes)") +
  theme_minimal(base_size = 12)

print(p_pred)


# =============================================================================
# 10. SUMMARY AND GUIDANCE FOR MODEL SELECTION
# =============================================================================

cat("
=================================================================
MODEL SELECTION GUIDANCE
=================================================================

Step 1: Check Tweedie power p from M1/M2
  p < 1.1  -> Poisson-like; zeros are not structural, just rare counts
              Consider nb() (negative binomial) family instead
  p in 1.2-1.8  -> Classic Tweedie regime; distribution is appropriate
  p > 1.9  -> Gamma-like; very few genuine zeros; maybe zeros are noise
              Consider using Gamma directly on log-transformed + epsilon

Step 2: Check AR(1) rho from M1 residuals
  rho < 0.1  -> AR correction is minor; M1 is adequate
  rho > 0.3  -> Use M2; the temporal autocorrelation is real and inflates edf
  rho > 0.7  -> Consider a full time series model (GAMM with corAR1 in nlme)

Step 3: Basis adequacy (k.check output)
  edf ~ k' and p < 0.05  -> Increase knot count
  Recommended starting point: N_KNOTS_EXPOSURE = 6, N_KNOTS_LAG = 6
  Maximum sensible: N_KNOTS_EXPOSURE = 8, N_KNOTS_LAG = 8
  (beyond this, computation dominates; use thin plate splines if needed)

Step 4: Tweedie vs hurdle
  If your scientific question is ONLY 'does glucose predict ketone level':
    -> Tweedie (M2) is sufficient
  If you want to separately quantify:
    'does glucose predict WHETHER ketones appear'
    AND 'does glucose predict HOW HIGH they go':
    -> Two-part hurdle (M3a + M3b) is more interpretable
  If these two questions have DIFFERENT optimal lags:
    -> Strong argument for M3 over M2

Step 5: Concurvity
  If concurvity > 0.8, the cross-basis is too flexible for your data.
  Reduce N_KNOTS_LAG by 1 and recheck.

Primary recommendation for your setup:
  M2 (Tweedie + AR(1) + subject RE) with N_KNOTS_EXPOSURE = 6, N_KNOTS_LAG = 5
  is the best single-model starting point before moving to mechanistic modelling.
=================================================================
")
