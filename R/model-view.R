# Flat "model view" read by the C++ backend (D-021).
#
# Everything the C++ rate kernel and solver need is flattened once per
# simulation into plain vectors with C++ conventions:
# * offsets and indices are 0-based;
# * matrices are column-major numeric vectors with their dimensions stated
#   in the field name, e.g. inv_w_f is [substage, genotype] (n_A x n_g);
# * logical flags are integers (0/1);
# * adult sub-stage widths are stored as inverse widths 1 / w (0 for classes a
#   genotype never reaches), so the kernel multiplies instead of dividing and
#   never sees Inf;
# * derived constants are precomputed (e.g. mu_L0 = -log(s_L)).
# The layout of every field is checked again on the C++ side
# (src/model_view.cpp).

.sg_model_view <- function(model) {
  if (!inherits(model, "sg_model")) {
    .sg_stop("`model` must be an sg_model object.")
  }
  s <- .sg_reference_setup(model)
  p <- s$p
  f <- s$forcing
  inv_width <- function(w) {
    iw <- 1 / w
    iw[!is.finite(w)] <- 0
    as.numeric(iw)
  }
  rel <- model$releases
  h <- model$hosts
  lay <- model$layout
  host_start <- function(b) if (is.null(h)) -1L else as.integer(lay$start[lay$block == b] - 1L)
  view <- list(
    # dimensions
    n_g = as.integer(s$n_g), n_E = as.integer(s$n_E), n_L = as.integer(s$n_L),
    n_P = as.integer(s$n_P), n_A = as.integer(s$n_A),
    n_states = as.integer(s$n_states), n_days = as.integer(f$n_days),
    n_counters = length(s$counter_names),
    # 0-based start of each state block, in the order E, L, P, M, U, F
    offsets = as.integer(lay$start[seq_len(6L)] - 1L),
    # host-wound model (D-030): flag, infestation sub-stages, 0-based starts of
    # HS, HW, HI (-1 without hosts), host count per node and per-host rates
    host_model = as.integer(!is.null(h)),
    n_I = as.integer(if (is.null(h)) 0L else h$n_I),
    host_offsets = c(host_start("HS"), host_start("HW"), host_start("HI")),
    N_H = if (is.null(h)) 0 else h$n,
    host_wounding = if (is.null(h)) 0 else h$wounding,
    host_healing = if (is.null(h)) 0 else h$healing,
    host_advance = if (is.null(h)) 0 else h$advance,
    host_treatment = if (is.null(h)) 0 else h$treatment,
    host_attractiveness = if (is.null(h)) 0 else h$attractiveness,
    # options
    mating_instant = as.integer(p$mating_mode == "instant"),
    search_all_adults = as.integer(p$wound_competitors == "all_adults"),
    larval_per_wound = as.integer(p$larval_crowding_mode == "per_wound"),
    larval_per_infested_host = as.integer(p$larval_crowding_mode == "per_infested_host"),
    oviposition_batches = as.integer(p$oviposition_mode == "batches"),
    # scalar parameters (counts per node where relevant)
    X_E = p$egg_duration_dd, tau_L = p$larva_duration_days,
    X_P = p$pupa_duration_dd,
    R = p$eggs_per_day_max,
    a_s = p$wound_search_rate,
    sigma0 = p$wound_finding_max, nu = p$mating_rate,
    mu_L0 = -log(p$larva_survival_daily),
    gamma_L = p$larval_crowding_coef,
    mu_W = p$pupa_drought_mortality, mu_b_F = p$female_background_mortality,
    mu_b_M = p$male_background_mortality,
    H = if (is.null(h)) s$counts$H else 0, h_M = s$counts$h_M,
    h_M_param = p$male_density_half_mating,
    clutch_size = if (is.na(p$clutch_size)) 0 else p$clutch_size,
    # adult chains: [substage, genotype], column-major
    inv_w_f = inv_width(s$w_f), inv_w_m = inv_width(s$w_m),
    last_f = as.integer(s$last_f), last_m = as.integer(s$last_m),
    laying = as.integer(s$laying), receptive = as.integer(s$receptive),
    mature = as.integer(s$mature),
    # cube: egg tensor [(g, m) with m fastest, k]; eta [g, m]; per genotype
    egg_tensor = as.numeric(s$egg_tensor), eta = as.numeric(s$eta),
    phi = as.numeric(s$phi), xiF = as.numeric(s$xiF), xiM = as.numeric(s$xiM),
    wild = as.integer(s$geno %in% model$cube$wild_type),
    # daily forcing
    dx = as.numeric(f$dx), dx_egg = as.numeric(f$dx_egg), mu_T = as.numeric(f$mu_T),
    mu_E = as.numeric(f$mu_E), phi_ov = as.numeric(f$phi_ov),
    moisture = as.numeric(f$moisture),
    # releases: day (1-based simulation day), state index (0-based), count
    rel_day = as.integer(rel$day), rel_index = as.integer(rel$index - 1L),
    rel_count = as.numeric(rel$count)
  )
  structure(view, class = "sg_model_view", counter_names = s$counter_names)
}
