# Scenario 22: sensitivity of the Texas transect (scenario 21)
#
# Runs scenario 21 for the transect and its base runs, then changes one thing
# at a time:
#   * the movement clock (design note 5.3, decision 2): calendar (the given
#     rates every day), capped degree-day scaling (the default,
#     movement_kappa_max = 1) and uncapped degree-day scaling
#     (movement_kappa_max = Inf). No study relates screwworm movement to
#     temperature, so the clock is a sensitivity dimension;
#   * the diffusion coefficient, 10 times larger;
#   * the time step of the chain-binomial solver, dt = 0.05 instead of 0.1:
#     newly emerged adults cannot move in the step they emerge in, which
#     slows the movement of unmated females by a share of order dt (D-032
#     phase 2 findings); does it change the front?
#   * the ODE without the node extinction rule: fractions of a fly ("atto-
#     flies") survive the winter in every cell and grow again in spring, so
#     cells come back without being recolonised (design note 5.7).
#
# Data as in scenario 21 (run from the repository root). Run time: a few
# minutes.

source(system.file("scenarios", "21_texas_transect.R", package = "ScrewDrivE"))

# 1. Variants ----------------------------------------------------------------------
# The functions below use objects made by scenario 21 (nodes, h, n_cells,
# d_coef, cube, land, all_nodes, transect_summary), which the linter cannot
# see through source().
# nolint start: object_usage_linter.
variant_model <- function(clock = "physiological", kappa_max = 1, d_scale = 1) {
  p_v <- sg_params(movement_kappa_max = kappa_max)
  layer_v <- sg_movement_diffusion(
    stats::setNames(nodes$area, nodes$id),
    data.frame(from = nodes$id[-nrow(nodes)], to = nodes$id[-1L], border = h, dist = h),
    D = lapply(d_coef, function(d) d * d_scale), outer_border = c(rep(0, n_cells), h), clock = clock
  )
  sg_model(p_v, cube, landscape = land, movement = layer_v)
}
variants <- list(
  "capped clock (default)" = list(),
  "calendar clock" = list(clock = "calendar"),
  "uncapped clock" = list(kappa_max = Inf),
  "D x 10" = list(d_scale = 10)
)

# The variants have the same states as the base model; its start is used as
# a plain vector (an sg_init() object belongs to the model it was built for).
x0_states <- as.numeric(x0)

summarise_runs <- function(label, m, dt = 0.1, n_rep = 10, seed = 22) {
  ode_v <- simulate_ode(m, x0_states, node_extinction = TRUE, rtol = 1e-6, atol = 1e-6)
  cb_v <- simulate_cb(m, x0_states, n_rep = n_rep, dt = dt, seed = seed)
  sp_o <- sg_spread(ode_v, nodes = all_nodes)
  sp_c <- sg_spread(cb_v, nodes = all_nodes)
  cb_rows <- do.call(rbind, lapply(seq_len(n_rep), function(r) {
    transect_summary(sp_c$time[sp_c$time$replicate == r, ])
  }))
  list(
    table = data.frame(
      variant = label, dt = dt,
      ode_front_winter_min_km = transect_summary(sp_o$time)$front_winter_min_km,
      cb_front_winter_min_km = mean(cb_rows$front_winter_min_km),
      ode_front_1962_08_01_km = transect_summary(sp_o$time)$front_1962_08_01_km,
      cb_front_1962_08_01_km = mean(cb_rows$front_1962_08_01_km),
      cb_front_speed_km_per_generation = mean(cb_rows$front_speed_km_per_generation),
      cb_front_speed_sd = stats::sd(cb_rows$front_speed_km_per_generation)
    ),
    fronts = rbind(
      data.frame(
        variant = label, date = sp_o$time$date, front = sp_o$time$front, run = "ODE", line = "ODE"
      ),
      data.frame(
        variant = label, date = sp_c$time$date, front = sp_c$time$front, run = "chain-binomial",
        line = paste("cb", sp_c$time$replicate)
      )
    )
  )
}

# nolint end
runs <- c(
  Map(function(label, args) {
    summarise_runs(label, do.call(variant_model, args))
  }, names(variants), variants),
  list(summarise_runs("capped clock, dt = 0.05", variant_model(), dt = 0.05))
)
sens <- do.call(rbind, lapply(runs, `[[`, "table"))
print(sens, row.names = FALSE, digits = 3)
utils::write.csv(sens, file.path(fig_dir, "transect_sensitivity.csv"), row.names = FALSE)

# 2. The ODE without the extinction rule ----------------------------------------
# Counting a cell as occupied from one fly, the winter looks the same; but
# the emptied cells keep fractions of a fly (here: more than 1e-9) and grow
# back in spring wherever the weather allows, without any immigration.
ode_atto <- simulate_ode(model, x0, rtol = 1e-6, atol = 1e-6)
sp_atto <- sg_spread(ode_atto, nodes = all_nodes)
traces <- sg_spread(ode_atto, nodes = all_nodes, threshold = 1e-9)
cat(sprintf(
  paste0(
    "\nODE without the extinction rule: front on 1 August 1962 at %g km (with the rule: %g km); ",
    "cells holding more than 1e-9 flies at the winter low: %d of %d.\n"
  ),
  transect_summary(sp_atto$time)$front_1962_08_01_km,
  transect_summary(sp_ode$time)$front_1962_08_01_km,
  min(traces$time$occupied[traces$time$date >= winter[1L] & traces$time$date <= winter[2L]]),
  nrow(nodes)
))

# 3. Figure ------------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  fronts <- do.call(rbind, lapply(runs, `[[`, "fronts"))
  fronts <- rbind(fronts, data.frame(
    variant = "ODE without extinction rule", date = sp_atto$time$date, front = sp_atto$time$front,
    run = "ODE", line = "ODE"
  ))
  fronts$variant <- factor(fronts$variant, levels = unique(fronts$variant))
  fig <- ggplot2::ggplot(fronts, ggplot2::aes(date, front, group = line, colour = run)) +
    ggplot2::geom_step(alpha = 0.7) +
    ggplot2::facet_wrap(ggplot2::vars(variant), ncol = 2) +
    ggplot2::scale_colour_manual(values = c(ODE = "#000000", `chain-binomial` = "#E69F00")) +
    ggplot2::labs(
      x = NULL, y = "km north of the source node", colour = NULL,
      title = "Northernmost cell with flies: sensitivity"
    ) +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(fig_dir, "fig3_front_sensitivity.png"), fig,
    width = 10, height = 9, dpi = 150
  )
}
