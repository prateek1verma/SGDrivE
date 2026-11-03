#' Read and validate an MGDrivE inheritance cube
#'
#' Checks an inheritance cube from the 'MGDrivE' package (for example
#' `MGDrivE::cubeMendelian()`) or a hand-built cube with the same structure,
#' converts its adult lifespan modifiers to lifespan fractions, and
#' precomputes the quantities the solvers need. The cube itself is used
#' unchanged.
#'
#' @param cube A list with the elements of an MGDrivE v1 cube:
#'   \describe{
#'     \item{`genotypesID`}{Character vector of the n genotype names.}
#'     \item{`ih`}{Array n x n x n: `ih[g, m, k]` is the fraction of offspring
#'       of genotype k from a female of genotype g mated to a male of
#'       genotype m. Each `[g, m, ]` slice sums to 1, or to 0 for crosses that
#'       cannot happen (e.g. "female" genotypes that are male in sex-linked
#'       cubes).}
#'     \item{`tau`}{Array n x n x n of offspring viability modifiers, >= 0.}
#'     \item{`eta`}{n x n matrix of mating fitness, `eta[g, m]` = preference
#'       of a female of genotype g for males of genotype m (a length-n vector
#'       over male genotypes is also accepted).}
#'     \item{`phi`}{Fraction of emerging adults that are female, per genotype.}
#'     \item{`xiF`, `xiM`}{Emergence success of females and males, per
#'       genotype, in \[0, 1\].}
#'     \item{`s`}{Fecundity multiplier of females, per genotype, >= 0. It
#'       scales the egg batches a female lays and, with them, her share of
#'       the wound search (D-031): `s = 0` for females that produce no eggs,
#'       such as irradiated sterile females, which then do not compete for
#'       wounds.}
#'     \item{`omega`}{MGDrivE v1 adult lifespan modifier, per genotype, >= 0
#'       (daily survival is multiplied by omega).}
#'     \item{`wildType`}{Name(s) of the wild-type genotype(s).}
#'     \item{`releaseType`}{Optional: name(s) of the release genotype(s).}
#'   }
#'   Per-genotype vectors may be named (they are reordered to `genotypesID`)
#'   or unnamed (assumed to be in `genotypesID` order). Other list elements
#'   are ignored.
#' @param reference_mortality Daily adult mortality (probability, in (0, 1)) that was used
#'   when the cube's `omega` values were computed with
#'   `MGDrivE::calcOmega(mu, L)`. Required when any `omega` differs from 1
#'   and `lifespan_fraction` is not given (D-012).
#' @param lifespan_fraction Optional lifespan fractions per genotype (> 0; 1 = wild-type
#'   lifespan, 0.5 = half), named or in `genotypesID` order. When given, the
#'   cube's `omega` is not used.
#'
#' @details
#' **Lifespan fractions.** ScrewDrivE represents a genotype's effect on adult
#' lifespan as a fraction l_g of the wild-type physiological lifespan `adult_lifespan_dd`
#' (D-002; model specification S5). Temperature mortality is the same for all
#' genotypes. MGDrivE v1 cubes instead store omega, a multiplier of daily
#' survival. For a cube built with `MGDrivE::calcOmega(mu, L)` (mu =
#' `reference_mortality`) the exact
#' inverse is
#' \deqn{l_g = mu / (1 - omega_g (1 - mu)),}
#' which recovers L. `omega = 0` (immediate death) maps to l_g = 0.001. See
#' [sg_omega_to_ell()].
#'
#' **Egg tensor.** For egg production (specification S9) the cube is reduced
#' to `egg_tensor[(g, m), k] = s_g * ih[g, m, k] * tau[g, m, k]`, an
#' n^2 x n matrix whose rows run over crosses with the female genotype g
#' varying slowest and the mate genotype m fastest (row `(g - 1) * n + m`),
#' matching the state layout of mated females. Its row sums are
#' `s_g * batch_viability`: a female of genotype g lays `s_g` times the
#' batches of a wild-type female (and searches for wounds in proportion),
#' each batch with `batch_viability[g|m]` viable eggs per egg (D-031).
#'
#' **Checks.** Dimensions and genotype names; `ih` slices summing to 1 or 0;
#' non-negative, finite entries; `phi`, `xiF`, `xiM` in \[0, 1\]; wild-type and
#' release genotypes present. A warning is given for `tau > 1` and for crosses
#' between a female-capable and a male-capable genotype that produce no
#' offspring.
#'
#' @return An object of class `sg_cube`: a list with `genotypes`,
#'   `n_genotypes`, `wild_type`, `release_type`, `ih`, `tau`, `eta` (n x n),
#'   per-genotype vectors `phi`, `xiF`, `xiM`, `s`, `omega` and `ell`,
#'   `ell_source` (`"none"`, `"omega"` or `"user"`), `mu_ref` (the reference mortality),
#'   `egg_tensor`
#'   (n^2 x n), `egg_total` (total viable-egg factor per cross) and
#'   `batch_viability` (viable eggs per egg laid in a batch, per cross: the
#'   egg total without the female's fecundity `s`, D-031).
#' @seealso [sg_omega_to_ell()]
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   cube <- sg_cube(MGDrivE::cubeMendelian())
#'   cube
#'
#'   # A cube with a lifespan cost built with calcOmega(0.055, 0.5)
#'   om <- MGDrivE::calcOmega(mu = 0.055, lifespanReduction = 0.5)
#'   cube2 <- sg_cube(MGDrivE::cubeMendelian(omega = c(aa = om)), reference_mortality = 0.055)
#'   cube2$ell
#' }
sg_cube <- function(cube, reference_mortality = NULL, lifespan_fraction = NULL) {
  if (!is.list(cube)) {
    .sg_stop("`cube` must be a list (an MGDrivE inheritance cube).")
  }
  required <- c("genotypesID", "ih", "tau", "eta", "phi", "xiF", "xiM", "s", "omega", "wildType")
  missing_el <- setdiff(required, names(cube))
  if (length(missing_el) > 0L) {
    .sg_stop(
      "The cube is missing element(s): ", paste(missing_el, collapse = ", "),
      ". MGDrivE cube functions (e.g. MGDrivE::cubeMendelian()) return all of them."
    )
  }

  geno <- .sg_cube_genotypes(cube)
  n <- length(geno)
  ih <- .sg_cube_array(cube$ih, "ih", geno)
  tau <- .sg_cube_array(cube$tau, "tau", geno)
  eta <- .sg_cube_eta(cube$eta, geno)
  phi <- .sg_cube_vector(cube$phi, "phi", geno, upper = 1)
  xi_f <- .sg_cube_vector(cube$xiF, "xiF", geno, upper = 1)
  xi_m <- .sg_cube_vector(cube$xiM, "xiM", geno, upper = 1)
  fec <- .sg_cube_vector(cube$s, "s", geno)
  omega <- .sg_cube_vector(cube$omega, "omega", geno)

  .sg_cube_check_ih(ih, phi, geno)
  if (any(tau > 1)) {
    .sg_warn(
      "Some tau (offspring viability) values exceed 1 (max ", .sg_fmt(max(tau)),
      "); viability modifiers are usually in [0, 1]."
    )
  }

  wild <- .sg_cube_types(cube$wildType, "wildType", geno, allow_empty = FALSE)
  rel <- .sg_cube_types(cube$releaseType, "releaseType", geno, allow_empty = TRUE)

  lifespan <- .sg_cube_ell(omega, reference_mortality, lifespan_fraction, geno)

  egg <- .sg_egg_tensor(ih, tau, fec, geno)

  structure(
    list(
      genotypes = geno,
      n_genotypes = n,
      wild_type = wild,
      release_type = rel,
      ih = ih,
      tau = tau,
      eta = eta,
      phi = phi,
      xiF = xi_f,
      xiM = xi_m,
      s = fec,
      omega = omega,
      ell = lifespan$ell,
      ell_source = lifespan$source,
      mu_ref = lifespan$mu_ref,
      egg_tensor = egg,
      egg_total = stats::setNames(rowSums(egg), rownames(egg)),
      batch_viability = .sg_batch_viability(ih, tau, geno)
    ),
    class = "sg_cube"
  )
}

#' Convert MGDrivE v1 lifespan modifiers to lifespan fractions
#'
#' Inverts `MGDrivE::calcOmega()`: given the omega values of a cube and the
#' daily adult mortality `reference_mortality` used to compute them, returns the fraction
#' of the wild-type adult lifespan retained by each genotype.
#'
#' @param omega Numeric vector of MGDrivE v1 omega values (>= 0; daily
#'   survival multipliers, so `(1 - reference_mortality) * omega` must be below 1).
#' @param reference_mortality Daily adult mortality probability in (0, 1) used when the
#'   omega values were computed.
#' @return Numeric vector of lifespan fractions l = mu /
#'   (1 - omega (1 - mu)) with mu = `reference_mortality`, with omega = 0 mapped to 0.001
#'   (D-002). Names
#'   are kept.
#' @seealso [sg_cube()]
#' @export
#' @examples
#' sg_omega_to_ell(c(1, 0.9418, 0), reference_mortality = 0.055)
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   # round trip: calcOmega(mu, L) -> L
#'   om <- MGDrivE::calcOmega(mu = 0.055, lifespanReduction = 0.5)
#'   sg_omega_to_ell(om, reference_mortality = 0.055)
#' }
sg_omega_to_ell <- function(omega, reference_mortality) {
  mu_ok <- is.numeric(reference_mortality) && length(reference_mortality) == 1L &&
    is.finite(reference_mortality) && reference_mortality > 0 && reference_mortality < 1
  if (!mu_ok) {
    .sg_stop(
      "`reference_mortality` must be a single daily mortality probability in (0, 1); got ",
      .sg_fmt(reference_mortality), "."
    )
  }
  if (!is.numeric(omega) || any(!is.finite(omega)) || any(omega < 0)) {
    .sg_stop("`omega` must be finite and >= 0.")
  }
  surv <- omega * (1 - reference_mortality)
  if (any(surv >= 1)) {
    bad <- which(surv >= 1)
    .sg_stop(
      "omega = ", .sg_fmt(omega[bad[1L]]), " with reference_mortality = ",
      .sg_fmt(reference_mortality),
      " gives daily survival (1 - reference_mortality) * omega >= 1, i.e. no adult ever dies. ",
      "Check reference_mortality: it must be the mortality used to compute omega."
    )
  }
  # 1 - omega (1 - mu_ref) written as mu_ref + (1 - mu_ref)(1 - omega): the same
  # value, but exact (l = 1) for omega = 1 in floating point.
  ell <- reference_mortality / (reference_mortality + (1 - reference_mortality) * (1 - omega))
  ell[omega == 0] <- 0.001
  ell
}

#' @export
print.sg_cube <- function(x, ...) {
  n <- x$n_genotypes
  cat("<sg_cube> ", n, " genotype", if (n != 1L) "s", "\n", sep = "")
  cat("  genotypes: ", .sg_wrap_names(x$genotypes), "\n", sep = "")
  cat("  wild type: ", paste(x$wild_type, collapse = ", "), "\n", sep = "")
  cat(
    "  release:   ",
    if (length(x$release_type)) paste(x$release_type, collapse = ", ") else "(none)",
    "\n",
    sep = ""
  )
  ell_label <- switch(x$ell_source,
    none = "all 1 (omega = 1)",
    omega = paste0("from omega with mu_ref = ", .sg_fmt(x$mu_ref)),
    user = "given by user"
  )
  cat("  lifespan fraction l_g: ", ell_label, "\n", sep = "")
  diffs <- list(
    phi = x$phi[x$phi != 0.5], xiF = x$xiF[x$xiF != 1], xiM = x$xiM[x$xiM != 1],
    s = x$s[x$s != 1], ell = x$ell[x$ell != 1]
  )
  diffs <- diffs[lengths(diffs) > 0L]
  if (length(diffs) > 0L) {
    cat("  non-default per-genotype values:\n")
    for (nm in names(diffs)) {
      v <- diffs[[nm]]
      shown <- utils::head(v, 6L)
      cat(sprintf(
        "    %-4s %s%s\n", nm,
        paste0(names(shown), "=", format(shown, digits = 4L), collapse = "  "),
        if (length(v) > 6L) paste0("  ... (", length(v), " genotypes)") else ""
      ))
    }
  }
  cat("  eta != 1: ", sum(x$eta != 1), " of ", n * n, " female x male pairs\n", sep = "")
  cat("  tau != 1: ", sum(x$tau != 1), " of ", n^3, " entries\n", sep = "")
  invisible(x)
}

# ---- internal helpers ---------------------------------------------------------

.sg_wrap_names <- function(x, max_show = 12L) {
  if (length(x) <= max_show) {
    return(paste(x, collapse = " "))
  }
  paste0(paste(x[seq_len(max_show)], collapse = " "), " ... (", length(x), " total)")
}

.sg_cube_genotypes <- function(cube) {
  geno <- cube$genotypesID
  if (!is.character(geno) || length(geno) == 0L || anyNA(geno) || any(!nzchar(geno))) {
    .sg_stop("cube$genotypesID must be a non-empty character vector of genotype names.")
  }
  if (anyDuplicated(geno) > 0L) {
    .sg_stop(
      "cube$genotypesID has duplicated names: ",
      paste(unique(geno[duplicated(geno)]), collapse = ", "), "."
    )
  }
  if (!is.null(cube$genotypesN) && !isTRUE(all.equal(as.numeric(cube$genotypesN), length(geno)))) {
    .sg_stop(
      "cube$genotypesN (", .sg_fmt(cube$genotypesN), ") does not match the number of ",
      "genotypes in cube$genotypesID (", length(geno), ")."
    )
  }
  geno
}

# n x n x n array with matching dimnames, non-negative and finite.
.sg_cube_array <- function(a, what, geno) {
  n <- length(geno)
  if (!is.numeric(a) || !identical(as.integer(dim(a)), rep(n, 3L))) {
    .sg_stop(
      "cube$", what, " must be a numeric array of dimension ", n, " x ", n, " x ", n,
      " (genotypes x genotypes x genotypes); got ",
      if (is.null(dim(a))) "no dimensions" else paste(dim(a), collapse = " x "), "."
    )
  }
  dn <- dimnames(a)
  if (!is.null(dn)) {
    for (d in 1:3) {
      if (!is.null(dn[[d]]) && !identical(as.character(dn[[d]]), geno)) {
        .sg_stop(
          "Dimension ", d, " names of cube$", what, " do not match cube$genotypesID ",
          "(same names in the same order are required)."
        )
      }
    }
  }
  if (any(!is.finite(a)) || any(a < 0)) {
    .sg_stop("cube$", what, " must contain finite values >= 0.")
  }
  a <- array(as.numeric(a), dim = rep(n, 3L), dimnames = list(geno, geno, geno))
  a
}

.sg_cube_eta <- function(eta, geno) {
  n <- length(geno)
  if (is.null(dim(eta)) && is.numeric(eta) && length(eta) == n) {
    eta <- .sg_cube_vector(eta, "eta", geno)
    eta <- matrix(eta, n, n, byrow = TRUE)
  }
  if (!is.numeric(eta) || !identical(as.integer(dim(eta)), c(n, n))) {
    .sg_stop(
      "cube$eta must be an ", n, " x ", n, " numeric matrix (female x male genotype) ",
      "or a length-", n, " vector over male genotypes."
    )
  }
  dn <- dimnames(eta)
  if (!is.null(dn)) {
    for (d in 1:2) {
      if (!is.null(dn[[d]]) && !identical(as.character(dn[[d]]), geno)) {
        .sg_stop("Dimension ", d, " names of cube$eta do not match cube$genotypesID.")
      }
    }
  }
  if (any(!is.finite(eta)) || any(eta < 0)) {
    .sg_stop("cube$eta must contain finite values >= 0.")
  }
  matrix(as.numeric(eta), n, n, dimnames = list(geno, geno))
}

# Per-genotype vector: length n, named (reordered) or unnamed (genotypesID order).
.sg_cube_vector <- function(v, what, geno, upper = Inf) {
  n <- length(geno)
  if (!is.numeric(v) || length(v) != n) {
    .sg_stop(
      "cube$", what, " must be a numeric vector with one value per genotype (", n,
      "); got length ", length(v), "."
    )
  }
  nm <- names(v)
  if (!is.null(nm) && all(nzchar(nm))) {
    if (!setequal(nm, geno) || anyDuplicated(nm) > 0L) {
      .sg_stop("Names of cube$", what, " do not match cube$genotypesID.")
    }
    v <- v[geno]
  }
  v <- stats::setNames(as.numeric(v), geno)
  if (any(!is.finite(v)) || any(v < 0) || any(v > upper)) {
    .sg_stop(
      "cube$", what, " must contain finite values in [0, ",
      if (is.finite(upper)) .sg_fmt(upper) else "Inf", "]."
    )
  }
  v
}

.sg_cube_types <- function(x, what, geno, allow_empty) {
  if (is.null(x) || length(x) == 0L) {
    if (allow_empty) {
      return(character(0))
    }
    .sg_stop("cube$", what, " must name at least one genotype.")
  }
  x <- as.character(x)
  bad <- setdiff(x, geno)
  if (length(bad) > 0L) {
    .sg_stop(
      "cube$", what, " contains genotype(s) not in cube$genotypesID: ",
      paste(bad, collapse = ", "), "."
    )
  }
  unique(x)
}

# Each [g, m, ] slice sums to 1 (possible cross) or 0 (impossible cross).
.sg_cube_check_ih <- function(ih, phi, geno) {
  tol <- 1e-8
  sums <- apply(ih, c(1L, 2L), sum)
  bad <- which(abs(sums - 1) > tol & abs(sums) > tol, arr.ind = TRUE)
  if (nrow(bad) > 0L) {
    g <- geno[bad[1L, 1L]]
    m <- geno[bad[1L, 2L]]
    .sg_stop(
      "cube$ih slices must sum to 1 over offspring genotypes (or to 0 for crosses ",
      "that cannot happen). ", nrow(bad), " slice(s) do not; e.g. ih[\"", g, "\", \"",
      m, "\", ] sums to ", .sg_fmt(sums[bad[1L, 1L], bad[1L, 2L]]), "."
    )
  }
  # A cross is possible when g can be female (phi > 0) and m can be male (phi < 1).
  possible <- outer(phi > 0, phi < 1)
  dead <- which(possible & abs(sums) <= tol, arr.ind = TRUE)
  if (nrow(dead) > 0L) {
    .sg_warn(
      nrow(dead), " possible cross(es) produce no offspring (ih slice sums to 0), e.g. ",
      "female ", geno[dead[1L, 1L]], " x male ", geno[dead[1L, 2L]], "."
    )
  }
  invisible(TRUE)
}

.sg_cube_ell <- function(omega, mu_ref, ell, geno) {
  if (!is.null(ell)) {
    ell <- .sg_cube_vector(ell, "lifespan_fraction", geno)
    names(ell) <- geno
    if (any(ell <= 0)) {
      .sg_stop("`lifespan_fraction` must be > 0.")
    }
    if (any(omega != 1)) {
      .sg_warn("`lifespan_fraction` was given, so the cube's omega values are not used.")
    }
    return(list(ell = ell, source = "user", mu_ref = if (is.null(mu_ref)) NA_real_ else mu_ref))
  }
  if (all(omega == 1)) {
    ones <- stats::setNames(rep(1, length(geno)), geno)
    return(list(ell = ones, source = "none", mu_ref = NA_real_))
  }
  if (is.null(mu_ref)) {
    .sg_stop(
      "The cube has omega values different from 1 (",
      paste(names(omega)[omega != 1], collapse = ", "), "). To convert them to lifespan ",
      "fractions, give `reference_mortality`, the daily adult mortality used when they were ",
      "computed with MGDrivE::calcOmega(), or give `lifespan_fraction` directly."
    )
  }
  list(ell = sg_omega_to_ell(omega, mu_ref), source = "omega", mu_ref = mu_ref)
}

# Viable eggs per egg in a batch for each cross, sum_k ih[g, m, k] *
# tau[g, m, k], named "g|m" with m fastest (D-031). The female's fecundity
# s_g sets how many batches she lays (and so how much she searches), not
# their viability.
.sg_batch_viability <- function(ih, tau, geno) {
  egg <- .sg_egg_tensor(ih, tau, rep(1, length(geno)), geno)
  stats::setNames(rowSums(egg), rownames(egg))
}

# egg_tensor[(g, m), k] = s_g * ih[g, m, k] * tau[g, m, k]; rows g-major.
.sg_egg_tensor <- function(ih, tau, fec, geno) {
  n <- length(geno)
  t3 <- ih * tau
  t3 <- sweep(t3, 1L, fec, `*`)
  # aperm puts the mate genotype first so that, after flattening column-major,
  # row index = m + (g - 1) * n (female genotype slowest).
  mat <- matrix(aperm(t3, c(2L, 1L, 3L)), n * n, n)
  rownames(mat) <- paste(rep(geno, each = n), rep(geno, times = n), sep = "|")
  colnames(mat) <- geno
  mat
}
