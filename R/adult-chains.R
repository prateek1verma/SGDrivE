# Segmented adult chains (model specification S6, decision D-018).
#
# Each adult chain (females: U and F, which share age classes; males: M) is
# split at its maturity thresholds into segments with exact mean durations.
# Sub-stage counts per segment are the same for all genotypes; genotype g
# with lifespan L_g = ell_g * X_A shortens (or, for ell_g > 1, lengthens) only
# the segment in which L_g falls, and segments starting at or after L_g are
# never reached.

# Allocate n sub-stages to segments in proportion to their lengths, at least
# one each, by largest remainder. Near-integer quotas (floating-point noise)
# are treated as integers.
.sg_allocate_substages <- function(n, lengths) {
  k <- length(lengths)
  if (n < k) {
    .sg_stop(
      "adult_substages = ", n, " is too small for an adult chain with ", k, " segments; ",
      "use adult_substages >= ", k, "."
    )
  }
  quota <- n * lengths / sum(lengths)
  alloc <- pmax(1L, as.integer(floor(quota + 1e-9)))
  while (sum(alloc) < n) {
    j <- which.max(quota - alloc)
    alloc[j] <- alloc[j] + 1L
  }
  while (sum(alloc) > n) {
    cand <- which(alloc > 1L)
    j <- cand[which.max((alloc - quota)[cand])]
    alloc[j] <- alloc[j] - 1L
  }
  alloc
}

# Build one adult chain.
#   thresholds: named threshold ages (dd), e.g. c(X_rec = 0, X_ovi = 61);
#               values of 0 add no boundary.
#   X_A: wild-type physiological lifespan (dd); n_A: sub-stages per sex;
#   ell: named lifespan fractions per genotype; substages: optional counts per
#   segment; sex: "female" or "male" (for messages).
# Returns a list with the segment structure, per-class flags (one logical
# vector per threshold: class segment starts at or after the threshold),
# per-genotype class widths (dd; NA for classes never reached) and the last
# class reached by each genotype.
# nolint start: object_name_linter.
.sg_adult_chain <- function(thresholds, X_A, n_A, ell, substages = NULL, sex) {
  # nolint end
  inner <- sort(unique(thresholds[thresholds > 0]))
  bounds <- c(0, unname(inner), X_A)
  k <- length(bounds) - 1L
  seg_start <- bounds[-(k + 1L)]
  seg_len <- diff(bounds)

  counts <- if (is.null(substages)) {
    .sg_allocate_substages(n_A, seg_len)
  } else {
    .sg_check_substages(substages, k, n_A, sex, bounds)
  }
  class_seg <- rep(seq_len(k), counts)

  flags <- lapply(thresholds, function(x) seg_start[class_seg] >= x)
  names(flags) <- names(thresholds)

  geno <- names(ell)
  width <- matrix(NA_real_, length(geno), n_A, dimnames = list(geno, NULL))
  last <- stats::setNames(integer(length(geno)), geno)
  for (g in geno) {
    life <- ell[[g]] * X_A
    seg_end <- pmin(bounds[-1L], life)
    seg_end[k] <- life
    reached <- seg_start < life
    j_last <- max(which(reached))
    w <- (seg_end - seg_start) / counts
    w[!reached] <- NA_real_
    width[g, ] <- w[class_seg]
    last[[g]] <- sum(counts[seq_len(j_last)])
  }

  list(
    sex = sex,
    thresholds = thresholds,
    boundaries = bounds,
    segment_counts = counts,
    class_segment = class_seg,
    flags = flags,
    width = width,
    last_class = last
  )
}

.sg_check_substages <- function(x, k, n_A, sex, bounds) { # nolint: object_name_linter.
  seg_txt <- paste0(
    "[", format(bounds[-length(bounds)], digits = 4L), ", ",
    format(bounds[-1L], digits = 4L), ") dd",
    collapse = "; "
  )
  ok <- is.numeric(x) && length(x) == k && all(is.finite(x)) &&
    all(abs(x - round(x)) < 1e-8) && all(x >= 1)
  if (!ok) {
    .sg_stop(
      "adult_segment_substages$", sex, " must give ", k, " whole number(s) >= 1, one per ",
      "segment of the ", sex, " chain (", seg_txt, ")."
    )
  }
  x <- as.integer(round(x))
  if (sum(x) != n_A) {
    .sg_stop(
      "adult_segment_substages$", sex, " must sum to adult_substages = ", n_A,
      " (it sums to ", sum(x), "). Change adult_substages with sg_params() if needed."
    )
  }
  x
}
