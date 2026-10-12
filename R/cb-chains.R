# Exact advancement along a chain of sub-stages over one time step, for the
# chain-binomial solver (D-022, model specification S16). Pure-R version: a
# matrix exponential by scaling and squaring, independent of the
# uniformization used in C++ (src/chains.h), so each checks the other.

#' Advancement probabilities along a chain over a time step
#'
#' @param q Rates of moving on from each class (per day).
#' @param last Number of classes that are reached; from class `last` the
#'   individuals leave the chain, classes after it are never reached.
#' @param dt Time step (days).
#' @return Matrix with `length(q)` rows and `length(q) + 1` columns: entry
#'   `[i, j]` is the probability of being in class `j` at the end of the step
#'   when starting in class `i`, the last column the probability of having
#'   left the chain. Rows after `last` are zero. Death is not included.
#' @noRd
.sg_chain_advance_r <- function(q, last, dt) {
  n <- length(q)
  out <- matrix(0, n, n + 1L)
  if (last < 1L) {
    return(out)
  }
  # generator on the reached classes plus the absorbing exit state
  m <- last + 1L
  gen <- matrix(0, m, m)
  for (i in seq_len(last)) {
    gen[i, i] <- -q[i]
    gen[i, i + 1L] <- q[i]
  }
  p <- .sg_expm(gen * dt)
  out[seq_len(last), seq_len(last)] <- p[seq_len(last), seq_len(last)]
  out[seq_len(last), n + 1L] <- p[seq_len(last), m]
  out
}

# Matrix exponential by scaling and squaring with a Taylor series (norm of
# the scaled matrix at most 1/2, 30 terms).
.sg_expm <- function(a) {
  nrm <- max(rowSums(abs(a)))
  s <- if (nrm > 0.5) ceiling(log2(nrm / 0.5)) else 0L
  a <- a / 2^s
  term <- diag(nrow(a))
  out <- term
  for (k in seq_len(30L)) {
    term <- term %*% a / k
    out <- out + term
  }
  for (i in seq_len(s)) out <- out %*% out
  out
}
