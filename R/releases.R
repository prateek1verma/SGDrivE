#' Release schedule
#'
#' Defines releases of a genotype at a life stage on a set of days (D-005).
#' Several schedules (different genotypes, stages or sexes) are combined with
#' `c()`. Released individuals enter the first sub-stage of their life stage
#' at the start of the release day (model specification section 9).
#'
#' @param genotype Name of the released genotype; must be one of the cube's
#'   genotypes (checked when the model is assembled).
#' @param stage Life stage released: `"egg"`, `"pupa"` or `"adult"`. Adults
#'   are released unmated: males enter M, females enter U.
#' @param sex For adult releases only: `"male"`, `"female"` or `"both"`.
#'   Must be left `NULL` for eggs and pupae, whose sex is assigned at
#'   emergence by the cube's `phi`.
#' @param times Release days: either day numbers of the simulation (whole
#'   numbers >= 1; day 1 is the first day of the forcing, D-015) or dates
#'   (class `Date` or `"YYYY-MM-DD"` text), converted to days when the model
#'   is assembled.
#' @param number Number released per event: one value for all events or one
#'   per element of `times`; >= 0.
#' @param per Unit of `number`: `"node"` (individuals per node) or `"km2"`
#'   (individuals per square kilometre, multiplied by the node area when the
#'   model is assembled; D-009).
#' @param split For `sex = "both"`: fraction of each release that is female
#'   (default 0.5, D-013).
#'
#' @return An object of class `sg_releases` holding a data frame of events
#'   (one row per event, stage and sex) with columns `time`, `time_type`
#'   (`"day"` or `"date"`), `genotype`, `stage`, `sex` (`NA` for eggs and
#'   pupae), `number` and `per`.
#' @seealso `c()` to combine schedules; `as.data.frame()` to inspect them.
#' @export
#' @examples
#' # 1000 sterile males per km2 every week for 26 weeks, from day 30
#' rel <- sg_releases("S",
#'   stage = "adult", sex = "male",
#'   times = seq(30, by = 7, length.out = 26),
#'   number = 1000, per = "km2"
#' )
#' rel
#'
#' # Releases on dates, males and females (50:50)
#' rel2 <- sg_releases("S",
#'   stage = "adult", sex = "both",
#'   times = as.Date("1971-06-01") + c(0, 14), number = 500
#' )
#'
#' rel2
#'
#' # Combine schedules given in the same way (day numbers here)
#' pupae <- sg_releases("S", stage = "pupa", times = 100, number = 50)
#' c(rel, pupae)
sg_releases <- function(genotype, stage = c("egg", "pupa", "adult"), sex = NULL,
                        times, number, per = c("node", "km2"), split = 0.5) {
  geno_ok <- is.character(genotype) && length(genotype) == 1L && !is.na(genotype) &&
    nzchar(genotype)
  if (!geno_ok) {
    .sg_stop("`genotype` must be a single genotype name, e.g. \"S\".")
  }
  stage <- .sg_match_choice(stage, c("egg", "pupa", "adult"), "stage")
  per <- .sg_match_choice(per, c("node", "km2"), "per")

  if (stage == "adult") {
    if (is.null(sex)) {
      .sg_stop("Adult releases need `sex`: \"male\", \"female\" or \"both\".")
    }
    sex <- .sg_match_choice(sex, c("male", "female", "both"), "sex")
  } else if (!is.null(sex)) {
    .sg_stop(
      "`sex` applies only to adult releases; leave it NULL for ", stage,
      " releases (their sex is assigned at emergence by the cube's phi)."
    )
  }
  split_ok <- is.numeric(split) && length(split) == 1L && is.finite(split) &&
    split >= 0 && split <= 1
  if (!split_ok) {
    .sg_stop("`split` (fraction female) must be a single number in [0, 1].")
  }
  if (!identical(sex, "both") && split != 0.5) {
    .sg_warn("`split` is used only when sex = \"both\"; it is ignored here.")
  }

  tm <- .sg_release_times(times)
  number <- .sg_release_number(number, length(tm$time))

  ev <- data.frame(
    time = tm$time, time_type = tm$type, genotype = genotype, stage = stage,
    sex = if (is.null(sex)) NA_character_ else sex, number = number, per = per,
    stringsAsFactors = FALSE
  )
  if (identical(sex, "both")) {
    males <- ev
    males$sex <- "male"
    males$number <- ev$number * (1 - split)
    females <- ev
    females$sex <- "female"
    females$number <- ev$number * split
    ev <- rbind(males, females)
  }
  .sg_new_releases(ev)
}

#' @export
c.sg_releases <- function(...) {
  parts <- list(...)
  parts <- parts[!vapply(parts, is.null, TRUE)]
  ok <- vapply(parts, inherits, TRUE, what = "sg_releases")
  if (!all(ok)) {
    .sg_stop("Only sg_releases objects can be combined with c().")
  }
  types <- unique(unlist(lapply(parts, function(x) x$events$time_type)))
  if (length(types) > 1L) {
    .sg_stop(
      "Cannot combine releases given as day numbers with releases given as dates; ",
      "use one or the other."
    )
  }
  .sg_new_releases(do.call(rbind, lapply(parts, function(x) x$events)))
}

#' @export
as.data.frame.sg_releases <- function(x, ...) {
  ev <- x$events
  if (identical(ev$time_type[1L], "date")) {
    ev$time <- as.Date(ev$time, origin = "1970-01-01")
  }
  ev
}

#' @export
print.sg_releases <- function(x, ...) {
  ev <- as.data.frame(x)
  n_events <- nrow(unique(ev[, c("time", "genotype", "stage")]))
  cat("<sg_releases> ", n_events, " release event", if (n_events != 1L) "s", "\n", sep = "")
  if (nrow(ev) == 0L) {
    return(invisible(x))
  }
  key <- paste(ev$genotype, ev$stage, ifelse(is.na(ev$sex), "", ev$sex), ev$per, sep = "\r")
  for (k in unique(key)) {
    d <- ev[key == k, , drop = FALSE]
    who <- paste0(
      d$genotype[1], " ", d$stage[1],
      if (!is.na(d$sex[1])) paste0(" (", d$sex[1], ")") else ""
    )
    is_date <- inherits(d$time, "Date")
    when <- if (nrow(d) == 1L) {
      paste("on", if (is_date) format(d$time) else paste("day", d$time))
    } else {
      paste0(
        if (is_date) "" else "days ", format(min(d$time)), " to ", format(max(d$time))
      )
    }
    cat(sprintf(
      "  %-24s %3d %s, %s, total %s per %s\n",
      who, nrow(d), if (nrow(d) == 1L) "event" else "events", when,
      format(sum(d$number), big.mark = ",", scientific = FALSE), d$per[1]
    ))
  }
  invisible(x)
}

# ---- internal ------------------------------------------------------------------

.sg_new_releases <- function(ev) {
  ev <- ev[order(ev$time, ev$genotype, ev$stage, ev$sex, na.last = TRUE), , drop = FALSE]
  rownames(ev) <- NULL
  structure(list(events = ev), class = "sg_releases")
}

.sg_match_choice <- function(x, choices, what) {
  if (identical(x, choices)) {
    return(choices[1L])
  }
  if (!is.character(x) || length(x) != 1L || !x %in% choices) {
    .sg_stop(
      "`", what, "` must be one of ", paste0("\"", choices, "\"", collapse = ", "),
      if (is.character(x) && length(x) == 1L) paste0("; got \"", x, "\".") else ".",
      if (is.character(x) && length(x) == 1L) .sg_suggest(x, choices) else ""
    )
  }
  x
}

# Day numbers (whole, >= 1) or dates, stored as numbers with their type.
.sg_release_times <- function(times) {
  if (missing(times) || length(times) == 0L) {
    .sg_stop("`times` must give at least one release day or date.")
  }
  if (inherits(times, "Date") || is.character(times)) {
    d <- tryCatch(as.Date(times), error = function(e) rep(NA, length(times)))
    if (anyNA(d)) {
      .sg_stop("Release dates in `times` must be valid dates (\"YYYY-MM-DD\").")
    }
    return(list(time = as.numeric(d), type = "date"))
  }
  if (!is.numeric(times) || anyNA(times) || any(!is.finite(times))) {
    .sg_stop("`times` must be day numbers (whole numbers >= 1) or dates.")
  }
  if (any(abs(times - round(times)) > 1e-8) || any(times < 1)) {
    .sg_stop(
      "Release day numbers must be whole numbers >= 1 (releases happen at the ",
      "start of a day; day 1 is the first day of the forcing)."
    )
  }
  list(time = round(as.numeric(times)), type = "day")
}

.sg_release_number <- function(number, n) {
  if (!is.numeric(number) || anyNA(number) || any(!is.finite(number)) || any(number < 0)) {
    .sg_stop("`number` must be finite and >= 0.")
  }
  if (length(number) == 1L) {
    return(rep(as.numeric(number), n))
  }
  if (length(number) != n) {
    .sg_stop(
      "`number` must have length 1 or the same length as `times` (", n,
      "); got length ", length(number), "."
    )
  }
  as.numeric(number)
}

# Resolve releases for a model (used by sg_model()): dates -> simulation days,
# per-km2 numbers -> counts per node, genotypes checked against the cube.
# Returns a data frame with columns day, genotype, stage, sex, count.
.sg_resolve_releases <- function(releases, dates, area, genotypes) {
  empty <- data.frame(
    day = integer(0), genotype = character(0), stage = character(0),
    sex = character(0), count = numeric(0), stringsAsFactors = FALSE
  )
  if (is.null(releases)) {
    return(empty)
  }
  if (!inherits(releases, "sg_releases")) {
    .sg_stop("`releases` must be an sg_releases object (see sg_releases()) or NULL.")
  }
  ev <- releases$events
  if (nrow(ev) == 0L) {
    return(empty)
  }
  bad <- setdiff(unique(ev$genotype), genotypes)
  if (length(bad) > 0L) {
    .sg_stop(
      "Released genotype(s) not in the cube: ", paste(bad, collapse = ", "), ".",
      .sg_suggest(bad[1L], genotypes), " Cube genotypes: ", .sg_wrap_names(genotypes), "."
    )
  }
  n_days <- length(dates)
  day <- if (ev$time_type[1L] == "date") {
    as.integer(ev$time - as.numeric(dates[1L])) + 1L
  } else {
    as.integer(ev$time)
  }
  outside <- day < 1L | day > n_days
  if (any(outside)) {
    i <- which(outside)[1L]
    when <- if (ev$time_type[1L] == "date") {
      format(as.Date(ev$time[i], origin = "1970-01-01"))
    } else {
      paste("day", ev$time[i])
    }
    .sg_stop(
      sum(outside), " release event(s) fall outside the simulation period (",
      format(dates[1L]), " to ", format(dates[n_days]), ", days 1 to ", n_days,
      "), e.g. ", when, "."
    )
  }
  count <- ifelse(ev$per == "km2", ev$number * area, ev$number)
  data.frame(
    day = day, genotype = ev$genotype, stage = ev$stage, sex = ev$sex,
    count = count, stringsAsFactors = FALSE
  )
}
