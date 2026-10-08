# Internal helpers for errors, warnings and messages (D-016: base R only).
#
# Messages say what is wrong and how to fix it. `call. = FALSE` keeps them
# readable: the user sees the message, not an internal function call.

.sg_stop <- function(...) {
  stop(paste0(...), call. = FALSE)
}

.sg_warn <- function(...) {
  warning(paste0(...), call. = FALSE)
}

# " Did you mean 'x'?" when a close match exists, otherwise "".
.sg_suggest <- function(x, choices) {
  if (length(choices) == 0L) {
    return("")
  }
  d <- utils::adist(x, choices, ignore.case = TRUE)[1L, ]
  best <- which.min(d)
  if (d[best] <= max(2L, nchar(x) %/% 3L)) {
    paste0(" Did you mean '", choices[best], "'?")
  } else {
    ""
  }
}

# Compact, readable representation of a value for messages and printing.
.sg_fmt <- function(x) {
  if (length(x) != 1L) {
    return(paste0("<length ", length(x), ">"))
  }
  if (is.na(x)) {
    return("NA")
  }
  if (is.character(x)) {
    return(paste0("\"", x, "\""))
  }
  if (is.numeric(x)) {
    return(format(x, digits = 6L))
  }
  as.character(x)
}
