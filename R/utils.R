is_numeric0 <- function(x) {
  is.numeric(x) && length(x) == 0
}

`%||%` <- function(x, y) if (is.null(x)) y else x
