#' @param k string of field (data type) to return headers for
#' 
#' @details
#' If \code{get_header} is called with no values, the entire available header list is returned.
#' @examples
#' \dontrun{
#' get_header("manifest")
#' }
#' @rdname headers
#' @export
#' 
get_header <- function(k = NULL) {
  checkmate::assert_string(k, null.ok = TRUE)
  avail_headers <- names(headers)
  if (is.null(k)) {
    headers
  } else if (k %in% avail_headers) {
    headers[[k]]
  } else {
    stop(sprintf("Unknown header. k must be NULL or %s", avail_headers))
  }
}


# list of headers
headers <- list(response_metrics = c(
  "x_mean",
  "x_AOC",
  "x_AOC_range",
  "xc50",
  "x_max",
  "ec50",
  "x_inf",
  "x_0",
  "h",
  "r2",
  "x_sd_avg",
  "fit_type"
))
