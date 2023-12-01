#' Search and filter datasets that align with the specified options.
#'
#' @name search_data
#' @aliases search_data
#' @description Show the different types of data sets in this package.
#' @usage search_data(...)
#' @param ... task, type, n, p, tag
#' @returns Return a list with the the specified options of dataset.
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @examples
#' search_data(n > 100, p < 10)
#' @export

search_data <- function(...) {
  data <- dataSDA::dataset
  filtered_data <- data %>% filter(...)
  return(filtered_data)
}

