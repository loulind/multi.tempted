#' @title multTEMPTED all
#' @description
#' Formats and performs full analysis from input dataset
#'
#' @param
#'
#' @returns
#' @export
#'
#' @examples
#'
multiTEMPTED <- function(featuretables, timepoints, subjID) {

  datlists <- format_multitempted(featuretables, timepoints, subjID)
  decomp <- multi_tempted_decomp(datlists)
  # odds ratio thing?
  # outputs list containing all pertinent output

}
