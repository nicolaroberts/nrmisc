#' Refactor factor variables in a data.frame
#'
#'
#' @param df data.frame
#'
#' @return A data.frame with each factor variable re-factored (old unused levels are dropped)
#'
#' @export
#'

refactor <- function(df){
    cat <- sapply(df, is.factor)
    df[cat] <- lapply(df[cat], factor)
    return(df)
}
