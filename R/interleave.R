#' Interleave values from two vectors
#'
#' Alternatingly interleaves values from two vectors of equal length,
#' returning the result in a vector of double the length.
#'
#' @param x First vector, will fill the odd elements of the result.
#'  May be of any type, e.g. numeric, character, logical, etc.
#' @param y Second vector, will fill the even elements of the result.
#'  Must match the type and length of \code{x}
#'
#' @return A vector of x and y interleaved.
#'
#' @export
#'
#' @examples
#' interleave(x=seq(from=1, by=2, length.out=5), y=seq(from=2, by=2, length.out=5))
#' interleave(x=rep(TRUE, 3), y=rep(FALSE, 3))

interleave <- function(x, y){
    if (mode(x) != mode(y)) stop("x and y must have the same mode")
    if (length(x) != length(y)) stop("x and y must have the same length")
    ans <- vector(mode=mode(x), length=length(x) + length(y))
    ans[seq(1, length(ans), 2)] <- x
    ans[seq(2, length(ans), 2)] <- y
    return(ans)
}
