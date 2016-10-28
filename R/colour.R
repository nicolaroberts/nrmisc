#' Colour palettes
#'
#' @param palette_name One of \code{"colblind"} (max 9), \code{"bigpal"} (max 29)
#' @param n Number of colours to return
#'
#' @return A vector of n colours
#' @importFrom grDevices rgb
#' @seealso \code{\link{col_show}}, \code{\link[grDevices]{colorRamp}}
#' @export
#'
#' @examples
#' mycol <- nr_col_pal("colblind", 9)
#' col_show(mycol)

nr_col_pal <- function(palette_name, n){

    nr_palettes <- list(
        colblind = rgb(c( 86, 204, 230,   0, 160, 240,   0, 213, 0),
                       c(180, 121, 159, 158, 160, 228, 114,  94, 0),
                       c(233, 167,   0, 115, 160,  66, 178,   0, 0),
                         maxColorValue=255),

        bigpal = c("gray10", "grey70",
                   "midnightblue", "cornflowerblue",
                   "purple3", "mediumpurple1",
                   "violetred2", "plum",
                   "darkorange2", "lightgoldenrod2",
                   "darkgreen", "olivedrab3",
                   "dodgerblue1", "skyblue2",
                   "magenta4", "lightpink3",
                   "firebrick", "gold1",
                   "grey40", "aquamarine3",
                   "darkslategrey", "cyan2",
                   "slateblue3", "indianred3",
                   "darkred", "mediumorchid2",
                   "chocolate4", "darkgoldenrod2",
                   "darkkhaki")
    )

    if (! palette_name %in% names(nr_palettes)){
        stop("palette_name is not valid. See ?nr_col_pal")
    }

    pal <- nr_palettes[[palette_name]]

    if (n < 1 | n %% 1 != 0) stop("n must be a postive integer")

    if(n > length(pal)) {
        stop("n is larger than the number of colours available in this palette.
             See ?nr_col_pick")
    }

    result <- pal[1:n]
    return(result)
}

#' Plot colour palette
#'
#' @param cols Vector of colours
#' @seealso \code{\link{nr_col_pal}}
#' @importFrom graphics plot rect
#' @export
#'
#' @examples
#' mycol <- nr_col_pal("colblind", 9)
#' col_show(mycol)
col_show <- function(cols){
    ncolor <- length(cols)
    side <- ceiling(sqrt(ncolor))
    plot(NULL, xlim=c(0,side), ylim=c(0,side), xlab="", ylab="",
         xaxt="n", yaxt="n", bty="n")
    rect(xleft=rep(seq(0, side - 1), length=ncolor),
         ybottom=rep(seq(side - 1, 0), each=side, length=ncolor),
         xright=rep(1:side, length=ncolor),
         ytop=rep(side:1, each=side, length=ncolor),
         col=cols)
}
