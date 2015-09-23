#' Collapse consecutive genomic ranges with same meta-data value
#'
#' @param gr A GRanges object with one meta-data column
#' @param na.rm Logical - if TRUE (default), ranges with NA in metadata column are removed from output.
#'  If FALSE, ranges with NA are retained AND gaps between ranges are returned (also with NA)
#' @return A GRanges object where consecutive ranges with the same meta-data
#'  value are collapsed together
#' @import GenomicRanges IRanges
#' @export
#' @examples
#' library(GenomicRanges)
#' my_gr <- GRanges(c("chr1", "chr1", "chr1", "chr1", "chr2"),
#'                  IRanges(start=c(10, 20, 40, 50, 60),
#'                          end=c(19, 29, 49, 59, 69)),
#'                  value=c(1, 2, 2, 2, 2))
#' collapse_granges_by_val(my_gr)

collapse_granges_by_val <- function(gr, na.rm=TRUE){
    if (class(gr) != "GRanges") stop("gr must be of class GRanges")
    if (ncol(mcols(gr)) != 1) {
        stop("gr must have (only) 1 metadata column")
    }

    colname <- names(mcols(gr))

    gr_gaps <- gaps(gr)
    mcols(gr_gaps)[[colname]] <- NA

    gr <- c(gr_gaps, gr)

    gr <- sort(gr)

    grlist <- split(gr, as.vector(seqnames(gr)))

    collapse <- function(x) {
        value_rle <- S4Vectors::Rle(mcols(x)[[colname]])

        runends <- cumsum(S4Vectors::runLength(value_rle))
        runstarts <- c(1, runends[-length(runends)] + 1)

        iranges <- ranges(x)

        ans <- GRanges(seqnames(x)[1],
                        IRanges(start(iranges)[runstarts],
                                end(iranges)[runends]))

        mcols(ans)[[colname]] <- S4Vectors::runValue(value_rle)

        return(ans)
    }

    ans <- unlist(GRangesList(lapply(grlist, collapse)))
    seqinfo(ans) <- seqinfo(gr)

    if (na.rm){
        ans <- ans[!is.na(mcols(ans)[[colname]])]
    }

    return(ans)
}
