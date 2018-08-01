#' Collapse consecutive genomic ranges with same meta-data value
#'
#' @param gr A GRanges object with one meta-data column
#' @param na.rm Logical - if TRUE (default), ranges with NA in metadata column are removed from output.
#'  If FALSE, ranges with NA are retained AND gaps between ranges are returned (also with NA)
#' @return A GRanges object where consecutive ranges with the same meta-data
#'  value are collapsed together
#' @import GenomicRanges IRanges GenomeInfoDb S4Vectors
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
    gr_gaps <- gr_gaps[-which(width(gr_gaps) %in% seqlengths(seqinfo(gr)))]

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


#' Sample random positions from a GRanges object
#'
#' @param gr A GRanges object
#' @param n Number of positions to sample
#' @return A GRanges object of n positions uniformly sampled from gr
#' @import GenomicRanges
#' @export
sample_GRanges <- function(gr, n){
    # randomly sample the ranges, weighted by their width
    ri <- sample(length(gr), n, replace=TRUE, prob=width(gr))
    rranges <- gr[ri]
    # randomly sample positions in those ranges
    rand <- sapply(width(rranges), sample, size=1)
    pos <- start(rranges)+rand-1
    # output
    ans <- GRanges(seqnames(rranges), IRanges(pos, width=1),
                   strand='*', seqinfo=seqinfo(gr), mcols(rranges))
    ans <- sort(ans)
    return(ans)
}




#' Plot one quantitative meta-data track from a GRanges object
#'
#' @param gr A GRanges object
#' @param track Name of quantitative meta-data track to plot
#' @param seqname seqname to plot for
#' @param start starting position of plot (start of seqname if left NULL)
#' @param end end position of plot (end of seqname if left NULL)
#' @import GenomicRanges ggplot2
#' @export
plot_gr_track <- function(gr, track, seqname, start=NULL, end=NULL){

    # set start and end if null
    if (is.null(start)) start <- start(range(gr[seqnames(gr)==seqname]))
    if (is.null(end)) end <- end(range(gr[seqnames(gr)==seqname]))

    # subset gr
    gr <- subsetByOverlaps(gr, GRanges(seqname, IRanges(start, end)))

    # get info needed
    mdata <- as.data.frame(gr)
    pdata <- data.frame(seqnames=rep(mdata$seqnames, each=2),
                        pos=interleave(mdata$start, mdata$end),
                        value=rep(mdata[[track]], each=2))
    # plot basics
    p <- ggplot(data=pdata, aes_string(x='pos', y='value')) +
        theme(panel.grid=element_blank(), text=element_text(size=16)) +
        scale_x_continuous(name=paste('position on', seqname)) +
        scale_y_continuous(name=track)

    # plot line if number of ranges less than 2000, otherwise stat_bin2d
    if (length(gr) <= 2000){
        ans <- p + geom_line()
    } else {
        ans <- p + stat_bin2d(binwidth=c(width(range(gr))/1000,
                                         diff(range(pdata$value))/100)) +
            scale_fill_continuous(guide=FALSE)
    }

    return(ans)
}
