#' heading
#'
#' @param gpos GRanges object of SV breakpoint positions (each 'range' is a 1bp position)
#' with at least two metadata columns: a) breakpoint_id - a factor of unique ids per breakpoint junction
#' so that the two positions in a junction can be matched together, and b) class - a factor
#' of SV class values that will correspond to plotting colour
#' @param plot_windows GRanges object of plotting windows to return
#' @param class_cols Vector of colours to match the factor levels of gpos$class
#' @param plot_title (Optional) Character vector of plot titles, one for each range in plot_windows
#' @return List of ggplot values, one for each plot window (can be added to afterwards) - WIP
#' @import GenomicRanges ggrepel ggplot2
#' @importFrom plyr ldply
#' @export

plot_sv_rainfall <- function(gpos, plot_windows, class_cols, plot_title=NULL){

    # check inputs
    if (class(gpos) != "GRanges" | class(plot_windows) != "GRanges") {
        stop("gpos and plot_windows must both be of class GRanges")
    }
    if (!"breakpoint_id" %in% colnames(mcols(gpos)) | ! "class" %in% colnames(mcols(gpos))) {
        stop("gpos must have metadata columns breakpoint_id and class")
    }
    if ( ! "factor" %in% class(gpos$breakpoint_id) | max(table(gpos$breakpoint_id)) != 2) {
        stop("gpos$breakpoint_id must be factor with no more than two positions per level")
    }
    if(! "factor" %in% class(gpos$class)) stop("gpos$class must be a factor")
    if(length(class_cols) != nlevels(gpos$class)) stop("class_cols must have length nlevels(gpos$class)")


    ans <- vector("list", length(plot_windows))


    for(i in seq_along(plot_windows)){

        window <- plot_windows[i]

        # overlapping genes
        near_ge <- suppressWarnings(subsetByOverlaps(nrmisc::txs, window))
        near_ge <- unlist(range(near_ge))
        near_ge <- near_ge[as.character(seqnames(near_ge))==seqnames(window)[1]]


        # immune loci in plotting window
        near_il <- suppressWarnings(subsetByOverlaps(nrmisc::immune, window))

        # expand plotting window so whole gene shown if onyl partially covered by original window
        window <- reduce(c(range(near_ge), range(near_il), range(window)),
                         ignore.strand=TRUE)

        # breakpoints in plotting window
        wpos <- subsetByOverlaps(gpos, window)
        wpos <- sort(wpos, ignore.strand=TRUE)
        # distance to next breakpoint (in any sample)
        mcols(wpos)[,'dist_to_next'] <- c(diff(start(wpos)), NA)
        wpos <- as.data.frame(wpos)
        # remove last point (no distance)
        wpos <- wpos[-nrow(wpos),]
        # refactor breakpoint_id
        wpos$breakpoint_id <- factor(wpos$breakpoint_id)

        # junctions with both ends in plotting window
        both <- split(wpos, wpos$breakpoint_id)
        both <- both[names(which(sapply(both, nrow)==2))]
        both <- plyr::ldply(both)

        if(nrow(both)>0){
            odds <- seq(1, nrow(both), 2)
            evens <- seq(2, nrow(both), 2)
            connex <- data.frame(ax=both$start[odds], ay=log10(both$dist_to_next[odds]+1),
                                 bx=both$start[evens], by=log10(both$dist_to_next[evens]+1),
                                 class=both$class[odds])
        }

        p.sv <- ggplot() +
            geom_jitter(data=wpos, aes_string(x='start', y='log10(dist_to_next+1)',
                                       colour='class')) +
            scale_colour_manual(values=class_cols, drop=FALSE) +
            scale_y_continuous(limits=c(-0.8, 4.8), breaks=0:4,
                               labels=c(0, 10, 100, 1000, 10000),
                               name="Distance to next breakpoint (any sample)") +
            scale_x_continuous(name=seqnames(window)[1], labels=scales::comma) +
            geom_segment(aes(x=(end(window)-1e4), y=4.8,
                             xend=end(window), yend=4.8)) +
            geom_text(aes(label="10kb", x=end(window)-5e3, y=4.65)) +
            theme(panel.grid.minor=element_blank(), text=element_text(size=14),
                  plot.title=element_text(size=12))


        if(!is.null(plot_title)) p.sv <- p.sv + ggtitle(plot_title[i])


        if(nrow(both)>0){
            p.sv <- p.sv + geom_curve(data=connex, size=0.2,
                                      aes_string(x='ax', xend='bx', y='ay', yend='by',
                            colour='factor(class, levels=levels(wpos$class))'))

        }

        if (length(near_ge) > 0) {

            gene_exons <- reduce(unlist(nrmisc::exs[names(near_ge)]))
            gene_exons <- gene_exons[as.character(seqnames(gene_exons))==seqnames(window)[1]]

            p.sv <- p.sv +
                geom_segment(aes(x=start(near_ge), y=-0.5,
                                 xend=end(near_ge), yend=-0.5), size=0.3) +
                geom_text_repel(aes(label=names(near_ge),
                                    x=mid(ranges(near_ge)), y=-0.7),
                                segment.size=0) +
                geom_segment(aes(x=start(gene_exons), y=-0.5,
                                 xend=end(gene_exons), yend=-0.5), size=5)
        }

        if (length(near_il) > 0) {
            near_il <- near_il[order(width(near_il), decreasing=TRUE)[1:5]]
            p.sv <- p.sv +
                geom_segment(aes(x=start(near_il), y=-0.5,
                                 xend=end(near_il), yend=-0.5), size=0.3) +
                geom_text_repel(aes(label=names(near_il),
                                    x=mid(ranges(near_il)), y=-0.5),
                                segment.size=0)
        }

        print(p.sv)
        ans[[i]] <- p.sv
    }

    return(ans)

}
