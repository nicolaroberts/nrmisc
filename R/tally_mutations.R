#' Tally mutations in the 96 base substitution classes defined by trinucleotide context
#'
#' @param gr A GRanges object with ranges of width 1 (variant positions)
#'  and at least three metadata columns for reference base, alternate base,
#'  and sample name.
#' @param genome A BSgenome object of the reference genome to use for
#'  trinucleotide context calculations (e.g. 'BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19')
#' @param ref Variable name for the reference base column in gr (default 'ref').
#' @param alt Variable name for the alternative base column in gr (default 'alt').
#' @param samp Variable name for the sample ID column in gr (default 'sampleID').
#' @return A matrix of mutation counts with one row per sample and one
#'  column per base substitution class (trinucleotide context)
#' @import GenomicRanges Biostrings GenomeInfoDb
#' @export
tally_mutations_96 <- function(gr, genome, ref='ref', alt='alt', samp='sampleID'){

    # check genome is BSgenome object
    if(class(genome)!='BSgenome') stop('genome must be BSgenome object')

    # check column names of gr
    if(any(!c(ref, alt, samp) %in% colnames(mcols(gr)))){
        stop('gr must have three metadata columns with names indicated by ref, alt and samp arguments')
    }

    # check widths all 1
    if (any(width(gr) != 1)) stop("all gr widths must be 1")

    # ensure gr seqlevels are in reference genome
    if (any(!seqlevels(gr) %in% seqlevels(genome))){
        stop("one or more seqlevels(gr) not present in reference genome provided")
    }

    # add seqinfo
    seqinfo(gr) <- seqinfo(genome)[seqlevels(gr)]

    # resort
    gr <- sort(gr)

    # find the trinuc context of each mutation
    bases <- c("A", "C", "G", "T")
    trinuc_levels <- paste0(rep(bases, each=16), rep(rep(bases, each=4), 4), rep(bases, 16))

    get_trinuc <- function(seqname){
        pos <- start(gr[seqnames(gr)==seqname])
        view <- Views(genome[[seqname]], start=pos-1, end=pos+1)
        ans <- factor(as.character(view), levels=trinuc_levels, labels=1:64)
        return(as.numeric(ans))
    }

    trinuc <- sapply(seqlevels(gr), get_trinuc)
    gr$trinuc <- factor(unlist(trinuc, use.names=FALSE), levels=1:64, labels=trinuc_levels)
    remove(trinuc)

    # are there any places where the reference doesn't match up?
    bad <- which(substr(gr$trinuc, 2, 2) != mcols(gr)[,ref])
    if (length(bad) > 0) {
        warning(sprintf("%d input positions had 'ref' base not matching provided reference - dropped %s", length(bad), paste(bad, collapse=', ')))

        gr <- gr[-bad]
    }
    if(length(gr)==0) stop('No gr entries match reference base in genome')

    # convert mutations wrt A or G ref base into reverse complement
    gr$REF <- mcols(gr)[,ref]
    gr$ALT <- mcols(gr)[,alt]
    gr$context <- gr$trinuc
    torc <- which(mcols(gr)[,ref] %in% c('A', 'G'))
    gr$REF[torc] <- as.character(reverseComplement(DNAStringSet(gr$REF[torc])))
    gr$ALT[torc] <- as.character(reverseComplement(DNAStringSet(gr$ALT[torc])))
    gr$context[torc] <- as.character(reverseComplement(DNAStringSet(gr$context[torc])))


    # define classes
    gr$class <- paste(gr$REF, gr$ALT, 'in', gr$context, sep='.')

    class_levels <- paste(rep(c("C", "T"), each=48),
                          rep(c("A", "G", "T", "A", "C", 'G'), each=16),
                          'in', paste0(rep(rep(bases, each=4), 6),
                                       rep(c("C", "T"), each=48),
                                       rep(bases, 24)), sep='.')

    gr$class <- factor(gr$class, levels=class_levels)

    # tally mutation classes in each tumour
    tally <- table(mcols(gr)[,samp], gr$class)
    class(tally) <- "matrix"

    return(tally)
}
