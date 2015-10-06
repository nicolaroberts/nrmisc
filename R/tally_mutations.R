#' Tally mutations in the 96 base substitution classes defined by trinucleotide context
#'
#' @param gr A GRanges object with ranges of width 1 (variant positions)
#'  and at least three metadata columns called 'ref' for reference base,
#'  'alt' for alternate base, and 'sampleID' for sample name.
#'  Only hg19 chroms chr1-chr22, chrX and chrY are kept.
#' @return A matrix of mutation counts with one row per sample and one
#'  column per base substitution class (trinucleotide context)
#' @import GenomicRanges IRanges S4Vectors GenomeInfoDb BSgenome.Hsapiens.UCSC.hg19 Biostrings
#' @export
#' @examples
#' data(ov_tcga, package='SomaticCancerAlterations')
#' ov_tcga <- ov_tcga[ov_tcga$Variant_Type=="SNP"]
#' ov_tcga <- ov_tcga[1:1000]
#' gr <- ov_tcga
#' S4Vectors::mcols(gr) <- NULL
#' gr$ref <- ov_tcga$Reference_Allele
#' gr$alt <- ov_tcga$Tumor_Seq_Allele2
#' gr$sampleID <- as.factor(ov_tcga$Sample_ID)
#' tally <- tally_mutations_96(gr)
tally_mutations_96 <- function(gr){

    # check widths all 1
    if (any(width(gr) != 1)) stop("all widths must be 1")


    # ensure chrom names start with chr
    if (all(substr(seqlevels(gr), 1, 3) != "chr")) {
        gr <- renameSeqlevels(gr, paste0('chr', seqlevels(gr)))
    } else if (any(substr(seqlevels(gr), 1, 3) != "chr")) {
        stop("Either ALL or NONE of seqlevels(gr) can start with 'chr'")
    }

    genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

    # ony keep chr1-chr22, chrX, chrY
    gr <- keepSeqlevels(gr, seqlevels(genome)[1:24])

    # add seqinfo
    seqinfo(gr) <- seqinfo(genome)

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
    bad <- which(substr(gr$trinuc, 2, 2) != gr$ref)
    if (length(bad) > 0) {
        warning(sprintf("%d input positions had 'ref' base not matching hg19 reference - dropped %s", length(bad), paste(bad, collapse=', ')))

        gr <- gr[-bad]
    }

    # convert mutations wrt A or G ref base into reverse complement
    gr$REF <- gr$ref
    gr$ALT <- gr$alt
    gr$context <- gr$trinuc
    torc <- which(gr$ref %in% c('A', 'G'))
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
    tally <- table(gr$sampleID, gr$class)
    class(tally) <- "matrix"

    return(tally)
}
