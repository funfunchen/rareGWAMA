#' Single variant meta-analysis
#'
#' @param score.stat.file files of score statistics
#' @param callrate.cutoff Cutoffs of call rate, lower than which will NOT be analyzed (labelled as missing)
#' @param hwe.cutoff Cutoffs of HWE p-values
#' @param refaltList A list that contains pos, ref and alt for calibrating variant sites to have the right ref and alt alleles;
#' @return a list consisting of results
#' @export
rareMETALS.cleanScore <- function(score.stat.file,refaltList,hwe.cutoff,callrate.cutoff,af.diff.max=0.7,checkAF=TRUE)
    {
        res <- read.table(gzfile(score.stat.file),as.is=TRUE,header=TRUE,fill=TRUE);
        pos.vec <- paste(res$CHROM,res$POS,sep=":");
        pos.both <- set.intersect(refaltList$pos,pos.vec);
        ix.match <- match(pos.both,pos.vec);
        res <- res[ix.match,];
        #for monomorphic variants:
        ix.mono <- which(rm.na(res$AF)==0 | rm.na(res$AF)==1);
        

    }
