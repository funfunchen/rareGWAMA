#' calculate GC by frequency bins;
#'
#' @param score.stat.file The file names for score statistics;
#' @param maf.bin The frequency bins. The GC for variants in each frequency bin will be calculated separately;
#' @return A list which consist of GCs for each study;
#' @export 
rareMETALS.calcGC <- function(score.stat.file,maf.bin) {
    if(maf.bin[1]!=0 | maf.bin[length(maf.bin)]!=1) stop('the start and end of the MAF bin have to be 0 and 1');
    res.out <- list();
    for(ii in 1:length(score.stat.file)) {
        cat('Calculating the genomic control value for ',score.stat.file[ii],' \n');
        
        file.ex <- substr(score.stat.file[ii],nchar(score.stat.file[ii])-2,nchar(score.stat.file[ii]));
        if(file.ex=='.gz') {
                                                                           }
        if(file.ex!='.gz') {
                                                                          }
        res <- as.data.frame(res);
        res.header <- colnames(res);
        ix.af <- which(res.header=="ALL_AF");
        if(length(ix.af)>0) {
            res.header[ix.af] <- "ALT_FREQ";
            colnames(res) <- gsub("#","",res.header);
        }
        maf <- as.numeric(res$ALT_FREQ);
        maf[which(maf>.5)] <- 1-maf[which(maf>.5)];
        lambda <- 0;
        
        for(jj in 1:(length(maf.bin)-1)) {
            ix.jj <- which(maf>=maf.bin[jj] & maf<=maf.bin[jj+1]);
            pval <- res$PVALUE[ix.jj];
            mid.pval <- median(pval,na.rm=TRUE);
            lambda[jj] <- qchisq(mid.pval,df=1,lower.tail=FALSE)/qchisq(.5,df=1,lower.tail=FALSE);
            
        }
        lambda[which(is.na(lambda))] <- 1;
        res.out[[ii]] <- cbind(maf.bin[1:(length(maf.bin)-1)],maf.bin[2:(length(maf.bin))],lambda);
        
    }
    return(res.out);
}
