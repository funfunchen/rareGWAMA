#' conditional analysis for single variant association test;
#'
#' @param score.stat.file the file names of score statistic files;
#' @param imp.qual.file the file names of imputation quality;
#' @param vcf.ref.file the file names of the reference panel file;
#' @param candidateVar the tabix range;
#' @param knownVar known variant;
#' @param alternative The alternative hypothesis. Default is two.sided;
#' @param col.impqual The column number for the imputation quality score;
#' @param impQual.lb The lower bound for the imputation quality. Variants with imputaiton quality less than impQual.lb will be labelled as missing;
#' @param impQualWeight Using imputation quality as weight
#' @param rmMultiAllelicSite Default is TRUE. Multi-allelic sites will be removed from the analyses if set TRUE, and a variable posMulti will be output; The variant site with multiple alleles can be analyzed using rareGWAMA.single.multiAllele function; 
#' @return A list of analysis results;
#' @export 
rareGWAMA.cond.single.msso <- function(score.stat.file,imp.qual.file=NULL,vcf.ref.file,candidateVar,knownVar,alternative="two.sided",...) {
    uniq.allele <- function(x) {x.tab <- table(x);return(paste(names(x.tab),sep=',',collapse=','))}
    extraPar <- list(...);
    sizePerBatch <- extraPar$sizePerBatch;
    if(is.null(sizePerBatch)) sizePerBatch <- 100;
    refGeno <- extraPar$refGeno;
    col.impqual <- extraPar$col.impqual;

    impQual.lb <- extraPar$impQual.lb;
    impQualWeight <- FALSE;
    rmMultiAllelicSite <- extraPar$rmMultiAllelicSite;
    if(is.null(col.impqual)) col.impqual <- 5;
    if(is.null(impQual.lb)) impQual.lb <- 0.7;
    if(is.null(rmMultiAllelicSite)) rmMultiAllelicSite <- TRUE;
    if(is.null(refGeno)) refGeno <- "DS";
    pseudoScore <- extraPar$pseudoScore;
    if(is.null(pseudoScore)) pseudoScore <- TRUE;
    beta.est <- 0;beta.se <- 0;statistic <- 0;p.value <- 0;ref.tab <- 0;alt.tab <- 0;pos.all <- 0;
    ii <- 0;batchEnd <- 0;batchStart <- 0;nSample <- 0;af <- 0;
    while(batchEnd<length(candidateVar)) {
        batchStart <- batchEnd+1;
        batchEnd <- batchStart+sizePerBatch;
        if(batchEnd>length(candidateVar)) batchEnd <- length(candidateVar);
        candidateVar.ii <- candidateVar[batchStart:batchEnd];
        tabix.range <- get.tabix.range(c(candidateVar.ii,knownVar));
        a <- Sys.time();
        capture.output(raw.data.all <- rvmeta.readDataByRange( score.stat.file, NULL, tabix.range,multiAllelic = TRUE));
        vcfIndv <- refGeno;
        annoType <- "";
        vcfColumn <- c("CHROM","POS","REF","ALT");
        vcfInfo <- NULL;      
        geno.list <- readVCFToListByRange(vcf.ref.file, tabix.range, "", vcfColumn, vcfInfo, vcfIndv)        
        raw.imp.qual <- NULL;
        if(!is.null(imp.qual.file))
            raw.imp.qual <- lapply(imp.qual.file,tabix.read.table,tabixRange=tabix.range);
        time.readData <- Sys.time()-a;
        b <- Sys.time();
        raw.data.all <- raw.data.all[[1]];
        cat('Read in',length(raw.data.all$ref[[1]]),'variants\n',sep=' ');
        dat <- GWAMA.formatData(raw.data.all,raw.imp.qual,impQualWeight,impQual.lb,col.impqual);
        if(rmMultiAllelicSite==TRUE) {
            tmp <- GWAMA.rmMulti(dat);
            dat <- tmp$dat;posMulti <- tmp$posMulti;
        }
        pos <- gsub("_.*","",dat$pos);
        if(refGeno=="DS") {
            gt <- geno.list$DS;
            gt <- matrix(as.numeric(gt),nrow=nrow(gt),ncol=ncol(gt));
        }
        if(refGeno=="GT") {
            gt.tmp <- geno.list$GT;
   
            gt <- matrix(NA,nrow=nrow(gt.tmp),ncol=ncol(gt.tmp));
            gt[which(gt.tmp=="0/0",arr.ind=T)] <- 0;
            gt[which(gt.tmp=="1/0",arr.ind=T)] <- 1;
            gt[which(gt.tmp=="0/1",arr.ind=T)] <- 1;
            gt[which(gt.tmp=="1/1",arr.ind=T)] <- 2
            gt[which(gt.tmp=="0|0",arr.ind=T)] <- 0;
            gt[which(gt.tmp=="1|0",arr.ind=T)] <- 1;
            gt[which(gt.tmp=="0|1",arr.ind=T)] <- 1;
            gt[which(gt.tmp=="1|1",arr.ind=T)] <- 2
            
        }
        
        r2.tmp <- cor(gt,use='pairwise.complete');
        r2.tmp <- rm.na(r2.tmp);
        pos.vcf <- paste(geno.list$CHROM,geno.list$POS,sep=":");
        r2 <- matrix(0,nrow=length(pos),ncol=length(pos));
        r2 <- as.matrix(r2.tmp[match(pos,pos.vcf),match(pos,pos.vcf)]);
        
        diag(r2) <- 1;
        for(kk in 1:length(candidateVar.ii)) {
            batchId <- (batchStart:batchEnd)[kk];
            ix.candidate <- match(intersect(pos,candidateVar.ii[kk]),pos);        
            ix.known <- match(intersect(pos,knownVar),pos);
            res.cond <- list();
            if(length(ix.candidate)>0 & length(ix.known)>0) {
                
                dat.replace0 <- dat;
                dat.replace0$ustat.mat <- rm.na(dat$ustat.mat);
                dat.replace0$vstat.mat <- rm.na(dat$vstat.mat);
                res.cond <- getCondUV(dat=dat.replace0,lambda=.1,ix.candidate=ix.candidate,ix.known=ix.known,r2=r2);
                
                statistic[batchId] <- (res.cond$conditional.ustat)^2/diag(res.cond$conditional.V);
                p.value[batchId] <- pchisq(statistic[batchId],df=1,lower.tail=F);
                beta.est[batchId] <- res.cond$conditional.beta.est;
                beta.se[batchId] <- sqrt(diag(res.cond$conditional.beta.var));
                
                nSample[batchId] <- res.cond$nSample[ix.candidate];
                ref.tab[batchId] <- apply(matrix(dat$ref.mat[ix.candidate,],nrow=1),1,uniq.allele);
                alt.tab[batchId] <- apply(matrix(dat$alt.mat[ix.candidate,],nrow=1),1,uniq.allele);
                pos.all[batchId] <- pos[ix.candidate];
                af.meta <- rowSums((dat$af.mat)*(dat$nSample.mat),na.rm=T)/rowSums((dat$nSample.mat),na.rm=T);
                af[batchId] <- af.meta[ix.candidate];
            }
        }
    }
    
    res.formatted <- cbind(pos.all,
                           ref.tab,
                           alt.tab,
                           format(af,digits=3),
                           format(statistic,digits=3),
                           format(p.value,digits=3),
                           format(beta.est,digits=3),
                           format(beta.se,digits=3),
                           nSample);
    colnames(res.formatted) <- c("POS","REF","ALT","AF","STAT","PVALUE","BETA","SD","N");
    return(list(res.formatted=res.formatted,
                dat=dat,
                raw.data.all=raw.data.all,
                res.cond=res.cond));
}
