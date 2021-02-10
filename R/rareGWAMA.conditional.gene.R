#' conditional analysis for single variant association test;
#'
#' @param score.stat.file the file names of score statistic files;
#' @param imp.qual.file the file names of imputation quality;
#' @param vcf.ref.file the file names of the reference panel file;
#' @param candidateVar.df a data frame with pos being the position (pos) of the variants, and groupName being the groups that should be analyzed together;
#' @param knownVar known variant;
#' @param alternative The alternative hypothesis. Default is two.sided;
#' @param col.impqual The column number for the imputation quality score;
#' @param impQual.lb The lower bound for the imputation quality. Variants with imputaiton quality less than impQual.lb will be labelled as missing;
#' @param impQualWeight Using imputation quality as weight
#' @param rmMultiAllelicSite Default is TRUE. Multi-allelic sites will be removed from the analyses if set TRUE, and a variable posMulti will be output; The variant site with multiple alleles can be analyzed using rareGWAMA.single.multiAllele function; 
#' @return A list of analysis results;
#' @export 
rareGWAMA.cond.gene <- function(score.stat.file,imp.qual.file=NULL,vcf.ref.file,candidateVar.df,knownVar,alternative="two.sided",...) {
    uniq.allele <- function(x) {x.tab <- table(x);return(paste(names(x.tab),sep=',',collapse=','))}
    extraPar <- list(...);
    sizePerBatch <- extraPar$sizePerBatch;
    if(is.null(sizePerBatch)) sizePerBatch <- 100;
    refGeno <- extraPar$refGeno;
    col.impqual <- extraPar$col.impqual;
    
    impQual.lb <- extraPar$impQual.lb;
    maf.cutoff <- extraPar$maf.cutoff;
    impQualWeight <- FALSE;
    rmMultiAllelicSite <- extraPar$rmMultiAllelicSite;
    if(is.null(col.impqual)) col.impqual <- 5;
    if(is.null(impQual.lb)) impQual.lb <- 0.7;
    if(is.null(rmMultiAllelicSite)) rmMultiAllelicSite <- TRUE;
    if(is.null(refGeno)) refGeno <- "GT";
    if(is.null(maf.cutoff)) maf.cutoff <- 0.05;
    missing <- extraPar$missing;
    if(is.null(missing)) missing <- 'pseudoscore';
    rvtest <- extraPar$rvtest;
    if(is.null(rvtest)) rvtest <- "burden";
    
    group.vec <- unique(candidateVar.df$groupName);
    statistic <- 0;p.value <- 0;no.site <- 0;pos.gene <- 0;variant.direction.effect <- 0;maf.cutoff.out <- 0;
    for(ii in 1:length(group.vec)) {
        candidateVar.ii <- candidateVar.df$pos[candidateVar.df$groupName==group.vec[ii]]
        candidateVar.ii <- candidateVar.ii[!candidateVar.ii%in%knownVar];
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
        cat('Analyzing',group.vec[ii],'\n',sep=' ');
        cat('Read in',length(raw.data.all$ref[[1]]),'variants\n',sep=' ');
        dat <- GWAMA.formatData(raw.data.all=raw.data.all,
                                raw.imp.qual=raw.imp.qual,
                                impQualWeight=impQualWeight,
                                impQual.lb=impQual.lb,
                                col.impqual=col.impqual,
                                maf.cutoff=maf.cutoff,
                                knownVar=knownVar);
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
            gt.tmp <- geno.list$GT
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
        ix.candidate <- match(intersect(pos,candidateVar.ii),pos);
        ix.known <- match(intersect(pos,knownVar),pos);
        res.cond <- list();
        if(length(ix.candidate)>0 & length(ix.known)>0) {
            if(missing=='pseudoscore') {
                res.cond <- getCondUV(dat=dat,lambda=.1,ix.candidate=ix.candidate,ix.known=ix.known,r2=r2);
            }
            if(missing=='dswm') {
                ustat.mat.known <- matrix(is.na(dat$ustat.mat[ix.known,]),nrow=length(ix.known));
                ix.missing <- which(colSums(ustat.mat.known)>0);
                if(length(ix.missing)<ncol(dat$ustat.mat)) {
                    if(length(ix.missing)>0) {
                        dat$ustat.mat <- as.matrix(dat$ustat.mat[,-ix.missing]);
                        dat$vstat.mat <- as.matrix(dat$vstat.mat[,-ix.missing]);
                        dat$nSample.mat <- as.matrix(dat$nSample.mat[,-ix.missing]);
                        dat$af.mat <- as.matrix(dat$af.mat[,-ix.missing]);
                    }
                    ustat.meta <- rowSums(dat$ustat.mat,na.rm=TRUE);
                    vstat.sq.meta <- rowSums((dat$vstat.mat)^2,na.rm=TRUE);
                    vstat.meta <- sqrt(vstat.sq.meta);
                    nSample.meta <- rowSums(dat$nSample.mat,na.rm=TRUE);
                    N <- max(nSample.meta,na.rm=TRUE);
                    X.T.times.X <- diag(vstat.meta)%*%r2%*%diag(vstat.meta);
                    res.cond <- get.conditional.score.stat(ustat.meta,X.T.times.X,N,ix.candidate,ix.known);
                }
            }
            
            if(missing=='msso') {
                warning("replace missing data with 0 may result in biased conditional analysis results");
                dat$ustat.mat <- rm.na(dat$ustat.mat);
                dat$vstat.mat <- rm.na(dat$vstat.mat);
                dat$nSample.mat <- rm.na(dat$nSample.mat);
                ustat.meta <- rowSums(dat$ustat.mat,na.rm=TRUE);
                vstat.sq.meta <- rowSums((dat$vstat.mat)^2,na.rm=TRUE);
                vstat.meta <- sqrt(vstat.sq.meta);
                nSample.meta <- rowSums(dat$nSample.mat,na.rm=TRUE);
                N <- max(nSample.meta,na.rm=TRUE);
                X.T.times.X <- diag(vstat.meta)%*%r2%*%diag(vstat.meta);
                res.cond <- get.conditional.score.stat(ustat.meta,X.T.times.X,N,ix.candidate,ix.known);           
            }
            res.cond$ustat.meta <- res.cond$conditional.ustat;
            variant.direction.effect.tmp <- rep("X",length(res.cond$conditional.ustat));
            variant.direction.effect.tmp[which(res.cond$conditional.ustat>0)] <- "+";
            variant.direction.effect.tmp[which(res.cond$conditional.ustat<0)] <- "-";
            variant.direction.effect.tmp[which(res.cond$conditional.ustat==0)] <- "=";
            variant.direction.effect[ii] <- paste(variant.direction.effect.tmp,sep='',collapse='');
            maf.cutoff.out[ii] <- maf.cutoff
            no.site[ii] <- length(res.cond$conditional.ustat);
            res.cond$V.meta <- res.cond$conditional.V;
            res.cond$maf.meta <- dat$maf.meta[ix.candidate];
            res.cond$nSample.meta <- rowSums(dat$nSample.mat,na.rm=TRUE)[ix.candidate];
            if(rvtest=='skat') res.assoc <- rareGWAMA.skat(dat=res.cond);
            if(rvtest=='burden') res.assoc <- rareGWAMA.burden(dat=res.cond);
            if(rvtest=='vt') {
                res.assoc <- rareGWAMA.vt(dat=res.cond);
                maf.cutoff.out[ii] <- res.assoc$maf.cutoff.vt;
                no.site[ii] <- res.assoc$no.site.VT;
            }
            if(rvtest=='all') {
                res.assoc.skat <- rareGWAMA.skat(dat=res.cond);
                res.assoc.burden <- rareGWAMA.burden(dat=res.cond);

                res.assoc.vt <- rareGWAMA.vt(dat=res.cond);
                res.assoc <- list(statistic=paste(format(res.assoc.burden$statistic,digits=2),format(res.assoc.skat$statistic,digits=2),format(res.assoc.vt$statistic,digits=2),sep=',',collapse=','),
                                  p.value=paste(format(res.assoc.burden$p.value,digits=2),format(res.assoc.skat$p.value,digits=2),format(res.assoc.vt$p.value,digits=2),sep=',',collapse=','))                                      
                maf.cutoff.out[ii] <- paste(maf.cutoff,format(res.assoc.vt$maf.cutoff.vt,digits=2),sep=',',collapse=',');
                no.site[ii] <- paste(no.site[ii],res.assoc.vt$no.site.VT,sep=',',collapse=',');
            }
            statistic[ii] <- format(res.assoc$statistic,digits=2);
            p.value[ii] <- format(res.assoc$p.value,digits=2);
            pos.gene[ii] <- paste(intersect(pos,candidateVar.ii),sep=",",collapse=",");
        }
    }
    res.formatted <- cbind(group.vec,
                          statistic,
                          p.value,
                          no.site,
                          maf.cutoff.out,
                          variant.direction.effect,
                          pos.gene);
    colnames(res.formatted) <- c("GENE","STATISTIC","PVALUE","N_SITE", "MAF_CUTOFF","SINGLE_VAR_EFFECT","VARIANTS");
    return(list(res.formatted=res.formatted));
}
