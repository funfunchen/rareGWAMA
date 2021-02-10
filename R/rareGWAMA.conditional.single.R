#' conditional analysis for single variant association test;
#'
#' @param score.stat.file the file names of score statistic files;
#' @param imp.qual.file the file names of imputation quality;
#' @param ref.file the file names of the reference panel file; For VCF references, use the VCF file name; for binary PLINK references, use the prefix; 
#' @param candidateVar the tabix range;
#' @param knownVar known variant;
#' @param alternative The alternative hypothesis. Default is two.sided;
#' @param col.impqual The column number for the imputation quality score;
#' @param impQual.lb The lower bound for the imputation quality. Variants with imputaiton quality less than impQual.lb will be labelled as missing;
#' @param impQualWeight Using imputation quality as weight
#' @param rmMultiAllelicSite Default is TRUE. Multi-allelic sites will be removed from the analyses if set TRUE, and a variable posMulti will be output; The variant site with multiple alleles can be analyzed using rareGWAMA.single.multiAllele function; 
#' @return A list of analysis results;
#' @export 
rareGWAMA.cond.single <- function(score.stat.file,imp.qual.file=NULL,ref.file,candidateVar,knownVar,alternative="two.sided",...) {
    candidateVar <- gsub("chr","",candidateVar);
    knownVar <- gsub("chr","",knownVar);
    uniq.allele <- function(x) {x.tab <- table(x);return(paste(names(x.tab),sep=',',collapse=','))}
    extraPar <- list(...);
    vcf.ref.file <- ref.file;
    if(is.null(extraPar$col.ref.plink)) extraPar$col.ref.plink <- 6;
    if(is.null(extraPar$col.alt.plink)) extraPar$col.alt.plink <- 5;
    if(is.null(extraPar$refFileFormat)) extraPar$refFileFormat <- "vcf";
    if(is.null(extraPar$trans.ethnic)) extraPar$trans.ethnic <- FALSE;
    if(extraPar$refFileFormat=="plink") extraPar$plinkPrefix <- ref.file;
    r2.cutoff <- extraPar$r2.cutoff;
    if(is.null(r2.cutoff)) r2.cutoff <- 0.95;
    nBoot <- extraPar$nBoot;
    if(is.null(nBoot)) nBoot <- 100;
    sizePerBatch <- extraPar$sizePerBatch;
    if(is.null(sizePerBatch)) sizePerBatch <- 100;
    refGeno <- extraPar$refGeno;
    col.impqual <- extraPar$col.impqual;
    if(is.null(extraPar$chrVcfPrefix)) extraPar$chrVcfPrefix <- "";
    if(is.null(extraPar$chrSumstatPrefix)) extraPar$chrSumstatPrefix <- "";
    impQual.lb <- extraPar$impQual.lb;
    impQualWeight <- FALSE;
    rmMultiAllelicSite <- extraPar$rmMultiAllelicSite;
    if(is.null(col.impqual)) col.impqual <- 5;
    if(is.null(impQual.lb)) impQual.lb <- 0.7;
    if(is.null(rmMultiAllelicSite)) rmMultiAllelicSite <- TRUE;
    if(is.null(refGeno)) refGeno <- "GT";
    pseudoScore <- extraPar$pseudoScore;
    method <- extraPar$method;
    if(is.null(method)) method <- "PCBS";
    
    if(is.null(pseudoScore)) pseudoScore <- TRUE;
    beta.est <- 0;beta.se <- 0;statistic <- 0;p.value <- 0;ref.tab <- 0;alt.tab <- 0;pos.all <- 0;marginal.statistic <- 0;marginal.p.value <- 0;pos.out <- 0;pval.memo <- 0;pval.re <- 0;
    ii <- 0;batchEnd <- 0;batchStart <- 0;nSample <- 0;af <- 0;numStability <- 0;refAvail <- 0;varId <- 0;
    chisq.mat <- matrix(nrow=length(candidateVar)*100,ncol=ncol(extraPar$af.pca));
    if(extraPar$refFileFormat=="plink") {
        fname.bim <- paste0(extraPar$plinkPrefix,".bim");
        fname.fam <- paste0(extraPar$plinkPrefix,".fam");
        
        bim <- as.data.frame(fread(fname.bim,header=FALSE));
        fam <- as.data.frame(fread(fname.fam,header=FALSE));
        ix.tmp <- match(fam[,1],extraPar$ref.ancestry[,1]);
        extraPar$ref.ancestry <- extraPar$ref.ancestry[ix.tmp,];
        indvIx <- 1:nrow(fam);
        chrpos.plink <- paste(bim[,1],bim[,4],sep=":");
        markerIx <- match(c(paste0(extraPar$chrSumstatPrefix,candidateVar),paste0(extraPar$chrSumstatPrefix,knownVar)),chrpos.plink);
        ix.rm <- which(is.na(markerIx));
        if(length(ix.rm)>0) markerIx <- markerIx[-ix.rm];
        geno.plink <- readPlinkToMatrixByIndex(extraPar$plinkPrefix, indvIx, markerIx);
        ref.plink <- unlist(bim[markerIx,5]);
        alt.plink <- unlist(bim[markerIx,6]);
        chr.plink <- unlist(bim[markerIx,1]);
        pos.plink <- unlist(bim[markerIx,4]);
        chrpos.plink <- paste(chr.plink,pos.plink,sep=":");
    }
    while(batchEnd<length(candidateVar)) {
        batchStart <- batchEnd+1;
        batchEnd <- batchStart+sizePerBatch;
        if(batchEnd>length(candidateVar)) batchEnd <- length(candidateVar);
        candidateVar.ii <- candidateVar[batchStart:batchEnd];
        tabix.range <- get.tabix.range(c(paste0(extraPar$chrSumstatPrefix,candidateVar.ii),paste0(extraPar$chrSumstatPrefix,knownVar)));
        a <- Sys.time();
        capture.output(raw.data.all <- rvmeta.readDataByRange( score.stat.file, NULL, tabix.range,multiAllelic = TRUE));
        raw.data.all[[1]]$pos <- gsub("chr","",raw.data.all[[1]]$pos);
        vcfIndv <- refGeno;
        annoType <- "";
        vcfColumn <- c("CHROM","POS","REF","ALT");
        vcfInfo <- NULL;
        chr.tmp <- gsub("chr","",gsub(":.*","",c(candidateVar.ii,knownVar)));
        
        pos.tmp <- gsub(".*:","",c(candidateVar.ii,knownVar));
        tabix.range.vcf <- paste(paste(extraPar$chrVcfPrefix,chr.tmp,sep=""),paste(pos.tmp,pos.tmp,sep="-"),sep=":",collapse=",");
        if(extraPar$refFileFormat=="vcf") {
            geno.list <- readVCFToListByRange(vcf.ref.file, tabix.range.vcf, "", vcfColumn, vcfInfo, vcfIndv)
            if(length(geno.list$POS)==0) {
                geno.list$GT <- matrix(nrow=length(geno.list$sampleId),ncol=0);
                geno.list$DS <- matrix(nrow=length(geno.list$sampleId),ncol=0);
            }
            rownames(geno.list$GT) <- geno.list$sampleId;
            rownames(geno.list$DS) <- geno.list$sampleId;
                
        }
        if(extraPar$refFileFormat=="vcf.vbi") {
            geno.list <- list();
            f0 <- file();
            sink(file=f0,type='output')
            sink(file=f0,type='message');

            geno.vbi <- (seqminer::readSingleChromosomeVCFToMatrixByRange(fileName=vcf.ref.file,range=tabix.range.vcf))[[1]];
            

            sink(type='output');
            sink(type='message');
            close(f0);
            chrompos<- gsub("chr","",gsub("_.*","",colnames(geno.vbi)));
            refalt <- gsub(".*_","",colnames(geno.vbi));
            geno.list$CHROM <- gsub(":.*","",chrompos);
            geno.list$POS <- gsub(".*:","",chrompos);
            geno.list$REF <- gsub("/.*","",refalt);
            geno.list$ALT <- gsub(".*/","",refalt);
            geno.list$GT <- geno.vbi;
            colnames(geno.list$GT) <- paste(geno.list$CHROM,geno.list$POS,sep=":");
            geno.list$sampleId <- rownames(geno.vbi);
            
        }
        raw.imp.qual <- NULL;
        if(!is.null(imp.qual.file))
            raw.imp.qual <- lapply(imp.qual.file,tabix.read.table,tabixRange=tabix.range);
        time.readData <- Sys.time()-a;
        b <- Sys.time();
        raw.data.all <- raw.data.all[[1]];
        cat('Read in',length(unique(gsub("_.*","",raw.data.all$pos))),'variants',sep=' ');
        dat <- GWAMA.formatData(raw.data.all,raw.imp.qual,impQualWeight,impQual.lb,col.impqual,rmMultiAllelicSite=rmMultiAllelicSite);
        if(rmMultiAllelicSite==TRUE) {
            tmp <- GWAMA.rmMulti(dat);
            dat <- tmp$dat;posMulti <- tmp$posMulti;
        }
           
        pos <- dat$pos;
        pos.tmp <- gsub("_.*","",dat$pos);
        if(extraPar$refFileFormat=='plink') {
            ix.match.plink <- match(pos.tmp,chrpos.plink)       
            gt <- geno.plink[,ix.match.plink];
            
            
            geno.list <- list(CHROM=chr.plink[ix.match.plink],POS=pos.plink[ix.match.plink],REF=ref.plink[ix.match.plink],ALT=alt.plink[ix.match.plink],sampleId=rownames(geno.plink));
            ix.flip <- which(dat$ref.tab!=geno.list$REF);
            geno.list$REF <- dat$ref.tab;
            geno.list$ALT <- dat$alt.tab;            
            if(length(ix.flip)>0) gt[,ix.flip] <- 2-gt[,ix.flip];
        }
        
        if(extraPar$refFileFormat=="vcf") {
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
        }
        if(extraPar$refFileFormat=="vcf.vbi") {
            gt <- geno.list$GT;
            ix.match.vbi <- match(pos.tmp,colnames(gt))
            gt <- Matrix(as.matrix(gt[,ix.match.vbi]));
            geno.list <- list(CHROM=geno.list$CHROM[ix.match.vbi],
                              POS=geno.list$POS[ix.match.vbi],
                              REF=geno.list$REF[ix.match.vbi],
                              ALT=geno.list$ALT[ix.match.vbi],
                              gt=gt,
                              sampleId=rownames(gt));
        }
        if(ncol(gt)>0) {
            r2.tmp <- corSparse(gt);
        }
        if(ncol(gt)==0) r2.tmp <- matrix(nrow=0,ncol=0);
        r2.list <- list();
        ancestry.grp <- NULL;
        if(extraPar$trans.ethnic==TRUE) {
            z2.simu <- matrix((rnorm(1e6*ncol(extraPar$af.pca)))^2,ncol=ncol(extraPar$af.pca));
            
            z2.simu.cumsum <- matrix(apply(z2.simu,1,cumsum),ncol=ncol(extraPar$af.pca),byrow=TRUE);
            p.simu.cumsum <- matrix(nrow=nrow(z2.simu.cumsum),ncol=ncol(z2.simu.cumsum));
            for(ii in 1:ncol(z2.simu.cumsum)) p.simu.cumsum[,ii] <- pchisq(z2.simu.cumsum[,ii],df=ii,lower.tail=FALSE)
            minp.vec <- apply(p.simu.cumsum,1,min);
            
            invZ.p.simu.cumsum <- qnorm(p.simu.cumsum,lower.tail=FALSE);
            cor.z2 <- cor(invZ.p.simu.cumsum);
            extraPar$cor.z2 <- cor.z2;
            
            ancestry.grp <- unique(unlist(strsplit(extraPar$ref.ancestry[,2],split=",")));
            for(ee in 1:length(ancestry.grp)) {
                id.ee <- extraPar$ref.ancestry[grep(paste0("\\b",ancestry.grp[ee],"\\b"),extraPar$ref.ancestry[,2]),1];
                ix.match <- match(id.ee,geno.list$sampleId);
                r2.list[[ee]] <- corSparse(Matrix(matrix(gt[ix.match,],ncol=ncol(gt))));
                
                colnames(r2.list[[ee]]) <- pos.tmp;
            }
            names(r2.list) <- ancestry.grp;
        }
        geno.list$CHROM <- gsub("chr","",geno.list$CHROM);
        pos.vcf <- paste(geno.list$CHROM,geno.list$POS,sep=":");
        refalt <- paste(geno.list$REF,geno.list$ALT,sep="/");
        pos.vcf <- paste(pos.vcf,refalt,sep="_");
        
        r2 <- as.matrix(r2.tmp[match(pos,pos.vcf),match(pos,pos.vcf)]);
        r2.raw <- r2;
        r2 <- rm.na(r2);
        colnames(r2) <- pos;
        diag(r2) <- 1;
        
        for(kk in 1:length(candidateVar.ii)) {
            batchId <- (batchStart:batchEnd)[kk];
            ix.candidate <- which(pos.tmp==intersect(pos.tmp,candidateVar.ii[kk]));
            ix.known <- (1:length(pos.tmp))[pos.tmp%in%intersect(pos.tmp,knownVar)];
            res.cond <- list();
            cond.ok <- 0;
            if(length(ix.candidate)>0 & length(ix.known)>0) {
                varId.start <- varId+1;
                varId.end <- varId+length(ix.candidate);
                ref.tab[varId.start:varId.end] <- apply(matrix(dat$ref.mat[ix.candidate,],nrow=1),1,uniq.allele);
                alt.tab[varId.start:varId.end] <- apply(matrix(dat$alt.mat[ix.candidate,],nrow=1),1,uniq.allele);
                pos.all[varId.start:varId.end] <- pos[ix.candidate];
                pos.out[varId.start:varId.end] <- gsub("_.*","",pos.all[varId.start:varId.end]);
                numStability[varId.start:varId.end] <- (length(which(abs(r2[ix.candidate,ix.known])>r2.cutoff))==0 && !is.na(sum(abs(r2[ix.candidate,ix.known]))))
                af.meta <- rowSums((dat$af.mat)*(dat$nSample.mat),na.rm=T)/rowSums((dat$nSample.mat),na.rm=T);
                af[varId.start:varId.end] <- af.meta[ix.candidate];
                refAvail[varId.start:varId.end] <- as.numeric(sum(is.na(r2.raw[ix.candidate,ix.known]))==0);
            }
            if(length(ix.candidate)>0 & length(ix.known)>0) {
                if(method=="PCBS") {
                    dat$ustat.mat.ori <- dat$ustat.mat;
                    dat$vstat.mat.ori <- dat$vstat.mat;
                    dat$z.mat <- dat$ustat.mat/dat$vstat.mat;
                    
                    dat$vstat.mat <- sqrt(2*dat$af.mat*(1-dat$af.mat)*dat$nSample.mat)*dat$w.mat;
                    dat$ustat.mat <- (dat$z.mat)*(dat$vstat.mat);
                    if(is.null(extraPar$sandwich)) {
                        extraPar$sandwich <- TRUE;
                    }
                    if(extraPar$sandwich==FALSE) {
                        res.cond <- getCondUV.wrapper(dat=dat,lambda=extraPar$lambda,ix.candidate=ix.candidate,ix.known=ix.known,r2=r2,sandwich=extraPar$sandwich,r2.list=r2.list,study.ancestry=extraPar$study.ancestry,trans.ethnic=extraPar$trans.ethnic,af.pca=extraPar$af.pca,cor.z2=extraPar$cor.z2);
                    }
                    if(extraPar$sandwich==TRUE) {
                        res.cond <- getCondUV.wrapper(dat=dat,lambda=extraPar$lambda,ix.candidate=ix.candidate,ix.known=ix.known,r2=r2.raw,sandwich=extraPar$sandwich,boot=TRUE,nBoot=nBoot,r2.list=r2.list,study.ancestry=extraPar$study.ancestry,trans.ethnic=extraPar$trans.ethnic,af.pca=extraPar$af.pca,cor.z2=extraPar$cor.z2);
                        
                    }
                    
                    cond.ok <- 1;
                    
                    if(any(res.cond$numStability==0)) {
                        numStability[varId.start:varId.end] <- 0;
                        cond.ok <- 0
                    }                    
                }
                
                if(method=="DISCARD") {
                    dat$ustat.mat.ori <- dat$ustat.mat;
                    dat$vstat.mat.ori <- dat$vstat.mat;
                    dat$z.mat <- dat$ustat.mat/dat$vstat.mat;
                    dat$vstat.mat <- sqrt(2*dat$af.mat*(1-dat$af.mat)*dat$nSample.mat)*dat$w.mat;
                    dat$ustat.mat <- (dat$z.mat)*(dat$vstat.mat);
                    tmp.ustat <- matrix(dat$ustat.mat[c(ix.candidate,ix.known),],nrow=length(ix.candidate)+length(ix.known));
                    ix.missing <- which(colSums(is.na(tmp.ustat))>0);
                    dat.rmMissing <- dat;
                    if(length(ix.missing)==0) {
                        res.cond <- getCondUV(dat=dat.rmMissing,lambda=extraPar$lambda,ix.candidate=ix.candidate,ix.known=ix.known,r2=r2);
                        cond.ok <- 1;
                    }
                    if(length(ix.missing)>0 & length(ix.missing)<ncol(dat$ustat.mat)) {
                        dat.rmMissing$ustat.mat <- matrix(dat$ustat.mat[,-ix.missing],ncol=ncol(dat$ustat.mat)-length(ix.missing));
                        dat.rmMissing$vstat.mat <- matrix(dat$vstat.mat[,-ix.missing],ncol=ncol(dat$ustat.mat)-length(ix.missing));
                        res.cond <- getCondUV(dat=dat.rmMissing,lambda=extraPar$lambda,ix.candidate=ix.candidate,ix.known=ix.known,r2=r2);
                        cond.ok <- 1;
                    }
                }

                if(method=="REPLACE0") {
                    dat$ustat.mat.ori <- dat$ustat.mat;
                    dat$vstat.mat.ori <- dat$vstat.mat;
                    dat$z.mat <- dat$ustat.mat/dat$vstat.mat;
                    dat$vstat.mat <- sqrt(2*dat$af.mat*(1-dat$af.mat)*dat$nSample.mat)*dat$w.mat;
                    dat$ustat.mat <- (dat$z.mat)*(dat$vstat.mat);
                    dat.rp0 <- dat;
                    dat.rp0$ustat.mat <- rm.na(dat$ustat.mat);
                    dat.rp0$vstat.mat <- rm.na(dat$vstat.mat);
                    dat.rp0$nSample.mat <- rm.na(dat$nSample.mat);
                    res.cond <- getCondUV(dat=dat.rp0,lambda=extraPar$lambda,ix.candidate=ix.candidate,ix.known=ix.known,r2=r2);
                    cond.ok <- 1;
                }
                if(method=="impZ") {
                    dat$ustat.mat.ori <- dat$ustat.mat;
                    dat$vstat.mat.ori <- dat$vstat.mat;
                    
                    dat$z.mat <- dat$ustat.mat/dat$vstat.mat;
                    dat$vstat.mat <- sqrt(2*dat$af.mat*(1-dat$af.mat)*dat$nSample.mat)*dat$w.mat;
                    dat$ustat.mat <- (dat$z.mat)*(dat$vstat.mat);
                    dat.impz <- dat;
                    af.vec <- rowSums(dat$af.mat*dat$nSample.mat,na.rm=TRUE)/rowSums(dat$nSample.mat,na.rm=TRUE);
                    n.vec <- apply(dat$nSample,2,median,na.rm=TRUE);
                    
                    dat.impz$z.mat <- dat$ustat.mat/dat$vstat.mat;
                    for(ss in 1:ncol(dat$vstat.mat)) {
                        
                        ix.missing.ss <- which(is.na(dat$ustat.mat[,ss]));
                        if(length(ix.missing.ss)>0) {
                            zz <- as.matrix(r2[ix.missing.ss,ix.missing.ss]);
                            zx <- matrix(r2[ix.missing.ss,-ix.missing.ss],nrow=length(ix.missing.ss),ncol=nrow(r2)-length(ix.missing.ss));
                            dat.impz$z.mat[ix.missing.ss,ss] <- ginv(zz)%*%zx%*%matrix(dat.impz$z.mat[-ix.missing.ss,ss],ncol=1);
                            dat.impz$vstat.mat[ix.missing.ss,ss] <- sqrt(2*af.vec[ix.missing.ss]*(1-af.vec[ix.missing.ss])*n.vec[ss]);
                        }
                    }
                    dat.impz$ustat.mat <- (dat.impz$z.mat)*(dat.impz$vstat.mat);
                    res.cond <- getCondUV(dat=dat.impz,lambda=extraPar$lambda,ix.candidate=ix.candidate,ix.known=ix.known,r2=r2);
                    cond.ok <- 1;
                  
                }
                if(method=="COJO") {
                    w.mat <- sqrt(dat$nSample*rm.na(dat$af.mat)*(1-rm.na(dat$af.mat)))*dat$w.mat;
                    z.mat <- dat$ustat.mat/dat$vstat.mat;
                    
                    statistic.meta <- rowSums(w.mat*z.mat,na.rm=TRUE)/sqrt(rowSums(w.mat^2,na.rm=TRUE));
                    n.vec <- rowSums(dat$nSample,na.rm=TRUE);
                    af.vec <- rowSums(dat$af.mat*dat$nSample.mat,na.rm=TRUE)/rowSums(dat$nSample.mat,na.rm=TRUE);
                    ustat.new <- rm.na(statistic.meta*sqrt(2*n.vec*af.vec*(1-af.vec)));
                    X.T.times.X <- matrix(0,nrow=length(n.vec),ncol=length(n.vec));
                    
                    for(mm1 in 1:nrow(dat$ustat.mat)) {
                        for(mm2 in 1:nrow(dat$ustat.mat)) {
                            X.T.times.X[mm1,mm2] <- rm.na(sqrt(2*af.vec[mm1]*(1-af.vec[mm1]))*sqrt(2*af.vec[mm2]*(1-af.vec[mm2]))*min(n.vec[mm1],n.vec[mm2],na.rm=TRUE)*r2[mm1,mm2]);
                        }
                    }
                    X.T.times.X <- regMat(X.T.times.X,extraPar$lambda);
                    res.cond <- get.conditional.score.stat(ustat.new,X.T.times.X,n.vec,ix.candidate,ix.known);
                    res.cond$conditional.V <- res.cond$conditional.V;
                    cond.ok <- 1;
                }
                if(method=="SYN+") {
                    dat$ustat.mat.ori <- dat$ustat.mat;
                    dat$vstat.mat.ori <- dat$vstat.mat;
                    dat$z.mat <- dat$ustat.mat/dat$vstat.mat;
                    dat$vstat.mat <- sqrt(2*dat$af.mat*(1-dat$af.mat)*dat$nSample.mat)*dat$w.mat;
                    dat$ustat.mat <- (dat$z.mat)*(dat$vstat.mat);
                    ix.tmp <- c(ix.candidate,ix.known);
                    ix.candidate.tmp <- 1:length(ix.candidate);
                    ix.known.tmp <- (length(ix.candidate)+1):length(ix.tmp);
                    ix.na <- which(rowSums(!is.na(dat$ustat.mat[ix.tmp,]))==0);
                    r2 <- regMat(r2,extraPar$lambda);
                    if(length(ix.na)==0) {
                        res.cond <- synthesis(matrix(dat$ustat.mat[ix.tmp,],nrow=length(ix.tmp)),matrix((dat$vstat.mat[ix.tmp,])^2,nrow=length(ix.tmp)),as.matrix(r2[ix.tmp,ix.tmp]));
                        res.cond$conditional.beta.est <- res.cond$beta.syn[ix.candidate.tmp];
                        res.cond$conditional.beta.var <- as.matrix(res.cond$var.syn[ix.candidate.tmp,ix.candidate.tmp]);
                        res.cond$conditional.ustat <- (res.cond$conditional.beta.est[ix.candidate.tmp])/(res.cond$conditional.beta.var[ix.candidate.tmp,ix.candidate.tmp]);
                        res.cond$conditional.V <- (as.matrix(1/(res.cond$conditional.beta.var[ix.candidate.tmp,ix.candidate.tmp])));
                        
                        res.cond$numStability <- 1;
                        res.cond$nSample <- rowSums(dat$nSample.mat,na.rm=TRUE);
                        cond.ok <- 1;
                    }
                }
                
                if(cond.ok==1) {
                    statistic[varId.start:varId.end] <- (res.cond$conditional.ustat)^2/diag(res.cond$conditional.V);
                    marginal.statistic[varId.start:varId.end] <- (sum(dat$ustat.mat[ix.candidate,],na.rm=TRUE))^2/sum((dat$vstat.mat[ix.candidate,])^2,na.rm=TRUE);
                    marginal.p.value[varId.start:varId.end] <- pchisq(marginal.statistic[varId.start:varId.end],df=1,lower.tail=FALSE)
                    p.value[varId.start:varId.end] <- pchisq(statistic[varId.start:varId.end],df=1,lower.tail=F);
                    beta.est[varId.start:varId.end] <- res.cond$conditional.beta.est;
                    beta.se[varId.start:varId.end] <- sqrt(diag(res.cond$conditional.beta.var));
                    
                    nSample[varId.start:varId.end] <- res.cond$nSample[ix.candidate];

                    for(ww in 1:length(ix.candidate)) {    
                        ref.tab[(varId.start:varId.end)[ww]] <- apply(matrix(dat$ref.mat[ix.candidate[ww],],nrow=1),1,uniq.allele);
                        alt.tab[(varId.start:varId.end)[ww]] <- apply(matrix(dat$alt.mat[ix.candidate[ww],],nrow=1),1,uniq.allele);
                    }
                    
                    pos.all[varId.start:varId.end] <- pos[ix.candidate];
                    pos.out[varId.start:varId.end] <- gsub("_.*","",pos.all[varId.start:varId.end]);
                    af.meta <- rowSums((dat$af.mat)*(dat$nSample.mat),na.rm=T)/rowSums((dat$nSample.mat),na.rm=T);
                    af[varId.start:varId.end] <- af.meta[ix.candidate];
                    pos.ma <- unique(pos.out[varId.start:varId.end]);
                    if(extraPar$trans.ethnic==TRUE) {
                        chisq.mat[varId.start:varId.end,] <- res.cond$chisq.mat;
                        pval.memo[varId.start:varId.end] <- res.cond$pval.memo;
                        pval.re[varId.start:varId.end] <- res.cond$pval.re
                    }
                    if(length(pos.ma)>1) {
                        tmp <- estimateMAeffect.meta(res.cond$conditional.ustat,sqrt(diag(res.cond$conditional.V)),boot,nBoot=100);
                        statistic[varId.start:varId.end] <- tmp$statistic;
                        p.value[varId.start:varId.end] <- tmp$p.value;
                        beta.est[varId.start:varId.end] <- tmp$beta.est;
                        beta.se[varId.start:varId.end] <- tmp$beta.se;
                    }
                }
                varId <- varId.end;
            }
        }
    }
    chisq.mat <- chisq.mat[1:varId.end,];
    if(extraPar$trans.ethnic==FALSE) {
        res.formatted <- cbind(pos.out,
                               ref.tab,
                               alt.tab,
                               format(af,digits=3),
                               format(statistic,digits=3),
                               format(p.value,digits=3),
                               format(beta.est,digits=3),
                               format(beta.se,digits=3),
                               nSample,
                               numStability,
                               refAvail);
        colnames(res.formatted) <- c("POS","REF","ALT","AF","STAT","PVALUE","BETA","SD","N","numStability",'refAvailability');
    }
    
    if(extraPar$trans.ethnic==TRUE) {
        res.formatted <- cbind(pos.out,
                               ref.tab,
                               alt.tab,
                               format(af,digits=3),
                               format(statistic,digits=3),
                               format(p.value,digits=3),
                               format(beta.est,digits=3),
                               format(beta.se,digits=3),
                               format(pval.memo,digits=3),
                               nSample,
                               numStability,
                               refAvail);
        colnames(res.formatted) <- c("POS","REF","ALT","AF","STAT","PVALUE","BETA","SD","PVALUE-MEMO","N","numStability",'refAvailability');
    }
    
    
    
    
    return(list(res.formatted=res.formatted,
                dat=dat,
                r2=r2,
                r2.list=r2.list,
                refSampleId=geno.list$sampleId,
                chisq.mat=chisq.mat,
                pval.memo=pval.memo,
                marginal.statistic=marginal.statistic,
                marginal.p.value=marginal.p.value,
                raw.data.all=raw.data.all,
                res.cond=res.cond));
}

#' fisher z transform;
#'
#' @param r2
#' @param n.r2
#' @return a list with z and its std;
#' @export
fisher.z <- function(r2,n.r2) {
    z <- 1/2*log((1+r2)/(1-r2));
    z.se <- 1/sqrt(n.r2-3);
    return(list(z=z,
                z.se=z.se));
}
#' inverse fisher z transform;
#'
#' @param z
#' @return a list with z and its std;
#' @export
inv.fisher.z <- function(z) {
    r2 <- (exp(2*z)-1)/(exp(2*z)+1);
    return(r2);
}

#' wrapper of getCondUV.sand for trans-ethnic and pooled analysis;
#'
#' @export 
getCondUV.wrapper <- function(dat,lambda,ix.candidate,ix.known,r2,sandwich,boot=TRUE,nBoot=100,r2.list,study.ancestry,trans.ethnic,af.pca=NULL,cor.z2=NULL,otho=TRUE) {
    nSample.meta <- rowSums(dat$nSample.mat,na.rm=TRUE);
    if(!is.null(af.pca)) {
        af.pca <- matrix(unlist(af.pca),nrow=nrow(af.pca),ncol=ncol(af.pca));
    }
    if(trans.ethnic==FALSE) {
        if(sandwich==TRUE) {
            return(getCondUV.sand(dat=dat,lambda=lambda,ix.candidate=ix.candidate,ix.known=ix.known,r2=r2,sandwich=TRUE,boot=TRUE,nBoot=nBoot));
            
        }
        if(sandwich==FALSE) {
            return(getCondUV(dat=dat,lambda=lambda,ix.candidate=ix.candidate,ix.known=ix.known,r2=r2));
        }
    }

    if(trans.ethnic==TRUE) {
        conditional.ustat <- 0;
        conditional.V <- 0;
        uniq.ancestry <- unique(study.ancestry);
        conditional.ustat.mat <- matrix(NA,nrow=length(ix.candidate),ncol=length(uniq.ancestry));
        conditional.V.mat <- conditional.ustat.mat;
        af.pca.avg <- matrix(nrow=length(uniq.ancestry),ncol=ncol(af.pca));
        
        for(ee in 1:length(uniq.ancestry)) {
            ix.ee <- which(names(r2.list)==uniq.ancestry[ee]);
            r2.ee <- r2.list[[ix.ee]];
            ix.study <- which(study.ancestry==uniq.ancestry[ee]);
 
            af.pca.avg[ee,] <- colMeans(matrix(af.pca[ix.study,],nrow=length(ix.study)),na.rm=TRUE);
            
            dat.ee <- list(ustat.mat=as.matrix(dat$ustat.mat[,ix.study]),
                           vstat.mat=as.matrix(dat$vstat.mat[,ix.study]),
                           nSample.mat=as.matrix(dat$nSample.mat[,ix.study]));
            refAvail.ee <- as.numeric(sum(is.na(r2.ee[ix.candidate,ix.known]))==0);

            if(sandwich==TRUE & refAvail.ee==0) {
                res.cond.ee <- getCondUV.sand(dat=dat.ee,lambda=lambda,ix.candidate=ix.candidate,ix.known=ix.known,r2=r2.ee,boot=TRUE,nBoot=nBoot);
            }
            if(sandwich==FALSE | refAvail.ee==1) {
                res.cond.ee <- getCondUV(dat=dat.ee,lambda=lambda,ix.candidate=ix.candidate,ix.known=ix.known,r2=r2.ee);
            }
            
            
            conditional.ustat <- conditional.ustat+rm.na(res.cond.ee$conditional.ustat);
            conditional.V <- conditional.V+rm.na(res.cond.ee$conditional.V);
            conditional.ustat.mat[,ee] <- rm.na(res.cond.ee$conditional.ustat);
            conditional.V.mat[,ee] <- rm.na(res.cond.ee$conditional.V);
                
        }
        conditional.beta.mat <- conditional.ustat.mat/conditional.V.mat;
        conditional.beta.se.mat <- sqrt(1/conditional.V.mat);
        res.trans.ethnic <- trans.ethnic.meta.core(conditional.beta.mat,conditional.beta.se.mat,af.pca.avg,cor.z2,otho,re=TRUE);
        conditional.beta.est <- conditional.ustat/conditional.V;
        conditional.beta.var <- 1/conditional.V;
        return(list(conditional.ustat=conditional.ustat,
                    conditional.V=conditional.V,
                    conditional.beta.est=conditional.beta.est,
                    nSample=nSample.meta,
                    chisq.mat=res.trans.ethnic$chisq.mat,
                    pval.re=res.trans.ethnic$pval.re,
                    pval.memo=res.trans.ethnic$p.value.memo,
                    conditional.beta.var=conditional.beta.var));
    }
        


 
}




#' get conditional U and V for sandwich estimator;
#'
#' @param ix.candidate candidate variant indices;
#' @param ix.known conditioned variant indices
#' @param dat the list with contributed summary stat;
#' @param r2 LD matrices;
#' @param n.r2 sample size for the ref panel; 
#' @return a list with condiitonal association analysis statistics;
#' @export 
getCondUV.sand <- function(...) {
    pars <- list(...);
    boot <- pars$boot;
    nBoot <- pars$nBoot;
    if(is.null(nBoot)) nBoot <- 100;
    if(is.null(boot)) boot <- FALSE;
    ix.candidate <- pars$ix.candidate;
    ix.known <- pars$ix.known;
    dat <- pars$dat;
    lambda <- pars$lambda;
    r2 <- pars$r2;
    
    n.ref <- pars$n.ref;
    method <- pars$method;
    if(is.null(method)) method <- "replace";
    if(is.null(lambda)) lambda <- 0.1;
    ustat.meta <- rowSums(dat$ustat.mat,na.rm=TRUE);
    vstat.sq.meta <- rowSums((dat$vstat.mat)^2,na.rm=TRUE);
    z.meta <- ustat.meta/sqrt(vstat.sq.meta);
    z.mat <- dat$ustat.mat/dat$vstat.mat;
    beta.meta <- ustat.meta/vstat.sq.meta;
    z.exp <- beta.meta*dat$vstat.mat;
    nSample.meta <- rowSums(dat$nSample.mat,na.rm=TRUE);
    
    conditional.ustat <- NA;
    conditional.V <- matrix(0,nrow=length(ix.candidate),ncol=length(ix.candidate));
    conditional.beta.var <- matrix(0,nrow=length(ix.candidate),ncol=length(ix.candidate));
    
    conditional.beta.est <- NA;
    statistic.joint.sand <- NA;
    pval.joint.sand <- NA;
    statistic.joint.std <- NA;
    pval.joint.std <- NA;
    pval.burden.std <- NA;
    statistic.burden.std <- NA;
    pval.burden.sand <- NA;
    pval.skat.sand <- NA;
    pval.skat.std <- NA;
    statistic.burden.sand <- NA;
    for(ii in 1:length(ix.candidate)) {
        ix.ii <- c(ix.candidate[ii],ix.known);
        r2.ii <- r2[ix.ii,ix.ii];
        V.ii <- matrix(0,nrow=length(ix.known)+1,ncol=length(ix.known)+1);
        V0.ii <- V.ii;
        z.centered <- rm.na(matrix((z.mat-z.exp)[ix.ii,],nrow=length(ix.known)+1));
        z.centered <- z.centered-rowMeans(z.centered);

        for(jj in 1:ncol(dat$ustat.mat)) {
            V.tmp <- r2.ii*rm.na(dat$vstat.mat[ix.ii,jj]);
            V.ii <- V.ii+t(V.tmp)*rm.na(dat$vstat.mat[ix.ii,jj]);
            V0.ii <- V0.ii+z.centered[,jj]%*%t(z.centered[,jj]);
            
        }
        V0.ii <- V0.ii/ncol(dat$ustat.mat);
        diag(V0.ii) <- 1;
        r2.sand <- V0.ii;
        
        if(method=="combined") {
            V0.ii <- r2.meta;
        }
        if(method=='sandwich') {
            V0.ii <- r2.sand;
            ix.rp <- which(!is.na(r2.sand),arr.ind=TRUE);
        }
            
        if(method=="ref") {
            V0.ii <- rm.na(r2);
        }
        if(method=='replace') {
            ix.na <- which(is.na(r2.ii),arr.ind=TRUE);
            r2.tmp <- r2.ii;
            r2.tmp[ix.na] <- V0.ii[ix.na];
            V0.ii <- r2.tmp;
            ix.rp <- ix.na;
        }

        var.boot <- 0;
        if(boot==TRUE) {
            dat.boot <- list(ustat.mat=matrix(dat$ustat.mat[ix.ii,],nrow=length(ix.ii)),
                             vstat.mat=matrix(dat$vstat.mat[ix.ii,],nrow=length(ix.ii)));
            var.boot <- bootVar(dat=dat.boot,r2=V0.ii,nBoot=nBoot,ix.candidate=1,ix.known=2:length(ix.ii),N=median(nSample.meta[which(nSample.meta>0)],na.rm=TRUE),ix.rp=ix.rp);
        }

        r2.final <- V0.ii;

        sd.mat <- matrix(0,nrow=length(ix.known)+1,ncol=length(ix.known)+1);
        diag(sd.mat) <- sqrt(vstat.sq.meta[ix.ii]);
        V0.ii <- sd.mat%*%V0.ii%*%sd.mat;
        if(boot==FALSE) var.boot <- V0.ii;
        beta.ii.joint <- ginv(V0.ii)%*%ustat.meta[ix.ii];
        beta.ii.sd <- sqrt(ginv(V0.ii)[1,1])
        
        res.cond <- try(get.conditional.score.stat(ustat.meta[ix.ii],V0.ii,median(nSample.meta[which(nSample.meta>0)]),1,2:length(ix.ii)),silent=TRUE);
        if(class(res.cond)!='try-error') {
            conditional.ustat[ii] <- res.cond$conditional.ustat;
            conditional.V[ii,ii] <- var.boot+as.numeric(res.cond$conditional.V);
            conditional.beta.est[ii] <- res.cond$conditional.beta.est;
            conditional.beta.var[ii] <- res.cond$conditional.beta.var;
        }
    }   
    return(list(conditional.ustat=rm.na(conditional.ustat),
                conditional.V=rm.na(conditional.V),
                z.mat=z.mat,
                z.centered=z.centered,
                z.exp=z.exp,
                r2.final=r2.final,
                r2.sand=r2.sand,
                numStability=rep(1,length(ix.candidate)),
                conditional.beta.est=rm.na(conditional.beta.est),
                conditional.beta.var=rm.na(conditional.beta.var),
                nSample=nSample.meta));
}

#' bootVar
#'
#' @param r2
#' @param dat the summary statistic calculated from the data;
#' @param nBoot default is 100;
#' @param ix.candidate the index of the candidate variants;
#' @param ix.known the index of the known variants; 
#' @return estimated conditional variance
#' @export
bootVar <- function(dat,r2,nBoot=1000,ix.candidate,ix.known,N,ix.rp) {
    nStudy <- ncol(dat$ustat.mat);
    nVar <- nrow(dat$ustat.mat);
    ustat.meta <- rowSums(dat$ustat.mat,na.rm=TRUE);
    vstat.sq.meta <- rowSums((dat$vstat.mat)^2,na.rm=TRUE);
    sd.mat <- matrix(0,nrow=nVar,ncol=nVar);
    diag(sd.mat) <- sqrt(vstat.sq.meta);
    z.simu <- rmvnorm(nStudy*nBoot,mean=rep(0,nVar),sigma=r2);
    conditional.ustat.boot <- 0;
    conditional.V.boot <- 0;
    for(ii in 1:nBoot) {
        z.ii <- t(matrix(z.simu[((ii-1)*nStudy+1):(ii*nStudy),],nrow=nStudy));
        ustat.meta.ii <- rowSums(z.ii*(dat$vstat.mat),na.rm=TRUE);
        
        cor.ii.tmp <- z.ii%*%t(z.ii)/nStudy;
        cor.ii <- r2;
        cor.ii[ix.rp] <- cor.ii.tmp[ix.rp];
        V.ii <- sd.mat%*%cor.ii%*%sd.mat;
        res.cond.ii <- try(get.conditional.score.stat(ustat.meta,V.ii,N,ix.candidate,ix.known),silent=TRUE);
        if(class(res.cond.ii)!='try-error') {
            conditional.ustat.boot[ii] <- res.cond.ii$conditional.ustat;
            conditional.V.boot[ii] <- as.numeric(res.cond.ii$conditional.V);
        }
    }
    return(var(conditional.ustat.boot,na.rm=TRUE));
}

#' get conditional U and V; implement our PCBS method for conditional association analyses;
#'
#' @param ix.canddidate candidate variant indices;
#' @param ix.known conditioned variant indices
#' @param dat the list with contributed summary stat;
#' @param r2 LD matrices;
#' @param sandwich use sandwich estimator or not; default is not to apply sandwich estimator;
#' @return a list with condiitonal association analysis statistics;
#' @export 
getCondUV <- function(...) {
    numStability <- 1;
    pars <- list(...);
    ix.candidate <- pars$ix.candidate;
    ix.known <- pars$ix.known;
    dat <- pars$dat;
    lambda <- pars$lambda;
    r2 <- pars$r2;
    sandwich <- pars$sandwich;
    if(is.null(sandwich)) sandwich <- FALSE;
    
    if(is.null(lambda)) lambda <- 0.1;
    ustat.meta <- rowSums(dat$ustat.mat,na.rm=TRUE);
    vstat.sq.meta <- rowSums((dat$vstat.mat)^2,na.rm=TRUE);
    nSample.meta <- rowSums(dat$nSample.mat,na.rm=TRUE);
    nSample.meta <- nSample.meta;
    
    V <- 0;
    N.V <- 0;
    V0 <- 0;
    for(jj in 1:ncol(dat$ustat.mat)) {
        
        vstat <- (rm.na(dat$vstat.mat[,jj]));
        V.ii <- r2*vstat;

        V <- V+t(V.ii)*vstat;
        V0 <- V0+rm.na(dat$ustat.mat[,jj]%*%t(dat$ustat.mat[,jj]))
        N.V <- N.V+rm.na(sqrt(dat$nSample.mat[,jj]))%*%t(rm.na(sqrt(dat$nSample.mat[,jj])))
    }
    covG <- V/N.V;
    covG <- rm.na(covG);
    covG <- regMat(covG,lambda);


    
    U.meta <- ustat.meta/nSample.meta;
    U.meta <- rm.na(U.meta);
    U.XY <- U.meta[ix.candidate];
    U.ZY <- U.meta[ix.known];
        
    V.XZ <- matrix(covG[ix.candidate,ix.known],nrow=length(ix.candidate),ncol=length(ix.known));
    V.ZZ <- matrix(covG[ix.known,ix.known],nrow=length(ix.known),ncol=length(ix.known));
    V.XX <- matrix(covG[ix.candidate,ix.candidate],nrow=length(ix.candidate),ncol=length(ix.candidate));
    conditional.ustat <- U.XY-V.XZ%*%ginv(V.ZZ)%*%U.ZY;
    
    beta.ZY <- ginv(V.ZZ)%*%U.ZY;
    scaleMat <- as.matrix(as.numeric(diag(as.matrix(N.V)))%*%t(as.numeric(diag(as.matrix(N.V)))))
    calc.varU <- function(V,ix.candidate,ix.known,scaleMat) {
        var.U.XY <- rm.na(V[ix.candidate,ix.candidate]/(scaleMat[ix.candidate,ix.candidate]));
        var.U.ZY <- rm.na(V[ix.known,ix.known]/(scaleMat[ix.known,ix.known]))
        cov.U.XY.U.ZY <- rm.na(V[ix.candidate,ix.known]/matrix(scaleMat[ix.candidate,ix.known],nrow=length(ix.candidate),ncol=length(ix.known)));
        return(list(var.U.XY=var.U.XY,
                    var.U.ZY=var.U.ZY,
                    cov.U.XY.U.ZY=cov.U.XY.U.ZY));
    }
    tmp <- calc.varU(V,ix.candidate,ix.known,scaleMat);
    var.U.XY <- tmp$var.U.XY;
    var.U.ZY <- tmp$var.U.ZY;
    cov.U.XY.U.ZY <- tmp$cov.U.XY.U.ZY;
    
    conditional.V <- var.U.XY+V.XZ%*%ginv(V.ZZ)%*%var.U.ZY%*%ginv(V.ZZ)%*%t(V.XZ)-cov.U.XY.U.ZY%*%t(V.XZ%*%ginv(V.ZZ))-(V.XZ%*%ginv(V.ZZ))%*%t(cov.U.XY.U.ZY);
    sigma.sq.est <- abs(1-(t(U.ZY)%*%ginv(V.ZZ)%*%U.ZY));
    conditional.V <- conditional.V*as.numeric(sigma.sq.est);

    if(det(conditional.V)<0) numStability <- 0;
    conditional.V <- regMat(conditional.V,lambda);
    conditional.beta.est <- ginv(V.XX-V.XZ%*%ginv(V.ZZ)%*%t(V.XZ))%*%conditional.ustat;
    conditional.beta.var <- ginv(V.XX-V.XZ%*%ginv(V.ZZ)%*%t(V.XZ))%*%conditional.V%*%ginv(V.XX-V.XZ%*%ginv(V.ZZ)%*%t(V.XZ));
    conditional.ustat <- ginv(conditional.beta.var)%*%conditional.beta.est;
    conditional.V <- ginv(conditional.beta.var);
    return(list(conditional.ustat=conditional.ustat,
                conditional.V=conditional.V,
                numStability=numStability,
                conditional.beta.est=conditional.beta.est,
                conditional.beta.var=conditional.beta.var,
                nSample=nSample.meta));
    
}
    
synthesis <- function(U,V,ref.cor) {
    no.var <- ncol(ref.cor)
    nstudy <- ncol(U)
    b.syn.list <- list(); cov.syn.list <- list(); w1.list <- list()
    w1 <- diag(no.var)
    v.dim.start <- 1;v.dim.end <- 1; ix.rm.list <- list()
    for(jj in 1:nstudy) {ix.rm.list[[jj]] <- which(is.na(U[,jj]))}
    for(jj in 1:nstudy){
        w1 <- diag(no.var)
        ix.rm <- ix.rm.list[[jj]]; u <- U[,jj]; v <- V[,jj]
        if(jj==1) v.dim.end[jj] <- no.var-length(ix.rm)
        if(jj>1) v.dim.end[jj] <- v.dim.end[jj-1]+no.var-length(ix.rm)
        v.dim.start[jj+1] <- v.dim.end[jj]+1
        if(length(ix.rm)==no.var) next
        cov.syn <- diag(sqrt(v)) %*% ref.cor %*% diag(sqrt(v))
        if(length(ix.rm)>0) {
            w1 <- matrix(w1[-ix.rm,], nrow=(no.var-(length(ix.rm))))
            XZ <- matrix(ref.cor[-ix.rm,ix.rm],nrow=no.var-length(ix.rm),ncol=length(ix.rm))
            XX <- as.matrix(ref.cor[-ix.rm,-ix.rm]);
            w1[,ix.rm] <- ginv(XX)%*%XZ
            cov.syn <- cov.syn[-ix.rm,-ix.rm]
        }
        cov.syn.list[[jj]] <- ginv(cov.syn)
        w1.list[[jj]] <- w1
        b.syn.list[[jj]] <- cov.syn.list[[jj]] %*% u[!is.na(u)]
    }
    w1 <- w1.list[[1]]
    v <- matrix(0,nrow=v.dim.end[nstudy],ncol=v.dim.end[nstudy])
    for(jj in 1:nstudy) {
        if(v.dim.start[jj]!=v.dim.start[jj+1]){
        if(jj>1) {w1 <- rbind(w1,w1.list[[jj]]);}
        v[(v.dim.start[jj]):(v.dim.end[jj]),(v.dim.start[jj]):(v.dim.end[jj])] <- cov.syn.list[[jj]];}
    }
    b <- unlist(b.syn.list)
    beta.syn <- (ginv(t(w1)%*%ginv(v)%*%(w1))%*%(t(w1)%*%ginv(v)%*%b))
    var.syn <- (ginv(t(w1)%*%ginv(v)%*%(w1)))
    return(list(beta.syn=beta.syn,var.syn=var.syn))
}


