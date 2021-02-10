#' get imputation quality from tables based upon imputation quality files;
#'
#' @param imp.qual the list with imputation qualities
#' @param pos A vector of positions where the imputatoin quality will be retrieved;
#' @param col.impqual The column where the imputation quality will be retrieved;
#' @return a matrix with imputaton qualities
#' @export
getImpQual <- function(imp.qual,pos,col.impqual,rmMultiAllelicSite=TRUE) {
    if(rmMultiAllelicSite==TRUE) {
        pos.imp <- try(paste(imp.qual[,1],imp.qual[,2],sep=":"),silent=TRUE);
        if(class(pos.imp)=='try-error') return(rep(NA,length(pos)))
        ix <- match(gsub("_.*","",pos),pos.imp);
    }
    if(rmMultiAllelicSite==FALSE) {
        pos.imp <- try(paste(paste(imp.qual[,1],imp.qual[,2],sep=":"),paste(imp.qual[,3],imp.qual[,4],sep="/"),sep="_"),silent=TRUE);
        if(class(pos.imp)=='try-error') return(rep(NA,length(pos)))
        ix <- match(pos,pos.imp);
    }

    
    return(as.numeric(imp.qual[ix,col.impqual]));
}

#' unique alternative alleles;
#'
#' @param x A vector of allele labels;
#' @return A vector of unique alleles;
#' @export
uniq.allele <- function(x) {x.tab <- table(x);return(paste(names(x.tab),sep=',',collapse=','))}

#' single variant meta-analysis integrating imputation quality;
#'
#' @param score.stat.file the file names of score statistic files;
#' @param imp.qual.file the file names of imputation quality;
#' @param tabix.range the tabix range. IT must be in quote and provided as a string;
#' @param alternative The alternative hypothesis. Default is two.sided;
#' @param col.impqual The column number for the imputation quality score;
#' @param impQual.lb The lower bound for the imputation quality. Variants with imputaiton quality less than impQual.lb will be labelled as missing;
#' @param impQualWeight Using imputation quality as weight
#' @param rmMultiAllelicSite Default is TRUE. Multi-allelic sites will be removed from the analyses if set TRUE, and a variable posMulti will be output; The variant site with multiple alleles can be analyzed using rareGWAMA.single.multiAllele function;
#' @param weight The weights used in meta-analysis; the default choice is Npq+impQ, which is a weighting scheme that takes into account both the imputation quality and allele freq difference between studies; 
#' @return A list of analysis results;
#' @export 
rareGWAMA.single <- function(score.stat.file,imp.qual.file=NULL,tabix.range,alternative="two.sided",col.impqual=5,impQual.lb=0.7,impQualWeight=FALSE,rmMultiAllelicSite=FALSE,gc=FALSE,...) {
    a <- Sys.time();
    extraPar <- list(...);
    trans.ethnic <- extraPar$trans.ethnic;
    if(is.null(extraPar$memo.recalibrate)) extraPar$memo.recalibrate <- FALSE;
    if(is.null(extraPar$re)) extraPar$re <- FALSE;
    
    if(is.null(trans.ethnic)) {
        trans.ethnic <- FALSE;
    }
    if(trans.ethnic==TRUE) {
        if(is.null(extraPar$af.pca)) {
            stop('study specific PCs (or MDS) are needed for trans-ethnic analysis');
        }
        z2.simu <- matrix((rnorm(1e6*ncol(extraPar$af.pca)))^2,ncol=ncol(extraPar$af.pca));
        
        z2.simu.cumsum <- matrix(apply(z2.simu,1,cumsum),ncol=ncol(extraPar$af.pca),byrow=TRUE);
        p.simu.cumsum <- matrix(nrow=nrow(z2.simu.cumsum),ncol=ncol(z2.simu.cumsum));
        for(ii in 1:ncol(z2.simu.cumsum)) p.simu.cumsum[,ii] <- pchisq(z2.simu.cumsum[,ii],df=ii,lower.tail=FALSE)
        minp.vec <- apply(p.simu.cumsum,1,min);

        invZ.p.simu.cumsum <- qnorm(p.simu.cumsum,lower.tail=FALSE);
        cor.z2 <- cor(invZ.p.simu.cumsum);
        cor.z2.list <- list();
        
        
    }
    weight <- extraPar$weight;
    if(is.null(weight)) weight <- "Npq+impQ";
    lambda <- matrix(rep(1,length(score.stat.file)),ncol=1);
    maf.bin <- matrix(c(0,1),ncol=2);
    var.y <- extraPar$var.y;
    if(is.null(var.y)) var.y <- 1;
    binaryTrait <- extraPar$binaryTrait;
    if(is.null(binaryTrait)) binaryTrait <- FALSE;
    if(gc==TRUE) {
        lambda <- extraPar$lambda;
        if(is.null(lambda)) {
            lambda <- matrix(1,nrow=length(score.stat.file),ncol=nrow(maf.bin));
        }
        maf.bin <- extraPar$maf.bin;
        if(is.null(maf.bin)) stop("maf bin must be provided");
    }
    capture.output(raw.data.all <- rvmeta.readDataByRange( score.stat.file, NULL, tabix.range, multiAllelic = TRUE));
    raw.imp.qual <- NULL;
    if(!is.null(imp.qual.file))
        raw.imp.qual <- lapply(imp.qual.file,tabix.read.table,tabixRange=tabix.range);
    time.readData <- Sys.time()-a;
    b <- Sys.time();
    raw.data.all <- raw.data.all[[1]];
    
    cat('Read in',length(unique(gsub("_.*","",raw.data.all$pos))),'variants',sep=' ');
    dat <- GWAMA.formatData(raw.data.all,raw.imp.qual,impQualWeight,impQual.lb,col.impqual,rmMultiAllelicSite=rmMultiAllelicSite);
    dat.withMulti <- dat;
    tmp <- GWAMA.rmMulti(dat);
    dat <- tmp$dat;posMulti <- tmp$posMulti;
    if(nrow(dat$ustat.mat)==0) stop('no variants');
    direction.meta <- apply(dat$direction.mat,1,paste,sep="",collapse="");   
    maf.meta <- rowSums((dat$af.mat)*(dat$nSample.mat),na.rm=TRUE)/rowSums(dat$nSample.mat,na.rm=TRUE);
    maf.meta[which(maf.meta>0.5)] <- 1-maf.meta[which(maf.meta>0.5)];
    for(bb in 1:nrow(maf.bin)) {
        ix.bb <- which(maf.meta>=maf.bin[bb,1] & maf.meta<maf.bin[bb,2]);
        if(length(ix.bb)>0) {
            dat$ustat.mat[ix.bb,] <- scale(matrix(dat$ustat.mat[ix.bb,],nrow=length(ix.bb)),center=FALSE,scale=sqrt(lambda[,bb]));
        }
    }
    ustat.meta <- rowSums((dat$ustat.mat)*(dat$w.mat),na.rm=TRUE);
    vstat.sq.meta <- rowSums((dat$w.mat)^2*(dat$vstat.mat)^2,na.rm=TRUE);
    
    beta.meta <- (ustat.meta)/(vstat.sq.meta);
    beta.sd.meta <- sqrt(1/vstat.sq.meta);
    if(weight=="N") {
        w.mat <- sqrt(dat$nSample);
        z.mat <- dat$ustat.mat/dat$vstat.mat;
        statistic.meta <- (rowSums(w.mat*z.mat,na.rm=TRUE))^2/rowSums(w.mat^2,na.rm=TRUE)
    }
    sum.weight <- rowSums((dat$nSample)*(dat$w.mat)^2,na.rm=TRUE);
    if(weight=='Npq+impQ') {
        w.mat <- sqrt(dat$nSample*rm.na(dat$af.mat)*(1-rm.na(dat$af.mat)))*dat$w.mat;
        z.mat <- dat$ustat.mat/dat$vstat.mat;
        statistic.meta <- (rowSums(w.mat*z.mat,na.rm=TRUE))^2/rowSums(w.mat^2,na.rm=TRUE)
        cat('bi-allelic variant association completed\n');
        cat('a total of ',nrow(length(statistic.meta)), ' bi-allelic variants analyzed\n')
    }
    
    if(weight=='N+impQ') {
        w.mat <- sqrt(dat$nSample)*dat$w.mat;
        z.mat <- dat$ustat.mat/dat$vstat.mat;
        statistic.meta <- (rowSums(w.mat*z.mat,na.rm=TRUE))^2/rowSums(w.mat^2,na.rm=TRUE)
    }
    if(trans.ethnic==TRUE) {
        res.trans <- trans.ethnic.meta(dat,extraPar$af.pca,cor.z2,re=extraPar$re,recalibrate=extraPar$memo.recalibrate);
        statistic.omnibus.mrvt <- res.trans$statistic.omnibus.mrvt;
        p.value.omnibus.mrvt<- res.trans$p.value.omnibus.mrvt;
        no.pc.omnibus.mrvt <- res.trans$no.pc.omnibus.mrvt;
    }
    I2 <- NULL;cochranQ.pVal.mixChisq <- NULL;cochranQ.pVal <- NULL;cochranQ.df <- NULL;cochranQ.stat <- NULL;beta.byStudy.mat <- NULL;beta.var.byStudy.mat <- NULL;
    if(!is.null(extraPar$hetStat)) {
        if(extraPar$hetStat==TRUE) {
            w0.mat <- sqrt(dat$nSample*rm.na(dat$af.mat)*(1-rm.na(dat$af.mat)))*dat$w.mat;
            z.mat <- dat$ustat.mat/dat$vstat.mat;
            beta.byStudy.mat <- z.mat/w0.mat;
            beta.var.byStudy.mat <- 1/w0.mat^2;
            
            for(ix.var in 1:nrow(beta.byStudy.mat)) {
                beta.byStudy <- beta.byStudy.mat[ix.var,];
                beta.var.byStudy <- beta.var.byStudy.mat[ix.var,];
              
                beta.meta[ix.var] <- as.numeric(sum(beta.byStudy*(1/beta.var.byStudy))/sum(1/beta.var.byStudy));
                beta.sd.meta[ix.var] <- sqrt(sum(1/beta.var.byStudy)/(sum(1/beta.var.byStudy))^2);
                w.mat <- matrix(0,nrow=length(beta.byStudy),ncol=length(beta.byStudy));
                diag(w.mat) <- 1;
                weight.byStudy <- (w0.mat[ix.var,])^2;
                weight.byStudy <- weight.byStudy/sum(weight.byStudy,na.rm=TRUE);
                w.mat <- w.mat+rm.na(matrix((-1)*rep(weight.byStudy,length(beta.byStudy)),nrow=length(beta.byStudy),ncol=length(beta.byStudy),byrow=TRUE));
                cochranQ.stat.mixChisq <- t(w.mat%*%rm.na(beta.byStudy))%*%(w.mat%*%rm.na(beta.byStudy))
                v.mat <- matrix(0,nrow=length(beta.byStudy),ncol=length(beta.byStudy));
                diag(v.mat) <- sqrt(rm.na(beta.var.byStudy));
                cochranQ.stat[ix.var] <- sum((beta.byStudy-beta.meta[ix.var])^2/(beta.var.byStudy+(beta.sd.meta[ix.var])^2-2*weight.byStudy*beta.var.byStudy),na.rm=TRUE);
                cochranQ.df[ix.var] <- length(which(!is.na(beta.byStudy-beta.meta[ix.var])^2))-1;
                if(cochranQ.df[ix.var]>0)
                {
                    cochranQ.pVal[ix.var] <- pchisq(cochranQ.stat[ix.var],df=cochranQ.df[ix.var],lower.tail=FALSE);
                    svd.mat <- try(svd(v.mat%*%w.mat%*%w.mat%*%v.mat),silent=TRUE);
                    cochranQ.pVal.mixChisq[ix.var] <- cochranQ.pVal[ix.var];
                    if(class(svd.mat)!="try-error") {
                        lambda <- svd.mat$d;
                        cochranQ.pVal.mixChisq[ix.var] <- try(liu(cochranQ.stat.mixChisq,lambda),silent=TRUE);
                    }
                }
                if(cochranQ.df[ix.var]<=0)
                {
                    cochranQ.pVal[ix.var] <- NA;
                    cochranQ.pVal.mixChisq[ix.var] <- NA;
                }
                
                I2[ix.var] <- (cochranQ.stat[ix.var]-cochranQ.df[ix.var])/cochranQ.stat[ix.var];
            }
            I2[which(I2<0)] <- 0;
            I2 <- paste(myFormat(I2*100,digits=2),"%",sep="");
        }
    }
    if(weight=="optim") statistic.meta <- ustat.meta^2/vstat.sq.meta;
    p.value.meta <- pchisq(statistic.meta,df=1,lower.tail=FALSE);
    nSample.meta <- rowSums(dat$nSample.mat,na.rm=TRUE);
    af.meta <- rowSums((dat$nSample.mat)*(dat$af.mat),na.rm=TRUE)/rowSums((dat$nSample.mat),na.rm=TRUE);
    if(binaryTrait==FALSE) {
        beta.meta <- sign(ustat.meta)*sqrt(statistic.meta)*sqrt(var.y)/sqrt(2*nSample.meta*af.meta*(1-af.meta));
        beta.sd.meta <- sqrt(var.y/(2*nSample.meta*af.meta*(1-af.meta)));
    }
    if(binaryTrait==TRUE) {
        beta.meta <- sign(ustat.meta)*sqrt(statistic.meta)/sqrt(2*nSample.meta*af.meta*(1-af.meta))/sqrt(var.y);
        beta.sd.meta <- sqrt(1/sqrt(2*nSample.meta*af.meta*(1-af.meta)*var.y));
    }
        
    time.compAssoc <- Sys.time()-b;
    res.out <- cbind(gsub("_.*","",dat$pos),dat$ref.tab,dat$alt.tab,statistic.meta,p.value.meta,beta.meta,beta.sd.meta,nSample.meta,direction.meta);
    
    res.formatted <- cbind(gsub("_.*","",dat$pos),
                           dat$ref.tab,
                           dat$alt.tab,
                           format(af.meta,digits=3),
                           format(statistic.meta,digits=3),
                           format(p.value.meta,digits=3),
                           format(beta.meta,digits=3),
                           format(beta.sd.meta,digits=3),
                           nSample.meta,
                           direction.meta,
                           as.integer(sum.weight),
                           rep(-1,length(af.meta)));
    ix.0 <- which(p.value.meta==0);
    if(length(ix.0)>0) {
        pval.precise <- get.precise.pval(statistic.meta[ix.0]);
        res.formatted[ix.0,6] <- pval.precise;
    }
    
    
    colnames(res.formatted) <- c("POS","REF","ALT","AF","STAT","PVALUE","BETA","SD","N","DIRECTION","EFFECTIVE_N","numStability");
    gamma.est.mat <- NULL;
    gamma.sd.mat <- NULL;
    if(trans.ethnic==TRUE) {
        gamma.est.mat <- res.trans$gamma.est.mat;
        gamma.sd.mat <- res.trans$gamma.sd.mat;

        if(extraPar$re==FALSE) {
            res.formatted <- cbind(gsub("_.*","",dat$pos),
                                   dat$ref.tab,
                                   dat$alt.tab,
                                   format(af.meta,digits=3),
                                   format(statistic.meta,digits=3),
                                   format(p.value.meta,digits=3),
                                   format(beta.meta,digits=3),
                                   format(beta.sd.meta,digits=3),
                                   nSample.meta,
                                   direction.meta,
                                   as.integer(sum.weight),
                                   rep(-1,length(af.meta)),
                                   format(statistic.omnibus.mrvt,digits=3),
                                   format(p.value.omnibus.mrvt,digits=3),
                                   no.pc.omnibus.mrvt,
                                   format(res.trans$p.value.original,digits=3),
                                   res.trans$chisq.mat);
            colnames(res.formatted) <- c("POS","REF","ALT","AF","STAT","PVALUE","BETA","SD","N","DIRECTION","EFFECTIVE_N","numStability","STAT_MEMO","PVAL_MEMO","NO_ANCESTRY_PC_MEMO","PVALUE_MR_MEGA",paste0("CHISQ-",0:(ncol(res.trans$chisq.mat)-1),"PC"));
            
        }
        if(extraPar$re==TRUE) {
            res.formatted <- cbind(gsub("_.*","",dat$pos),
                                   dat$ref.tab,
                                   dat$alt.tab,
                                   format(af.meta,digits=3),
                                   format(statistic.meta,digits=3),
                                   format(p.value.meta,digits=3),
                                   format(beta.meta,digits=3),
                                   format(beta.sd.meta,digits=3),
                                   nSample.meta,
                                   direction.meta,
                                   as.integer(sum.weight),
                                   rep(-1,length(af.meta)),
                                   format(statistic.omnibus.mrvt,digits=3),
                                   format(p.value.omnibus.mrvt,digits=3),
                                   no.pc.omnibus.mrvt,
                                   format(res.trans$p.value.original,digits=3),
                                   res.trans$chisq.mat,
                                   res.trans$statistic.re,
                                   res.trans$pval.re);
            colnames(res.formatted) <- c("POS","REF","ALT","AF","STAT","PVALUE","BETA","SD","N","DIRECTION","EFFECTIVE_N","numStability","STAT_MEMO","PVAL_MEMO","NO_ANCESTRY_PC_MEMO","PVALUE_MR_MEGA",paste0("CHISQ-",0:(ncol(res.trans$chisq.mat)-1),"PC"),"STAT_RE","PVALUE_RE");
            
        }
        
    }
    pos.only <- gsub("_.*","",dat.withMulti$pos);
    res.formatted.multi <- matrix(nrow=0,ncol=ncol(res.formatted));
    res.collapse <- matrix(nrow=0,ncol=2);
    if(rmMultiAllelicSite==FALSE) {
        cat('start multi-allelic analysis\n');
        cat('a total of  ',length(dat$posMulti),' will be analyzed\n');
        for(ii in 1:length(dat$posMulti))  {
            if(ii/100==as.integer(ii/100)) {
                cat(ii,' multi-allelic variants analyzed\n');
                
            }
            ix.ii <- which(pos.only==dat$posMulti[ii]);
            if(length(ix.ii)>0) {
                cor.multi <- list();
                cor.multi$chrpos <- dat$posMulti[ii]
                
                cor.multi$ref <- paste(apply(matrix(dat.withMulti$ref.mat[ix.ii,],nrow=length(ix.ii)),1,uniq.allele),sep=',',collapse=',');
                
                cor.multi$alt <- paste(apply(matrix(dat.withMulti$alt.mat[ix.ii,],nrow=length(ix.ii)),1,uniq.allele),sep=',',collapse=',');
                
                af.meta <- rowSums((matrix(dat.withMulti$af.mat[ix.ii,],nrow=length(ix.ii)))*(matrix(dat.withMulti$nSample.mat[ix.ii,],nrow=length(ix.ii))),na.rm=TRUE)/rowSums(matrix(dat.withMulti$nSample.mat[ix.ii,],nrow=length(ix.ii)),na.rm=TRUE);
                cor.mat <- (-2)*as.vector(af.meta)%*%t(as.vector(af.meta));
                
                
                diag(cor.mat) <- 2*as.vector(af.meta)*(1-as.vector(af.meta));
                cor.mat <- rm.na(cov2cor(cor.mat));                
                diag(cor.mat) <- 1;
                cor.multi$cor <- paste(as.vector(cor.mat),sep=',',collapse=',')
                
                res.multi.ii <- (multiAlleleAssoc(dat$posMulti[ii],dat.withMulti,cor.multi,sandwich=extraPar$sandwich,boot=extraPar$boot));
                
                res.formatted <- rbind(res.formatted,res.multi.ii$res.formatted);
                res.collapse <- rbind(res.collapse,res.multi.ii$res.collapse);
            }
        }
    }
    
    return(list(res.out=res.out,
                res.formatted=res.formatted,
                res.collapse=res.collapse,
                I2=I2,
                gamma.est.mat=gamma.est.mat,
                gamma.sd.mat=gamma.sd.mat,
                cochranQ.pVal.mixChisq=cochranQ.pVal.mixChisq,
                cochranQ.pVal=cochranQ.pVal,
                cochranQ.df=cochranQ.df,
                cochranQ.stat=cochranQ.stat,
                beta.byStudy.mat=beta.byStudy.mat,
                beta.var.byStudy.mat=beta.var.byStudy.mat,
                formattedData=dat,
                formattedData.withMulti=dat.withMulti,
                posMulti=dat$posMulti,
                raw.data.all=raw.data.all,
                raw.imp.qual=raw.imp.qual,
                time.compAssoc=time.compAssoc,
                time.readData=time.readData));
}

#' estimate multi-allelic effect from summary stat;
#'
#' @param ustat ustat;
#' @param vstat ustat;
#' @param af af;
#' @param pos pos;
#' @export
estimateMAeffect <- function(ustat,vstat,af,pos,boot=FALSE,nBoot=100) {
    pos.ma <- unique(pos[duplicated(pos)]);
    statistic <- ustat/vstat;
    beta.est <- ustat/vstat^2;
    beta.se <- 1/vstat;
    p.value <- pchisq(statistic,df=1,lower.tail=FALSE);
    for(ii in 1:length(pos.ma)) {
        ix.ii <- which(pos==pos.ma[ii]);
        af.ii <- af[ix.ii];
        sd.mat <- matrix(0,nrow=length(ix.ii),ncol=length(ix.ii));
        diag(sd.mat) <- vstat[ix.ii];
        cor0.ii <- cov2cor(ustat[ix.ii]%*%t(ustat[ix.ii]));
        cov0.ii <- sd.mat%*%cor0.ii%*%sd.mat;
        beta.est[ix.ii] <- ginv(cov0.ii)%*%ustat[ix.ii];
        beta.se[ix.ii] <- sqrt(as.vector(diag(ginv(cov0.ii))));
        p.value[ix.ii] <- pchisq((beta.est[ix.ii]/beta.se[ix.ii])^2,df=1,lower.tail=FALSE);
        statistic[ix.ii] <- (beta.est[ix.ii]/beta.se[ix.ii])^2
        
    }
    return(list(statistic=statistic,
                p.value=p.value,
                beta.est=beta.est,
                r2=cor0.ii,
                beta.se=beta.se));

}

#' estimate multi-allelic effect from summary stat;
#'
#' @param ustat.mat ustat;
#' @param vstat.mat ustat;
#' @export
estimateMAeffect.meta <- function(ustat.mat,vstat.mat,boot=FALSE,nBoot=100) {
    nVar <- nrow(ustat.mat);
    nStudy <- ncol(ustat.mat);
    ustat.meta <- rowSums(ustat.mat,na.rm=TRUE);
    vstat.sq.meta <- rowSums(vstat.mat^2,na.rm=TRUE);
    vstat.meta <- sqrt(vstat.sq.meta);
    z.meta <- ustat.meta/sqrt(vstat.sq.meta);
    z.mat <- ustat.mat/vstat.mat;
    beta.meta <- ustat.meta/vstat.sq.meta;
    z.exp <- beta.meta*vstat.mat;
    
    z.centered <- rm.na(matrix((z.mat-z.exp),nrow=nVar));
    z.centered <- z.centered-rowMeans(z.centered);
    V0.ii <- 0;
    for(jj in 1:ncol(ustat.mat)) {
        V0.ii <- V0.ii+z.centered[,jj]%*%t(z.centered[,jj]);
    }
    V0.ii <- V0.ii/ncol(ustat.mat);
    diag(V0.ii) <- 1;
    r2.sand <- V0.ii;
    sd.mat <- matrix(0,nrow=nVar,ncol=nVar);
    diag(sd.mat) <- vstat.meta;
    V0.ii <- rm.na(sd.mat%*%V0.ii%*%sd.mat);
    ix.inf <- which(V0.ii==Inf | V0.ii==-Inf, arr.ind=TRUE);
    if(length(ix.inf)>0) V0.ii[ix.inf] <- 0;
    beta.ii.joint <- ginv(V0.ii)%*%ustat.meta;
    boot.var <- rep(0,length(beta.ii.joint));

    if(boot==TRUE) {
        boot.var <- bootVar.multiallelic(ustat.meta,vstat.meta,r2.sand,nStudy,nBoot=100);
    }
    beta.ii.sd <- sqrt(diag(ginv(V0.ii))+boot.var);
    statistic <- (beta.ii.joint/beta.ii.sd)^2;
    p.value <- pchisq(statistic,df=1,lower.tail=FALSE);
    
    return(list(beta.est=beta.ii.joint,
                beta.sd=beta.ii.sd,
                statistic=statistic,
                p.value=p.value,
                r2=r2.sand));
}


#' Trans-ethnic meta-analysis with simple interface;
#'
#' @param beta.mat genetic effect estimates (one row per variant, and one column per study);
#' @param se.mat matrix of standard errors for the beta.mat
#' @param af.pca PCA based upon allele frequencies
#' @param cor.z2 correlation of z2;
#' @param otho if to othorgonalize the af pca in case studies contain missing data;
#' @param re if to include random effect to capture nonancestral heterogeniety
#' @export
trans.ethnic.meta.core <- function(beta.mat,se.mat,af.pca,cor.z2,otho=TRUE,re=FALSE,recalibrate=FALSE) {
    if(re==TRUE) recalibrate <- TRUE;
    cor.z2 <- as.matrix(cor.z2);
    beta.mat[which(beta.mat==Inf)] <- NA;
    beta.mat[which(beta.mat==-Inf)] <- NA;
    
    se.mat[which(se.mat==Inf)] <- NA;
    se.mat[which(se.mat==-Inf)] <- NA;
    w.mat <- 1/se.mat;
    z.mat <- beta.mat/se.mat;
    statistic <- rep(NA,nrow(z.mat));p.value <- rep(NA,nrow(z.mat));
    statistic.omnibus.mrvt <- statistic;
    p.value.omnibus.mrvt <- statistic;
    no.pc <- statistic;
    no.pc.aic <- 0;
    no.pc.bic <- 0;
    p.value.original <- statistic;
    chisq.mat <- matrix(nrow=length(statistic),ncol=nrow(cor.z2));
    gamma.est.mat <- chisq.mat;
    gamma.sd.mat <- chisq.mat;
    pval.chisq.mat <- chisq.mat;
    gamma.est.optim.mat <- chisq.mat;
    gamma.est.aic.mat <- chisq.mat;
    gamma.est.bic.mat <- chisq.mat;
    aic.mat <- chisq.mat;
    bic.mat <- aic.mat;
    gamma.est.list <- list();
    pval.re <- rep(NA,nrow(z.mat));
    statistic.re <- rep(NA,nrow(z.mat));
    for(jj in 1:nrow(z.mat)) {
        
        if(jj/1000==as.integer(jj/1000)) cat(jj, ' variants finished\n');
        if(!all(is.na(beta.mat[jj,]))) {
            lm.weight.all <- matrix(0,nrow=length(w.mat[jj,]),ncol=length(w.mat[jj,]));
            diag(lm.weight.all) <- (w.mat[jj,])^2;
            ix.rm <- which(is.na(beta.mat[jj,]) | is.na(w.mat[jj,]));
            beta.jj <- beta.mat[jj,];
            se.jj <- se.mat[jj,];
            af.pca.jj <- af.pca;
            
            if(length(ix.rm)>0) {
                lm.weight.all <- matrix(lm.weight.all[-ix.rm,-ix.rm],nrow=nrow(lm.weight.all)-length(ix.rm));
                beta.jj <- beta.mat[jj,-ix.rm];
                se.jj <- se.mat[jj,-ix.rm];
                af.pca.jj <- matrix(af.pca[-ix.rm,],ncol=ncol(af.pca));
            }
            if(otho==TRUE) af.pca.jj <- gram.schmitt(af.pca.jj,ip.w,w=lm.weight.all);
            z2.cum <- 0;
            aic <- 0;
            bic <- 0;
            for(ll in 1:ncol(af.pca.jj)) {
                af.pca.tmp <- matrix(af.pca.jj[,1:ll],ncol=ll);
                gamma.est <- ginv(t(af.pca.tmp)%*%lm.weight.all%*%af.pca.tmp)%*%(t(af.pca.tmp)%*%lm.weight.all%*%beta.jj);
                gamma.est.list[[ll]] <- rep(0,ncol(af.pca.jj));
                gamma.est.list[[ll]][1:ll] <- gamma.est;
                gamma.var <- ginv(t(af.pca.tmp)%*%lm.weight.all%*%af.pca.tmp);
                z2.cum[ll] <- t(gamma.est)%*%ginv(gamma.var)%*%gamma.est;
                res.beta <- beta.jj-(af.pca.tmp)%*%as.matrix(gamma.est)
                
                aic[ll] <-  t(res.beta)%*%lm.weight.all%*%res.beta+2*ll;
                bic[ll] <- t(res.beta)%*%lm.weight.all%*%res.beta+ll*log(nrow(af.pca.jj));
            }
            if(recalibrate==TRUE & re==FALSE) {
                mu.z <- rm.na(sqrt(z2.cum-c(0,z2.cum[-1])));
            
                z.simu <- matrix(rnorm(1e3*length(mu.z),mean=mu.z),ncol=length(mu.z),byrow=TRUE);
                z2.simu <- z.simu^2;
                z2.simu.cumsum <- matrix(apply(z2.simu,1,cumsum),ncol=4,byrow=TRUE);
                cor.z2 <- cor(z2.simu.cumsum);
            }
            
            if(re==TRUE) {
                res.re <- re2.mr(beta.jj,se.jj,af.pca.jj);
                pval.re.jj <- res.re$p.value;
                statistic.re.jj <- res.re$statistic;
                z2.cum.re <- c(z2.cum,statistic.re.jj);
                
                mu.z <- rm.na(sqrt(z2.cum-c(0,z2.cum[-1])));
                
                z.simu <- matrix(rnorm(1e3*length(mu.z),mean=mu.z),ncol=length(mu.z),byrow=TRUE);
                z2.simu <- z.simu^2;
                z2.simu.cumsum <- matrix(apply(z2.simu,1,cumsum),ncol=4,byrow=TRUE);
                ind <- rbinom(1e3,1,.5);
                z2.simu.cumsum.re <- cbind(z2.simu.cumsum, ind*z2.simu.cumsum[,ncol(z2.simu.cumsum)]+ (1-ind)*(z2.simu.cumsum[,ncol(z2.simu.cumsum)]+rchisq(1e3,df=1)))
                cor.z2.re <- cor(z2.simu.cumsum.re);
            
            }
            aic.mat[jj,] <- aic;
            bic.mat[jj,] <- bic;
            gamma.est.mat[jj,] <- gamma.est;
            gamma.sd.mat[jj,] <- sqrt(diag(gamma.var));
            p.z2 <- pchisq(z2.cum,df=1:length(z2.cum),lower.tail=FALSE);
            
            pval.chisq.mat[jj,] <- p.z2;
            invZ.minp <- qnorm(min(p.z2),lower.tail=FALSE);
            p.value.omnibus.mrvt[jj] <- pvt(invZ.minp,mu=rep(0,nrow(cor.z2)),sigma=cor.z2,alternative="greater");
            statistic.omnibus.mrvt[jj] <- z2.cum[which.min(p.z2)];
            gamma.est.optim.mat[jj,] <- gamma.est.list[[which.min(p.z2)]];
            gamma.est.aic.mat[jj,] <- gamma.est.list[[which.min(aic)]];
            gamma.est.bic.mat[jj,] <- gamma.est.list[[which.min(bic)]];
            
            no.pc[jj] <- which.min(p.z2)-1;
            no.pc.aic[jj] <- which.min(aic)-1;
            no.pc.bic[jj] <- which.min(bic)-1;
            chisq.mat[jj,] <- z2.cum;
            
            if(re==TRUE) {
                invZ.minp <- qnorm(min(c(p.z2,pval.re.jj)),lower.tail=FALSE)
                pval.re[jj] <- pval.re.jj;
                statistic.re[jj] <- statistic.re.jj;
                p.value.omnibus.mrvt[jj] <- pvt(invZ.minp,mu=rep(0,nrow(cor.z2.re)),sigma=cor.z2.re,alternative="greater");
            }
            
            p.value.original[jj] <- pchisq(t(gamma.est)%*%ginv(gamma.var)%*%gamma.est,lower.tail=FALSE,df=length(gamma.est));
            
        }
    }
    return(list(statistic.omnibus.mrvt=statistic.omnibus.mrvt,
                p.value.omnibus.mrvt=p.value.omnibus.mrvt,
                p.value.original=p.value.original,
                chisq.mat=chisq.mat,
                p.value.memo=p.value.omnibus.mrvt,
                p.value.mrmega=p.value.original,
                pval.re=pval.re,
                statistic.re=statistic.re,
                pval.chisq.mat=pval.chisq.mat,
                gamma.est.mat=gamma.est.mat,
                gamma.sd.mat=gamma.sd.mat,
                gamma.est.optim.mat=gamma.est.optim.mat,
                gamma.est.aic.mat=gamma.est.aic.mat,
                gamma.est.bic.mat=gamma.est.bic.mat,
                no.pc.omnibus.mrvt=no.pc,
                no.pc.aic=no.pc.aic,
                no.pc.bic=no.pc.bic,
                no.pc=no.pc));
}

#' trans-ethnic meta-analysis
#'
#' @param dat the formatted data;
#' @param af.pca the PCs or MDS coordinate for the allele frequencies
#' @return trans-ethnic results; 
#' @export 
trans.ethnic.meta <- function(dat,af.pca,cor.z2,otho=TRUE,re=FALSE,recalibrate=FALSE) {
    z.mat <- matrix(dat$ustat.mat/dat$vstat.mat,nrow=nrow(dat$ustat.mat),ncol=ncol(dat$ustat.mat));
    w.mat <- matrix(rm.na(sqrt(rm.na(dat$nSample)*rm.na(dat$af.mat)*(1-rm.na(dat$af.mat)))*dat$w.mat),nrow=nrow(dat$ustat.mat),ncol=ncol(dat$ustat.mat));
    beta.mat <- z.mat/w.mat;
    se.mat <- 1/w.mat;

    return(trans.ethnic.meta.core(beta.mat,se.mat,af.pca,cor.z2,otho,re,recalibrate=recalibrate));
}


#' estimate the empirical correlation between statistics;
#'
#' @param statistic.mat matrix of statistics 
#' @param p.value.mat matrix of p-values that corresponds to the statistic
#' @return p.value.minp and the estimated covariance matrix;
#' @export
emp.copula <- function(statistic.mat,p.value.mat) {
    minp <- apply(p.value.mat,1,min,na.rm=TRUE);
    ix.sig <- which(minp<.05);
    ix.insig <- which(minp>.05);
    cor.mat.sig <- cor(statistic.mat[ix.sig,]);
    cor.mat.insig <- cor(statistic.mat[ix.insig,]);

    p.value.minp <- 0;
    for(ii in 1:length(minp)) {
        if(minp[ii]<0.05) {
            p.value.minp[ii] <- pvt(qnorm(minp[ii],lower.tail=FALSE),mu=rep(0,nrow(cor.mat.sig)),sigma=cor.mat.sig,alternative="greater");
        }
        if(minp[ii]>=0.05) {
            p.value.minp[ii] <- pvt(qnorm(minp[ii],lower.tail=FALSE),mu=rep(0,nrow(cor.mat.insig)),sigma=cor.mat.insig,alternative="greater");
        }
    }
    return(list(p.value.minp=p.value.minp,
                cor.mat.sig=cor.mat.sig,
                cor.mat.insig=cor.mat.insig));
}

#' Perform multi-allelic association tests;
#'
#' @param pos The position of multi-allelic site;
#' @param dat The formated data;
#' @param corMultiAllele.mat the correlation matrix generated from multiple alleles;
#' @return a matrix with formatted result;
#' @export 
multiAlleleAssoc <- function(pos,dat,corMultiAllele.mat,...) {
    pars <- list(...);
    boot <- pars$boot;
    if(is.null(boot)) boot <- FALSE;
    nBoot <- 100;
    if(is.null(pars$sandwich)) {pars$sandwich <- TRUE;}
    res.formatted <- matrix(rep(NA,11),nrow=1);
    res.collapse <- matrix(rep(NA,2),nrow=1);
    ix <- which(corMultiAllele.mat$chrpos==pos);
    if(length(ix)==1) {
        ref.cor <- unlist(strsplit(corMultiAllele.mat$ref[ix],split=','));
        alt.cor <- unlist(strsplit(corMultiAllele.mat$alt[ix],split=','));
        if(length(ref.cor)<length(alt.cor))  ref.cor[(length(ref.cor)+1):length(alt.cor)] <- ref.cor[1];
        if(length(ref.cor)>length(alt.cor))  alt.cor[(length(alt.cor)+1):length(ref.cor)] <- alt.cor[1];
        
        cor.mat <- rm.na(matrix(as.numeric(unlist(strsplit(corMultiAllele.mat$cor[ix],split=','))),nrow=length(ref.cor),ncol=length(ref.cor)));
        ref.alt.cor <- paste(ref.cor,alt.cor,sep=',');
        pos.dat <- gsub("_.*","",dat$pos);
        ix <- which(pos.dat==pos);
        dat$pos <- pos.dat[ix];
        dat$ref.mat <- matrix(dat$ref.mat[ix,],nrow=length(ix))
        dat$alt.mat <- matrix(dat$alt.mat[ix,],nrow=length(ix));
        dat$nref.mat <- matrix(dat$nref.mat[ix,],nrow=length(ix))
        dat$nalt.mat <- matrix(dat$nalt.mat[ix,],nrow=length(ix));
        dat$nhet.mat <- matrix(dat$nhet.mat[ix,],nrow=length(ix));
        dat$ustat.mat <- matrix(dat$ustat.mat[ix,],nrow=length(ix));
        dat$vstat.mat <- matrix(dat$vstat.mat[ix,],nrow=length(ix));
        dat$af.mat <- matrix(dat$af.mat[ix,],nrow=length(ix));
        dat$w.mat <- matrix(dat$w.mat[ix,],nrow=length(ix));
        dat$nSample.mat <- matrix(dat$nSample.mat[ix,],nrow=length(ix));
        dat$ref.alt.mat <- matrix(paste(dat$ref.mat,dat$alt.mat,sep=','),nrow=length(ix));
        ref.alt.dat <- unique(as.character(paste(dat$ref.mat,dat$alt.mat,sep=',')));
        ref.alt.both <- intersect(ref.alt.dat,ref.alt.cor);
        
        if(length(ref.alt.both)>0) {            
            ref.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$ref.mat));
            alt.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$alt.mat));
            ref.alt.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$ref.alt.mat));
            
            nref.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$nref.mat));
            nalt.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$nalt.mat));
            nhet.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$nhet.mat));
            ustat.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$ustat.mat));
            vstat.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$vstat.mat));
            af.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$af.mat));
            w.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$w.mat));
            nSample.mat <- matrix(nrow=length(ref.alt.both),ncol=ncol(dat$nSample.mat));
            dat$pos <- pos.dat;
            
            for(ii in 1:ncol(ref.mat)) {
                ix.match <- match(ref.alt.both,dat$ref.alt.mat[,ii]);
                ref.mat[,ii] <- dat$ref.mat[ix.match,ii];
                alt.mat[,ii] <- dat$alt.mat[ix.match,ii];
                nref.mat[,ii] <- dat$nref.mat[ix.match,ii];
                nalt.mat[,ii] <- dat$nalt.mat[ix.match,ii];
                nhet.mat[,ii] <- dat$nhet.mat[ix.match,ii];
                ustat.mat[,ii] <- dat$ustat.mat[ix.match,ii];
                vstat.mat[,ii] <- dat$vstat.mat[ix.match,ii];
                af.mat[,ii] <- dat$af.mat[ix.match,ii];
                w.mat[,ii] <- dat$w.mat[ix.match,ii];
                nSample.mat[,ii] <- dat$nSample.mat[ix.match,ii];                
            }
            dat$ref.mat <- ref.mat;
            dat$alt.mat <- alt.mat;
            dat$nref.mat <- nref.mat;
            dat$nalt.mat <- nalt.mat;
            dat$nhet.mat <- nhet.mat;
            dat$ustat.mat <- ustat.mat;
            dat$vstat.mat <- vstat.mat;
            dat$af.mat <- af.mat;
            dat$w.mat <- w.mat;
            dat$nSample.mat <- nSample.mat;
            dat$direction.mat <- dat$ustat.mat;
            dat$direction.mat[which(dat$ustat.mat>0,arr.ind=T)] <- "+";
            dat$direction.mat[which(dat$ustat.mat<0,arr.ind=T)] <- "-";
            dat$direction.mat[which(dat$ustat.mat==0,arr.ind=T)] <- "=";
            dat$direction.mat[which(is.na(dat$ustat.mat),arr.ind=T)] <- "X";
            direction.meta <- apply(dat$direction.mat,1,paste,sep="",collapse="");
            
            af.meta <- rowSums((dat$af.mat)*(dat$nSample.mat),na.rm=T)/rowSums((dat$nSample.mat),na.rm=T);
            
            ix <- match(ref.alt.both,ref.alt.cor);
            if(length(ix)>0) {
                cor.mat <- matrix(as.numeric(cor.mat[ix,ix]),nrow=length(ix),ncol=length(ix));
                ustat.meta.tmp <- rowSums((dat$ustat.mat)*(dat$w.mat),na.rm=TRUE);
                vstat.meta.tmp <- sqrt(rowSums((dat$vstat.mat)^2*(dat$w.mat)^2,na.rm=TRUE));
                V <- matrix(0,nrow=nrow(cor.mat),ncol=ncol(cor.mat));
                V0 <- V;
                for(ii in 1:ncol(dat$ustat.mat)) {
                    S <- matrix(0,nrow=nrow(cor.mat),ncol=ncol(cor.mat));
                    diag(S) <- rm.na(sqrt((dat$w.mat[,ii])^2*(dat$vstat.mat[,ii])^2))
                    V0 <- V0+(as.vector(rm.na(dat$w.mat[,ii]*dat$ustat.mat[,ii]))%*%t(as.vector(rm.na(dat$w.mat[,ii]*dat$ustat.mat[,ii]))));
                    V <- V+S%*%cor.mat%*%S;
                }
                V <- regMat(V,0.1);
                V0 <- regMat(V0,.1);
                cor.V <- rm.na(cov2cor(V));
                diag(cor.V) <- 1;
                
                cor.V0 <- rm.na(cov2cor(V0));
                diag(cor.V0) <- 1;
                sd.mat <- matrix(0,nrow=nrow(V),ncol=ncol(V));
                diag(sd.mat) <- vstat.meta.tmp;
                V0 <- sd.mat%*%cor.V0%*%sd.mat;
                numStability <- rep(1,length(ustat.meta.tmp));
                beta.est.meta <- 0;
                beta.sd.meta <- 0;
                beta.sd.meta.std <- 0;
                nVar <- nrow(dat$ustat.mat);
                if(boot==TRUE) {
                    nStudy <- ncol(dat$ustat.mat);
                    boot.var <- bootVar.multiallelic(ustat.meta.tmp,vstat.meta.tmp,cor.V0,nStudy);

                }
                if(boot==FALSE) boot.var <- rep(0,nVar);
                
                for(jj in 1:length(ustat.meta.tmp)) {
                    ix.keep <- 1:nrow(V);
                    V.tmp <- as.matrix(V[ix.keep,ix.keep]);
                    ustat.meta.tmp1 <- ustat.meta.tmp[ix.keep];
          
                    beta.est.tmp <- NA;beta.sd.tmp <- NA;beta.sd.std <- NA;
                    beta.est.tmp[ix.keep] <- ginv(V0)%*%ustat.meta.tmp1;
                    beta.sd.tmp[ix.keep] <- sqrt(as.vector(diag(ginv(V0)))+boot.var);
                    beta.sd.std[ix.keep] <- sqrt(as.vector(diag(ginv(V.tmp))));
                    beta.est.meta[jj] <- beta.est.tmp[jj];
                    beta.sd.meta[jj] <- beta.sd.tmp[jj];
                    beta.sd.meta.std[jj] <- beta.sd.std[jj];
                    numStability[jj] <- max(abs(cor.V[jj,-jj]),na.rm=TRUE);
                    if(length(ix.keep)<length(ustat.meta.tmp)) numStability[jj] <- 0;
                }
                statistic.meta <- (beta.est.meta/beta.sd.meta)^2;
                statistic.meta.std <- (beta.est.meta/beta.sd.meta.std)^2;
                p.value.meta <- pchisq(statistic.meta,df=1,lower.tail=FALSE);
                p.value.meta.std <- pchisq(statistic.meta.std,df=1,lower.tail=FALSE);

                nSample.meta <- rowSums(dat$nSample.mat,na.rm=T);
                sum.weight <- rowSums((dat$nSample)*(dat$w.mat)^2,na.rm=TRUE);

                ref.tab <- apply(dat$ref.mat,1,uniq.allele);
                alt.tab <- apply(dat$alt.mat,1,uniq.allele);
                res.formatted <- cbind(pos,
                                       ref.tab,
                                       alt.tab,
                                       format(af.meta,digits=3),
                                       format(statistic.meta,digits=3),
                                       format(p.value.meta,digits=3),
                                       format(beta.est.meta,digits=3),
                                       format(beta.sd.meta,digits=3),
                                       nSample.meta,
                                       direction.meta,
                                       sum.weight,
                                       numStability);
                ix.0 <- which(p.value.meta==0);
                if(length(ix.0)>0) {
                    pval.precise <- get.precise.pval(statistic.meta[ix.0]);
                    res.formatted[ix.0,6] <- pval.precise;
                }
                
                res.formatted.std <- cbind(pos,
                                           ref.tab,
                                           alt.tab,
                                           format(af.meta,digits=3),
                                           format(statistic.meta.std,digits=3),
                                           format(p.value.meta.std,digits=3),
                                           format(beta.est.meta,digits=3),
                                           format(beta.sd.meta.std,digits=3),
                                           nSample.meta,
                                           direction.meta,
                                           nSample.meta,
                                           numStability);
                if(pars$sandwich==FALSE) {
                    res.formatted <- res.formatted.std;
                }
                stat.collapse <- (sum(ustat.meta.tmp))^2/sum(V);
                p.value.collapse <- pchisq(stat.collapse,df=1,lower.tail=FALSE);
                res.collapse <- matrix(c(unique(pos),p.value.collapse),nrow=1);
            }
        }
    }
    return(list(res.formatted=res.formatted,
                res.formatted.std=res.formatted.std,
                res.collapse=res.collapse)); 
}
#' get variance from the bootstrap for multi-allelic site;
#'
#' @param ustat.meta score statistic from meta-analysis
#' @param vstat.meta standard deviation of the meta-analysis score statistics
#' @param r2 the LD matrix;
#' @param nStudy the number of studies; 
#' @export 
bootVar.multiallelic <- function(ustat.meta,vstat.meta,r2,nStudy,nBoot=100) {
    nVar <- length(ustat.meta);
    vstat.sq.meta <- vstat.meta^2;
    sd.mat <- matrix(0,nrow=nVar,ncol=nVar);
    diag(sd.mat) <- sqrt(vstat.sq.meta);
    z.simu <- rmvnorm(nStudy*nBoot,mean=rep(0,nVar),sigma=r2);
    beta.boot <- 0;
    beta.var.boot <- 0;
    beta.boot <- matrix(0,nrow=nBoot,ncol=nVar);
    for(ii in 1:nBoot) {
        z.ii <- t(matrix(z.simu[((ii-1)*nStudy+1):(ii*nStudy),],nrow=nStudy));
        cor.ii <- z.ii%*%t(z.ii)/nStudy;
        
        V.ii <- sd.mat%*%cor.ii%*%sd.mat;
        beta.boot[ii,] <- ginv(V.ii)%*%ustat.meta;
    }
    boot.var <- apply(beta.boot,2,var,na.rm=TRUE);
    return(boot.var);
}

#' format data into matrices;
#' @param raw.data.all The read in from summary assoicaiton statistic files;
#' @param raw.imp.qual The raw imputation quality files;
#' @param impQualWeight If using imputation quality as weight
#' @return a list of converted data for downstream meta-analysis;
#' @export
GWAMA.formatData <- function(raw.data.all,raw.imp.qual,impQualWeight,impQual.lb,col.impqual,...) {
    pos <- raw.data.all$pos;
    extraPar <- list(...);
    rmMultiAllelicSite <- extraPar$rmMultiAllelicSite;
    
    maf.cutoff <- extraPar$maf.cutoff;
    if(is.null(maf.cutoff)) maf.cutoff <- 1.01;
    knownVar <- extraPar$knownVar;
    ref.mat <- matrix(unlist(raw.data.all$ref),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    alt.mat <- matrix(unlist(raw.data.all$alt),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    
    ref.tab <- apply(ref.mat,1,uniq.allele);
    alt.tab <- apply(alt.mat,1,uniq.allele);
    nSample.mat <- matrix(unlist(raw.data.all$nSample),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    nref.mat <- matrix(unlist(raw.data.all$nref),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    nalt.mat <- matrix(unlist(raw.data.all$nalt),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    nhet.mat <- matrix(unlist(raw.data.all$nhet),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    af.mat <- matrix(unlist(raw.data.all$af),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    ustat.mat <- matrix(unlist(raw.data.all$ustat),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    vstat.mat <- matrix(unlist(raw.data.all$vstat),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    ref.mat <- matrix(unlist(raw.data.all$ref),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    alt.mat <- matrix(unlist(raw.data.all$alt),ncol=length(raw.data.all$ref),nrow=length(raw.data.all$ref[[1]]));
    ustat.mat[which(ustat.mat==Inf | ustat.mat==-Inf,arr.ind=TRUE)] <- NA;
    vstat.mat[which(vstat.mat==Inf | vstat.mat==-Inf,arr.ind=TRUE)] <- NA;
    
    if(impQualWeight==FALSE) 
        w.mat <- matrix(1,nrow=nrow(vstat.mat),ncol=ncol(vstat.mat));
    if(!is.null(raw.imp.qual)) {
        imp.qual <- matrix(unlist(lapply(raw.imp.qual,getImpQual,rmMultiAllelicSite=rmMultiAllelicSite,pos=pos,col.impqual=col.impqual)),ncol=ncol(ref.mat),nrow=nrow(ref.mat));
        ix.lowQual <- which(imp.qual<impQual.lb,arr.ind=TRUE);
        ustat.mat[ix.lowQual] <- NA;
        vstat.mat[ix.lowQual] <- NA;
        w.mat <- sqrt(imp.qual);
        w.mat[which(is.na(w.mat),arr.ind=TRUE)] <- 1;
        nSample.mat[ix.lowQual] <- NA;
    }
    direction.mat <- ustat.mat;
    direction.mat[which(ustat.mat>0,arr.ind=T)] <- "+";
    direction.mat[which(ustat.mat<0,arr.ind=T)] <- "-";
    direction.mat[which(ustat.mat==0,arr.ind=T)] <- "=";
    direction.mat[which(is.na(ustat.mat),arr.ind=T)] <- "X";

    af.vec <- rowSums(af.mat*nSample.mat,na.rm=TRUE)/rowSums(nSample.mat,na.rm=TRUE);
    maf.vec <- rm.na(af.vec);
    maf.vec[which(maf.vec>.5)] <- 1-maf.vec[which(maf.vec>.5)];
    pos.norefalt <- gsub("_.*","",pos);
    ix.rare <- unique(c(which(maf.vec<=maf.cutoff & maf.vec>=0),which(pos.norefalt%in%knownVar)));
    nref.mat <- matrix(nref.mat[ix.rare,],nrow=length(ix.rare));
    nalt.mat <- matrix(nalt.mat[ix.rare,],nrow=length(ix.rare));
    ref.tab <- ref.tab[ix.rare];
    alt.tab <- alt.tab[ix.rare];
    af.mat <- matrix(af.mat[ix.rare,],nrow=length(ix.rare));
    ref.mat <- matrix(ref.mat[ix.rare,],nrow=length(ix.rare));
    alt.mat <- matrix(alt.mat[ix.rare,],nrow=length(ix.rare));
    nhet.mat <- matrix(nhet.mat[ix.rare,],nrow=length(ix.rare));
    ustat.mat <- matrix(ustat.mat[ix.rare,],nrow=length(ix.rare));
    w.mat <- matrix(w.mat[ix.rare,],nrow=length(ix.rare));
    pos <- pos[ix.rare]
    nSample.mat <- matrix(nSample.mat[ix.rare,],nrow=length(ix.rare));
    direction.mat <- matrix(direction.mat[ix.rare,],nrow=length(ix.rare));
    vstat.mat <- matrix(vstat.mat[ix.rare,],nrow=length(ix.rare));
    
    return(list(nref.mat=nref.mat,
                nalt.mat=nalt.mat,
                ref.tab=ref.tab,
                alt.tab=alt.tab,
                af.mat=af.mat,
                ref.mat=ref.mat,
                alt.mat=alt.mat,
                nhet.mat=nhet.mat,
                ustat.mat=ustat.mat,
                w.mat=w.mat,
                pos=pos,
                maf.meta=maf.vec[ix.rare],
                nSample.mat=nSample.mat,
                direction.mat=direction.mat,
                vstat.mat=vstat.mat));
}
#' remove multi-allelic sites;
#'
#' @param dat The dataset;
#' @return the data set without multi-allelic sites;
#' @export
GWAMA.rmMulti <- function(dat) {
    pos.tmp <- gsub("_.*","",dat$pos);
    pos.tab <- table(pos.tmp);
    pos.multi <- names(pos.tab)[which(pos.tab>1)];
    if(length(pos.multi)>0) ix.multi <- (1:length(pos.tmp))[pos.tmp %in% pos.multi];
    if(length(pos.multi)==0) ix.multi <- integer(0);
    if(length(ix.multi)>0) {
        dat$pos <- dat$pos[-ix.multi];
        dat$ref.mat <- matrix(dat$ref.mat[-ix.multi,],nrow=nrow(dat$ref.mat)-length(ix.multi),ncol=ncol(dat$ref.mat));
        dat$alt.mat <- matrix(dat$alt.mat[-ix.multi,],nrow=nrow(dat$alt.mat)-length(ix.multi),ncol=ncol(dat$alt.mat));
        dat$ref.tab <- dat$ref.tab[-ix.multi];
        dat$alt.tab <- dat$alt.tab[-ix.multi];
        dat$w.mat <- matrix(dat$w.mat[-ix.multi,],nrow=nrow(dat$w.mat)-length(ix.multi),ncol=ncol(dat$w.mat));
        dat$nSample.mat <- matrix(dat$nSample.mat[-ix.multi,],nrow=nrow(dat$nSample.mat)-length(ix.multi),ncol=ncol(dat$nSample.mat));
        dat$nref.mat <- matrix(dat$nref.mat[-ix.multi,],nrow=nrow(dat$nref.mat)-length(ix.multi),ncol=ncol(dat$nref.mat));
        dat$nalt.mat <- matrix(dat$nalt.mat[-ix.multi,],nrow=nrow(dat$nalt.mat)-length(ix.multi),ncol=ncol(dat$nalt.mat));
        dat$nhet.mat <- matrix(dat$nhet.mat[-ix.multi,],nrow=nrow(dat$nhet.mat)-length(ix.multi),ncol=ncol(dat$nhet.mat));
        dat$af.mat <- matrix(dat$af.mat[-ix.multi,],nrow=nrow(dat$af.mat)-length(ix.multi),ncol=ncol(dat$af.mat));
        dat$ustat.mat <- matrix(dat$ustat.mat[-ix.multi,],nrow=nrow(dat$ustat.mat)-length(ix.multi),ncol=ncol(dat$ustat.mat));
        dat$vstat.mat <- matrix(dat$vstat.mat[-ix.multi,],nrow=nrow(dat$vstat.mat)-length(ix.multi),ncol=ncol(dat$vstat.mat));
        dat$direction.mat <- matrix(dat$direction.mat[-ix.multi,],nrow=nrow(dat$direction.mat)-length(ix.multi),ncol=ncol(dat$direction.mat));
        dat$posMulti <- pos.multi;
    }
    return(list(dat=dat,
                posMulti=pos.multi));
}

#' gateway function for multi-allelic analysis in rareGWAMA
#'
#' @param score stat.file The score statistics file names;
#' @param imp.qual.file The imputation quality file names;
#' @param cor.multiallele.file The correlation information for multi-allelic sites;
#' @param tabix.range The tabix range for the variants to be analyzed;
#' @param alternative The alternative hypothesis; only the default choice two.sided is implemented now;
#' @param col.impqual The column number for the imputation quality files; The default choice is 5.
#' @param impQual.lb The lower bound for the imputation quality cutoff; Variants with imputation quality less than the cutoffs are labelled as missing and removed from the meta-analysis
#' @param impQualWeight Whether to apply the imputation quality based optimal weighting;
#' @export
#' @return A list of formatted output; 
rareGWAMA.single.multiAllele <- function(score.stat.file,imp.qual.file=NULL,tabix.range,alternative="two.sided",col.impqual=5,impQual.lb=0.7,impQualWeight=FALSE) {
    a <- Sys.time();
    capture.output(raw.data.all <- rvmeta.readDataByRange( score.stat.file, NULL, tabix.range,multiAllelic = TRUE));
    raw.imp.qual <- NULL;
    if(!is.null(imp.qual.file))
        raw.imp.qual <- lapply(imp.qual.file,tabix.read.table,tabixRange=tabix.range);
    
    time.readData <- Sys.time()-a;
    b <- Sys.time();
    raw.data.all <- raw.data.all[[1]];
    cat('Read in',length(raw.data.all$ref[[1]]),'variants',sep=' ');
    dat <- GWAMA.formatData(raw.data.all,raw.imp.qual,impQualWeight,impQual.lb,col.impqual,rmMultiAllelicSite=rmMultiAllelicSite);
    pos <- unique(gsub("_.*","",dat$pos));
    res.assoc.tmp <- sapply(pos,multiAlleleAssoc,dat=dat,corMultiAllele.mat=corMultiAllele.mat);
    res.formatted <- do.call(rbind,res.assoc.tmp);
    colnames(res.formatted) <- c("POS","REF","ALT","AF","STAT","PVALUE","BETA","SD","N","DIRECTION");
    return(list(res.formatted=res.formatted,
                formattedData=dat,
                raw.data.all=raw.data.all,
                raw.imp.qual=raw.imp.qual));
}
#' convert chisq statistic to beta for binary trait assuming the variance for y;
#' @param statistic chisq stat;
#' @param beta.ori original beta; for getting the sign;
#' @param var.y the variance for Y. for binary trait, the variance of y is frac.case*(1-frac.case); for continuous trait, the variance for y is typically set to 1; 
#' @param af allele freq;
#' @param N sample size;
#' @export
convertChisq2Beta <- function(statistic,beta.ori,var.y,af,N,binaryTrait=FALSE) {

    if(binaryTrait==FALSE) {
        beta.out <- sign(beta.ori)*sqrt(statistic)*sqrt(var.y)/sqrt(2*N*af*(1-af));
        beta.sd <- sqrt(var.y/(2*N*af*(1-af)));
    }
    if(binaryTrait==TRUE) {
        beta.out <- sign(beta.ori)*sqrt(statistic)/sqrt(2*N*af*(1-af))/sqrt(var.y);
        beta.sd <- (1/sqrt(2*N*af*(1-af)*var.y));
    }
    return(list(beta.est=beta.out,
                beta.sd=beta.sd));
}

#' meta-regression likelihood with variance heterogeneity;
#' @param pars the parameters for the gamma and the log variance 
#' @param b.vec genetic effect estimates
#' @param se.vec the standard deviation of the genetic effect estimates
#' @param af.pca the covariates to be included
#' @return the likelihood values; 
#' @export
l.mr <- function(pars,b.vec,se.vec,af.pca) {
    af.pca <- as.matrix(af.pca);
    gamma <- pars[1:ncol(af.pca)];
    log.tau <- pars[ncol(af.pca)+1];
    tau2 <- exp(2*log.tau);
    l <- (-1)*sum(dnorm(b.vec,af.pca%*%gamma,sqrt(se.vec^2+tau2),log=T));
    return(l);

}
#' restricted likelihood with no variance heterogeneity;
#' @param pars.0 the parameters for the gamma
#' @param b.vec genetic effect estimates
#' @param se.vec the standard deviation of the genetic effect estimates
#' @param af.pca the covariates to be included
#' @return the likelihood values; 
#' @export
l.mr.0 <- function(pars.0,b.vec,se.vec,af.pca) {
    pars <- c(pars.0,-Inf);
    return(l.mr(pars,b.vec,se.vec,af.pca));
}

#' mixed effect meta-regression
#'
#' @param b.vec genetic effect estimates
#' @param se.vec standard deviation for the genetic effect estimates
#' @param af.pca the covariates to be included in the meta-regression
#' @return a list consists of p.value, statistic, gamma.est and tau2;
#' @export
re2.mr <- function(b.vec,se.vec,af.pca) {
    af.pca <- as.matrix(af.pca);
    pars.init <- c(rep(0,ncol(af.pca)),log(.02));
    res.est <- nlminb(pars.init,l.mr,b.vec=b.vec,se.vec=se.vec,af.pca=af.pca);
    
    pars.0.init <- rep(0,ncol(af.pca));
    res.est.0 <- nlminb(pars.0.init,l.mr.0,b.vec=b.vec,se.vec=se.vec,af.pca=af.pca);
    pars.est <- c(res.est.0$par,-Inf);
    if(res.est$objective < res.est.0$objective) pars.est <- res.est$par;
    gamma.est <- pars.est[1:ncol(af.pca)];
    tau2 <- exp(2*pars.est[ncol(af.pca)+1]);
    statistic <- sum(log(se.vec^2/(se.vec^2+tau2))+b.vec^2/se.vec^2-(b.vec-af.pca%*%gamma.est)^2/(se.vec^2+tau2));


    p.value <- pmixchisq(statistic,lambda=c(.5,.5),h=c(length(gamma.est),length(gamma.est)+1));
    return(list(p.value=p.value,
           statistic=statistic,
           gamma.est=gamma.est,
           tau2.est=tau2))
    
}

#' calculate the mixture chisq p-value;
#'
#' @param statistic the test statistic;
#' @param lambda mixture proportion
#' @param h multiplicity, and the default is rep(1,length(lambda))
#' @return p values for mixture chisq distributed random variables;
#' @export
pmixchisq <- function(statistic,lambda,h=NULL) {
    if(is.null(h)) h <- rep(1,length(lambda));
    p.value <- try(davies(statistic,lambda=lambda,h=h)$Qq,silent=TRUE);
    p.value.liu <- try(liu(statistic,lambda=lambda,h=h),silent=TRUE);

    if(class(p.value)=='try-error' | class(p.value.liu)=='try-error') return(NA);
    if(class(p.value)!='try-error') {
        if(p.value<=0 | p.value>=1) p.value <- p.value.liu
    }
    if(is.null(p.value)) p.value <- NA;
    return(p.value);
}
