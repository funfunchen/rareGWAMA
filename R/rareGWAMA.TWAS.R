#' trans-ethnic TWAS method;
#'
#' @param dat the data that contains everything needed;
#' @return a bunch of results, such as p-values;
#' @export
rareGWAMA.twas.predEff <- function(dat,af.pca,af.pca.eqtl) {

    z.mat <- matrix(dat$ustat.mat/dat$vstat.mat,nrow=nrow(dat$ustat.mat),ncol=ncol(dat$ustat.mat));

    w.mat <- matrix(rm.na(sqrt(rm.na(dat$nSample.mat)*rm.na(dat$af.mat)*(1-rm.na(dat$af.mat)))*dat$w.mat),nrow=nrow(dat$ustat.mat),ncol=ncol(dat$ustat.mat));
    z.mat <- matrix(z.mat[dat$ix.rare,],nrow=length(dat$ix.rare));
    w.mat <- matrix(w.mat[dat$ix.rare,],nrow=length(dat$ix.rare));
    beta.mat <- z.mat/w.mat;
    se.mat <- 1/w.mat;
    cov.beta.pred <- matrix(0,nrow=nrow(z.mat),ncol=nrow(z.mat));
    beta.pred <- rep(0,nrow(z.mat));
    af.pca <- as.matrix(af.pca);
    beta.pred.list <- list();
    cov.beta.pred.list <- list();
    gcov.beta.pred <- matrix(0,nrow=nrow(z.mat)*ncol(af.pca),ncol=nrow(z.mat)*ncol(af.pca));
    for(kk in 1:ncol(af.pca)) {
        for(jj1 in 1:nrow(z.mat)) {
            for(jj2 in 1:nrow(z.mat)) {
                if(jj1>=jj2) {
                    if(!all(is.na(beta.mat[jj1,])) & !all(is.na(beta.mat[jj2,]))) {
                        lm.weight.all1 <- matrix(0,nrow=length(w.mat[jj1,]),ncol=length(w.mat[jj1,]));
                        diag(lm.weight.all1) <- (w.mat[jj1,])^2;
                        lm.weight.all2 <- matrix(0,nrow=length(w.mat[jj2,]),ncol=length(w.mat[jj2,]));
                        diag(lm.weight.all2) <- (w.mat[jj2,])^2;                    
                        ix.rm <- which(is.na(beta.mat[jj1,]) | is.na(w.mat[jj1,]) | is.na(beta.mat[jj2,]) | is.na(w.mat[jj2,]));
                        beta.jj1 <- beta.mat[jj1,];
                        se.jj1 <- se.mat[jj1,];
                        beta.jj2 <- beta.mat[jj2,];
                        se.jj2 <- se.mat[jj2,];
                        af.pca.jj <- af.pca[,1:kk];
                        
                        if(length(ix.rm)>0) {
                            lm.weight.all1 <- matrix(lm.weight.all1[-ix.rm,-ix.rm],nrow=nrow(lm.weight.all1)-length(ix.rm));
                            lm.weight.all2 <- matrix(lm.weight.all2[-ix.rm,-ix.rm],nrow=nrow(lm.weight.all2)-length(ix.rm));
                            
                            beta.jj1 <- beta.mat[jj1,-ix.rm];
                            se.jj1 <- se.mat[jj1,-ix.rm];
                            beta.jj2 <- beta.mat[jj2,-ix.rm];
                            se.jj2 <- se.mat[jj2,-ix.rm];
                            
                            af.pca.jj <- matrix(af.pca[-ix.rm,1:kk],ncol=kk);
                        }
                        if(length(ix.rm)<ncol(beta.mat)) {
                            gamma.est1 <- ginv(t(af.pca.jj)%*%lm.weight.all1%*%af.pca.jj)%*%(t(af.pca.jj)%*%lm.weight.all1%*%beta.jj1);
                            beta.pred[jj1] <- as.numeric(t(af.pca.eqtl[1:kk])%*%gamma.est1);
                            A.L <- ginv(t(af.pca.jj)%*%lm.weight.all1%*%af.pca.jj)%*%(t(af.pca.jj)%*%lm.weight.all1);
                            A.R <- ginv(t(af.pca.jj)%*%lm.weight.all2%*%af.pca.jj)%*%(t(af.pca.jj)%*%lm.weight.all2);
                            cov.beta <- matrix(0,nrow=ncol(beta.mat)-length(ix.rm),ncol=ncol(beta.mat)-length(ix.rm));
                            for(bb in 1:nrow(cov.beta)) {
                                    cov.beta[bb,bb] <- se.jj2[bb]*se.jj1[bb]*dat$r2.by.study[[bb]][jj1,jj2];
                                }
                            cov.gamma.1.2 <- A.L%*%cov.beta%*%t(A.R);
                            
                            cov.beta.pred[jj1,jj2] <- t(af.pca.eqtl[1:kk])%*%cov.gamma.1.2%*%(af.pca.eqtl[1:kk]);
                            cov.beta.pred[jj2,jj1] <- cov.beta.pred[jj1,jj2];
                            
                            if(kk==1) {
                                for(ll1 in 1:ncol(af.pca)) {
                                    for(ll2 in 1:ncol(af.pca)) {
                                        af.pca.ll1 <- as.matrix(af.pca[,1:ll1]);
                                        af.pca.ll2 <- as.matrix(af.pca[,1:ll2]);
                                        if(length(ix.rm)>0) {
                                            af.pca.ll1 <- matrix(af.pca.ll1[-ix.rm,],ncol=ll1);
                                            af.pca.ll2 <- matrix(af.pca.ll2[-ix.rm,],ncol=ll2);
                                        }
                                            
                                        A.L1 <- ginv(t(af.pca.ll1)%*%lm.weight.all1%*%af.pca.ll1)%*%(t(af.pca.ll1)%*%lm.weight.all1);
                                        A.R2 <- ginv(t(af.pca.ll2)%*%lm.weight.all2%*%af.pca.ll2)%*%(t(af.pca.ll2)%*%lm.weight.all2);
                                        cov.gamma.1.2.tmp <- A.L1%*%cov.beta%*%t(A.R2);
                                        
                                        gcov.beta.pred[jj1+(ll1-1)*nrow(z.mat),jj2+(ll2-1)*nrow(z.mat)] <- rm.na(t(af.pca.eqtl[1:ll1])%*%cov.gamma.1.2.tmp%*%(af.pca.eqtl[1:ll2]));
                                        
                                        gcov.beta.pred[jj2+(ll1-1)*nrow(z.mat),jj1+(ll2-1)*nrow(z.mat)] <- gcov.beta.pred[jj1+(ll1-1)*nrow(z.mat),jj2+(ll2-1)*nrow(z.mat)];
                                        
                                    }
                                }
                            }
                        }
                    }
                }
                
            }
        }
        beta.pred.list[[kk]] <- beta.pred;
        cov.beta.pred.list[[kk]] <- cov.beta.pred;
        
    }
    
    
    return(list(beta.pred.list=beta.pred.list,
                cov.beta.pred.list=cov.beta.pred.list,
                gcov.beta.pred=gcov.beta.pred));
}

#' rareGWAMA.twas function
#'
#' @param dat the formatted data from indv studies;
#' @param af.pca the pca of allele frequencies from participating studies;
#' @param af.pca.eqtl
#' @export
rareGWAMA.twas <- function(dat,af.pca,pca.eqtl) {

    if(nrow(dat$ustat.mat)==0 | length(dat$ix.rare)==0)
        return(list(statistic=NA,
                    p.value=NA,
                    p.value.minp=NA,
                    maf.cutoff=1))
    res.in <- rareGWAMA.twas.predEff(dat,af.pca,pca.eqtl);
    
    beta.pred.list <- res.in$beta.pred.list;
    cov.beta.pred.list <- res.in$cov.beta.pred.list;

    statistic <- NA;p.value <- NA;
    z.stat <- NA;numer <- NA;
    for(kk in 1:length(beta.pred.list)) {
        beta.pred <- rm.na(beta.pred.list[[kk]]);
        
        cov.beta.pred <- rm.na(as.matrix(cov.beta.pred.list[[kk]]));
        reg.cov.beta.pred <- cov.beta.pred+beta.pred%*%t(beta.pred);
        cor.beta.pred <- rm.na(cov2cor(cov.beta.pred));
        reg.cor.beta.pred <- rm.na(cov2cor(reg.cov.beta.pred));
        diag(cor.beta.pred) <- 1;
        diag(reg.cor.beta.pred) <- 1;
        statistic.tmp <- sum(dat$twas.weight*(beta.pred/sqrt(diag(cov.beta.pred))),na.rm=TRUE);
        
        ix.na <- which(is.na(dat$twas.weight))
        if(length(ix.na)>0){
          statistic.var <- t(dat$twas.weight[-ix.na])%*%cor.beta.pred[-ix.na, -ix.na]%*%(dat$twas.weight[-ix.na]);
        } else {
          statistic.var <- t(dat$twas.weight)%*%cor.beta.pred%*%(dat$twas.weight);
        }
        
        statistic[kk] <- as.numeric(statistic.tmp^2/statistic.var);
        z.stat[kk] <- statistic.tmp/sqrt(statistic.var);
        numer[kk] <- statistic.tmp;
        p.value[kk] <- pchisq(statistic[kk],df=1,lower.tail=FALSE);
    }
    cov.stat <- matrix(0,nrow=length(beta.pred.list),ncol=length(beta.pred.list));
    for(kk1 in 1:length(beta.pred.list)) {
        for(kk2 in 1:length(beta.pred.list)) {
            cov.tmp <- res.in$gcov.beta.pred[((kk1-1)*length(dat$twas.weight)+1):((kk1)*length(dat$twas.weight)),((kk2-1)*length(dat$twas.weight)+1):((kk2)*length(dat$twas.weight))];
            
            ix.na <- which(is.na(dat$twas.weight))
            if (length(ix.na)>0){
              cov.stat[kk1,kk2] <- t(dat$twas.weight[-ix.na])%*%cov.tmp[-ix.na, -ix.na]%*%(dat$twas.weight[-ix.na]);  
            } else{
              cov.stat[kk1,kk2] <- t(dat$twas.weight)%*%cov.tmp%*%(dat$twas.weight);  
            }
        }
    }
    cov.stat <- make.pos.def(cov.stat);
    cor.stat <- rm.na(cov2cor(cov.stat));
    diag(cor.stat) <- 1;
    
    cor.stat <- 1/2*(cor.stat+t(cor.stat));
    
    if(all(is.na(statistic))){  ## debug added: in case all the value in statistics are NaN NaN NaN NaN
      return(list(statistic=NA,
                  p.value=NA,
                  p.value.minp=NA,
                  maf.cutoff=1))
    }

    pvalue.minp <- pvt(max(statistic),mu=rep(0,length(beta.pred.list)),sigma=cor.stat);
    return(list(statistic=statistic,
                p.value=p.value,
                p.value.minp=pvalue.minp,
                maf.cutoff=1))
}
