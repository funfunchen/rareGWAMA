#' estimate correlations between multiple twas statistics;
#'
#' @param dat dataset that contains of eQTL weight, score statistics, and covariance matrix
#' @param varList variant list as well as corresponding transcript; 
#' @param af.pca allele frequency PCs
#' @param af.pca.eqtl PCs for the eqtl dataset
#' @export
estimatePredictedBeta.cor <- function(dat,varList,af.pca,af.pca.eqtl) {
    z.mat <- matrix(dat$ustat.mat/dat$vstat.mat,nrow=nrow(dat$ustat.mat),ncol=ncol(dat$ustat.mat));
    w.mat <- matrix(rm.na(sqrt(rm.na(dat$nSample.mat)*rm.na(dat$af.mat)*(1-rm.na(dat$af.mat)))*dat$w.mat),nrow=nrow(dat$ustat.mat),ncol=ncol(dat$ustat.mat));


    z.mat <- matrix(z.mat[dat$ix.rare,],nrow=length(dat$ix.rare));
    w.mat <- matrix(w.mat[dat$ix.rare,],nrow=length(dat$ix.rare));
    beta.mat <- z.mat/w.mat;
    se.mat <- 1/w.mat;
    cov.beta.pred <- matrix(0,nrow=nrow(z.mat),ncol=nrow(z.mat));
    
    af.pca <- as.matrix(af.pca);
    gene.vec <- unique(varList[,1]);
    no.gene <- length(gene.vec);
    no.pc <- ncol(af.pca);
    beta.pred <- rep(0,nrow(z.mat)*no.pc);
    beta.pred.list <- list();
    cov.beta.pred.list <- list();

    gcov.beta.pred <- matrix(0,nrow=nrow(z.mat)*ncol(af.pca),ncol=nrow(z.mat)*ncol(af.pca));
    for(jj1 in 1:nrow(varList)) {
        for(jj2 in jj1:nrow(varList)) {
                
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
                
                if(length(ix.rm)>0) {
                    lm.weight.all1 <- matrix(lm.weight.all1[-ix.rm,-ix.rm],nrow=nrow(lm.weight.all1)-length(ix.rm));
                    lm.weight.all2 <- matrix(lm.weight.all2[-ix.rm,-ix.rm],nrow=nrow(lm.weight.all2)-length(ix.rm));
                    
                    beta.jj1 <- beta.mat[jj1,-ix.rm];
                    se.jj1 <- se.mat[jj1,-ix.rm];
                    beta.jj2 <- beta.mat[jj2,-ix.rm];
                    se.jj2 <- se.mat[jj2,-ix.rm];
                    
                }
                if(length(ix.rm)<ncol(beta.mat)) {
                    cov.beta <- matrix(0,nrow=ncol(beta.mat)-length(ix.rm),ncol=ncol(beta.mat)-length(ix.rm));
                    for(bb in 1:nrow(cov.beta)) {
                        cov.beta[bb,bb] <- se.jj2[bb]*se.jj1[bb]*dat$r2.by.study[[bb]][jj1,jj2];
                    }
                    
                
                    for(ll1 in 1:ncol(af.pca)) {
                        for(ll2 in ll1:ncol(af.pca)) {
                            af.pca.ll1 <- as.matrix(af.pca[,1:ll1]);
                            af.pca.ll2 <- as.matrix(af.pca[,1:ll2]);
                            if(length(ix.rm)>0) {
                                af.pca.ll1 <- matrix(af.pca.ll1[-ix.rm,],ncol=ll1);
                                af.pca.ll2 <- matrix(af.pca.ll2[-ix.rm,],ncol=ll2);
                            }
                            
                            A.L1 <- ginv(t(af.pca.ll1)%*%lm.weight.all1%*%af.pca.ll1)%*%(t(af.pca.ll1)%*%lm.weight.all1);
                            A.R2 <- ginv(t(af.pca.ll2)%*%lm.weight.all2%*%af.pca.ll2)%*%(t(af.pca.ll2)%*%lm.weight.all2);
                            cov.gamma.1.2.tmp <- A.L1%*%cov.beta%*%t(A.R2);
                            gcov.beta.pred[jj1+(ll1-1)*nrow(z.mat),jj2+(ll2-1)*nrow(z.mat)] <- t(af.pca.eqtl[1:ll1])%*%cov.gamma.1.2.tmp%*%(af.pca.eqtl[1:ll2]);
                            gcov.beta.pred[jj2+(ll1-1)*nrow(z.mat),jj1+(ll2-1)*nrow(z.mat)] <- gcov.beta.pred[jj1+(ll1-1)*nrow(z.mat),jj2+(ll2-1)*nrow(z.mat)];
                            
                            ll1.0 <- ll2;ll2.0 <- ll1;
                            
                            af.pca.ll1 <- as.matrix(af.pca[,1:ll1.0]);
                            af.pca.ll2 <- as.matrix(af.pca[,1:ll2.0]);
                            if(length(ix.rm)>0) {
                                af.pca.ll1 <- matrix(af.pca.ll1[-ix.rm,],ncol=ll1.0);
                                af.pca.ll2 <- matrix(af.pca.ll2[-ix.rm,],ncol=ll2.0);
                            }
                            
                            A.L1 <- ginv(t(af.pca.ll1)%*%lm.weight.all1%*%af.pca.ll1)%*%(t(af.pca.ll1)%*%lm.weight.all1);
                            A.R2 <- ginv(t(af.pca.ll2)%*%lm.weight.all2%*%af.pca.ll2)%*%(t(af.pca.ll2)%*%lm.weight.all2);
                            cov.gamma.1.2.tmp <- A.L1%*%cov.beta%*%t(A.R2);
                            gcov.beta.pred[jj1+(ll1.0-1)*nrow(z.mat),jj2+(ll2.0-1)*nrow(z.mat)] <- t(af.pca.eqtl[1:ll1.0])%*%cov.gamma.1.2.tmp%*%(af.pca.eqtl[1:ll2.0]);
                            gcov.beta.pred[jj2+(ll1.0-1)*nrow(z.mat),jj1+(ll2.0-1)*nrow(z.mat)] <- gcov.beta.pred[jj1+(ll1.0-1)*nrow(z.mat),jj2+(ll2.0-1)*nrow(z.mat)];
                                                        
                            
                        }
                    }
                }
            }
        }
        
    }
    
    
    return(list(z.mat=z.mat,
                beta.mat=beta.mat,
                se.mat=se.mat,
                gcov.beta.pred=gcov.beta.pred));
    
}

#' estimate TWAS Statistic Correlation 
#'
#' @param dat dataset;
#' @param varList variant list with genes being the first column and the variant positions (in the format of 1:123) in the second column;
#' @param af.pca MDS component for the GWAS dataset;
#' @param af.pca.eqtl MDS component for the eqtl dataset
#' @export
estimateTWAS.cor <- function(dat,varList,af.pca,af.pca.eqtl,maxNumVar=500) {
    if(nrow(varList)>maxNumVar) res.beta.pred <- estimatePredictedBeta.cor.fast(dat,varList,af.pca,af.pca.eqtl);
    if(nrow(varList)<maxNumVar) res.beta.pred <- estimatePredictedBeta.cor(dat,varList,af.pca,af.pca.eqtl);
    
    af.pca <- as.matrix(af.pca);
    gene.vec <- unique(varList[,1]);
    gcov.twas <- matrix(nrow=length(gene.vec)*ncol(af.pca),ncol=length(gene.vec)*ncol(af.pca));
    for(gg1 in 1:length(gene.vec)) {
        for(gg2 in 1:length(gene.vec)) {
            ix.var1 <- which(varList[,1]==gene.vec[gg1]);
            ix.var2 <- which(varList[,1]==gene.vec[gg2]);
            for(kk1 in 1:ncol(af.pca)) {
                for(kk2 in 1:ncol(af.pca)) {
                    cov.tmp <- res.beta.pred$gcov.beta.pred[(kk1-1)*nrow(res.beta.pred$z.mat)+ix.var1,(kk2-1)*nrow(res.beta.pred$z.mat)+ix.var2];
                    gcov.twas[(kk1-1)*length(gene.vec)+gg1,(kk2-1)*length(gene.vec)+gg2] <- t(dat$twas.weight[ix.var1])%*%cov.tmp%*%(dat$twas.weight[ix.var2]);
                }
            }
            
        }
    }
    return(list(gcov.twas=gcov.twas));
}

#' estimate the correlation between minimal p-value statistic
#'
#' @param gcov.twas
#' @param no.pc no.pc;default is 3;
#' @export
estimateTESLA.cor <- function(gcov.twas,no.pc=3) {
    no.gene <- nrow(gcov.twas)/(no.pc+1);
    gcor.twas <- rm.na(cov2cor(gcov.twas));
    diag(gcor.twas) <- 1;
    gcor.twas <- (gcor.twas+t(gcor.twas))/2;
    n.simu <- 20000;
    z.tmp <- rmvnorm(n.simu,sigma=gcor.twas);
    z.absmax <- matrix(nrow=nrow(z.tmp),ncol=ncol(z.tmp)/4);
    for(ii in 1:(ncol(gcor.twas)/4)) {
        z.absmax[,ii] <- apply(z.tmp[,((ii-1)*4+1):(ii*4)],1,function(x) return(max(abs(x))));
    }
    cor.tesla <- cor(z.absmax)
    return(cor.tesla);
    
    
}


#' a faster way (with some approximations) to estimate the correlations between multiple TWAS statistics
#' @param dat dataset that contains of eQTL weight, score statistics, and covariance matrix
#' @param varList variant list as well as corresponding transcript; 
#' @param af.pca allele frequency PCs
#' @param af.pca.eqtl PCs for the eqtl dataset
#' @export
estimatePredictedBeta.cor.fast <- function(dat,varList,af.pca,af.pca.eqtl) {
    z.mat <- matrix(dat$ustat.mat/dat$vstat.mat,nrow=nrow(dat$ustat.mat),ncol=ncol(dat$ustat.mat));
    w.mat <- matrix(rm.na(sqrt(rm.na(dat$nSample.mat)*rm.na(dat$af.mat)*(1-rm.na(dat$af.mat)))*dat$w.mat),nrow=nrow(dat$ustat.mat),ncol=ncol(dat$ustat.mat));
    z.mat <- matrix(z.mat[dat$ix.rare,],nrow=length(dat$ix.rare));
    w.mat <- matrix(w.mat[dat$ix.rare,],nrow=length(dat$ix.rare));
    beta.mat <- z.mat/w.mat;
    se.mat <- 1/w.mat;
    cov.beta.pred <- matrix(0,nrow=nrow(z.mat),ncol=nrow(z.mat));
    
    af.pca <- as.matrix(af.pca);
    gene.vec <- unique(varList[,1]);
    no.gene <- length(gene.vec);
    no.pc <- ncol(af.pca);
    beta.pred <- rep(0,nrow(z.mat)*no.pc);
    beta.pred.list <- list();
    cov.beta.pred.list <- list();

    gcov.beta.pred <- matrix(0,nrow=nrow(z.mat)*ncol(af.pca),ncol=nrow(z.mat)*ncol(af.pca));
    af.pca.eqtl.mat <- matrix(0,nrow=length(af.pca.eqtl),ncol=length(af.pca.eqtl));
    for(ii in 1:length(af.pca.eqtl))
        af.pca.eqtl.mat[1:ii,ii] <- af.pca.eqtl[1:ii];
    for(jj1 in 1:nrow(varList)) {
        for(jj2 in jj1:nrow(varList)) {
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
                
                if(length(ix.rm)>0) {
                    lm.weight.all1 <- matrix(lm.weight.all1[-ix.rm,-ix.rm],nrow=nrow(lm.weight.all1)-length(ix.rm));
                    lm.weight.all2 <- matrix(lm.weight.all2[-ix.rm,-ix.rm],nrow=nrow(lm.weight.all2)-length(ix.rm));
                    
                    beta.jj1 <- beta.mat[jj1,-ix.rm];
                    se.jj1 <- se.mat[jj1,-ix.rm];
                    beta.jj2 <- beta.mat[jj2,-ix.rm];
                    se.jj2 <- se.mat[jj2,-ix.rm];
                    
                }
                if(length(ix.rm)<ncol(beta.mat)) {
                    cov.beta <- matrix(0,nrow=ncol(beta.mat)-length(ix.rm),ncol=ncol(beta.mat)-length(ix.rm));
                    for(bb in 1:nrow(cov.beta)) {
                        cov.beta[bb,bb] <- se.jj2[bb]*se.jj1[bb]*dat$r2.by.study[[bb]][jj1,jj2];
                    }
                    af.pca.1 <- as.matrix(af.pca);
                    af.pca.2 <- as.matrix(af.pca);
                    if(length(ix.rm)>0) {
                        af.pca.1 <- matrix(af.pca.1[-ix.rm,],ncol=ncol(af.pca));
                        af.pca.2 <- matrix(af.pca.2[-ix.rm,],ncol=ncol(af.pca));
                    }
                    A.L1 <- ginv(t(af.pca.1)%*%lm.weight.all1%*%af.pca.1)%*%(t(af.pca.1)%*%lm.weight.all1);
                    A.R2 <- ginv(t(af.pca.2)%*%lm.weight.all2%*%af.pca.2)%*%(t(af.pca.2)%*%lm.weight.all2);
                    cov.gamma.1.2.tmp <- A.L1%*%cov.beta%*%t(A.R2);
                    gcov.beta.pred[jj1+((1:length(af.pca.eqtl))-1)*nrow(z.mat),jj2+((1:length(af.pca.eqtl))-1)*nrow(z.mat)] <- t(af.pca.eqtl.mat)%*%cov.gamma.1.2.tmp%*%(af.pca.eqtl.mat)
                    gcov.beta.pred[jj2+((1:length(af.pca.eqtl))-1)*nrow(z.mat),jj1+((1:length(af.pca.eqtl))-1)*nrow(z.mat)] <- gcov.beta.pred[jj1+((1:length(af.pca.eqtl))-1)*nrow(z.mat),jj2+((1:length(af.pca.eqtl))-1)*nrow(z.mat)];
                
            
                }
            }
        }
     
    }
    return(list(z.mat=z.mat,
                beta.mat=beta.mat,
                se.mat=se.mat,
                gcov.beta.pred=gcov.beta.pred))
}
