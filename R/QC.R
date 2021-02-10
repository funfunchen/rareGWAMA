
#' Impute missing summary association statistics assuming
#'
#' @param ustat.list the score statistics;
#' @param vstat.list the vstat list;
#' @param cov.mat.list the list of the covariance matrix
#' @param N.mat the matrix of sample sizes;each row for a study and each column for a variant site;
#' @export
imputeMeta <- function(ustat.list,vstat.list,cov.mat.list,N.mat,beta.vec=NULL,ix.known,lambda=0.01) {
    U.imp <- 0;nSample.U <- 0;
    covG <- matrix(0,nrow=nrow(cov.mat.list[[1]]),ncol=ncol(cov.mat.list[[1]]));
    nSample.covG <- covG;
    N.mat.imp <- N.mat;
    U.meta <- 0;
    for(ii in 1:length(ustat.list))
    {
        N.mat.imp[ii,] <- max(rm.na(N.mat[ii,]));
        U.meta <- U.meta+rm.na(ustat.list[[ii]]);
        
        nSample.U <- nSample.U+rm.na(N.mat[ii,]);
        for(jj in 1:length(ustat.list[[1]]))
        {
            for(kk in 1:jj)
            {
                covG[jj,kk] <- covG[jj,kk]+rm.na(sqrt(N.mat[ii,jj]*N.mat[ii,kk])*cov.mat.list[[ii]][jj,kk]);
                covG[kk,jj] <- covG[jj,kk];
                nSample.covG[jj,kk] <- nSample.covG[jj,kk]+sqrt(rm.na(N.mat[ii,jj])*rm.na(N.mat[ii,kk]));
                nSample.covG[kk,jj] <- nSample.covG[jj,kk];
            }
        }
    }
    
    covG.ori <- covG;    
    covG <- (covG/nSample.covG);
    corG <- cov2cor(covG);
    Id <- matrix(0,nrow=nrow(covG),ncol=ncol(covG));
    diag(Id) <- 1;
    
    N.meta.ori <- apply(N.mat,2,sum,na.rm=T);
    N.meta <- rep(max(N.meta.ori),length(N.meta.ori));
    U.meta.imp <- U.meta*(rm.na(N.meta/N.meta.ori));
    V.tmp <- diag(sqrt(N.meta))%*%covG%*%diag(sqrt(N.meta))
    beta.imp <- ginv(V.tmp)%*%U.meta.imp;
    scalar <- matrix(0,nrow=length(ustat.list[[1]]),ncol=length(ustat.list[[1]]));
    diag(scalar) <- (rm.na(N.meta/N.meta.ori));    
    cov.U.meta.imp <- scalar%*%covG.ori%*%scalar;
    cov.beta.imp <- ginv(V.tmp)%*%cov.U.meta.imp%*%ginv(V.tmp);
    V.meta.imp <- ginv(cov.beta.imp);
    
    U.meta.imp <- V.meta.imp%*%beta.imp;
    
    return(list(covG=covG,
                nSample.covG=nSample.covG,
                N.mat.imp=N.mat.imp,
                N.meta.ori=N.meta.ori,
                N.meta=N.meta,
                scalar.diag=diag(scalar),
                U.meta.imp=U.meta.imp,
                V.meta.imp=V.meta.imp,
                N.meta.imp=N.meta));
}



#' Impute missing summary association statistics assuming
#'
#' @param ustat.list the score statistics;
#' @param vstat.list the vstat list;
#' @param cov.mat.list the list of the covariance matrix
#' @param N.mat the matrix of sample sizes;each row for a study and each column for a variant site;
#' @export
imputeConditional <- function(ustat.list,vstat.list,cov.mat.list,N.mat,beta.vec=NULL,ix.candidate,ix.known) {
    U.imp <- 0;nSample.U <- 0;
    covG <- matrix(0,nrow=nrow(cov.mat.list[[1]]),ncol=ncol(cov.mat.list[[1]]));
    nSample.covG <- covG;
    N.mat.imp <- N.mat;
    U.meta <- 0;
    for(ii in 1:length(ustat.list))
    {
        N.mat.imp[ii,] <- max(rm.na(N.mat[ii,]));
        U.meta <- U.meta+rm.na(ustat.list[[ii]]);
        
        nSample.U <- nSample.U+rm.na(N.mat[ii,]);
        for(jj in 1:length(ustat.list[[1]]))
        {
            for(kk in 1:jj)
            {
                covG[jj,kk] <- covG[jj,kk]+rm.na(sqrt(N.mat[ii,jj]*N.mat[ii,kk])*cov.mat.list[[ii]][jj,kk]);
                covG[kk,jj] <- covG[jj,kk];
                nSample.covG[jj,kk] <- nSample.covG[jj,kk]+sqrt(rm.na(N.mat[ii,jj])*rm.na(N.mat[ii,kk]));
                nSample.covG[kk,jj] <- nSample.covG[jj,kk];
            }
        }
    }
    U.meta <- U.meta/nSample.U;
    U.XY <- U.meta[ix.candidate];
    U.ZY <- U.meta[ix.known];
    covG <- covG/nSample.covG;
    V.XZ <- matrix(covG[ix.candidate,ix.known],nrow=length(ix.candidate),ncol=length(ix.known));
    V.ZZ <- matrix(covG[ix.known,ix.known],nrow=length(ix.known),ncol=length(ix.known));
    V.XX <- matrix(covG[ix.candidate,ix.candidate],nrow=length(ix.candidate),ncol=length(ix.candidate));
    conditional.ustat <- U.XY-V.XZ%*%ginv(V.ZZ)%*%U.ZY;
    
    var.U.XY <- V.XX/(nSample.covG[ix.candidate,ix.candidate]);
    var.U.ZY <- V.ZZ/(nSample.covG[ix.known,ix.known]);
    cov.U.XY.U.ZY <- V.XZ/matrix(nSample.covG[ix.candidate,ix.known],nrow=length(ix.candidate),ncol=length(ix.known));
    conditional.V <- var.U.XY+V.XZ%*%ginv(V.ZZ)%*%var.U.ZY%*%ginv(V.ZZ)%*%t(V.XZ)-cov.U.XY.U.ZY%*%t(V.XZ%*%ginv(V.ZZ))-(V.XZ%*%ginv(V.ZZ))%*%t(cov.U.XY.U.ZY);
    conditional.V <- regMat(conditional.V,0.1);

    
    return(list(conditional.ustat=conditional.ustat,
                conditional.V=conditional.V));
}
#' regularize matrix;
#' @param M matrix
#' @param lambda regularization parameter
#' @export
regMat <- function(M,lambda) {
    cor.tmp <- rm.na(cov2cor(M));
    diag(cor.tmp) <- 1;
   
    sd.mat <- matrix(0,nrow=nrow(M),ncol=ncol(M));
    id.mat <- matrix(0,nrow=nrow(M),ncol=ncol(M));
    diag(id.mat) <- 1;
    diag(sd.mat) <- sqrt(abs(diag(M)));
    cor.tmp <- cov2cor(cor.tmp+lambda*id.mat);
    M.reg <- sd.mat%*%(cor.tmp)%*%sd.mat;
    return(M.reg);
}
