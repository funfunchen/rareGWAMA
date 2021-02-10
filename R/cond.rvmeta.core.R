cond.rvmeta.core <- function(X.T.times.Y.all,X.T.times.X.all,maf.vec,N,var.Y,ix.X1,ix.X2,rv.test=c('GRANVIL',"WSS","VT","SKAT"),alternative,no.boot,alpha,extra.par=list())
  {
    #W~X1 Z~X2 X1|X2
    X.T.times.Y <- X.T.times.Y.all;
    X.T.times.X <- X.T.times.X.all;
    X2.T.times.Y <- X.T.times.Y[ix.X2];
    X1.T.times.Y <- X.T.times.Y[ix.X1];
    X2.T.times.X2 <- matrix((X.T.times.X[ix.X2,ix.X2]),nrow=length(ix.X2));
    X2.T.times.X1 <- matrix(X.T.times.X[ix.X2,ix.X1],nrow=length(ix.X2));
    X1.T.times.X2 <- matrix(X.T.times.X[ix.X1,ix.X2],nrow=length(ix.X1));
    X1.T.times.X1 <- matrix(X.T.times.X[ix.X1,ix.X1],nrow=length(ix.X1));
    gamma.est <- as.vector(ginv(X2.T.times.X2)%*%X2.T.times.Y);
    sigma.sq.est <- var.Y-2/N*t(gamma.est)%*%X2.T.times.Y+1/N*t(gamma.est)%*%X2.T.times.X2%*%(gamma.est);
    U <- X1.T.times.Y-X1.T.times.X2%*%gamma.est;
    V <- X1.T.times.X1-X1.T.times.X2%*%ginv(X2.T.times.X2)%*%X2.T.times.X1;
    N.list <- list();var.Y.list <- list();N.list <- list();X.T.times.Y.list <- list();maf.vec.list <- list();X.T.times.X.list <- list();cov.mat.list <- list();
    N.list[[1]] <- N;
    var.Y.list[[1]] <- as.numeric(sigma.sq.est);
    X.T.times.Y.list[[1]] <- U;
    X.T.times.X.list[[1]] <- V;
    cov.mat.list[[1]] <- V/N;
    maf.vec.new <- maf.vec[ix.X1];
    maf.vec.list[[1]] <- maf.vec.new;
    if(rv.test=='generic')
      {
          return(list(X.T.times.Y=U,
                      X.T.times.X=V*as.numeric(sigma.sq.est),
                      cov.mat=V*as.numeric(sigma.sq.est)/N,
                      maf.vec=maf.vec.new));
      }
    return(res);
}
#' This is the function to obtain conditional score statistics;
#'
#' @param ustat The non-standardized score statistics
#' @param X.T.times.X The covariance matrix for ustat;
#' @param ix.candidate the indics of candidate variant;
#' @param ix.known the indices of known variants that are conditioned on;
#' @return the score statistic and their covariance matrix;
#' @export
get.conditional.score.stat <- function(ustat,X.T.times.X,N,ix.candidate,ix.known,impState=NULL)
{
    X.T.times.Y <- ustat;
    ix.X1 <- ix.candidate;
    ix.X2 <- ix.known;
    var.Y <- 1;
    N <- as.integer(mean(N,na.rm=T));
    if(length(which(is.na(ustat)))+length(which(is.na(X.T.times.X)))>0)
    {
        return(list(conditional.ustat=as.matrix(NA),
                    conditional.V=as.matrix(NA)));
    }
    X2.T.times.Y <- X.T.times.Y[ix.X2];
    X1.T.times.Y <- X.T.times.Y[ix.X1];
    ix.type <- ix.X2;
    ix.X2.type <- intersect(ix.type,ix.X2);
    X2.T.times.X1 <- matrix(X.T.times.X[ix.X2,ix.X1],nrow=length(ix.X2),ncol=length(ix.X1));
    X1.T.times.X2 <- matrix(X.T.times.X[ix.X1,ix.X2],nrow=length(ix.X1),ncol=length(ix.X2));
    X1.T.times.X1 <- matrix(X.T.times.X[ix.X1,ix.X1],nrow=length(ix.X1),ncol=length(ix.X1));
    X2.T.times.X2 <- matrix((X.T.times.X[ix.X2,ix.X2]),nrow=length(ix.X2),ncol=length(ix.X2));

    if(length(ix.X2.type)>0) {
        X2.T.times.X2.type <- matrix((X.T.times.X[ix.X2.type,ix.X2.type]),nrow=length(ix.X2.type),ncol=length(ix.X2.type));
        X2.T.times.X1.type <- matrix(X.T.times.X[ix.X2.type,ix.X1],nrow=length(ix.X2.type),ncol=length(ix.X1));
        X1.T.times.X2.type <- matrix(X.T.times.X[ix.X1,ix.X2.type],nrow=length(ix.X1),ncol=length(ix.X2.type));
        
        X2.T.times.Y.type <- X.T.times.Y[ix.X2.type];
    }
    gamma.est <- as.vector(ginv(X2.T.times.X2)%*%X2.T.times.Y);
    sigma.sq.est <- var.Y-(t(X2.T.times.Y)%*%ginv(X2.T.times.X2)%*%X2.T.times.Y)/N;
    U <- X1.T.times.Y-X1.T.times.X2%*%gamma.est;
    if(length(ix.X2.type)>0) {
        V <- X1.T.times.X1-X1.T.times.X2.type%*%ginv(X2.T.times.X2.type)%*%X2.T.times.X1.type;
    }
    if(length(ix.X2.type)==0) {
        V <- X1.T.times.X1;        
    }
    V <- V*as.numeric(sigma.sq.est);

    return(list(conditional.ustat=as.matrix(U),
                sigma.sq.est=sigma.sq.est,
                conditional.beta.est=ginv(V)%*%U,
                conditional.beta.var=ginv(V),
                nSample=rep(N,length(ustat)),
                conditional.V=as.matrix(V)));
    
}

#' This is the function to obtain conditional score statistics when conditional on the genes that are detected by SKAT test;
#'
#' @param ustat The non-standardized score statistics
#' @param X.T.times.X The covariance matrix for ustat;
#' @param ix.candidate the indics of candidate variant;
#' @param ix.known the indices of known variants that are conditioned on;
#' @return the score statistic and their covariance matrix;
#' @export
get.conditional.score.stat.RE <- function(ustat,X.T.times.X,N,ix.candidate,ix.known) {
    ustat.candidate <- ustat[ix.candidiate];
    ustat.known <- ustat[ix.known];
    conditional.ustat <- ustat.candidate-X.T.times.X[ix.candidate,ix.known]%*%ginv(X.T.times.X[ix.known,ix.known])%*%ustat.known;
   
}

