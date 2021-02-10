pvt.core <- function(x,mu,sigma,alternative=c('two.sided','greater'))
  {
    sum.all <- 0;
    if(alternative=='greater')
      {
        dimn <- length(mu);
        for(ii in 1:dimn)
          {
            mat.ii <- combn(1:dimn,ii);
            sum.ii <- 0;
            for(jj in 1:ncol(mat.ii))
              {
                  sigma.jj <- sigma[mat.ii[,jj],mat.ii[,jj]];
                sum.ii <- sum.ii+as.numeric(pmvnorm(upper=rep(-x,ii),lower=rep(-Inf,ii),mean=mu[mat.ii[,jj]],sigma=sigma.jj));
              }
            if(ii%%2==1)
              sum.all <- sum.all+sum.ii;
            if(ii%%2==0)
              sum.all <- sum.all-sum.ii;
          }
        return(sum.all);
      }
    
    if(alternative=='two.sided')
      {
        dimn <- length(mu);
        for(ii in 1:dimn)
          {
            mat.ii <- combn(1:dimn,ii);
            sum.ii <- 0;
            for(jj in 1:ncol(mat.ii))
              {
                  sigma.jj <- sigma[mat.ii[,jj],mat.ii[,jj]];
                  sum.ii <- sum.ii+2*as.numeric(pmvnorm(upper=rep(-sqrt(x),ii),lower=rep(-Inf,ii),mean=mu[mat.ii[,jj]],sigma=sigma.jj));
            }
            if(ii%%2==1)
              sum.all <- sum.all+sum.ii;
            if(ii%%2==0)
              sum.all <- sum.all-sum.ii;
          }
        return(sum.all);
      }
    
  }
#' function to calculate VT p-values;
#'
#' @param statistic the z-score statistic that corresponds to the minimal p-values;
#' @param mu the mean value of the MVT distribution;
#' @param sigm the covariance matrix
#' @export 
pvt <- function(statistic,mu,sigma,alternative=c('greater','less','two.sided'))
  {
      sigma <- make.pos.def(sigma);
     
      
    p.value <- NA;
    if(length(alternative)>1) alternative <- "two.sided";
    if(alternative=='two.sided')
      {
        p.tmp <- 1-pmaxnormsq(statistic,mu,sigma);
        if(!is.na(p.tmp)){
          if(p.tmp>0) p.value <- p.tmp;
          
          if(p.tmp==0)
            {              
              p.value <- pvt.core(statistic,mu,sigma,'two.sided')
            }
        }
      }
    if(alternative=='less')
      {
        p.tmp <- 1-pmaxnormsq(-statistic,mu,sigma);
        if(p.tmp>0) p.value <- p.tmp;
        if(p.tmp==0)
          {
            p.value <- pvt.core(statistic,mu,sigma,'greater')
          }        
      }
    if(alternative=='greater')
      {
        p.tmp <- 1-pmaxnorm(statistic,mu,sigma);
        if(p.tmp>0) p.value <- p.tmp;
        if(p.tmp==0)
          {
            p.value <- pvt.core(statistic,mu,sigma,'greater')
          }
        
      }
    return(p.value);
  }


#' p-value for minimal stat for MVN;
#'
#' @param x is the square of the max statistic for MVN
#' @param mu mean value
#' @param sigma the covariance matrix;
#' @export
pmaxnormsq <- function(x,mu,sigma)
  {
    if(x<0) return(0);
    tmp <- try(pmvnorm(upper=rep(sqrt(x),length(mu)),lower=rep(-sqrt(x),length(mu)),mean=mu,sigma=sigma),silent=TRUE);

    if(class(tmp)=='try-error') return(NA);
    return(as.numeric(tmp));
  }

#' make a matrix being positive definite;
#'
#' @param x.mat input matrix;
#' @export 
make.pos.def <- function(x.mat) {
    x.mat <- (x.mat+t(x.mat))/2;

    x.square <- x.mat%*%t(x.mat);

    res.svd <- svd(x.square);
    d <- matrix(0,nrow=nrow(x.mat),ncol=ncol(x.mat));
    diag(d) <- sqrt(res.svd$d);
    x.mat <- (res.svd$u)%*%d%*%t(res.svd$u);
    return(x.mat);
}
    
