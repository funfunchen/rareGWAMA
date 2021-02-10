#' lassosum with essential predictors;
#'
#' @param b.vec genetic effect estimates;
#' @param s.vec standard error of genetic effect estimates;
#' @param r2.mat residual errors;
#' @param n sample size;
#' @param lambda l1 penalty parameter
#' @param alpha l2 penalty parameter
#' @export 
lasso.sum.ess <- function(b.vec,s.vec,r2.mat,n,group,lambda,alpha) {
    z.vec <- b.vec/s.vec;
    u.vec <- b.vec/s.vec^2;
    v.vec <- 1/s.vec;
    cov.mat <- cor2cov(r2.mat,v.vec);
    #beta.vec <- ginv(cov.mat)%*%u.vec;
    beta.vec <- b.vec;
    beta0.vec <- rep(0,length(z.vec));
    while(sum(abs(beta.vec-beta0.vec))>1e-5) {
        beta0.vec <- beta.vec;
        for(jj in 1:length(z.vec)) {
            x.j.times.r <- u.vec[jj]-sum(cov.mat[jj,-jj]*beta.vec[-jj]);
            x.j.times.x <- cov.mat[jj,jj];
            if(x.j.times.r>n*alpha[group[jj]] )
                beta.vec[jj] <- (x.j.times.r-n*(alpha[group[jj]]))/(x.j.times.x + 2*n*lambda[group[jj]]);
            if(x.j.times.r < -n*alpha[group[jj]] )
                beta.vec[jj] <- (x.j.times.r+n*(alpha[group[jj]]))/(x.j.times.x + 2*n*lambda[group[jj]]);
            if(x.j.times.r < n*alpha[group[jj]] & x.j.times.r > -n*alpha[group[jj]])
                beta.vec[jj] <- 0;
            
        }
    }
    return(beta.vec);
}
