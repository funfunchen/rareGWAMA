#' Gram-Schmitt othogonalization
#'
#' @param A a matrix; 
#' @param ip a bivariate function that serves as inner product;
#' @return 
#' @export
gram.schmitt <- function(A,ip,...) {
    p <- list(...);
    B <- A;
    for(ii in 2:ncol(A)) {
        proj.ii <- 0;
        for(jj in 1:(ii-1)) {
            if(ip(B[,jj],B[,jj],p$w)>1e-10) 
                proj.ii <- proj.ii+rm.na(ip(A[,ii],B[,jj],p$w)/ip(B[,jj],B[,jj],p$w)*B[,jj]);
        }
        B[,ii] <- B[,ii]-proj.ii;
    }
    B <- rm.na(B);
    return(B);
}

#' inner product;
#'
#' @param v1 vector 1;
#' @param v2 vector 2
#' @param w weight matrix;
#' @return
#' @export
ip.w <- function(v1,v2,w) {
    v1 <- as.matrix(v1);
    v2 <- as.matrix(v2);
    w <- as.matrix(w);
    return(as.numeric(t(v1)%*%w%*%v2));
}
