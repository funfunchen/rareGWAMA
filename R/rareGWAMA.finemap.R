#' convert z-score to posterior inclusion probability;
#'
#' @param z a vector of z-scores for the locus;
#' @param cov.z the covariance between z-scores; order has to be the same as in z; 
#' @param gene.vec a vector of gene names that matches z; 
#' @param cs.level credible set level default is 90%;
#' @param tau prior on the effect sizes. default is 10000 which is equivalent to a non-informative prior;
#' @export
zscore2pip <- function(z,cov.z,gene.vec,cs.level=.90,tau=10000) {
    ppi <- zscore.abf(z,tau);
    ix.order <- order(ppi,decreasing=TRUE);
    ix.cs95 <- min(which(cumsum(ppi[ix.order])>cs.level));
    numStability <- 1;
    if(is.infinite(ix.cs95)) {
        ix.cs95 <- 1:length(ix.order);
        numStability <- 0;
    }
    ix.gene <- ix.order[1:ix.cs95];
    gene.cs95 <- gene.vec[ix.gene];
    ix.top <- ix.gene[1];
    pip <- ppi;
    ppi <- as.data.frame(cbind(gene.vec,pip));
    return(list(pip=ppi,
                numStability=numStability,
                gene.cs=ppi[which(ppi[,1]%in%gene.cs95),],
                gene.top=gene.vec[ix.top]));
}

#' fine mapping main function;
#'
#' @param pval a vector p-values;
#' @param cov.z the covariance matrix between converted z-scores;
#' @param gene.vec a vector of gene names;
#' @param max.iter the max number of iteration conducted; default is 3, which means a maximum of 3 causal genes;
#' @param cs.level credible set level; default is set to 90%;
#' @param alpha significant threshold used to determine indpendently associated genes; default is 5e-8;
#' @param tau prior effect sizes; default is 10000, equivalent to non-informative prior;
#' @export
zscore.finemap <- function(pval,cov.z,gene.vec,max.iter=3,cs.level=.90,alpha=5e-8,tau=10000) {
    ix.0 <- which(pval==0);
    if(length(ix.0)>0) pval[ix.0] <- 1e-320;
    ix.rm <- which(pval==1);
    if(length(ix.rm)>0) {
        pval <- pval[-ix.rm];
        cov.z <- as.matrix(cov.z[-ix.rm,-ix.rm]);
        gene.vec <- gene.vec[-ix.rm];
    }
    if(length(gene.vec)==0) stop("no gene");
   z <- qnorm(pval,lower.tail=FALSE);
   no.iter <- 1;
   res <- list();
   res[[1]] <- zscore2pip(z,cov.z,gene.vec,tau=tau);
   rr <- 2;
   z.in <- z;
   p.z.in <- pnorm(z.in,lower.tail=FALSE);
   p.z.in[p.z.in == "0"] = .Machine$double.xmin
   cov.z.in <- cov.z;
   gene.vec.in <- gene.vec;
   while(rr<=max.iter & min(p.z.in)<alpha) {
 
       ix.top <- which.max(abs(z.in));
       z.in <- z.in[-ix.top]-matrix(cov.z.in[-ix.top,ix.top],ncol=1)%*%ginv(cov.z.in[ix.top,ix.top])%*%z.in[ix.top];
       cov.z.in <- as.matrix(cov.z.in[-ix.top,-ix.top])-matrix(cov.z.in[-ix.top,ix.top],ncol=length(ix.top))%*%ginv(cov.z.in[ix.top,ix.top])%*%matrix(cov.z.in[ix.top,-ix.top],nrow=length(ix.top));
       p.z.in <- pnorm(z.in/sqrt(abs(diag(cov.z.in))),lower.tail=FALSE);
 
       gene.vec.in <- gene.vec.in[-ix.top];
       if(min(p.z.in)<alpha) {
           res[[rr]] <- zscore2pip(z.in,cov.z.in,gene.vec.in,cs.level);
       }
       rr <- rr+1;
   }
   res.tmp <- res;
   res <- list();
   kk <- 1;
   for(ii in 1:length(res.tmp)) {
       if(res.tmp[[ii]]$numStability==1) {
           res[[kk]] <- res.tmp[[ii]];
           kk <- kk+1;
       }
   }

   return(res);
}

#' compute ABF from summary stat;
#' 
#' @param beta genetic effect estimates;
#' @param se the SE of the genetic effect estimates;
#' @return
#'  @export
zscore.abf <- function(z,tau=10000) {
    beta <- z;se <- 1;
    lambda <- sqrt((se)^2/((se)^2+tau^2)) * exp(tau^2*(beta)^2/(2*(se)^2*((se)^2+tau^2)))
    lambda[is.infinite(lambda)] = .Machine$double.xmax
    Lambda <- sum(lambda)
    ppi <- lambda/Lambda;
    return(ppi);
}

