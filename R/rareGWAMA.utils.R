#' format according to publication standard;
#'
#' @param n The number to be formatted;
#' @param digits how many digits to be retained;
#' @return formatted result a vector;
#' @export 
myFormat <- function(n,digits) {
    res <- rep(NA,length(n));
    if(length(which(abs(n)>1))>0)
        res[which(abs(n)>1)] <- trimws(as.character(round(n[which(abs(n)>1)],digits=digits)));
    if(length(which(abs(n)>0.0001 & abs(n) <=1)))
        res[which(abs(n)>0.0001 & abs(n) <=1 )] <- trimws(as.character(apply(as.matrix(n[which(abs(n)>0.0001 & abs(n) <=1 )]),1,format,digits=digits,scientific=FALSE)))
    if(length(which(abs(n)<=0.0001))>0)
        res[which(abs(n)<=0.0001)] <- trimws(as.character(apply(as.matrix(n[which(abs(n)<=0.0001)]),1,format,digits=digits,scientific=TRUE)));
    return(res);
}

#' compute ABF from summary stat;
#'
#' @param beta genetic effect estimates;
#' @param se the SE of the genetic effect estimates;
#' @return
#' @export
rareGWAMA.abf <- function(beta,se,tau=0.04) {
    lambda <- sqrt((se)^2/((se)^2+tau^2)) * exp(tau^2*(beta)^2/(2*(se)^2*((se)^2+tau^2)))
    Lambda <- sum(lambda)
    ppa <- lambda/Lambda;
    return(ppa);
}

#' get precise value for highly significant p-values;
#'
#' @param statistic the chisq statistic;
#' @return a vector of character-based p-values;
#' @export
get.precise.pval <- function(statistic,df=1) {
    log.pval <- pchisq(statistic,df=df,lower.tail=FALSE,log=TRUE)/log(10);
    log.pval.int <- as.integer(log.pval)-1;
    log.pval.deci <- log.pval-log.pval.int;
    pval.precise <- paste(format(10^log.pval.deci,digit=2),'e',log.pval.int,sep='');
    return(pval.precise);
    
}


#' combine p-values from Cauchy distribution;
#'
#' @param pval a vector of p-values
#' @param weight a vector of weights, the default is equal weight, and the weights add up to 1
#' @return combined p-values;
#' @export
cauchy.p <- function(pval,weight=NULL) {
    if(is.null(weight)) weight <- rep(1/length(pval),length(pval));
    cauchy.stat <- qcauchy(unlist(pval),lower.tail=FALSE);
    cauchy.combined <- sum(cauchy.stat*weight);
    p.value <- pcauchy(cauchy.combined,lower.tail=FALSE);
    return(list(p.value=p.value,
                statistic=cauchy.combined));
}


#' check if the variable exist before reading large files using fread; 
#'
#' @param varName variable name;
#' @param header if the file contain header default is NO;
#' @param file the file name to be read;
#' @export
fread.big <- function(...) {
    extraPar <- list(...);
    varName <- extraPar$varName;
    if(is.null(extraPar$varName)) extraPar$varName <- "res";
    if(is.null(extraPar$header)) extraPar$header <- FALSE;
    if(is.null(extraPar$sep)) extraPar$sep <- "\t";
    if(length(which(ls(envir =.GlobalEnv)==varName))==0) {
        res <- as.data.frame(fread(extraPar$file,header=extraPar$header,sep=extraPar$sep));
        assign(extraPar$varName,res,envir =.GlobalEnv)
    }

}

#' harmonic mean p combination
#'
#' @param pval.mat matrix of p-values; each row being one set of p-values;
#' @return a vector of combined p-values using hmp; 
#' @export
rareGWAMA.hmp <- function(pval.mat) {
    my.hmp <- function(x) {
        ix.rm <- which(is.na(x));
        if(length(ix.rm)>0) x <- x[-ix.rm];
        if(length(x)==0) return(NA);
        if(length(x)>0) {
            x[which(x==0)] <- 1e-321;
            L <- length(x);
            w <- rep(1/L,L);
            tmp <- try(p.hmp(x,w=w,L=L));
            if(class(tmp)=='try-error') tmp <- NA;
            return(tmp)
                        
        }
    }
    pval.hmp <- apply(pval.mat,1,my.hmp);
    return(pval.hmp);
    
}

#' a faster version to convert correlation matrix to covariance matrix;
#'
#' @param cor.mat correlation matrix;
#' @param sd.vec vector of standard deviation;
#' @export
cor2cov <- function(cor.mat,sd.vec) {
    return(t(t(cor.mat*sd.vec)*sd.vec));
}
