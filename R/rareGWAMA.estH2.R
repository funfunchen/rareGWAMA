#' Estimate the variance explained;
#'
#' @param pos a vector of positions;
#' @param beta.vec A vector of estimated genetic effects;
#' @param beta.se A vector of estimated genetic effects se;
#' @param vcf.ref.file the file name of the VCF file
#' @return A list of analysis results;
#' @export 
rareGWAMA.estH2 <- function(pos,beta.vec,beta.se,vcf.ref.file,...) {
    extraPar <- list(...);
    lambda <- extraPar$lambda;
    if(length(lambda)==0) lambda <- 0;
    refGeno <- extraPar$refGeno;
    if(length(refGeno)==0) refGeno <- "GT";
    tabix.range <- get.tabix.range(pos);
    vcfIndv <- refGeno;
    annoType <- "";
    vcfColumn <- c("CHROM","POS","REF","ALT");
    vcfInfo <- NULL;      
    geno.list <- readVCFToListByRange(vcf.ref.file, tabix.range, "", vcfColumn, vcfInfo, vcfIndv)        
    pos.vcf <- paste(geno.list$CHROM,geno.list$POS,sep=":");

    if(refGeno=="DS") {
        gt <- geno.list$DS;
        gt <- matrix(as.numeric(gt),nrow=nrow(gt),ncol=ncol(gt));
    }
    if(refGeno=="GT") {
        gt.tmp <- geno.list$GT;
        gt <- matrix(NA,nrow=nrow(gt.tmp),ncol=ncol(gt.tmp));
        gt[which(gt.tmp=="0/0",arr.ind=T)] <- 0;
        gt[which(gt.tmp=="1/0",arr.ind=T)] <- 1;
        gt[which(gt.tmp=="0/1",arr.ind=T)] <- 1;
        gt[which(gt.tmp=="1/1",arr.ind=T)] <- 2
        gt[which(gt.tmp=="0|0",arr.ind=T)] <- 0;
        gt[which(gt.tmp=="1|0",arr.ind=T)] <- 1;
        gt[which(gt.tmp=="0|1",arr.ind=T)] <- 1;
        gt[which(gt.tmp=="1|1",arr.ind=T)] <- 2
    }
    
    r2 <- rm.na(cor(gt,use='pairwise.complete'));
    r2 <- rm.na(as.matrix(r2[match(pos,pos.vcf),match(pos,pos.vcf)]));
    covG <- rm.na(cov(gt,use='pairwise.complete'));
    covG <- rm.na(as.matrix(covG[match(pos,pos.vcf),match(pos,pos.vcf)]));

    ustat <- beta.vec/beta.se^2;
    vstat <- 1/(beta.se);
    diag.vstat <- matrix(0,nrow=length(vstat),ncol=length(vstat));
    diag(diag.vstat) <- vstat;
    V <- diag.vstat%*%r2%*%diag.vstat;
    V <- regMat(V,lambda);
    beta.joint <- as.vector(ginv(V)%*%ustat);
    
    h2.est <- t(as.matrix(beta.joint))%*%covG%*%as.matrix(beta.joint);
    return(h2.est);

}
