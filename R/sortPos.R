#' Sort genomic positions;
#'
#' @param pos Genomic positions in the form of 1:12345;
#' @return Return sorted genomic positions;
#' @export
sortPos <- function(pos) {   
    pos.in <- pos;
    chrpos.mat <- matrix(unlist(strsplit(pos,split=":")),ncol=2,byrow=T);
    chr <- as.numeric(chrpos.mat[,1])
    pos <- as.numeric(chrpos.mat[,2]);
    ix.order <- order(chr,pos);
    chr <- chr[ix.order];
    pos <- pos[ix.order];
    return(list(pos=paste(chr,pos,sep=":"),
                ix=ix.order));
}
#' Get genomic distance between two variant positions;
#'
#' @param pos1 Position for variant 1;
#' @param pos2 Position for variant 2;
#' @return Return variant positions;
#' @export
genomic.dist <- function(pos1,pos2) {   
    pos1.in <- pos1;
    pos2.in <- pos2;
    a <- strsplit(pos1.in,split=":");
    b <- strsplit(pos2.in,split=":");
    chr1 <- a[[1]][1];pos1 <- a[[1]][2];
    chr2 <- b[[1]][1];pos2 <- b[[1]][2];
    if(chr1!=chr2) return(Inf);
    if(chr1==chr2) return(abs(as.numeric(pos1)-as.numeric(pos2)));
}
