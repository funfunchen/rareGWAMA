#' Perform conditional analysis for single variant tests
#'
#' @param candidate.variant name of the range to be analyzed (for example, it can be a gene name e.g. APOE)
#' @param score.stat.file files of score statistics
#' @param cov.file covariance matrix files
#' @param known.variant.vec range of candidate variant, expressed in a vector, e.g. c("1:12345","1:234567");
#' @param known.variant.range range of known variant, expressed in tabix format, e.g. c("1:123456","1:2345678");
#' @param test test of rare variant tests
#' @param maf.cutoff Cutoffs of MAF used for determining rare variants
#' @param alternative Alternative hypothesis to be tested
#' @param out.digits The number of digits used in the output
#' @param callrate.cutoff Cutoff of call rates. Sites with callrates lower than the cutoff will be labeled as missing
#' @param hwe.cutoff Cutoff of HWE p-values. Sites with HWE pvalues lower than the cutoff will be labeled as missing
#' @return return a list of meta-analysis results
#' @export
conditional.rareMETALS.single <- function(candidate.variant.vec,score.stat.file,cov.file,known.variant.vec,maf.cutoff,no.boot=10000,alternative=c('two.sided','greater','less'),ix.gold=1,out.digits=4,callrate.cutoff=0,hwe.cutoff=0,p.value.known.variant.vec=NA,anno.known.variant.vec=NA,anno.candidate.variant.vec=NA,knownCoding='identity')
  {
    alpha <- 0.05;maf.cutoff <- 1;
    if(sum(is.na(p.value.known.variant.vec))>0) p.value.known.variant.vec <- rep("N/A",length(p.value.known.variant.vec));
    if(sum(is.na(anno.known.variant.vec))>0) anno.known.variant.vec <- rep("N/A",length(known.variant.vec));
    if(sum(is.na(anno.candidate.variant.vec))>0) anno.candidate.variant.vec <- rep("N/A",length(candidate.variant.vec));

    gene.file <- "refFlat_hg19.txt.gz";
    
    test <- "SINGLE";
    gene.name <- anno.candidate.variant.vec;
    extra.pars <- list(alpha=alpha,ix.gold=ix.gold,out.digits=out.digits,QC.par=list(callrate.cutoff=callrate.cutoff,hwe.cutoff=hwe.cutoff));
    raw.data.all <- list();
    for(ii in 1:length(candidate.variant.vec))
      {
        tabix.range <- get.tabix.range(c(candidate.variant.vec[ii],known.variant.vec));
        capture.output(raw.data.all[[ii]] <- rvmeta.readDataByRange( score.stat.file, cov.file,tabix.range));
      }
    
    res.null <- list(gene.name=NA,
                     p.value=NA,
                     statistic=NA,
                     no.var=NA,
                     no.sample=NA,
                     anno=NA,
                     ref=NA,
                     alt=NA,
                     maf.cutoff=NA,
                     direction.burden.by.study=NA,
                     direction.meta.single.var=NA,
                     beta1.est=NA,
                     no.site=NA,
                     beta1.sd=NA,
                     hsq.est=NA,
                     pos=NA);
    if(length(raw.data.all)==0)
      {
        return(list(res.null));        
      }
    res <- list(list(p.value=NA,statistic=NA,no.var=NA,no.sample=NA));
    p.value <- double(0);
    ref <- character(0);
    alt <- ref;
    direction.by.study <- character(0);
    ix.gold <- extra.pars$ix.gold;
    if(length(extra.pars$ix.gold)==0) {
      ix.gold <- 1;
    }
    out.digits <- extra.pars$out.digits;
    if(length(out.digits)==0) out.digits <- 4;
    statistic <- double(0);pos <- integer(0);anno <- character(0);direction <- integer(0);res.maf.vec <- double(0);beta1.est.vec <- double(0);beta1.sd.vec <- double(0);
    gene.name.out <- 0;p.value.out <- 0;statistic.out <- 0;no.site.out <- 0;beta1.est.out <- 0;beta1.sd.out <- 0;maf.cutoff.out <- 0;direction.burden.by.study.out <- 0;direction.meta.single.var.out <- 0;pos.ref.alt.out <- 0;pos.ref.alt.known.out <- 0;p.value.known.out <- 0;anno.known.out <- 0;
    pos.single.out <- character(0);p.value.single.out <- double(0);ref.single.out <- character(0);alt.single.out <- character(0);anno.single.out <- character(0);maf.single.out <- double(0);beta1.est.single.out <- double(0);beta1.sd.single.out <- double(0);pos.ref.alt.known.single.out <- character(0);direction.single.out <- character(0);p.value.known.single.out <- character(0);anno.known.single.out <- character(0);
    for(kk in 1:length(raw.data.all))
      {

        candidate.variant <- candidate.variant.vec[kk];
        anno.candidate.variant <- anno.candidate.variant.vec[kk];
        raw.data <- raw.data.all[[kk]][[1]];
        if(length(extra.pars$QC.par)>0) 
          raw.data <- QC(raw.data,extra.pars$QC.par,cov=1);        
        ix.var <- 1:length(raw.data$ustat[[1]]);
        ix.var <- c(ix.var,match(known.variant.vec,raw.data$pos));
        ix.var <- sort(unique(ix.var));
        ix.tmp <- match(known.variant.vec,(raw.data$pos)[ix.var]);
        ix.var <- c(ix.var[-ix.tmp],ix.var[ix.tmp]);
        if(length(ix.var)==length(ix.tmp))
          {
            res.null$gene.name <- gene.name[kk];
            
            res[[kk]] <- res.null;
            if(test=="SINGLE")
              {
                ix.candidate <- which(raw.data$pos==candidate.variant);
                p.value.single.out[kk] <- NA;
                pos.single.out[kk] <- c(paste(raw.data$pos[ix.candidate],sep='',collapse=","));
                direction.single.out[kk] <- NA;
                ref.single.out[kk] <- paste(raw.data$ref[[ix.gold]][ix.candidate],collapse=",");
                alt.single.out[kk] <- paste(raw.data$alt[[ix.gold]][ix.candidate],collapse=",");
                beta1.est.single.out[kk] <- NA;
                beta1.sd.single.out[kk] <- NA;
                maf.single.out[kk] <- NA;
                anno.single.out[kk] <- paste(raw.data$anno[ix.candidate],collapse=",");
                p.value.known.single.out[kk] <- NA;
                pos.ref.alt.known.single.out[kk] <- paste(raw.data$pos,raw.data$ref[[ix.gold]],raw.data$alt[[ix.gold]],sep="/",collapse=",");
              }
          }
        if(length(ix.var)>length(ix.tmp)) {
          
          ix.pop <- 1:length(raw.data$nSample);

          score.stat.vec.list <- list();mac.vec.list <- list();maf.vec.list <- list();cov.mat.list <- list();var.Y.list <- list();N.list <- list();mean.Y.list <- list();pos.list <- list();anno.list <- list();ac.vec.list <- list();af.vec.list <- list();
          ref.list <- list();
          alt.list <- list();
          for(ii in 1:length(ix.pop))
              {
                  if(length(raw.data$covXZ[[ii]])>0) {
                      warning(paste0("Study  ",ii," is analyzed as binary trait. It is advised to use rareMETALS2 for meta-analysis"))
                  }
                  
                  N.list[[ii]] <- rm.na(as.integer(mean(raw.data$nSample[[ii]],na.rm=TRUE)));
                  
                  U.stat <- rm.na(raw.data$ustat[[ii]][ix.var]);
              V.stat <- rm.na(raw.data$vstat[[ii]][ix.var]);
              score.stat.vec.list[[ii]] <- rm.na(U.stat/V.stat);
              cov.mat.list[[ii]] <- as.matrix(rm.na(as.matrix(raw.data$cov[[ii]])[ix.var,ix.var]));
              var.Y.list[[ii]] <- 1;
              mean.Y.list[[ii]] <- 0;
              af.vec.list[[ii]] <- rm.na((raw.data$af[[ii]])[ix.var]);
              ac.vec.list[[ii]] <- rm.na((raw.data$ac[[ii]])[ix.var]);
              pos.list[[ii]] <- (raw.data$pos)[ix.var];
              ref.list[[ii]] <- (raw.data$ref)[[ii]][ix.var];
              alt.list[[ii]] <- (raw.data$alt)[[ii]][ix.var];
              anno.list[[ii]] <- (raw.data$anno)[ix.var];
            }          
          if(length(ix.pop)>1)
            {
              for(ii in 1:length(ix.var))
                {
                  for(jj in (1:length(ix.pop))[-ix.gold])
                    {
                      if(is.na(ref.list[[ix.gold]][ii]) & !is.na(ref.list[[jj]][ii]))
                        {
                          ref.list[[ix.gold]][ii] <- ref.list[[jj]][ii];
                          alt.list[[ix.gold]][ii] <- alt.list[[jj]][ii];
                        }
                      if((is.na(alt.list[[ix.gold]][ii]) | (alt.list[[ix.gold]][ii])=='0') & !is.na(alt.list[[jj]][ii]))
                        {
                          alt.list[[ix.gold]][ii] <- alt.list[[jj]][ii];
                        }
                      if(!is.na(ref.list[[jj]][ii]) & !is.na(ref.list[[ix.gold]][ii]))
                        {
                          if(ref.list[[jj]][ii]==alt.list[[ix.gold]][ii] & (ref.list[[jj]][ii])!=(ref.list[[ix.gold]][ii]))
                            {
                              tmp <- ref.list[[jj]][ii];
                              ref.list[[jj]][ii] <- alt.list[[jj]][ii];
                              alt.list[[jj]][ii] <- tmp;                             
                              score.stat.vec.list[[jj]][ii] <- (-1)*(score.stat.vec.list[[jj]][ii]);
                              cov.mat.list[[jj]][ii,] <- (-1)*(cov.mat.list[[jj]][ii,]);
                              cov.mat.list[[jj]][,ii] <- (-1)*(cov.mat.list[[jj]][,ii]);
                              af.vec.list[[jj]][ii] <- 1-af.vec.list[[jj]][ii];
                              ac.vec.list[[jj]][ii] <- 2*N.list[[jj]]-ac.vec.list[[jj]][ii];
                            }
                        }
                    }
                }
            }
          maf.vec <- rep(0,length(af.vec.list[[1]]));
          af.vec <- maf.vec;
          mac.vec <- 0;ac.vec <- 0;no.sample <- 0;
          for(ii in 1:length(ix.pop))
            {
              ix.1 <- which(af.vec.list[[ii]]==1);
              if(length(ix.1)>0) {
                af.vec.list[[ii]][ix.1] <- 0;
                ac.vec.list[[ii]][ix.1] <- 0;
              }
              af.vec <- af.vec+(af.vec.list[[ii]])*(2*N.list[[ii]]);
              ac.vec <- ac.vec+ac.vec.list[[ii]];
              no.sample <- no.sample+N.list[[ii]];
            }
          af.vec <- af.vec/sum(2*unlist(N.list));
          maf.vec <- af.vec;
          mac.vec <- ac.vec;
          ix.major <- which(af.vec>0.5);
          if(length(ix.major)>0)
            {
              maf.vec[ix.major] <- 1-maf.vec[ix.major];
              mac.vec[ix.major] <- no.sample-ac.vec[ix.major];
            }
          ix.rare <- which(maf.vec<maf.cutoff);
          ix.tmp <- match(known.variant.vec,pos.list[[ix.gold]]);
          ix.rare <- unique(c(ix.rare,ix.tmp));
          ix.tmp <- match(known.variant.vec,pos.list[[ix.gold]][ix.rare]);
          ix.rare <- c(ix.rare[-ix.tmp],ix.rare[ix.tmp]);
          maf.vec.rare <- maf.vec[ix.rare];
          mac.vec.rare <- mac.vec[ix.rare];
          if(length(ix.rare)==length(ix.tmp))
            {
              res.null$gene.name <- gene.name[kk];
              res[[kk]] <- res.null;
            }
          if(length(ix.rare)>length(ix.tmp))
            {
              for(ii in 1:length(ix.pop))
                {
                  score.stat.vec.list[[ii]] <- score.stat.vec.list[[ii]][ix.rare];
                  cov.mat.list[[ii]] <- as.matrix(cov.mat.list[[ii]][ix.rare,ix.rare]);
                  af.vec.list[[ii]] <- af.vec.list[[ii]][ix.rare];
                  ac.vec.list[[ii]] <- ac.vec.list[[ii]][ix.rare];
                  anno.list[[ii]] <- anno.list[[ix.gold]][ix.rare];
                  pos.list[[ii]] <- pos.list[[ix.gold]][ix.rare];
                  ref.list[[ii]] <- ref.list[[ix.gold]][ix.rare];
                  alt.list[[ii]] <- alt.list[[ix.gold]][ix.rare];
                }
              res.extra <- list(anno=anno.list[[ix.gold]],
                                pos=pos.list[[ix.gold]],
                                ref=ref.list[[ix.gold]],
                                alt=alt.list[[ix.gold]],
                                af.vec.list=af.vec.list,
                                score.stat.vec.list=score.stat.vec.list,
                                cov.mat.list=cov.mat.list,
                                pos.list=pos.list,
                                alt.list=alt.list,
                                ref.list=ref.list,
                                raw.data=raw.data,
                                gene.name=gene.name[kk]);
              ix.X1 <- 1:(length(ix.rare)-length(ix.tmp));
              ix.X2 <- (length(ix.rare)-length(ix.tmp)+1):length(ix.rare);
              
              
              res.kk <- cond.rvmeta(score.stat.vec.list,af.vec.list,cov.mat.list,var.Y.list,N.list,alternative,no.boot,alpha,rv.test='SINGLE',extra.pars=list(ix.X1=ix.X1,ix.X2=ix.X2,maf.vec=maf.vec.rare,mac.vec=mac.vec.rare),knownCoding);
              
              res[[kk]] <- c(res.kk,res.extra);
              res[[kk]]$maf.cutoff <- maf.cutoff;
              
              gene.name.out[kk] <- res[[kk]]$gene.name;
              p.value.out[kk] <- format(res[[kk]]$p.value,digits=out.digits);
              statistic.out[kk] <- format(res[[kk]]$statistic,digits=out.digits);
              no.site.out[kk] <- res[[kk]]$no.site;
              beta1.est.out[kk] <- format(res[[kk]]$beta1.est,digits=out.digits);
              beta1.sd.out[kk] <- format(res[[kk]]$beta1.sd,digits=out.digits);
              maf.cutoff.out[kk] <- format(res[[kk]]$maf.cutoff,digits=out.digits);
              direction.burden.by.study.out[kk] <- res[[kk]]$direction.burden.by.study;
              direction.meta.single.var.out[kk] <- res[[kk]]$direction.meta.single.var;          
              pos.ref.alt.out[kk] <- paste(res[[kk]]$pos[ix.X1],res[[kk]]$ref[ix.X1],res[[kk]]$alt[ix.X1],res[[kk]]$anno[ix.X1],sep='/',collapse=',');
              pos.ref.alt.known.out[kk] <- paste(res[[kk]]$pos[ix.X2],res[[kk]]$ref[ix.X2],res[[kk]]$alt[ix.X2],res[[kk]]$anno[ix.X2],sep='/',collapse=',');
              p.value.known.out[kk] <- p.value.known.variant.vec[kk];
              anno.known.out[kk] <- anno.known.variant.vec[kk];
              
              
              p.value.single.out[kk] <- c((format(res[[kk]]$p.value.single,digit=out.digits)));
              pos.single.out[kk] <- c(paste(res[[kk]]$pos[ix.X1],sep=''));
              direction.single.out[kk] <- c(res[[kk]]$direction.single.vec);
              ref.single.out[kk] <- c(res[[kk]]$ref[ix.X1]);
              alt.single.out[kk] <- c(res[[kk]]$alt[ix.X1]);
              beta1.est.single.out[kk] <- c(format(res[[kk]]$beta1.est.single[ix.X1],digits=out.digits));
              beta1.sd.single.out[kk] <- c(format(res[[kk]]$beta1.sd.single[ix.X1],digits=out.digits));
              maf.single.out[kk] <- c(format(res[[kk]]$maf.vec[ix.X1],digits=out.digits));
              anno.single.out[kk] <- anno.candidate.variant;
              p.value.known.single.out[kk] <- c(rep(paste(format(p.value.known.variant.vec[kk],digits=out.digits),sep=',',collapse=','),length(ix.X1)));
              anno.known.single.out[kk] <- c(rep(paste(anno.known.variant.vec[kk],sep=',',collapse=','),length(ix.X1)));
              pos.ref.alt.known.single.out[kk] <- rep(format(pos.ref.alt.known.out[kk],digits=out.digits),length(ix.X1));
              
          }
      }
    }

    res.single.out <- cbind(pos.single.out,ref.single.out,alt.single.out,p.value.single.out,maf.single.out,beta1.est.single.out,beta1.sd.single.out,direction.single.out,anno.single.out,pos.ref.alt.known.single.out);
    colnames(res.single.out) <- c("POS","REF","ALT","PVALUE","AF","BETA_EST","BETA_SD","DIRECTION_BY_STUDY","ANNO","POS_REF_ALT_ANNO_KNOWN");
    
    return(list(res.list=res,
                res.out=res.single.out));
    
    
}
r2cov.mat <- function(r2.mat,maf.vec)
    {
        var.vec <- sqrt(maf.vec*(1-maf.vec)*2);
        var.mat <- (var.vec%*%t(var.vec));
        cov.mat <- r2.mat*var.mat;
        return(cov.mat);
    }
