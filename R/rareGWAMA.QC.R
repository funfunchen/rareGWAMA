



#' #' QC
#' @export
QC <- function(raw.data,QC.par,cov=1)
  {
    
      hwe.cutoff <- QC.par$hwe.cutoff;
      callrate.cutoff <- QC.par$callrate.cutoff;
      if(length(hwe.cutoff)==0) hwe.cutoff <- 0;
      if(length(callrate.cutoff)==0) callrate.cutoff <- 0;
      log.mat <- matrix("Untyped",nrow=length(raw.data$ref[[1]]),ncol=length(raw.data$ref));
      for(ii in 1:length(raw.data$hwe))
          {
              hwe.ii <- raw.data$hwe[[ii]];
              callrate.ii <- raw.data$callrate[[ii]];
              af.ii <- raw.data$af[[ii]];
              ix.rm <- which((callrate.ii<callrate.cutoff | hwe.ii<hwe.cutoff ) & af.ii!=0 & af.ii!=1);
              ix.hwe <- which(hwe.ii<hwe.cutoff & af.ii!=0 & af.ii!=1);
              ix.callrate <- which(callrate.ii<callrate.cutoff & af.ii!=0 & af.ii!=1);
              if(length(ix.hwe)>0) {
                  log.mat[ix.hwe,ii] <- "HWE";

              }
              if(length(ix.callrate)>0) {
                  log.mat[ix.callrate,ii] <- "CallRate";
              }
              ix.bug <- integer(0);
              if(cov==1)
                  {

                      res.diag <- rm.na(diag(raw.data$cov[[ii]])*raw.data$nSample[[ii]]);
                      res.vstat <- rm.na((raw.data$vstat[[ii]])^2)
                      diag.diff <- abs(res.diag-res.vstat);
                      ix.bug <- which(diag.diff>0.5);
                      msg <- paste(c('study',ii,'different missingness among variants'),sep=" ",collapse=' ');
                      
                      if(length(ix.bug)>0)
                          {
                              log.mat[ix.bug,ii] <- "bug";
                              warning(msg);
                          }
                      
                  }
              ix.rm <- unique(ix.rm);              
              if(length(ix.rm)>0)
                  {
                      raw.data$ustat[[ii]][ix.rm] <- NA;
                      raw.data$vstat[[ii]][ix.rm] <- NA;
                      if(cov==1)
                          {
                              raw.data$cov[[ii]][ix.rm,] <- NA;
                              raw.data$cov[[ii]][,ix.rm] <- NA;
                          }
                      raw.data$ref[[ii]][ix.rm] <- NA;
                      raw.data$alt[[ii]][ix.rm] <- NA;
                      raw.data$nSample[[ii]][ix.rm] <- NA;
                      raw.data$af[[ii]][ix.rm] <- NA;
                      raw.data$ac[[ii]][ix.rm] <- NA;
                      raw.data$nref[[ii]][ix.rm] <- NA;
                      raw.data$nhet[[ii]][ix.rm] <- NA;
                      raw.data$nalt[[ii]][ix.rm] <- NA;
                      raw.data$effect[[ii]][ix.rm] <- NA;
                      raw.data$pVal[[ii]][ix.rm] <- NA; 
                  }
          }
      raw.data$log.mat <- log.mat;
      return(raw.data);
  }

#' flip allele and fix the files
#' @param fnameIn input file name;
#' @param fnameOut output file name;
#' @param refaltList the list of group files;
#' @export
rareGWAMA.flip <- function(fnameIn,fnameOut,refaltList) {
    tabix.all <- paste(1:22,"1-1000000000",sep=":",collapse=",");
        
    raw.data <- rvmeta.readDataByRange(scoreTestFiles=fnameIn , covFiles=NULL,        ranges= tabix.all,         multiAllelic=F);
    raw.data <- raw.data[[1]];
    raw.data.ori <- raw.data;
    ii <- 1;ix.pop <- 1;ix.include <- 0;log.mat.var <- 0;correctFlip=TRUE;analyzeRefAltListOnly=TRUE
    for(ix.var in 1:length(raw.data$af[[1]])) {
        if(as.integer(ix.var/1000)==ix.var/1000) {
        }
        if(is.na(raw.data$af[[ii]][ix.var])) {
            raw.data$af[[ii]][ix.var] <- (raw.data$nalt[[ii]][ix.var]+1/2*raw.data$nhet[[ii]][ix.var])/raw.data$nSample[[ii]][ix.var];
        }
        if((is.na(refaltList$ref[ix.var]) | is.na(refaltList$alt[ix.var])) & analyzeRefAltListOnly ) {
            
            ii <- ix.pop;
            
            ix.include <- rep(0,length(raw.data$ustat));
            
                
            log.mat.var[ii] <- "NotInRefAltList";
            raw.data$nSample[[ii]][ix.var] <- NA;
            raw.data$af[[ii]][ix.var] <- NA;
            raw.data$ac[[ii]][ix.var] <- NA;
            raw.data$ustat[[ii]][ix.var] <- NA;
            raw.data$vstat[[ii]][ix.var] <- NA;
            raw.data$nref[[ii]][ix.var] <- NA;
            raw.data$nhet[[ii]][ix.var] <- NA;
            raw.data$nalt[[ii]][ix.var] <- NA;                
        }
        
        if((is.na(refaltList$ref[ix.var]) | is.na(refaltList$alt[ix.var])) & !analyzeRefAltListOnly ) {
            
            ii <- ix.pop;                
            ix.include <- rep(0,length(raw.data$ustat));
            log.mat.var[ii] <- "NotInRefAltList";
        }
    
    
    
        ii <- ix.pop;
        ref.gold <- refaltList$ref;alt.gold <- refaltList$alt;af.gold <- refaltList$af;af.diff.max <- refaltList$af.diff.max;checkAF <- refaltList$checkAF;
        if(length(checkAF)==0) checkAF <- FALSE;
        
        af.diff <- abs(raw.data$af[[ii]][ix.var]-af.gold[ix.var]);
        ix.include <- rep(0,length(raw.data$ustat));
        if(is.na(af.diff)) af.diff <- 0;
        if(!is.na(raw.data.ori$ref[[ii]][ix.var]) | !is.na(raw.data.ori$alt[[ii]][ix.var])) {
            
        if(rm.na(raw.data$af[[ii]][ix.var])==0 | rm.na(raw.data$af[[ii]][ix.var])==1) {
            
                
            if(rm.na(af.gold[ix.var])<.5) {
                
                ix.include[ii] <- 1;
                raw.data$ustat[[ii]][ix.var] <- 0;
                raw.data$af[[ii]][ix.var] <- 0;
                raw.data$vstat[[ii]][ix.var] <- 0;
                log.mat.var[ii] <- "Monomorphic";
            }
            if(rm.na(af.gold[ix.var])>=.5) {
                
                raw.data$af[[ii]][ix.var] <- 1;
                    
                raw.data$ustat[[ii]][ix.var] <- 0;
                raw.data$vstat[[ii]][ix.var] <- 0;
                nref.tmp <- (raw.data$nalt[[ii]][ix.var]);
                nalt.tmp <- (raw.data$nref[[ii]][ix.var]);
                nhet.tmp <- (raw.data$nhet[[ii]][ix.var]);
                raw.data$nref[[ii]][ix.var] <- nref.tmp;
                raw.data$nalt[[ii]][ix.var] <- nalt.tmp;
                raw.data$nhet[[ii]][ix.var] <- nhet.tmp;
                ix.include[ii] <- 1;
                log.mat.var[ii] <- "Monomorphic";
            }
        }           
    }
        
        if(!is.na(raw.data$ref[[ii]][ix.var]) & !is.na(raw.data$alt[[ii]][ix.var]) & !(rm.na(raw.data$af[[ii]][ix.var])==0 | rm.na(raw.data$af[[ii]][ix.var])==1)) {
            
            strandAmbiguous <- (((ref.gold[ix.var]=="A") & (alt.gold[ix.var]=="T")) | ((ref.gold[ix.var]=="T") & (alt.gold[ix.var]=="A")) | ((ref.gold[ix.var]=="C") & (alt.gold[ix.var]=="G")) | ((ref.gold[ix.var]=="G") & (alt.gold[ix.var]=="C")));
            flip.ref.alt <- (raw.data$ref[[ii]][ix.var]==refaltList$alt[ix.var] & raw.data$alt[[ii]][ix.var]==refaltList$ref[ix.var]);
            match.ref.alt <- (refaltList$ref[ix.var]==(raw.data$ref[[ii]][ix.var]) & refaltList$alt[ix.var]==raw.data$alt[[ii]][ix.var]);
            mono <- (((raw.data$ref[[ii]][ix.var]==".") | (raw.data$ref[[ii]][ix.var]==0) | (raw.data$alt[[ii]][ix.var]==".") | (raw.data$alt[[ii]][ix.var]==0) ) & (rm.na(raw.data$af[[ii]][ix.var])==0 | rm.na(raw.data$af[[ii]][ix.var])==1) & !is.na(raw.data.ori$nSample[[ii]][ix.var]))
            if(!match.ref.alt & !correctFlip &!mono) {
                
                
                ix.include[ii] <- 1;
                    
                log.mat.var[ii] <- "MismatchRemove";
                raw.data$nSample[[ii]][ix.var] <- NA;
                raw.data$af[[ii]][ix.var] <- NA;
                raw.data$ac[[ii]][ix.var] <- NA;
                raw.data$ustat[[ii]][ix.var] <- NA;
                raw.data$vstat[[ii]][ix.var] <- NA;
                raw.data$nref[[ii]][ix.var] <- NA;
                raw.data$nhet[[ii]][ix.var] <- NA;
                raw.data$nalt[[ii]][ix.var] <- NA;
            }
            
            if(match.ref.alt & (!strandAmbiguous)) {
                
                ix.include[ii] <- 1;
                log.mat.var[ii] <- "Match";
            }
            
            if(flip.ref.alt & (!strandAmbiguous)) {
                
                raw.data$af[[ii]][ix.var] <- 1-raw.data$af[[ii]][ix.var];
                raw.data$ustat[[ii]][ix.var] <- (-1)*(raw.data$ustat[[ii]][ix.var]);
                    
                
                nref.tmp <- (raw.data$nalt[[ii]][ix.var]);
                nalt.tmp <- (raw.data$nref[[ii]][ix.var]);
                nhet.tmp <- (raw.data$nhet[[ii]][ix.var]);
                raw.data$nref[[ii]][ix.var] <- nref.tmp;
                raw.data$nalt[[ii]][ix.var] <- nalt.tmp;
                raw.data$nhet[[ii]][ix.var] <- nhet.tmp;
                ix.include[ii] <- 1;
                
                log.mat.var[ii] <- "FlipRefAlt";
                
            }
            af.diff.min <- 0.05;
            if(flip.ref.alt & strandAmbiguous) {
                
                if(af.diff<=af.diff.min) {
                
                    ix.include[ii] <- 1;
                    
                    log.mat.var[ii] <- "FlipStrand";
                    
                    
                }
                if(af.diff>af.diff.max) {
                    
                    raw.data$af[[ii]][ix.var] <- 1-raw.data$af[[ii]][ix.var];
                    raw.data$ustat[[ii]][ix.var] <- (-1)*(raw.data$ustat[[ii]][ix.var]);
                        
                        
                        
                    nref.tmp <- (raw.data$nalt[[ii]][ix.var]);
                    nalt.tmp <- (raw.data$nref[[ii]][ix.var]);
                    nhet.tmp <- (raw.data$nhet[[ii]][ix.var]);
                    raw.data$nref[[ii]][ix.var] <- nref.tmp;
                    raw.data$nalt[[ii]][ix.var] <- nalt.tmp;
                    raw.data$nhet[[ii]][ix.var] <- nhet.tmp;
                    ix.include[ii] <- 1;
                    
                    log.mat.var[ii] <- "FlipRefAlt";
                }
            }
            
            if(match.ref.alt & strandAmbiguous) {
                
                if(af.diff<=af.diff.min) {
                    
                    ix.include[ii] <- 1;
                    
                    log.mat.var[ii] <- "Match";
                }
                if(af.diff>af.diff.max) {
                    
                    log.mat.var[ii] <- "FlipStrand";
                        
                    
                    raw.data$af[[ii]][ix.var] <- 1-raw.data$af[[ii]][ix.var]
                    raw.data$ac[[ii]][ix.var] <- 2*raw.data$nSample[[ii]][ix.var]-raw.data$ac[[ii]][ix.var];
                    raw.data$ustat[[ii]][ix.var] <- raw.data$ustat[[ii]][ix.var]*(-1);
                    tmp <- raw.data$nref[[ii]][ix.var];
                    raw.data$nref[[ii]][ix.var] <- 2*raw.data$nSample[[ii]][ix.var]-raw.data$nref[[ii]][ix.var];
                    raw.data$nalt[[ii]][ix.var] <- tmp;
                }
            }                
            if(mono) {
                
                if(rm.na(af.gold[ix.var])<.5) {
                    
                    ix.include[ii] <- 1;
                    log.mat.var[ii] <- "Monomorphic";
                    raw.data$ustat[[ii]][ix.var] <- 0;
                    raw.data$af[[ii]][ix.var] <- 0;
                }
                if(rm.na(af.gold[ix.var])>=.5) {
                    
                    raw.data$af[[ii]][ix.var] <- 1-raw.data$af[[ii]][ix.var];
                    raw.data$ustat[[ii]][ix.var] <- 0;
                        
                
                    nref.tmp <- (raw.data$nalt[[ii]][ix.var]);
                    nalt.tmp <- (raw.data$nref[[ii]][ix.var]);
                    nhet.tmp <- (raw.data$nhet[[ii]][ix.var]);
                    raw.data$nalt[[ii]][ix.var] <- 2*raw.data$nSample[[ii]][ix.var];
                    raw.data$nref[[ii]][ix.var] <- 0;
                    raw.data$nhet[[ii]][ix.var] <- 0;
                    ix.include[ii] <- 1;
                    log.mat.var[ii] <- "Monomorphic";
                    
                }                        
            }
            if(!match.ref.alt & !flip.ref.alt & !mono) {
                
                ix.include[ii] <- 1;
                
                
                log.mat.var[ii] <- "Unmatched";
                raw.data$nSample[[ii]][ix.var] <- NA;
                raw.data$af[[ii]][ix.var] <- NA;
                raw.data$ac[[ii]][ix.var] <- NA;
                raw.data$ustat[[ii]][ix.var] <- NA;
                raw.data$vstat[[ii]][ix.var] <- NA;
                raw.data$nref[[ii]][ix.var] <- NA;
                raw.data$nhet[[ii]][ix.var] <- NA;
                raw.data$nalt[[ii]][ix.var] <- NA;
            }
            if(checkAF==TRUE) {
                
                af.diff.new <- abs(raw.data$af[[ii]][ix.var]-af.gold[ix.var]);
                if(is.na(af.diff.new)) af.diff.new <- 0;
                if(af.diff.new>af.diff.max) {
                    
                    ix.include[ii] <- 1;
                    log.mat.var[ii] <- "DiffAF";
                    raw.data$nSample[[ii]][ix.var] <- NA;
                    raw.data$af[[ii]][ix.var] <- NA;
                    raw.data$ac[[ii]][ix.var] <- NA;
                    raw.data$ustat[[ii]][ix.var] <- NA;
                    raw.data$vstat[[ii]][ix.var] <- NA;
                    
                    raw.data$nref[[ii]][ix.var] <- NA;
                    raw.data$nhet[[ii]][ix.var] <- NA;
                    raw.data$nalt[[ii]][ix.var] <- NA;                                
                }
            }
            
        }
        
    }
    res.out <- list(raw.data);
    names(res.out) <- 'Range1';
    seqminer::rvmeta.writeScoreData(rvmetaData=res.out,   outName=fnameOut,      createIndex=FALSE);
         
    
}
    


#' This is the function for flipping alleles
#'
#' @param raw.data The input datasets to be considered flipped
#' @param raw.data.ori The input datasets to be considered flipped
#' @param refaltList The list consists of ref, alt, pos, af and af.diff.max, as well as the option of whether throw away sites with large af.differences checkAF;
#' @param ix.pop index of the population
#' @param ix.var index of the variant;
#' @param log.mat.var The log for QC procedure;
#' @param correctFlip Correct for score and covariance matrices for flipped alleles;
#' @return A list consist of modified raw.data, ix.include and log.mat.var
#' @export
flipAllele <- function(raw.data,raw.data.ori,refaltList,ix.pop,ix.var,log.mat.var,correctFlip=TRUE,analyzeRefAltListOnly=TRUE)
    {
        ii <- ix.pop;
        if(is.na(raw.data$af[[ii]][ix.var])) {
            raw.data$af[[ii]][ix.var] <- (raw.data$nalt[[ii]][ix.var]+1/2*raw.data$nhet[[ii]][ix.var])/raw.data$nSample[[ii]][ix.var];
        }
        if((is.na(refaltList$ref[ix.var]) | is.na(refaltList$alt[ix.var])) & analyzeRefAltListOnly )
            {
                ii <- ix.pop;
                
                ix.include <- rep(0,length(raw.data$ustat));
                
                if(length(raw.data$cov)>0)
                    {
                        raw.data$cov[[ii]][ix.var,] <- NA;
                        raw.data$cov[[ii]][,ix.var] <- NA;
                    }
                log.mat.var[ii] <- "NotInRefAltList";
                raw.data$nSample[[ii]][ix.var] <- NA;
                raw.data$af[[ii]][ix.var] <- NA;
                raw.data$ac[[ii]][ix.var] <- NA;
                raw.data$ustat[[ii]][ix.var] <- NA;
                raw.data$vstat[[ii]][ix.var] <- NA;
                raw.data$nref[[ii]][ix.var] <- NA;
                raw.data$nhet[[ii]][ix.var] <- NA;
                raw.data$nalt[[ii]][ix.var] <- NA;                
                return(list(raw.data=raw.data,
                            log.mat.var=log.mat.var,
                            ix.include=ix.include));
            }

        if((is.na(refaltList$ref[ix.var]) | is.na(refaltList$alt[ix.var])) & !analyzeRefAltListOnly )
            {
                ii <- ix.pop;                
                ix.include <- rep(0,length(raw.data$ustat));
                log.mat.var[ii] <- "NotInRefAltList";
                return(list(raw.data=raw.data,
                            log.mat.var=log.mat.var,
                            ix.include=ix.include));
            }


        
        ii <- ix.pop;
        ref.gold <- refaltList$ref;alt.gold <- refaltList$alt;af.gold <- refaltList$af;af.diff.max <- refaltList$af.diff.max;checkAF <- refaltList$checkAF;
        if(length(checkAF)==0) checkAF <- FALSE;
        
        af.diff <- abs(raw.data$af[[ii]][ix.var]-af.gold[ix.var]);
        ix.include <- rep(0,length(raw.data$ustat));
        if(is.na(af.diff)) af.diff <- 0;
        if(!is.na(raw.data.ori$ref[[ii]][ix.var]) | !is.na(raw.data.ori$alt[[ii]][ix.var]))
            {
                if(rm.na(raw.data$af[[ii]][ix.var])==0 | rm.na(raw.data$af[[ii]][ix.var])==1)
                    {

                        if(rm.na(af.gold[ix.var])<.5)
                            {
                                ix.include[ii] <- 1;
                                raw.data$ustat[[ii]][ix.var] <- 0;
                                raw.data$af[[ii]][ix.var] <- 0;
                                raw.data$vstat[[ii]][ix.var] <- 0;
                                log.mat.var[ii] <- "Monomorphic";
                            }
                        if(rm.na(af.gold[ix.var])>=.5)
                            {
                                raw.data$af[[ii]][ix.var] <- 1;
                                if(length(raw.data$cov)>0)
                                    {
                                        raw.data$cov[[ii]][ix.var,] <- (-1)*raw.data$cov[[ii]][ix.var,];
                                        raw.data$cov[[ii]][,ix.var] <- (-1)*raw.data$cov[[ii]][,ix.var]
                                    }
                                raw.data$ustat[[ii]][ix.var] <- 0;
                                raw.data$vstat[[ii]][ix.var] <- 0;
                                nref.tmp <- (raw.data$nalt[[ii]][ix.var]);
                                nalt.tmp <- (raw.data$nref[[ii]][ix.var]);
                                nhet.tmp <- (raw.data$nhet[[ii]][ix.var]);
                                raw.data$nref[[ii]][ix.var] <- nref.tmp;
                                raw.data$nalt[[ii]][ix.var] <- nalt.tmp;
                                raw.data$nhet[[ii]][ix.var] <- nhet.tmp;
                                ix.include[ii] <- 1;
                                log.mat.var[ii] <- "Monomorphic";
                            }
                    }           
            }
        
        if(!is.na(raw.data$ref[[ii]][ix.var]) & !is.na(raw.data$alt[[ii]][ix.var]) & !(rm.na(raw.data$af[[ii]][ix.var])==0 | rm.na(raw.data$af[[ii]][ix.var])==1))
            {
                strandAmbiguous <- (((ref.gold[ix.var]=="A") & (alt.gold[ix.var]=="T")) | ((ref.gold[ix.var]=="T") & (alt.gold[ix.var]=="A")) | ((ref.gold[ix.var]=="C") & (alt.gold[ix.var]=="G")) | ((ref.gold[ix.var]=="G") & (alt.gold[ix.var]=="C")));
                flip.ref.alt <- (raw.data$ref[[ii]][ix.var]==refaltList$alt[ix.var] & raw.data$alt[[ii]][ix.var]==refaltList$ref[ix.var]);
                match.ref.alt <- (refaltList$ref[ix.var]==(raw.data$ref[[ii]][ix.var]) & refaltList$alt[ix.var]==raw.data$alt[[ii]][ix.var]);
                mono <- (((raw.data$ref[[ii]][ix.var]==".") | (raw.data$ref[[ii]][ix.var]==0) | (raw.data$alt[[ii]][ix.var]==".") | (raw.data$alt[[ii]][ix.var]==0) ) & (rm.na(raw.data$af[[ii]][ix.var])==0 | rm.na(raw.data$af[[ii]][ix.var])==1) & !is.na(raw.data.ori$nSample[[ii]][ix.var]))
                if(!match.ref.alt & !correctFlip &!mono)
                    {
                        
                        ix.include[ii] <- 1;
                        if(length(raw.data$cov)>0)
                            {
                                raw.data$cov[[ii]][ix.var,] <- NA;
                                raw.data$cov[[ii]][,ix.var] <- NA;
                            }
                        log.mat.var[ii] <- "MismatchRemove";
                        raw.data$nSample[[ii]][ix.var] <- NA;
                        raw.data$af[[ii]][ix.var] <- NA;
                        raw.data$ac[[ii]][ix.var] <- NA;
                        raw.data$ustat[[ii]][ix.var] <- NA;
                        raw.data$vstat[[ii]][ix.var] <- NA;
                        raw.data$nref[[ii]][ix.var] <- NA;
                        raw.data$nhet[[ii]][ix.var] <- NA;
                        raw.data$nalt[[ii]][ix.var] <- NA;
                    }

                if(match.ref.alt & (!strandAmbiguous))
                    {
                        ix.include[ii] <- 1;
                        log.mat.var[ii] <- "Match";
                    }

                if(flip.ref.alt & (!strandAmbiguous))
                {
                    raw.data$af[[ii]][ix.var] <- 1-raw.data$af[[ii]][ix.var];
                    raw.data$ustat[[ii]][ix.var] <- (-1)*(raw.data$ustat[[ii]][ix.var]);
                    if(length(raw.data$cov)>0)
                    {
                        raw.data$cov[[ii]][ix.var,] <- (-1)*raw.data$cov[[ii]][ix.var,];
                        raw.data$cov[[ii]][,ix.var] <- (-1)*raw.data$cov[[ii]][,ix.var]
                    }
                    
                    nref.tmp <- (raw.data$nalt[[ii]][ix.var]);
                    nalt.tmp <- (raw.data$nref[[ii]][ix.var]);
                    nhet.tmp <- (raw.data$nhet[[ii]][ix.var]);
                    raw.data$nref[[ii]][ix.var] <- nref.tmp;
                    raw.data$nalt[[ii]][ix.var] <- nalt.tmp;
                    raw.data$nhet[[ii]][ix.var] <- nhet.tmp;
                    ix.include[ii] <- 1;
                    
                    log.mat.var[ii] <- "FlipRefAlt";
                    
                }
                af.diff.min <- 0.05;
                if(flip.ref.alt & strandAmbiguous)
                    {                        
                        if(af.diff<=af.diff.min)
                            {
                                ix.include[ii] <- 1;
                                
                                log.mat.var[ii] <- "FlipStrand";
                                
                                
                            }
                        if(af.diff>af.diff.max)
                            {
                                raw.data$af[[ii]][ix.var] <- 1-raw.data$af[[ii]][ix.var];
                                raw.data$ustat[[ii]][ix.var] <- (-1)*(raw.data$ustat[[ii]][ix.var]);
                                if(length(raw.data$cov)>0)
                                    {
                                        
                                        raw.data$cov[[ii]][ix.var,] <- (-1)*raw.data$cov[[ii]][ix.var,];
                                        raw.data$cov[[ii]][,ix.var] <- (-1)*raw.data$cov[[ii]][,ix.var]

                                    }
                                nref.tmp <- (raw.data$nalt[[ii]][ix.var]);
                                nalt.tmp <- (raw.data$nref[[ii]][ix.var]);
                                nhet.tmp <- (raw.data$nhet[[ii]][ix.var]);
                                raw.data$nref[[ii]][ix.var] <- nref.tmp;
                                raw.data$nalt[[ii]][ix.var] <- nalt.tmp;
                                raw.data$nhet[[ii]][ix.var] <- nhet.tmp;
                                ix.include[ii] <- 1;
                                
                                log.mat.var[ii] <- "FlipRefAlt";
                            }
                    }

                if(match.ref.alt & strandAmbiguous)
                    {
                        if(af.diff<=af.diff.min)
                            {
                                ix.include[ii] <- 1;
                                
                                log.mat.var[ii] <- "Match";
                            }
                        if(af.diff>af.diff.max)
                            {
                                log.mat.var[ii] <- "FlipStrand";
                                if(length(raw.data$cov)>0)
                                    {
                                        raw.data$cov[[ii]][ix.var,] <- NA;
                                        raw.data$cov[[ii]][,ix.var] <- NA;
                                    }
                                
                                raw.data$af[[ii]][ix.var] <- 1-raw.data$af[[ii]][ix.var]
                                raw.data$ac[[ii]][ix.var] <- 2*raw.data$nSample[[ii]][ix.var]-raw.data$ac[[ii]][ix.var];
                                raw.data$ustat[[ii]][ix.var] <- raw.data$ustat[[ii]][ix.var]*(-1);
                                tmp <- raw.data$nref[[ii]][ix.var];
                                raw.data$nref[[ii]][ix.var] <- 2*raw.data$nSample[[ii]][ix.var]-raw.data$nref[[ii]][ix.var];
                                raw.data$nalt[[ii]][ix.var] <- tmp;
                            }
                    }                
                if(mono)
                    {
                        if(rm.na(af.gold[ix.var])<.5)
                            {
                                ix.include[ii] <- 1;
                                log.mat.var[ii] <- "Monomorphic";
                                raw.data$ustat[[ii]][ix.var] <- 0;
                                raw.data$af[[ii]][ix.var] <- 0;
                            }
                        if(rm.na(af.gold[ix.var])>=.5)
                            {
                                raw.data$af[[ii]][ix.var] <- 1-raw.data$af[[ii]][ix.var];
                                raw.data$ustat[[ii]][ix.var] <- 0;
                                if(length(raw.data$cov)>0)
                                    {
                                        raw.data$cov[[ii]][ix.var,] <- (-1)*raw.data$cov[[ii]][ix.var,];
                                        raw.data$cov[[ii]][,ix.var] <- (-1)*raw.data$cov[[ii]][,ix.var]
                                    }
                                
                                nref.tmp <- (raw.data$nalt[[ii]][ix.var]);
                                nalt.tmp <- (raw.data$nref[[ii]][ix.var]);
                                nhet.tmp <- (raw.data$nhet[[ii]][ix.var]);
                                raw.data$nalt[[ii]][ix.var] <- 2*raw.data$nSample[[ii]][ix.var];
                                raw.data$nref[[ii]][ix.var] <- 0;
                                raw.data$nhet[[ii]][ix.var] <- 0;
                                ix.include[ii] <- 1;
                                log.mat.var[ii] <- "Monomorphic";

                            }                        
                    }
                if(!match.ref.alt & !flip.ref.alt & !mono)
                    {
                        ix.include[ii] <- 1;
                        if(length(raw.data$cov)>0)
                            {
                                raw.data$cov[[ii]][ix.var,] <- NA;
                                raw.data$cov[[ii]][,ix.var] <- NA;
                            }
                        
                        log.mat.var[ii] <- "Unmatched";
                        raw.data$nSample[[ii]][ix.var] <- NA;
                        raw.data$af[[ii]][ix.var] <- NA;
                        raw.data$ac[[ii]][ix.var] <- NA;
                        raw.data$ustat[[ii]][ix.var] <- NA;
                        raw.data$vstat[[ii]][ix.var] <- NA;
                        raw.data$nref[[ii]][ix.var] <- NA;
                        raw.data$nhet[[ii]][ix.var] <- NA;
                        raw.data$nalt[[ii]][ix.var] <- NA;
                    }
                if(checkAF==TRUE)
                    {
                        af.diff.new <- abs(raw.data$af[[ii]][ix.var]-af.gold[ix.var]);
                        if(is.na(af.diff.new)) af.diff.new <- 0;
                        if(af.diff.new>af.diff.max)
                            {
                                ix.include[ii] <- 1;
                                log.mat.var[ii] <- "DiffAF";
                                raw.data$nSample[[ii]][ix.var] <- NA;
                                raw.data$af[[ii]][ix.var] <- NA;
                                raw.data$ac[[ii]][ix.var] <- NA;
                                raw.data$ustat[[ii]][ix.var] <- NA;
                                raw.data$vstat[[ii]][ix.var] <- NA;
                                
                                raw.data$nref[[ii]][ix.var] <- NA;
                                raw.data$nhet[[ii]][ix.var] <- NA;
                                raw.data$nalt[[ii]][ix.var] <- NA;                                
                            }
                    }
     
            }
        return(list(raw.data=raw.data,
                    log.mat.var=log.mat.var,
                    ix.include=ix.include));
    }

