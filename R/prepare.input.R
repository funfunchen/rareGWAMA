#combines ped files into a single list for mega analyses
prepare.mega.input <- function(dat.ped.list,sumstat.list,ix.pop)
  {
    maf.vec.list <- list();
    covar.mat.list <- list();
    if(!is.na(ix.pop))
      {
        dat.ped <- dat.ped.list[[ix.pop]];
        maf.vec <- (sumstat.list[[ix.pop]])$maf.vec;
        dat.ped.list <- list(dat.ped);
        maf.vec.list <- list(maf.vec);
        covar.mat.list <- list(integer(0));
      }
    if(is.na(ix.pop))
      {
        for(ii in 1:length(dat.ped.list))
          {          
            maf.vec.list[[ii]] <- sumstat.list[[ii]]$maf.vec;
            covar.mat.list[[ii]] <- integer(0);
          }
      }
    return(list(dat.ped.list=dat.ped.list,
                maf.vec.list=maf.vec.list,
                covar.mat.list=covar.mat.list));
  }

prepare.pool.input <- function(dat.ped.list,sumstat.list,ix.pop)
  {
    maf.vec.list <- list();
    covar.mat.list <- list();
    if(!is.na(ix.pop))
      {
        dat.ped <- dat.ped.list[[ix.pop]];
        maf.vec <- (sumstat.list[[ix.pop]])$maf.vec;
        dat.ped.list <- list(dat.ped);
        maf.vec.list <- list(maf.vec);
        covar.mat.list <- list(integer(0));
      }
    if(is.na(ix.pop))
      {
        dat.ped.new <- list();        
        genotype.new <- dat.ped.list[[1]]$genotype
        N <- nrow(genotype.new);
        maf.vec.new <- (sumstat.list[[1]]$maf.vec)*N;
        phenotype.new <- dat.ped.list[[1]]$phenotype;
        phenotype.ori.new <- dat.ped.list[[1]]$phenotype.ori;
        if(length(dat.ped.list)>1)
          {
            for(ii in 2:length(dat.ped.list))
              {
                genotype.new <- rbind(genotype.new,dat.ped.list[[ii]]$genotype);
                phenotype.new <- c(phenotype.new,dat.ped.list[[ii]]$phenotype);
                phenotype.ori.new <- c(phenotype.ori.new,dat.ped.list[[ii]]$phenotype.ori);
                N <- nrow(dat.ped.list[[ii]]$genotype);
                maf.vec.new <- maf.vec.new+(sumstat.list[[ii]]$maf.vec)*N;
                
              }
          }
        maf.vec.new <- maf.vec.new/nrow(genotype.new);
        dat.ped.new <- list(genotype=genotype.new,
                            phenotype.ori=phenotype.ori.new,
                            phenotype=phenotype.new);
        
        dat.ped.list <- list(dat.ped.new);
        maf.vec.list <- list(maf.vec.new);
        covar.mat.list <- list(integer(0));
      }
    return(list(dat.ped.list=dat.ped.list,
                maf.vec.list=maf.vec.list,
                covar.mat.list=covar.mat.list));
  }

prepare.meta.input <- function(sumstat.list,ix.pop)
  {
    
    score.stat.vec.list <- list();maf.vec.list <- list();cov.mat.list <- list();var.Y.list <- list();N.list <- list();
    mean.Y.list <- list();
    if(is.na(ix.pop)) ix.pop <- 1:length(sumstat.list);
    ix.na <- list();
    for(ii in 1:length(ix.pop))
      {
        score.stat.vec.list[[ii]] <- sqrt((sumstat.list[[ix.pop[ii] ]])$score.stat.sq.vec)*((sumstat.list[[ix.pop[ii] ]])$direction.vec);
        maf.vec.list[[ii]] <- ((sumstat.list[[ix.pop[ii] ]])$maf.vec);
        cov.mat.list[[ii]] <- ((sumstat.list[[ix.pop[ii] ]])$cov.mat);
        var.Y.list[[ii]] <- ((sumstat.list[[ix.pop[ii] ]])$var.Y);
        mean.Y.list[[ii]] <- ((sumstat.list[[ix.pop[ii] ]])$mean.Y);
        N.list[[ii]] <- ((sumstat.list[[ix.pop[ii] ]])$N);
      }
    return(list(score.stat.vec.list=score.stat.vec.list,
                maf.vec.list=maf.vec.list,
                cov.mat.list=cov.mat.list,
                var.Y.list=var.Y.list,
                mean.Y.list=mean.Y.list,
                N.list=N.list))
  }
#' Remove nas
#'
#' @param x input;
#' @export
rm.na <- function(x) {
    
    ix.0 <- which(as.logical(x==Inf) | as.logical(x==-Inf) | as.logical(is.na(x)));
    if(length(ix.0)>0) x[ix.0] <- 0;
    return(x);
}


