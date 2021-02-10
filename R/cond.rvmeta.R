cond.rvmeta <- function(score.stat.vec.list,maf.vec.list,cov.mat.list,var.Y.list,N.list,alternative=c('two.sided','greater','less'),no.boot,alpha=0.05,rv.test,extra.pars=list(),knownCoding='identity')
  {
    if(length(alternative)>1) alternative <- "two.sided";
    res.list <- list();
    ix.X1 <- extra.pars$ix.X1;
    ix.X2 <- extra.pars$ix.X2
    X.T.times.Y.centered.list <- list();
    X.T.times.X.list <- list();
    X.T.times.Y.centered.uncond.list <- list();
    X.T.times.X.uncond.list <- list();
    res <- list();

    af.vec <- rep(0,length(ix.X1));
    maf.vec <- af.vec;
    ac.vec <- af.vec;
    ac.vec.list <- extra.pars$mac.vec.list;
    X.T.times.Y <- rep(0,length(ix.X1));
    X.T.times.X <- matrix(0,nrow=length(ix.X1),ncol=length(ix.X1));
    direction.burden.by.study.vec <- rep('',length(score.stat.vec.list));
    direction.meta.single.var.vec <- rep('',length(ix.X1));
    af.vec.cond.list <- list();
    ac.vec.cond.list <- list();
    burden.stat.by.study <- 0;
    direction.single.mat <- matrix(nrow=length(score.stat.vec.list),ncol=length(ix.X1))
    
    direction.code <- c("+","=","-");
    for(ii in 1:length(score.stat.vec.list))
      {
        score.stat.vec.list[[ii]] <- rm.na(score.stat.vec.list[[ii]]);

        
        maf.vec.list[[ii]] <- rm.na(maf.vec.list[[ii]]);
        ac.vec.list[[ii]] <- rm.na(ac.vec.list[[ii]]);
        af.vec.cond.list[[ii]] <- maf.vec.list[[ii]][ix.X1];
        ac.vec.cond.list[[ii]] <- ac.vec.list[[ii]][ix.X1];
        cov.mat.list[[ii]] <- as.matrix(rm.na(cov.mat.list[[ii]]));
        var.Y.list[[ii]] <- 1;
        X.T.times.X.uncond.list[[ii]] <- N.list[[ii]]*(cov.mat.list[[ii]])*(var.Y.list[[ii]]);
        
        X.T.times.Y.centered.uncond.list[[ii]] <- (sqrt(N.list[[ii]]))*(score.stat.vec.list[[ii]])*sqrt(diag(cov.mat.list[[ii]]))*sqrt(var.Y.list[[ii]]);
        maf.vec <- maf.vec.list[[ii]][ix.X1]*(N.list[[ii]])+maf.vec
        ac.vec <- ac.vec.list[[ii]][ix.X1]+ac.vec;
        if(knownCoding=="identity") {
            res.cond.ii <- cond.rvmeta.core(X.T.times.Y.centered.uncond.list[[ii]],X.T.times.X.uncond.list[[ii]],extra.pars$maf.vec,N.list[[ii]],var.Y.list[[ii]],ix.X1,ix.X2,"generic",alternative,no.boot,list());
        }
        if(knownCoding=="burden") {
            X.T.times.Y.ii <- c(X.T.times.Y.centered.uncond.list[[ii]][ix.X1],sum(X.T.times.Y.centered.uncond.list[[ii]][ix.X2]));
            ix.candidate <- ix.X1;ix.known <- ix.X2;
            X.T.times.X.ii <- rbind(cbind(matrix(X.T.times.X.uncond.list[[ii]][ix.X1,ix.X1],nrow=length(ix.candidate),ncol=length(ix.candidate)),
                                          matrix(rowSums(matrix(X.T.times.X.uncond.list[[ii]][ix.candidate,ix.known],nrow=length(ix.candidate),ncol=length(ix.known))),nrow=length(ix.candidate),ncol=1)),
                                    cbind(matrix(rowSums(matrix(X.T.times.X.uncond.list[[ii]][ix.candidate,ix.known],nrow=length(ix.candidate),ncol=length(ix.known))),ncol=length(ix.candidate),nrow=1),
                                      matrix(sum(X.T.times.X.uncond.list[[ii]][ix.known,ix.known]),nrow=1,ncol=1)));            
            res.cond.ii <- cond.rvmeta.core(X.T.times.Y.ii,X.T.times.X.ii,extra.pars$maf.vec,N.list[[ii]],var.Y.list[[ii]],ix.X1,length(ix.X1)+1,"generic",alternative,no.boot,list());
        }
        
        X.T.times.Y.centered.list[[ii]] <- rm.na(as.vector(res.cond.ii$X.T.times.Y));
        direction.single.mat[ii,] <- direction.code[sign(X.T.times.Y.centered.list[[ii]])+2];
        X.T.times.X.list[[ii]] <- rm.na(res.cond.ii$X.T.times.X);
        U.ii <- sum(X.T.times.Y.centered.list[[ii]]);
        burden.stat.by.study[ii] <- U.ii;
        cov.mat.list[[ii]] <- rm.na(X.T.times.X.list[[ii]]/N.list[[ii]]);
        if(U.ii>0) direction.burden.by.study.vec[ii] <- "+";
        if(U.ii<0) direction.burden.by.study.vec[ii] <- "-";
        if(U.ii==0) direction.burden.by.study.vec[ii] <- "?";                
        X.T.times.Y <- X.T.times.Y+X.T.times.Y.centered.list[[ii]];
        X.T.times.X <- X.T.times.X+X.T.times.X.list[[ii]];
      }
    direction.single.vec <- apply(direction.single.mat,2,paste,sep='',collapse='');
    maf.vec.cond <- (extra.pars$maf.vec)[ix.X1];
    ac.vec.cond <- ac.vec[ix.X1];

    for(ii in 1:length(ix.X1))
      {
        U.ii <- X.T.times.Y[ii];
        if(U.ii>0) direction.meta.single.var.vec[ii] <- "+";
        if(U.ii<0) direction.meta.single.var.vec[ii] <- "-";
        if(U.ii==0) direction.meta.single.var.vec[ii] <- "?"; 
      }
    direction.meta.single.var <- paste(direction.meta.single.var.vec,sep='',collapse='');
    direction.burden.by.study <- paste(direction.burden.by.study.vec,sep='',collapse='');
    N <- sum(unlist(N.list));
    maf.vec <- maf.vec/N;
    singlevar.af.vec <- maf.vec;
    cov.mat <- X.T.times.X/N;
    ustat.vec <- X.T.times.Y;
    X.T.times.X <- matrix(X.T.times.X,nrow=length(X.T.times.Y),ncol=length(X.T.times.Y));
    vstat.vec <- sqrt(diag(X.T.times.X));
    singlevar.stat.vec <- rm.na((ustat.vec/vstat.vec)^2);
    singlevar.pval.vec <- pchisq(singlevar.stat.vec,df=1,lower.tail=FALSE);
    ix.best <- which.min(singlevar.pval.vec);
    if(length(X.T.times.Y)==1)
      {
        var.X <- cov.mat;
        beta1.est <- ustat.vec/vstat.vec^2;
        beta1.sd <- 1/vstat.vec^2
        hsq.est <- beta1.est*beta1.est*var.X;
        statistic <- singlevar.stat.vec;
        p.value <- singlevar.pval.vec;
        return(list(statistic=statistic,
                    p.value=p.value,
                    p.value.single=p.value,
                    no.site=1,
                    beta1.est=beta1.est,
                    beta1.sd=beta1.sd,
                    beta1.est.single=beta1.est,
                    beta1.sd.single=beta1.sd,
                    hsq.est=hsq.est,
                    maf.vec=maf.vec,
                    ix.best=ix.best,
                    singlevar.af.vec=singlevar.af.vec,
                    singlevar.stat.vec=singlevar.stat.vec,
                    singlevar.pval.vec=singlevar.pval.vec,
                    direction.single.vec=direction.single.vec,
                    direction.meta.single.var=direction.meta.single.var,
                    direction.burden.by.study=direction.burden.by.study));
      }
    
    U.stat <- X.T.times.Y;
    V.stat.sq <- as.numeric(diag(X.T.times.X));
    if(alternative=='two.sided')
      {
        statistic.single <- U.stat^2/V.stat.sq;
        p.value.single <- pchisq(statistic.single,df=1,lower.tail=FALSE);
      }
    if(alternative=='greater')
      {
        statistic.single <- U.stat/sqrt(V.stat.sq);
        p.value.single <- pnorm(statistic.single,lower.tail=FALSE);
      }
    if(alternative=='less')
      {
        statistic.single <- U.stat/sqrt(V.stat.sq);
        p.value.single <- pnorm(statistic,lower.tail=TRUE);
      }
    beta1.est.single <- U.stat/V.stat.sq;
    beta1.sd.single <- sqrt(1/V.stat.sq);
    
    if(rv.test=='WSS')
      {
        weight <- extra.pars$weight;
        if(length(weight)!=1) {weight <- 'MZ';extra.pars$weight <- weight}
        res <- rvmeta.CMH.wss(X.T.times.Y.centered.list,X.T.times.X.list,maf.vec.cond,cov.mat.list,var.Y.list,N.list,alternative,no.boot,alpha,weight);
      }

    if(rv.test=='VT')
      {
        res <- rvmeta.CMH.vt(X.T.times.Y.centered.list,X.T.times.X.list,ac.vec.cond,cov.mat.list,var.Y.list,N.list,alternative,no.boot,alpha,extra.pars=list(max.TH=extra.pars$max.TH));
      }

    if(rv.test=='SKAT')
      {
        kernel <- extra.pars$kernel;
        if(length(kernel)!=1) kernel <- "beta";
        res <- rvmeta.CMH.skat(X.T.times.Y.centered.list,X.T.times.X.list,maf.vec.cond,cov.mat.list,var.Y.list,N.list,alternative,no.boot,alpha,kernel);
      }
    ix.var <- 1:length(maf.vec);

    if(rv.test=='VT')
      {
        ix.var <- which(ac.vec<=res$mac.cutoff);
      }
    ixVar.VT <- NULL;
    if(rv.test=="VT") ixVar.VT <- ix.var;
    X.T.times.X <- X.T.times.X[ix.var,ix.var];
    w <- rep(1,length(X.T.times.Y));
    if(rv.test=='WSS')
      {
        if(extra.pars$weight=='MB')
          {
            q <- ((maf.vec)+1/(2*N))*(2*N)/(2*N+2);
            w <- 1/sqrt(N*q*(1-q));
          }
      }
    beta1.est <- sum((w[ix.var])*(X.T.times.Y[ix.var]))/as.numeric(t(w[ix.var])%*%X.T.times.X%*%(w[ix.var]));
    beta1.sd <- sqrt(1/as.numeric(t(w[ix.var])%*%X.T.times.X%*%(w[ix.var])));
    macf.vec <- 2*(maf.vec[ix.var])*(1-maf.vec[ix.var])+(maf.vec[ix.var])*(maf.vec[ix.var]);
    hsq.est <- beta1.est*beta1.est*as.numeric(t(w[ix.var])%*%X.T.times.X%*%(w[ix.var]));
    beta1.conf.lower <- beta1.est-1.96*beta1.sd;
    beta1.conf.upper <- beta1.est+1.96*beta1.sd;
    return(c(res,list(direction.meta.single.var=direction.meta.single.var,
                      direction.burden.by.study=direction.burden.by.study,
                      direction.single.vec=direction.single.vec,
                      beta1.est=beta1.est,
                      beta1.sd=beta1.sd,
                      beta1.est.single=beta1.est.single,
                      beta1.sd.single=beta1.sd.single,
                      no.site=length(ix.var),
                      hsq.est=hsq.est,
                      maf.vec=maf.vec,
                      X.T.times.Y=X.T.times.Y,
                      p.value.single=p.value.single,
                      statistic.single=statistic.single,
                      singlevar.af.vec=singlevar.af.vec,
                      ix.best=ix.best,
                      ixVar.VT=ix.var,
                      singlevar.stat.vec=singlevar.stat.vec,
                      singlevar.pval.vec=singlevar.pval.vec,
                      burden.stat.by.study=burden.stat.by.study,
                      beta1.conf.lower=beta1.conf.lower,
                      beta1.conf.upper=beta1.conf.upper)));
  }
