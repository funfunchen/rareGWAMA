#rvmeta.CMH.vt(X.T.times.Y.centered.list,X.T.times.X.list,maf.vec.list,cov.mat.list,var.Y.list,N.list,alternative,no.boot);
rvmeta.CMH.vt <- function(X.T.times.Y.centered.list,X.T.times.X.list,mac.vec,cov.mat.list,var.Y.list=list(1),N.list,alternative,no.boot,alpha=0.05,extra.pars=list())
  {
    X.T.times.X <- 0*diag(length(mac.vec));#this is actually centered;
    X.T.times.Y <- rep(0,length(mac.vec));
    for(ii in 1:length(X.T.times.X.list))
      {
        X.T.times.Y <- X.T.times.Y+X.T.times.Y.centered.list[[ii]];
        X.T.times.X <- (N.list[[ii]])*(cov.mat.list[[ii]])*(var.Y.list[[ii]])+X.T.times.X;
      }
    
    N <- sum(unlist(N.list));
    maf.vec <- mac.vec/N/2;
    X.T.times.X <- as.matrix(X.T.times.X);
    
    X.T.times.X.NA <- X.T.times.X;
    ix.rm.NA <- which(diag(X.T.times.X.NA)==0);
    if(length(ix.rm.NA)>0)
      {

        X.T.times.X.NA <- X.T.times.X.NA[-ix.rm.NA,-ix.rm.NA];
        X.T.times.X.NA <- X.T.times.X.NA/sqrt(diag( X.T.times.X.NA)%*%t(diag(X.T.times.X.NA)))
      }
    diag(X.T.times.X.NA) <- NA;
    maf.TH <- sort(unique(maf.vec));
    max.TH <- extra.pars$max.TH;
    if(length(max.TH)>0)
      {
        ix.tmp <- as.integer(seq(1,length(maf.TH),length=max.TH));
        maf.TH <- maf.TH[ix.tmp];
        maf.TH <- unique(sort(maf.TH));
      }
    maf.TH.old <- maf.TH;
    if(maf.TH[1]==0) {maf.TH <- maf.TH[-1];}
    if(length(maf.TH)==0)
        {
            return(list(p.value=NA,statistic=NA,maf.cutoff=0,mac.cutoff=0,no.site.VT=length(which(mac.vec==0))));
        }
    if(length(maf.TH)==1)
      {
        ix.ii <- which(maf.vec<=maf.TH);    
        U.stat <- sum(X.T.times.Y[ix.ii]);
        ind.ii <- rep(0,length(maf.vec));
        ind.ii[ix.ii] <- 1;
        V.stat.sq <- as.numeric(t(ind.ii)%*%X.T.times.X%*%(ind.ii));
        if(alternative=='two.sided')
          {
            statistic <- U.stat^2/V.stat.sq;
            p.value <- pchisq(statistic,df=1,lower.tail=FALSE);
          }
        if(alternative=='greater')
          {
            statistic <- U.stat/sqrt(V.stat.sq);
            p.value <- pnorm(statistic,lower.tail=FALSE);
          }
        if(alternative=='less')
          {
            statistic <- U.stat/sqrt(V.stat.sq);
            p.value <- pnorm(statistic,lower.tail=TRUE);
          }
        return(list(p.value=p.value,
                    statistic=statistic,
                    maf.cutoff=maf.TH,
                    mac.cutoff=maf.TH*N*2,
                    no.site.VT=length(ix.ii)));
      }
    err.msg <- vector(length=0);
    ix.list <- list();
    X.T.times.Y.VT <- 0;
    VT.mat <- matrix(0,nrow=length(maf.TH),ncol=ncol(X.T.times.X));
    no.TH <- length(maf.TH);
    for(ii in 1:length(maf.TH))
      {
        ix.ii <- which(maf.vec<=maf.TH[ii]);
        ix.list[[ii]] <- ix.ii;
        VT.mat[ii,ix.ii] <- 1;
        X.T.times.Y.VT[ii] <- sum(X.T.times.Y[ix.ii]);
      }
    cov.X.VT <- matrix(0,nrow=length(ix.list),ncol=length(ix.list));
    cov.mat <- X.T.times.X;
    a=Sys.time();
    for(ii in 1:length(maf.TH))
      {
        for(jj in 1:length(maf.TH))
          {
            ix.ii <- rep(0,nrow(cov.mat));
            ix.jj <- rep(0,nrow(cov.mat));
            ix.ii[ix.list[[ii]]] <- 1;
            ix.jj[ix.list[[jj]]] <- 1;
            cov.X.VT[ii,jj] <- as.numeric(t(ix.ii)%*%cov.mat%*%(ix.jj));
          }
      }

    cor.X.VT <- cov.X.VT/sqrt(diag(cov.X.VT)%*%t(diag(cov.X.VT)));
    if(alternative=='two.sided')
      {
        vt.stat.vec <- (X.T.times.Y.VT^2/(diag(cov.X.VT)));
        ix.rm <- which(is.na(vt.stat.vec));
        if(length(ix.rm)>0) {
          vt.stat.vec <- vt.stat.vec[-ix.rm];
          cor.X.VT <- as.matrix(cor.X.VT[-ix.rm,-ix.rm]);          
          maf.TH <- maf.TH[-ix.rm];
        }
        ix.max <- which.max(vt.stat.vec);
        vt.max.stat <- vt.stat.vec[ix.max];
        maf.cutoff <- maf.TH[ix.max];
        if(length(ix.max)==0) {
            vt.max.stat <- NA;
            maf.cutoff <- maf.TH.old[length(maf.TH.old)];
        }
      }
    if(alternative=='greater')
      {
          vt.stat.vec <-(X.T.times.Y.VT/sqrt(diag(cov.X.VT)));
          ix.rm <- which(is.na(vt.stat.vec));
          if(length(ix.rm)>0) {
              vt.stat.vec <- vt.stat.vec[-ix.rm];
              cor.X.VT <- as.matrix(cor.X.VT[-ix.rm,-ix.rm]);
              maf.TH <- maf.TH[-ix.rm];
          }
          ix.max <- which.max(vt.stat.vec);
          vt.max.stat <- vt.stat.vec[ix.max];
          maf.cutoff <- maf.TH[ix.max];
      }
    if(alternative=='less')
      {
        vt.stat.vec <- (X.T.times.Y.VT/sqrt(diag(cov.X.VT)));
        ix.rm <- which(is.na(vt.stat.vec));
        if(length(ix.rm)>0) {
          vt.stat.vec <- vt.stat.vec[-ix.rm];
          cor.X.VT <- as.matrix(cor.X.VT[-ix.rm,-ix.rm]);
          maf.TH <- maf.TH[-ix.rm];
        }
        ix.max <- which.min(vt.stat.vec);
        vt.max.stat <- vt.stat.vec[ix.max];
        maf.cutoff <- maf.TH[ix.max];
      }
    p.value <- NA;p.value.analytic <- NA;
    
    if(length(vt.max.stat)>0) {
      if(1)
        {
          a <- Sys.time();
          if(!is.na(vt.max.stat)){
              
            p.value.analytic <- pvt(vt.max.stat,mu=rep(0,nrow(cor.X.VT)),sigma=cor.X.VT,alternative);
        }
        }
      if(no.boot>0)
        {
          res.svd <- svd(X.T.times.X)
          U <- res.svd$u;
          V <- res.svd$v;
          D <- res.svd$d;
          res.C <- .C("genericBoot",
                      U=as.double(U),
                      V=as.double(V),
                      D=as.double(D),
                      noMarkerPt=as.integer(nrow(cov.mat)),
                      X_T_times_XPt=as.double(X.T.times.X),
                      NPt=as.integer(N),
                      noBootPt=as.integer(no.boot),
                      varYPt=as.double(1),
                      statisticPt=as.double(rep(0,1000)),
                      alternativePt=as.integer(0),
                      testPt=as.integer(3),
                      extraPar=as.double(c(vt.max.stat,alpha,0.0,no.TH,as.vector(VT.mat))));
          p.value <- res.C$extraPar[3];
          statistic <- res.C$statisticPt;
        }
  }
    if(length(vt.max.stat)==0) {
        
        vt.max.stat <- NA;
        p.value.analytic <- NA;
        maf.cutoff <- maf.TH.old[length(maf.TH.old)];
    }
    p.value.boot <- p.value;
    p.value <- p.value.analytic;
    
    return(list(statistic=vt.max.stat,
                p.value=p.value,
                p.value.boot=p.value.boot,
                p.value.analytic=p.value.analytic,
                err.msg=err.msg,
                X.T.times.X=X.T.times.X,
                X.T.times.Y=X.T.times.Y,
                X.T.times.Y.VT=X.T.times.Y.VT,
                no.TH=nrow(cor.X.VT),
                maf.TH=maf.TH,
                mac.TH=maf.TH*N*2,
                cor.X.VT=cor.X.VT,
                maf.cutoff=maf.cutoff,
                mac.cutoff=maf.cutoff*N*2,
                no.site.VT=length(which(maf.vec<=maf.cutoff)),
                cor.X.VT=cor.X.VT));
}
pmaxnorm <- function(x,mu,sigma)
  {
    tmp <- try(as.numeric(pmvnorm(upper=rep(x,length(mu)),lower=rep(-Inf,length(mu)),mu,sigma)),silent=TRUE);
    if(length(attr(tmp,'class'))>0) return(NA);
    return(tmp);
  }
pminnorm <- function(x,mu,sigma)
  {
    tmp <- try(as.numeric(pmvnorm(upper=rep(Inf,length(mu)),lower=rep(x,length(mu)),mean=mu,sigma=sigma)),silent=TRUE);
    if(length(attr(tmp,'class'))>0) return(NA);
    return(1-tmp);
  }

prune.cor <- function(cor.mat,r2.cutoff=1-1e-5)
  {
    ix.rm <- integer(0);
    for(ii in 1:(nrow(cor.mat)-1))
      {
        for(jj in (ii+1):ncol(cor.mat))
          {
            if(!is.na(cor.mat[ii,jj]))
              {
                if(cor.mat[ii,jj]>r2.cutoff) {
                  if(length(which(ix.rm==ii))==0 & length(which(ix.rm==jj))==0)
                    ix.rm <- c(ix.rm,ii);
                }
              }
          }
      }
    ix.rm <- unique(ix.rm);
    ix.keep <- 1:nrow(cor.mat);
    if(length(ix.rm)>0) ix.keep <- ix.keep[-ix.rm];
    return(list(ix.keep=ix.keep,ix.rm=ix.rm));
  }
