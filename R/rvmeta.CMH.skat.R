rvmeta.CMH.skat <- function(X.T.times.Y.centered.list,X.T.times.X.list,maf.vec,cov.mat.list,var.Y.list,N.list,alternative,no.boot,alpha,kernel,...)
  {
      extraPar <- list(...);
    err.msg <- vector(length=0);
    if(alternative!="two.sided")
      {
        alternative <- 'two.sided';
        err.msg <- c(err.msg,'Alternative reset to two.sided');
      }
    if(length(kernel)>1) kernel <- "beta";
    X.T.times.X <- 0*diag(length(maf.vec));
    X.T.times.Y <- rep(0,length(maf.vec));
    if(length(maf.vec)<=1)
      {
          return(list(statistic=NA,
                      p.value=NA,
                      p.value.liu=NA,
                      p.value.imhof=NA,
                      p.value.davies=NA,
                      X.T.times.Y=X.T.times.Y,
                      X.T.times.X=X.T.times.X,
                      err.msg=err.msg));
      }
    if(kernel!="optimal-beta" & kernel!="optimal-linear")
      {
          if(kernel=='twas') {
              W <- diag(extraPar$twas.weight);
          }
        if(kernel=='linear')
          {
            W <- diag(length(maf.vec));
          }
        if(kernel=='beta')
          {
            W <- diag((dbeta(maf.vec,1,25))^2);
            if(length(maf.vec)==1) W <- as.matrix((dbeta(maf.vec,1,25))^2);
          }
        Q <- 0;lambda <- double(0);
        U <- rep(0,length(X.T.times.Y.centered.list[[1]]));
        V <- matrix(0,nrow=length(X.T.times.Y.centered.list[[1]]),ncol=length(X.T.times.Y.centered.list[[1]]));
        for(ii in 1:length(X.T.times.Y.centered.list))
          {
            U <- U+(X.T.times.Y.centered.list[[ii]]);
            V <- V+((N.list[[ii]])*(cov.mat.list[[ii]])*(var.Y.list[[ii]]));
          }

        Q <- as.numeric(t(U)%*%W%*%U);
        svd.V <- svd(V);
        lambda.V <- abs(svd.V$d);

        L <- (svd.V$u)%*%diag(sqrt(lambda.V))%*%t(svd.V$v);
        lambda <- try(get.eigen(W,L,t(L)),silent=TRUE);
        if(class(lambda)=='try-error')
          {
            return(list(statistic=NA,
                        p.value=NA,
                        p.value.liu=NA,
                        p.value.imhof=NA,
                        p.value.davies=NA,
                        err.msg=err.msg));
          }
        p.value.liu <- NA;
        p.value.davies <- NA;
        p.value.imhof <- NA;
        p.value <- NA;
        if(no.boot==0)
          {
            p.value.davies <- try(davies(Q,lambda=lambda)$Qq,silent=TRUE);
            p.value.liu <- try(liu(Q,lambda=lambda),silent=TRUE);
            p.value.imhof <- try(imhof(Q,lambda=lambda)$Qq,silent=TRUE);
            if(length(attr(p.value.davies,'class'))+length(attr(p.value.liu,'class'))+length(attr(p.value.imhof,'class'))>0)
              return(list(statistic=NA,
                          p.value=NA,
                          p.value.liu=NA,
                          p.value.imhof=NA,
                          p.value.davies=NA,
                          X.T.times.X=V,
                          err.msg=err.msg));
            p.value <- p.value.davies;
            if(p.value<=0 | p.value>=1) p.value <- p.value.liu;
          }
        if(no.boot>0)
          {

            Q.perm <- 0;
            for(ii in 1:length(X.T.times.Y.centered.list))
              {
                X.T.times.Y.perm <- rmvnorm(no.boot,mean=rep(0,nrow(X.T.times.X)),sigma=(N.list[[ii]])*(var.Y.list[[ii]])*(cov.mat.list[[ii]]));
                Q.perm <- Q.perm+diag(X.T.times.Y.perm%*%W%*%t(X.T.times.Y.perm));
              }
            if(alternative=='two.sided')
              p.value <- mean(Q.perm>Q)/2+mean(Q.perm>=Q)/2;
          }
      }
    if(kernel=="optimal-beta" | kernel=='optimal-linear')
      {
        if(length(X.T.times.Y.centered.list)>1)
          {
            warning("approx.mega test should be used");
            statistic <- NA;
            p.value <- NA;
          }
        X.T.times.X <- (N.list[[1]]*cov.mat.list[[1]]);
        J <- matrix(1,nrow=length(maf.vec),ncol=length(maf.vec));
        if(length(X.T.times.Y.centered.list)==1)
          {
            if(kernel=='optimal-linear')
              {
                W0 <- diag(length(maf.vec));
                if(length(maf.vec)==1) W0 <- as.matrix(1);
              }
            if(kernel=='optimal-beta')
              {
                W0 <- diag((dbeta(maf.vec,1,25))^2);
                if(length(maf.vec)==1) W0 <- as.matrix((dbeta(maf.vec,1,25))^2);
              }
            lambda.mat <- matrix(nrow=11,ncol=nrow(X.T.times.X));

            rho <- seq(0,1,by=0.1);
            for(ii in 1:11)
              {
               W <- (rho[ii]*W0+(1-rho[ii])*J)*var.Y.list[[1]];
               lambda.mat[ii,] <- svd(W%*%X.T.times.X)$d;
              }
            res.C.data <- .C('genericStat',
                             X_T_times_Y_permPt=as.double(X.T.times.Y.centered.list[[1]]),
                             X_T_times_XPt=as.double(X.T.times.X),
                             noMarkerPt=as.integer(nrow(X.T.times.X)),
                             Npt=as.integer(N.list[[1]]),
                             noBootPt=as.integer(1),
                             varYPt=as.double(var.Y.list[[1]]),
                             statisticPt=as.double(1.5),
                             alternativePt=as.integer(0),
                             testPt=as.integer(2),
                             extraPar=as.double(c(as.double(diag(W0)),
                               as.double(lambda.mat))));
            statData <- (res.C.data$statisticPt)[1];
            res.svd <- svd(X.T.times.X*var.Y.list[[1]])
            U <- res.svd$u;
            V <- res.svd$v;
            D <- res.svd$d;
            res.C.boot <- .C("genericBoot",
                             U <- as.double(U),
                             V <- as.double(V),
                             D <- as.double(D),
                             noMarkerPt=as.integer(nrow(X.T.times.X)),
                             X_T_times_XPt=as.double(X.T.times.X),
                             Npt=as.integer(N.list[[1]]),
                             noBootPt=as.integer(no.boot),
                             varYPt=as.double(var.Y.list[[1]]),
                             statisticPt=as.double(rep(1.5,1000)),
                             alternativePt=as.integer(0),
                             testPt=as.integer(2),
                             extraPar=as.double(c(statData,alpha,0,diag(W0),as.vector(lambda.mat))));

            p.value <- (res.C.boot$extraPar)[3];
            Q <- statData;
          }
      }
    return(list(statistic=Q,
                p.value=p.value,
                p.value.liu=p.value.liu,
                p.value.imhof=p.value.imhof,
                p.value.davies=p.value.davies,
                no.site=length(maf.vec),
                err.msg=err.msg));
  }
SKAT.perm.stat <- function(X.T.times.Y.perm,W,X.T.times.X,var.Y)
  {
    Q.perm <- diag(X.T.times.Y.perm%*%W%*%t(X.T.times.Y.perm));
    lambda <- svd(W%*%X.T.times.X*var.Y)$d;
    p.value <- davies.vec(Q.perm,lambda=lambda,acc=5e-5)$Qq;
    return(p.value)
  }

get.eigen <- function(kernel,res.svd.left,res.svd.right)
  {
    lambda1 <- (svd(res.svd.left%*%kernel%*%res.svd.right)$d);

    IDX1<-which(lambda1 >= 0)
    IDX2<-which(lambda1 > mean(lambda1[IDX1])/100000)
    if(length(IDX2) == 0){
      stop("No Eigenvalue is bigger than 0!!")
    }
    lambda<-lambda1[IDX2]
    return(lambda);
  }
