#' @param score.stat.file the file names of score statistic files;
#' @param imp.qual.file the file names of imputation quality;
#' @param vcf.ref.file the file names of the reference panel file;
#' @param tabix.range Tabix range for variants;
#' @return formatted data and assocation testing results; 
#' @export
rareGWAMA.conditional.gene <- function(score.stat.file,imp.qual.file=NULL,vcf.ref.file,anno,gc=FALSE,...) {
  uniq.allele <- function(x) {x.tab <- table(x);return(paste(names(x.tab),sep=',',collapse=','))}
  extraPar <- list(...);
  maf.cutoff <- extraPar$maf.cutoff;
  if(is.null(maf.cutoff)) maf.cutoff <- 1;
  r2.cutoff <- extraPar$r2.cutoff;
  if(is.null(r2.cutoff)) r2.cutoff <- 0.95;
  sizePerBatch <- extraPar$sizePerBatch;
  if(is.null(sizePerBatch)) sizePerBatch <- 100;
  if(is.null(extraPar$trans.ethnic)) {
    extraPar$trans.ethnic <- FALSE;
  }
  if(is.null(extraPar$maxNumVar)) extraPar$maxNumVar <- 500;
  if(is.null(extraPar$twas)) extraPar$twas <- FALSE;
  refGeno <- extraPar$refGeno;
  col.impqual <- extraPar$col.impqual;
  if(is.null(extraPar$regCovTWAS)) extraPar$regCovTWAS <- TRUE;
  if(is.null(extraPar$rvtest)) extraPar$rvtest <- "BURDEN";
  if(extraPar$rvtest=="TE-TWAS" | extraPar$rvtest=="TESLA") extraPar$twas <- TRUE;
  
  if(is.null(extraPar$refFileFormat)) extraPar$refFileFormat <- "vcf";
  gc.lambda <- matrix(rep(1,length(score.stat.file)),ncol=1);
  maf.bin <- matrix(c(0,1),ncol=2);
  if(gc==TRUE) {
    gc.lambda <- extraPar$gc.lambda;
    if(is.null(gc.lambda)) {
      gc.lambda <- matrix(1,nrow=length(score.stat.file),ncol=nrow(maf.bin));
    }
    maf.bin <- extraPar$maf.bin;
    if(is.null(maf.bin)) stop("maf bin must be provided");
  }
  impQual.lb <- extraPar$impQual.lb;
  impQualWeight <- FALSE;
  rmMultiAllelicSite <- extraPar$rmMultiAllelicSite;
  if(is.null(col.impqual)) col.impqual <- 5;
  if(is.null(impQual.lb)) impQual.lb <- 0.7;
  if(is.null(rmMultiAllelicSite)) rmMultiAllelicSite <- TRUE;
  if(is.null(refGeno)) refGeno <- "GT";
  pseudoScore <- extraPar$pseudoScore;
  if(is.null(pseudoScore)) pseudoScore <- TRUE;
  anno$chrom <- gsub("chr","",anno$chrom);
  
  beta.est <- 0;beta.se <- 0;statistic <- 0;p.value <- 0;ref.tab <- 0;alt.tab <- 0;pos.all <- 0;marginal.statistic <- 0;marginal.p.value <- 0;
  ii <- 0;batchEnd <- 0;batchStart <- 0;nSample <- 0;af <- 0;numStability <- 0;
  colnames(anno) <- tolower(colnames(anno));
  
  ix.intergenic <- which(gene.vec=="Intergenic");
  if(length(ix.intergenic)>0) gene.vec <- gene.vec[-ix.intergenic];
  numBatch <- length(extraPar$geneGroup);
  res.out <- matrix(nrow=length(extraPar$geneGroup),ncol=20);
  res.out[,1] <- extraPar$geneGroup;
  
  for(ii in 1:numBatch) {
    batchStart <- ii;
    batchEnd <- ii;
    #if(batchEnd>length(gene.vec)) batchEnd <- length(gene.vec);
    gene.ii <- unlist(strsplit(extraPar$geneGroup[ii],split=","));
    ix.var <- which(anno$gene%in% gene.ii & anno$anno%in%extraPar$annoType);
    if(length(ix.var)>0) {
      anno.ii <- anno[ix.var,];
      pos.ii <- paste0(anno$chrom[ix.var],":",anno$pos[ix.var],"_",anno$ref[ix.var],"/",anno$alt[ix.var]);
      # tabix.ii <- get.tabix.range(paste0(anno$chrom[ix.var],":",anno$pos[ix.var]));
      # tabix.ii.multichr <- tabix.ii;
      tabix.ii <- get.tabix.range(paste0(paste0(extraPar$chrSumstatPrefix, anno$chrom[ix.var]),":",anno$pos[ix.var]));
      tabix.ii.multichr <- tabix.ii;
      
      bb <- 1;
      vcf.tabix.ii.multichr <- list();
      for(cc in 1:length(vcf.ref.file)) {
        vcf.tabix.ii.multichr[[cc]] <- get.tabix.range(paste0(paste0(extraPar$chrVcfPrefix[cc],anno$chrom[ix.var]),":",anno$pos[ix.var])); 
      }
      if(length(unique(anno$chrom[ix.var]))>1) {
        for(aa in unique(anno$chrom[ix.var])) {
          ix.aa <- which(anno$chrom[ix.var]==aa);
          tabix.ii.multichr[bb] <- get.tabix.range(paste0(paste0(extraPar$chrSumstatPrefix,anno$chrom[ix.var[ix.aa]]),":",anno$pos[ix.var[ix.aa]]));
          for(cc in 1:length(vcf.ref.file)) {
            vcf.tabix.ii.multichr[[cc]][bb] <- get.tabix.range(paste0(paste0(extraPar$chrVcfPrefix[cc],anno$chrom[ix.var[ix.aa]]),":",anno$pos[ix.var[ix.aa]]));
          }
          bb <- bb+1;
        }
      }
      a <- Sys.time();
      cat('reading summary stat\n');
      f0 <- file();
      sink(file=f0,type='output')
      sink(file=f0,type='message');
      capture.output(raw.data.all <- rvmeta.readDataByRange( score.stat.file, NULL, tabix.ii,multiAllelic = TRUE));
      sink(type='output');
      sink(type='message');
      close(f0);
      ix.vcf <- unique(anno$chrom[ix.var]);
      
      ix.vcf[which(ix.vcf=="X")] <- 23;
      ix.vcf <- as.numeric(ix.vcf);
      cat('reading ref panel\n');
      if(extraPar$refFileFormat=="vcf") {
        vcfIndv <- refGeno;
        vcfColumn <- c("CHROM","POS","REF","ALT");
        vcfInfo <- NULL;
        f0 <- file();
        sink(file=f0,type='output')
        sink(file=f0,type='message');
        
        geno.list <- readVCFToListByRange(vcf.ref.file[ix.vcf], tabix.ii, "", vcfColumn, vcfInfo, vcfIndv)
        sink(type='output');
        sink(type='message');
        close(f0);
        pos.vcf <- paste0(geno.list$CHROM,":",geno.list$POS,"_",geno.list$REF,"/",geno.list$ALT);
        pos.vcf <- gsub("chr",'',pos.vcf);
      }
      
      if(extraPar$refFileFormat=="vcf.vbi") {
        geno.tmp <- list();
        var.count <- 0;
        pos.tmp <- character(0);
        gt.list <- list();
        ix.vcf.ori <- ix.vcf;
        
        geno.list <- list();
        pos.vcf.list <- list();
        
        for(bb in 1:length(vcf.ref.file)) {
          geno.tmp[[bb]] <- list();
          ix.vcf <- ix.vcf.ori;
          for(aa in 1:length(ix.vcf)) {
            f0 <- file();
            sink(file=f0,type='output')
            sink(file=f0,type='message');
            tmp <- try(seqminer::readSingleChromosomeVCFToMatrixByRange(fileName=vcf.ref.file[[bb]][ix.vcf[aa]],range=vcf.tabix.ii.multichr[[bb]][aa]),silent=TRUE);
            sink(type='output');
            sink(type='message');
            close(f0);
            
            geno.tmp[[bb]][[aa]] <- tmp[[1]];
            pos.tmp <- c(pos.tmp,colnames(tmp[[1]]));
            
            var.count[aa] <- ncol(tmp[[1]]);
            
          }
          if(aa > 1){
            ix.0 <- which(var.count==0);
            if(length(ix.0)>0) {
              var.count <- var.count[-ix.0];
              ix.vcf <- ix.vcf[-ix.0];
            }
          }
          gt.list[[bb]] <- Matrix(nrow=nrow(tmp[[1]]),ncol=length(pos.tmp));
          if(length(var.count)>0) {
            var.count <- cumsum(c(0,var.count));
            
            for(aa in 1:length(ix.vcf)) {
              if(var.count[aa+1]>var.count[aa] & var.count[aa]>=0)
                gt.list[[bb]][,(var.count[aa]+1):(var.count[aa+1])] <- geno.tmp[[bb]][[aa]];
            }
          }
          colnames(gt.list[[bb]]) <- pos.tmp;
          rownames(gt.list[[bb]]) <- rownames(tmp[[1]]);
          pos.vcf.list[[bb]] <- pos.tmp;
          pos.vcf.list[[bb]] <- gsub("chr",'',pos.vcf.list[[bb]]);
          geno.list[[bb]] <- list(sampleId = rownames(gt.list[[bb]]));
          
        }
      }
      
      pos.vcf <- unique(unlist(pos.vcf.list));
      
      raw.imp.qual <- NULL;
      if(!is.null(imp.qual.file))
        raw.imp.qual <- lapply(imp.qual.file,tabix.read.table,tabixRange=tabix.ii);
      time.readData <- Sys.time()-a;
      b <- Sys.time();
      raw.data.all <- raw.data.all[[1]];
      cat('Read in',length(raw.data.all$ref[[1]]),'variants in',format(b-a),'\n',sep=' ');
      
      a <- Sys.time();
      dat <- GWAMA.formatData(raw.data.all,raw.imp.qual,impQualWeight,impQual.lb,col.impqual,rmMultiAllelicSite=rmMultiAllelicSite);
      if(rmMultiAllelicSite==TRUE) {
        tmp <- GWAMA.rmMulti(dat);
        dat <- tmp$dat;posMulti <- tmp$posMulti;
      }
      maf.meta <- rowSums((dat$af.mat)*(dat$nSample.mat),na.rm=TRUE)/rowSums(dat$nSample.mat,na.rm=TRUE);
      maf.meta[which(maf.meta>0.5)] <- 1-maf.meta[which(maf.meta>0.5)];
      for(bb in 1:nrow(maf.bin)) {
        ix.bb <- which(maf.meta>=maf.bin[bb,1] & maf.meta<maf.bin[bb,2]);
        if(length(ix.bb)>0) {
          dat$ustat.mat[ix.bb,] <- scale(matrix(dat$ustat.mat[ix.bb,],nrow=length(ix.bb)),center=FALSE,scale=sqrt(gc.lambda[,bb]));
        }
      }
      
      pos <- gsub("_.*","",dat$pos);
      if(refGeno=="DS" & extraPar$refFileFormat=="vcf") {
        gt <- geno.list$DS;
        gt <- Matrix(as.numeric(gt),nrow=nrow(gt),ncol=ncol(gt));
        colnames(gt) <- pos.vcf;
        gt <- rm.na(gt);
      }
      if(refGeno=="GT" & extraPar$refFileFormat=="vcf") {
        gt.tmp <- geno.list$GT;
        gt <- Matrix(0,nrow=nrow(gt.tmp),ncol=ncol(gt.tmp));
        gt[which(gt.tmp=="0/0",arr.ind=T)] <- 0;
        gt[which(gt.tmp=="1/0",arr.ind=T)] <- 1;
        gt[which(gt.tmp=="0/1",arr.ind=T)] <- 1;
        gt[which(gt.tmp=="1/1",arr.ind=T)] <- 2
        gt[which(gt.tmp=="0|0",arr.ind=T)] <- 0;
        gt[which(gt.tmp=="1|0",arr.ind=T)] <- 1;
        gt[which(gt.tmp=="0|1",arr.ind=T)] <- 1;
        gt[which(gt.tmp=="1|1",arr.ind=T)] <- 2
        colnames(gt) <- pos.vcf;
        gt <- rm.na(gt);
        
      }
      jj <- 1;
      ix.gene <- (jj-1)+batchStart;
      pos.jj.tmp <- intersect(intersect(gsub("chr","",dat$pos), gsub("chr", "", pos.vcf)),
                              gsub("chr", "", pos.ii));
      
      ix.match.tmp <- match(pos.jj.tmp, gsub("chr", "", pos.vcf));
      anno.jj.tmp <- anno.ii[match(pos.jj.tmp, gsub("chr", "", pos.ii)),];
      all.gene.jj <- unique(anno.jj.tmp$gene);
      cor.tesla <- matrix(0,nrow=length(all.gene.jj),ncol=length(all.gene.jj));
      rownames(cor.tesla) <- all.gene.jj;
      colnames(cor.tesla) <- all.gene.jj;
      
      downsample.rep <- 1;
      downsample.rep <- 1;
      ds <- 1; {
        anno.jj <- anno.jj.tmp;
        pos.jj <- paste(paste(anno.jj$chrom,anno.jj$pos,sep=":"),paste(anno.jj$ref,anno.jj$alt,sep="/"),sep="_");
        # pos.jj <- intersect(intersect(gsub("chr","",dat$pos), 
        #                               gsub("chr", "", pos.vcf)),
        #                     gsub("chr", "", pos.jj));
        # ix.match <- match(gsub("chr", "", pos.jj), gsub("chr", "", pos.vcf));
        # ix.out <- which(ix.match > ncol(gt.list[[1]])) ## out of bound index, i.e. ix.match > ncol
        # ix.match <- ix.match[-ix.out]
        # pos.jj <- pos.jj[-ix.out]
        pos.jj <- intersect(intersect(gsub("chr","",dat$pos),
                                      gsub("chr", "", pos.vcf)),
                            gsub("chr", "", pos.jj));
        ix.match <- match(pos.jj, pos.vcf)
        cat('calculating LD\n');
        if(length(ix.match)>0) {
          r2 <- corSparse(Matrix(gt.list[[1]][,ix.match], ncol=length(ix.match)));
          r2[which(is.na(r2) | r2==Inf | r2==-Inf)] <- 0;
          diag(r2) <- 1;
          colnames(r2) <- pos.jj;
          rownames(r2) <- pos.jj;
          r2.list <- list();
          if(extraPar$trans.ethnic==TRUE) {
            for(ff in 1:length(vcf.ref.file)) {
              r2.list[[ff]] <- list();
              ancestry.grp <- unique(unlist(strsplit(extraPar$ref.ancestry[[ff]][,2],split=",")));
              for(ee in 1:length(ancestry.grp)) {
                id.ee <- extraPar$ref.ancestry[[ff]][grep(paste0("\\b",ancestry.grp[ee],"\\b"),extraPar$ref.ancestry[[ff]][,2]),1];
                ix.match.sample <- match(id.ee,geno.list[[ff]]$sampleId);
                r2.list[[ff]][[ee]] <- corSparse(Matrix(gt.list[[ff]][ix.match.sample,ix.match],ncol=length(ix.match)));
                
                r2.list[[ff]][[ee]][which(is.na(r2.list[[ff]][[ee]]) | r2.list[[ff]][[ee]]==Inf | r2.list[[ff]][[ee]]==-Inf)] <- 0;
                diag(r2.list[[ff]][[ee]]) <- 1;
                colnames(r2.list[[ff]][[ee]]) <- pos.jj;
              }
              names(r2.list[[ff]]) <- ancestry.grp;
            }
          }
          b <- Sys.time();
          cat('calc LD in ', format(b-a),'\n');
          
          names(r2.list) <- names(vcf.ref.file);
          diag(r2) <- 1;
          r2 <- rm.na(r2);
          dat$r2 <- r2;
          if(extraPar$twas==FALSE) {
            
            dat.jj <- rareGWAMA.formatGene(dat,maf.cutoff=maf.cutoff,lambda=extraPar$lambda,pos=pos.jj,r2.list=r2.list,trans.ethnic=extraPar$trans.ethnic,study.ancestry=extraPar$study.ancestry,study.ref.panel=extraPar$study.ref.panel);
            
          }
          if(extraPar$twas==TRUE) {
            cat('formating data\n');
            a <- Sys.time();
            dat.jj <- rareGWAMA.formatGene(dat,maf.cutoff=maf.cutoff,lambda=extraPar$lambda,pos=pos.jj,r2.list=r2.list,trans.ethnic=extraPar$trans.ethnic,study.ancestry=extraPar$study.ancestry,study.ref.panel=extraPar$study.ref.panel,anno=anno.jj,twas=extraPar$twas);
            b <- Sys.time();
            cat('format Data in ', format(b-a),'\n');
          }
          res.test <- list();
          stat.single <- (dat.jj$ustat.meta)^2/diag(dat.jj$V.meta);
          
          pval.single <- pchisq(stat.single,df=1,lower.tail=FALSE);
          beta.single <- (dat.jj$ustat.meta)/diag(dat.jj$V.meta);
          sd.single <- sqrt(1/diag(dat.jj$V.meta));
          ix.min <- which.min(pval.single);
          pos.min <- dat.jj$pos.rare[ix.min];
          
          if(extraPar$rvtest=="BURDEN") {
            res.test <- rareGWAMA.burden(dat.jj);
            res.test$maf.cutoff <- maf.cutoff;
          }
          if(extraPar$rvtest=="SKAT") {
            res.test <- rareGWAMA.skat(dat.jj);
            res.test$maf.cutoff <- maf.cutoff;
          }
          if(extraPar$rvtest=="VT") {
            res.test <- rareGWAMA.vt(dat.jj,max.TH=10);
          }
          if(extraPar$rvtest=='TE-TWAS' | extraPar$rvtest=="TESLA") {
            cat('calculating cov\n');
            a <- Sys.time();
            varList.tmp <- cbind(anno.jj$gene,paste(paste(anno.jj$chrom,anno.jj$pos,sep=":"),paste(anno.jj$ref,anno.jj$alt,sep="/"),sep="_"));
            ix.varList <- match(intersect(intersect(gsub("chr", "", varList.tmp[,2]), gsub("chr", "", dat.jj$pos)),
                                          gsub("chr", "", pos.jj)),
                                gsub("chr", "", varList.tmp[,2])); # remove out of bound var
            ix.varList <- match(intersect(gsub("chr", "", varList.tmp[,2]), 
                                          gsub("chr", "", dat.jj$pos)),
                                gsub("chr", "", varList.tmp[,2]));
            varList <- matrix(varList.tmp[ix.varList,],ncol=2);
            gcov.twas <- estimateTWAS.cor(dat.jj,varList,extraPar$af.pca,extraPar$af.pca.eqtl,extraPar$maxNumVar)$gcov.twas;
            cor.tesla.ds <- estimateTESLA.cor(gcov.twas);
            colnames(cor.tesla.ds) <- unique(varList[,1]);
            rownames(cor.tesla.ds) <- colnames(cor.tesla.ds);
            b <- Sys.time();
            cat('calculate cov in ',format(b-a),'\n');
            cor.tesla.tmp <- matrix(0,nrow=nrow(cor.tesla),ncol=ncol(cor.tesla));
            cor.tesla.tmp[match(rownames(cor.tesla.ds),all.gene.jj),match(rownames(cor.tesla.ds),all.gene.jj)] <- cor.tesla.ds;
            cor.tesla <- cor.tesla+cor.tesla.tmp;
          }
          
        }
      }
      if(extraPar$rvtest=='TE-TWAS' | extraPar$rvtest=="TESLA") {
        jj.out <- c(paste(colnames(cor.tesla),sep=',',collapse=','),paste(as.vector(cor.tesla/downsample.rep),sep=',',collapse=','));
        res.out[ii,1:length(jj.out)] <- jj.out;
      }
    }
  }
  
  if(extraPar$rvtest!="TE-TWAS") {
    colnames(res.out) <- paste('col',1:ncol(res.out),sep="_");
    
  }
  if(extraPar$rvtest=="TE-TWAS" | extraPar$rvtest=="TESLA") {
    colnames(res.out) <- paste('col',1:ncol(res.out),sep="_");
    colnames(res.out)[1:(length(jj.out))] <- c("GENE","COV");
  }   
  return(list(res.formatted=res.out));    
}
