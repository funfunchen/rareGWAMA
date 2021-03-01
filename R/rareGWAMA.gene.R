#' organize the formatted stat into analyzable format;
#'
#' @param dat formatted data;
#' @param maf.cutoff The minor allele frequency cutoff
#' @return dat which consist of formatted data for gene-level tests
#' @export
rareGWAMA.formatGene <- function(dat,...) {
    extraPar <- list(...);
    maf.cutoff <- extraPar$maf.cutoff;
    method <- extraPar$method;
    if(is.null(extraPar$twas)) extraPar$twas <- FALSE;
    twas <- extraPar$twas;
    pos <- extraPar$pos;
    if(is.null(extraPar$regMat.lambda)) extraPar$regMat.lambda <- 0;
    if(is.null(method)) method <- "PCBS";
    if(is.null(maf.cutoff)) maf.cutoff <- 0.05;
    if(is.null(extraPar$estJointEff)) extraPar$estJointEff <- FALSE;
    lambda <- extraPar$lambda;
    if(is.null(lambda)) lambda <- 0.1;

    af.meta <- rowSums(dat$af.mat*dat$nSample.mat,na.rm=TRUE)/rowSums(dat$nSample.mat,na.rm=TRUE)
    maf.meta <- rm.na(af.meta);
    
    maf.meta[which(maf.meta>.5)] <- 1-maf.meta[which(maf.meta>.5)];
    pos.rare <- intersect(gsub('chr','',dat$pos[which(maf.meta<maf.cutoff)]),gsub('chr','',pos));
    
    ix.rare <- match(gsub('chr','',pos.rare),gsub('chr','',dat$pos));
    ix.rare.r2 <- match(gsub('chr','',pos.rare),gsub('chr','',colnames(dat$r2)));
    twas.weight <- NULL;
    if(twas==TRUE) {
        anno <- extraPar$anno;
        pos.anno <- paste(paste(anno$chrom,anno$pos,sep=":"),paste(anno$ref,anno$alt,sep="/"),sep="_");
       
        ix.rare.pos <- match(gsub('chr','',pos.rare),gsub('chr','',pos.anno))
        twas.weight <- anno$weight[ix.rare.pos];
    }
    if(length(ix.rare)>0) {
        
        r2.rare <- as.matrix(dat$r2[ix.rare.r2,ix.rare.r2]);
        ustat.meta <- rowSums(matrix(dat$ustat.mat[ix.rare,],nrow=length(ix.rare)),na.rm=TRUE);
        N.ustat.meta <- rowSums(matrix(dat$nSample.mat[ix.rare,],nrow=length(ix.rare)),na.rm=TRUE);
        
        V.meta <- matrix(0,nrow=length(ix.rare),ncol=length(ix.rare));
        V.meta.reg <- V.meta;
        N.V.meta <- matrix(0,nrow=length(ix.rare),ncol=length(ix.rare));
        V.by.study <- list();
        r2.by.study <- list();
        for(ii in 1:ncol(dat$ustat.mat)) {

            sd.vec <- rm.na(dat$vstat.mat[ix.rare,ii]);
            if(extraPar$trans.ethnic==TRUE) {
                
                ix.ref.panel <- which(names(extraPar$r2.list)==extraPar$study.ref.panel[ii]);
            
                ix.study <- which(names(extraPar$r2.list[[ix.ref.panel]])==extraPar$study.ancestry[ii]);
                ix.rare.r2.ii <- match(gsub('chr','',pos.rare),gsub('chr','',colnames(extraPar$r2.list[[ix.ref.panel]][[ix.study]])));
                
                V.meta <- V.meta+cor2cov(rm.na(extraPar$r2.list[[ix.ref.panel]][[ix.study]][ix.rare.r2.ii,ix.rare.r2.ii]),sd.vec);
                V.by.study[[ii]] <- regMat(cor2cov((extraPar$r2.list[[ix.ref.panel]][[ix.study]][ix.rare.r2.ii,ix.rare.r2.ii]),sd.vec),extraPar$regMat.lambda);
                r2.by.study[[ii]] <- regMat(as.matrix((extraPar$r2.list[[ix.ref.panel]][[ix.study]][ix.rare.r2.ii,ix.rare.r2.ii])),extraPar$regMat.lambda);
            }
            if(extraPar$trans.ethnic==FALSE) {
                V.meta <- V.meta+cor2cov(r2.rare,sd.vec);
                V.by.study[[ii]] <- cor2cov(r2.rare,sd.vec);
                r2.by.study[[ii]] <- r2.rare;
            }
            V.meta.reg <- V.meta.reg+cor2cov(cov2cor(r2.rare+lambda*diag(nrow(r2.rare))),sd.vec);
            sqrt.nSample <- as.matrix(sqrt(rm.na(dat$nSample.mat[ix.rare,ii])));
            N.V.meta <- N.V.meta+sqrt.nSample%*%t(sqrt.nSample);
        }

        beta.est <- NULL;beta.var <- NULL;beta.est.reg <- NULL;beta.var.reg <- NULL;
        if(extraPar$estJointEff) {
            beta.est <- ginv(V.meta)%*%ustat.meta;
            beta.var <- ginv(V.meta);
            
            beta.est.reg <- ginv(V.meta.reg)%*%ustat.meta;
            beta.var.reg <- ginv(V.meta.reg);
        }
        V.meta.weighted <- rm.na(V.meta/N.V.meta);
        V.meta.weighted.reg <- rm.na(V.meta.reg/N.V.meta);
        
        ustat.meta.weighted <- ustat.meta/N.ustat.meta;
        
        weight <- matrix(0,nrow=length(ustat.meta.weighted),ncol=length(ustat.meta.weighted));
        diag(weight) <- 1/N.ustat.meta;
        var.ustat.meta.weighted <- weight%*%V.meta%*%weight;

        direction.singlevar.assoc <- rep("X",length(ustat.meta));
        direction.singlevar.assoc[which(ustat.meta>0)] <- "+";
        direction.singlevar.assoc[which(ustat.meta<0)] <- "-";
        
        dat$ustat.meta <- ustat.meta;
        dat$V.meta <- V.meta;
        dat$direction.singlevar.assoc <- direction.singlevar.assoc;
        dat$pos.rare <- pos.rare;
        dat$maf.meta <- maf.meta;
        dat$nSample.meta <- N.ustat.meta;
        dat$beta.est <- beta.est;
        dat$beta.var <- beta.var;
        dat$ustat.meta.weighted <- ustat.meta.weighted;
        dat$N.ustat.meta <- N.ustat.meta;
        dat$N.V.meta <- N.V.meta;
        dat$V.meta.weighted <- V.meta.weighted;
        dat$r2.rare <- r2.rare;
        dat$maf.rare <- maf.meta[ix.rare];
        dat$maf.meta <- dat$maf.rare;
        dat$pos.rare <- pos.rare;
        dat$ix.rare <- ix.rare;
        dat$beta.est.reg <- beta.est.reg;
        dat$beta.var.reg <- beta.var.reg;
        dat$twas.weight <- twas.weight;
        sdG <- matrix(0,nrow=length(dat$maf.rare),ncol=length(dat$maf.rare));
        diag(sdG) <- sqrt(2*(dat$maf.rare)*(1-dat$maf.rare));
        covG <- sdG%*%dat$r2.rare%*%sdG;
        dat$r2.by.study <- r2.by.study;
        dat$V.by.study <- V.by.study;
        
    }
    return(dat);
    
}

#' Conduct approximate gene-level tests and estimate variance explained; 
#' @param score.stat.file the file names of score statistic files;
#' @param imp.qual.file the file names of imputation quality;
#' @param vcf.ref.file the file names of the reference panel file;
#' @param tabix.range Tabix range for variants;
#' @return formatted data and assocation testing results; 
#' @export
rareGWAMA.gene <- function(score.stat.file,imp.qual.file=NULL,vcf.ref.file,anno,gc=FALSE,...) {
    uniq.allele <- function(x) {x.tab <- table(x);return(paste(names(x.tab),sep=',',collapse=','))}
    extraPar <- list(...);
    anno$chrom <- gsub("chr","",anno$chrom); # format the chr column
    maf.cutoff <- extraPar$maf.cutoff;
    if(is.null(maf.cutoff)) maf.cutoff <- 1;
    r2.cutoff <- extraPar$r2.cutoff;
    if(is.null(r2.cutoff)) r2.cutoff <- 0.95;
    sizePerBatch <- extraPar$sizePerBatch;
    if(is.null(sizePerBatch)) sizePerBatch <- 100;
    if(is.null(extraPar$trans.ethnic)) {
        extraPar$trans.ethnic <- FALSE;
    }
    if(is.null(extraPar$pc.no)) {
        extraPar$pc.no <- 3;
    }
    if(is.null(extraPar$twas)) extraPar$twas <- FALSE;
    refGeno <- extraPar$refGeno;
    col.impqual <- extraPar$col.impqual;
    if(is.null(extraPar$regCovTWAS)) extraPar$regCovTWAS <- FALSE;
    if(is.null(extraPar$rvtest)) extraPar$rvtest <- "OMNI";
    
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
    beta.est <- 0;beta.se <- 0;statistic <- 0;p.value <- 0;ref.tab <- 0;alt.tab <- 0;pos.all <- 0;marginal.statistic <- 0;marginal.p.value <- 0;
    ii <- 0;batchEnd <- 0;batchStart <- 0;nSample <- 0;af <- 0;numStability <- 0;
    colnames(anno) <- tolower(colnames(anno));
    gene.vec <- unique(anno$gene);
    ix.intergenic <- which(gene.vec=="Intergenic");
    if(length(ix.intergenic)>0) gene.vec <- gene.vec[-ix.intergenic];
    numBatch <- as.integer(length(gene.vec)/sizePerBatch);
    if(numBatch*sizePerBatch<length(gene.vec)) numBatch <- numBatch+1;
    col.no <-  14 + (2 * extraPar$pc.no) # setup results column numbers 
    res.out <- matrix(nrow=length(gene.vec),ncol=20);
    res.out[,1] <- gene.vec;
    
    for(ii in 1:numBatch) {
        batchStart <- (ii-1)*sizePerBatch+1;
        batchEnd <- ii*sizePerBatch;
        if(batchEnd>length(gene.vec)) batchEnd <- length(gene.vec);
        gene.ii <- gene.vec[batchStart:batchEnd];
        ix.var <- which(anno$gene%in% gene.ii & anno$anno%in%extraPar$annoType);
        if(length(ix.var)>0) {
            anno.ii <- anno[ix.var,];
            pos.ii <- paste0(anno$chrom[ix.var],":",anno$pos[ix.var],"_",anno$ref[ix.var],"/",anno$alt[ix.var]);
            tabix.ii <- get.tabix.range(paste0(paste0(extraPar$chrSumstatPrefix,anno$chrom[ix.var]),":",anno$pos[ix.var]));
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
            capture.output(raw.data.all <- rvmeta.readDataByRange( score.stat.file, NULL, tabix.ii,multiAllelic = TRUE));
            ix.vcf <- unique(anno$chrom[ix.var]);
            
            ix.vcf[which(ix.vcf=="X")] <- 23;
            ix.vcf <- as.numeric(ix.vcf);

            if(extraPar$refFileFormat=="vcf") {
                vcfIndv <- refGeno;
                
                vcfColumn <- c("CHROM","POS","REF","ALT");
                vcfInfo <- NULL;      
                geno.list <- readVCFToListByRange(vcf.ref.file[ix.vcf], tabix.ii, "", vcfColumn, vcfInfo, vcfIndv)
                pos.vcf <- paste0(geno.list$CHROM,":",geno.list$POS,"_",geno.list$REF,"/",geno.list$ALT);
                pos.vcf <- gsub("chr",'',pos.vcf);
            }

            if(extraPar$refFileFormat=="vcf.vbi") {
                geno.tmp <- list();
                var.count <- 0;
                pos.tmp <- character(0);
                gt.list <- list();
                geno.list <- list();
                ix.vcf.ori <- ix.vcf;
                pos.vcf.list <- list();
                for(bb in 1:length(vcf.ref.file)) {
                    geno.tmp[[bb]] <- list();
                    ix.vcf <- ix.vcf.ori;
                    var.count <- 0;
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
                    ix.0 <- which(var.count==0);
                    if(length(ix.0)>0) {
                        var.count <- var.count[-ix.0];
                        ix.vcf <- ix.vcf[-ix.0];
                    }
                    gt.list[[bb]] <- matrix(nrow=nrow(tmp[[1]]),ncol=length(pos.tmp));
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
            cat('Read in',length(raw.data.all$ref[[1]]),'variants from summary statistics\n',sep=' ');
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
            for(jj in 1:length(gene.ii)) {
                ix.gene <- (jj-1)+batchStart;
                ix.jj <- which(anno$gene[ix.var] == gene.ii[jj]);
                pos.jj <- intersect(intersect(gsub("chr","",dat$pos),gsub("chr","",pos.vcf.list[[1]])),gsub("chr","",pos.ii[ix.jj]));
                ix.match <- match(gsub("chr","",pos.jj),gsub("chr","",pos.vcf));
                anno.jj <- anno.ii[match(gsub("chr","",pos.jj),gsub("chr","",pos.ii)),]

               
                if(length(ix.match)>0) {
                    r2 <- corSparse(Matrix(gt.list[[1]][,ix.match],ncol=length(ix.match)));
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
                                pos.jj <- intersect(intersect(gsub("chr","",dat$pos),gsub("chr","",pos.vcf.list[[ff]])),gsub("chr","",pos.ii[ix.jj]));
                                ix.match <- match(gsub("chr","",pos.jj),gsub("chr","",pos.vcf.list[[ff]]));
                                r2.list[[ff]][[ee]] <- corSparse(rm.na.colMeans(Matrix(gt.list[[ff]][ix.match.sample,ix.match],ncol=length(ix.match))));
                                r2.list[[ff]][[ee]][which(is.na(r2.list[[ff]][[ee]]) | r2.list[[ff]][[ee]]==Inf | r2.list[[ff]][[ee]]==-Inf)] <- 0;
                                diag(r2.list[[ff]][[ee]]) <- 1;
                                colnames(r2.list[[ff]][[ee]]) <- pos.jj;
                            }
                            names(r2.list[[ff]]) <- ancestry.grp;
                        }
                    }
                    names(r2.list) <- names(vcf.ref.file);
                    diag(r2) <- 1;
                    r2 <- rm.na(r2);
                    dat$r2 <- r2;
                    if(extraPar$twas==FALSE) {
                        dat.jj <- rareGWAMA.formatGene(dat,maf.cutoff=maf.cutoff,lambda=extraPar$lambda,pos=pos.jj,r2.list=r2.list,trans.ethnic=extraPar$trans.ethnic,study.ancestry=extraPar$study.ancestry,study.ref.panel=extraPar$study.ref.panel);
                    }
                    if(extraPar$twas==TRUE) {
                        dat.jj <- rareGWAMA.formatGene(dat,maf.cutoff=maf.cutoff,lambda=extraPar$lambda,pos=pos.jj,r2.list=r2.list,trans.ethnic=extraPar$trans.ethnic,study.ancestry=extraPar$study.ancestry,study.ref.panel=extraPar$study.ref.panel,anno=anno.jj,twas=extraPar$twas);

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
                    if(extraPar$rvtest=="OMNI") {
                        maf.vec <- dat.jj$maf.rare;
                        res.test <- list(statistic=NA,p.value=NA,maf.cutoff=extraPar$maf.cutoff);

                        if(length(maf.vec)>0) {
                            maf.TH <- sort(unique(maf.vec));
                            max.TH <- extraPar$max.TH;
                            if(length(max.TH)>0) {
                                if(length(maf.TH)>max.TH) {
                                    ix.tmp <- as.integer(seq(1,length(maf.TH),length=max.TH));
                                    maf.TH <- maf.TH[ix.tmp];
                                    maf.TH <- unique(sort(maf.TH));
                                }
                            }
                            ustat.meta <- dat.jj$ustat.meta
                            alternative <- "two.sided";
                            pval.burden <- 0;pval.skat <- 0;
                            pval.both <- 0;
                            for(tt in 1:length(maf.TH)) {
                                dat.jj.maf <- dat.jj;
                                ix.var <- which(dat.jj$maf.meta<=maf.TH[tt]);
                                dat.jj.maf$ustat.meta <- dat.jj.maf$ustat.meta[ix.var];
                                
                                dat.jj.maf$V.meta <- as.matrix(dat.jj.maf$V.meta[ix.var,ix.var]);
                                
                                dat.jj.maf$maf.meta <- dat.jj.maf$maf.meta[ix.var];

                                res.test.burden <- rareGWAMA.burden(dat.jj.maf);
                                res.test.skat <- rareGWAMA.skat(dat.jj.maf);
                                pval.burden[tt] <- res.test.burden$p.value;
                                pval.skat[tt] <- res.test.skat$p.value;
                                pval.both[tt] <- cauchy.p(as.numeric(c(pval.burden[tt],pval.skat[tt])));
                            }
                        
                            res.test <- cauchy.p(pval.both)
                            maf.cutoff <- maf.TH[which.min(pval.both)];
                            res.test$maf.cutoff <- maf.cutoff;
                        }
                    }
                    if(extraPar$rvtest=='TE-TWAS') {
                        res.test <- rareGWAMA.twas(dat.jj,extraPar$af.pca,extraPar$af.pca.eqtl);
                    }
                    
                    gene.jj.pos.range <- paste(range(gsub(".*:","",gsub("_.*","",pos.jj))),sep="-",collapse="-");
                    gene.jj.pos <- paste0(unique(gsub(":.*","",pos.jj)),":",gene.jj.pos.range);
                    if(extraPar$rvtest=='TE-TWAS') {
                        
                        if(extraPar$regCovTWAS==TRUE) {
                            jj.out <- c(gene.jj.pos,format(res.test$statistic,digits=3),format(res.test$p.value,digits=3),format(res.test$reg.p.value.minp,digits=3),res.test$maf.cutoff,length(dat.jj$pos.rare),format(sum(dat.jj$maf.rare),digits=3),paste(dat.jj$pos.rare,sep=",",collapse=","),as.integer(mean(dat.jj$nSample.meta,na.rm=TRUE)),pos.min, beta.single[ix.min],sd.single[ix.min]);
                        }
                        if(extraPar$regCovTWAS==FALSE) {
                            jj.out <- c(gene.jj.pos,format(res.test$statistic,digits=3),format(res.test$p.value,digits=3),format(res.test$p.value.minp,digits=3),res.test$maf.cutoff,length(dat.jj$pos.rare),format(sum(dat.jj$maf.rare),digits=3),paste(dat.jj$pos.rare,sep=",",collapse=","),as.integer(mean(dat.jj$nSample.meta,na.rm=TRUE)),pos.min, beta.single[ix.min],sd.single[ix.min]);
                        }
                    }
                    if(extraPar$rvtest!='TE-TWAS') { 
                        jj.out <- c(gene.jj.pos,format(res.test$statistic,digits=3),format(res.test$p.value,digits=3),res.test$maf.cutoff,length(dat.jj$pos.rare),format(sum(dat.jj$maf.rare),digits=3),paste(dat.jj$pos.rare,sep=",",collapse=","),as.integer(mean(dat.jj$nSample.meta,na.rm=TRUE)),pos.min, beta.single[ix.min],sd.single[ix.min]);
                    }
                    res.out[ix.gene,2:(length(jj.out)+1)] <- jj.out;
                }
            }    
        }   
    }
    if(extraPar$rvtest!="TE-TWAS") {
        colnames(res.out) <- paste('col',1:ncol(res.out),sep="_");
        colname.tmp <- c("GENE","RANGE","STAT","P-VALUE","MAF_CUTOFF","NUM_VAR","TOTAL_MAF","POS_VAR","N","POS_SINGLE_MINP","BETA_SINGLE_MINP","SD_SINGLE_MINP")
        colnames(res.out)[1:(length(colname.tmp))] <- colname.tmp;
    }
    if(extraPar$rvtest=="TE-TWAS") {
        colnames(res.out) <- paste('col',1:ncol(res.out),sep="_");
        colname.tmp <- c("GENE","RANGE",paste("STAT_PC",0:(ncol(extraPar$af.pca)-1),sep=''),paste("PVALUE_PC",0:(ncol(extraPar$af.pca)-1),sep=''),"PVALUE_TETWAS","MAF_CUTOFF","NUM_VAR","TOTAL_MAF","POS_VAR","N","POS_SINGLE_MINP","BETA_SINGLE_MINP","SD_SINGLE_MINP");
        colnames(res.out)[1:length(colname.tmp)] <- colname.tmp;
    }

    return(list(res.formatted=res.out,
                beta.single=beta.single,
                sd.single=sd.single,
                r2.list=r2.list,
                pval.single=pchisq((beta.single/sd.single)^2,df=1,lower.tail=FALSE)));
                
}
#' Estimate the variance explained by the variants in a locus;
#'
#' @param r2.rare the r2 matrix for rare variants;
#' @param maf.rare the maf vector for rare variants;
#' @return r2 estimates;
#' @export
estimateH2 <- function(...) {
    dat <- list(...);
    r2 <- rm.na(dat$r2.rare)
    maf <- rm.na(dat$maf.rare)
   
    beta.est <- rm.na(as.vector(dat$beta.est))
    beta.var <- rm.na(dat$beta.var);
   
    h2 <- 0;h2.var <- 0;
    
    if(!is.null(r2)) {
        sd.mat <- matrix(0,nrow=length(maf),ncol=length(maf));
        diag(sd.mat) <- sqrt(2*maf*(1-maf));
        v.mat <- sd.mat%*%r2%*%sd.mat
        h2 <- t(beta.est)%*%v.mat%*%beta.est;
        svd.v <- svd(v.mat);
        svd.omega <- svd(beta.var%*%v.mat%*%beta.var);
        h2.var <- sum(2*(svd.v$d)^2*svd.omega$d);
    }
    return(list(h2.est=as.numeric(h2),
                h2.var=h2.var));
    
}

#' gene-level test
#'
#' @param dat
#' @return a list with statistics and p-values;
#' @export
rareGWAMA.burden <- function(dat,...) {
    
    if(is.null(dat$ustat.meta)) {
        return(list(statistic=NULL, p.value=NULL));
    }

    ustat.burden <- sum(dat$ustat.meta,na.rm=TRUE);
    vstat.sq.burden <- (sum(dat$V.meta,na.rm=TRUE));
    statistic <- ustat.burden^2/vstat.sq.burden;
    p.value <- pchisq(statistic,df=1,lower.tail=FALSE);
          
    return(list(statistic=statistic,
                p.value=p.value));
}



#' t test;
#'
#' @param dat;
#' @return a list with statistic and p-value;
#' @export
rareGWAMA.t <- function(dat,...) {
    statistic <- t(dat$ustat.meta)%*%ginv(dat$V.meta)%*%dat$ustat.meta;
    p.value <- pchisq(statistic,df=length(dat$ustat.meta),lower.tail=FALSE);
    return(list(statistic=statistic,
                p.value=p.value));

}


#' SKAT test
#' @param dat
#' @return a list with statistics and p-values;
#' @export
rareGWAMA.skat <- function(dat,...) {
    extraPar <- list(...);
    weight <- extraPar$weight;
    if(is.null(weight)) weight <- 'linear';
    
    if(is.null(dat$ustat.meta)) return(list(statistic=NULL, p.value=NULL));
    W <- matrix(0,nrow=length(dat$ustat.meta),ncol=length(dat$ustat.meta));
    
    if(weight=='beta')  diag(W) <- dbeta(dat$maf.meta,1,25)^2;
    if(weight=='linear') diag(W) <- 1;
    if(weight=='twas') diag(W) <- (extraPar$twas.weight)^2;
    Q <- sum((dat$ustat.meta)^2*diag(W));

    svd.V <- svd(dat$V.meta);
    lambda.V <- matrix(0,nrow=length(abs(svd.V$d)),ncol=length(abs(svd.V$d)));
    diag(lambda.V) <- abs(svd.V$d);
    
    L <- (svd.V$u)%*%(sqrt(lambda.V))%*%t(svd.V$v);
    lambda <- try(get.eigen(W,L,t(L)),silent=TRUE);
    if(class(lambda)=='try-error')
    {
        return(list(statistic=NA,
                    p.value=NA));
    }
    p.value.liu <- NA;
    p.value.davies <- NA;
    p.value.imhof <- NA;
    p.value <- NA;
    p.value.davies <- try(davies(Q,lambda=lambda)$Qq,silent=TRUE);
    p.value.liu <- try(liu(Q,lambda=lambda),silent=TRUE);
    p.value.imhof <- try(imhof(Q,lambda=lambda)$Qq,silent=TRUE);
    if(length(attr(p.value.davies,'class'))+length(attr(p.value.liu,'class'))+length(attr(p.value.imhof,'class'))>0)
        return(list(statistic=NA,
                    p.value=NA));
    p.value <- p.value.davies;
    if(p.value<=0 | p.value>=1) p.value <- p.value.liu;
    return(list(statistic=Q,
                p.value=p.value));

}

#' vt test;
#'
#' @param dat
#' @return a list with statistics and p-values;
#' @export
rareGWAMA.vt <- function(dat,...) {
    if(is.null(dat$ustat.meta)) return(list(statistic=NULL, p.value=NULL));
    extraPars <- list(...);
    maf.vec <- dat$maf.meta;
    maf.TH <- sort(unique(maf.vec));
    max.TH <- extraPars$max.TH;
    ustat.meta <- dat$ustat.meta
    alternative <- "two.sided";
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
        return(list(p.value=NA,
                    statistic=NA,
                    maf.cutoff.vt=0,
                    no.site.VT=0))
    }
    if(length(maf.TH)==1)
    {
        ix.ii <- which(maf.vec<=maf.TH);    
        dat$ustat.meta <- sum(ustat.meta[ix.ii]);
        ind.ii <- rep(0,length(maf.vec));
        ind.ii[ix.ii] <- 1;
        V.stat.sq <- as.numeric(t(ind.ii)%*%dat$V.meta%*%(ind.ii));
 
        statistic <- dat$ustat.meta^2/V.stat.sq;
        p.value <- pchisq(statistic,df=1,lower.tail=FALSE);
        
        return(list(p.value=p.value,
                    statistic=statistic,
                    maf.cutoff.vt=maf.TH,
                    no.site.VT=length(ix.ii)));
    }
    err.msg <- vector(length=0);
    ix.list <- list();
    ustat.meta.VT <- 0;
    VT.mat <- matrix(0,nrow=length(maf.TH),ncol=ncol(dat$V.meta));
    no.TH <- length(maf.TH);
    for(ii in 1:length(maf.TH))
    {
        
        ix.ii <- which(maf.vec<=maf.TH[ii]);

        ix.list[[ii]] <- ix.ii;
        VT.mat[ii,ix.ii] <- 1;
        ustat.meta.VT[ii] <- sum(ustat.meta[ix.ii]);
    }
    cov.X.VT <- matrix(0,nrow=length(ix.list),ncol=length(ix.list));
    cov.mat <- dat$V.meta;
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
    
    vt.stat.vec <- (ustat.meta.VT^2/(diag(cov.X.VT)));
    pval.vt.stat.vec <- pchisq(vt.stat.vec,df=1,lower.tail=FALSE);
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
    p.value <- NA;
    if(!is.na(vt.max.stat)){
        p.value <- pvt(vt.max.stat,mu=rep(0,nrow(cor.X.VT)),sigma=cor.X.VT,alternative);
    }
    return(list(p.value=p.value,
                statistic=vt.max.stat,
                maf.cutoff.vt=maf.cutoff,
                pval.vt.stat.vec=pval.vt.stat.vec,
                no.site.VT=length(ix.ii)));
    
}


rm.na.colMeans <- function(mat) {
    
    aa <- colMeans(as.matrix(mat),na.rm=TRUE);
    for(ii in 1:ncol(mat)) {
        ix.na <- which(is.na(mat[,ii]));
        if(length(ix.na)>0) 
            mat[ix.na,ii] <- rep(aa[ii],length(ix.na));
    }
    return(mat);
}
