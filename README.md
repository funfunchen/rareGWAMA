# A rareGWAMA R package
## and its Trans-ethnic TWAS application

**Table of Contents**

- [Introduction](#introduction)
- [Citation](#citation)
- [Installation](#Installing-the-rareGWAMA-R-package)
- [Quick tutorial](#quick-tutorial)
- [Feedback/Contact](#Feedback/Contact)


## Introduction

**rareGWAMA** is a flexible, swiss army knife software package for imputation based GWAS meta-analysis.   
It is developed and maintained by [Dajiang Liu's Group](https://dajiangliu.blog/).

This repository is for `TESLA` (*T*rans-*E*thnic transcriptome-wide association *S*tudy approach using an optimal *L*inear combination of *A*ssociation statistics), an improved TWAS method that optimally integrates trans-ethnic GWAS data with eQTL datasets. 

TESLA is implemented in the `rareGWAMA` package, if you want to use all the `Rare-Variant Association Analysis` related functions in raraGWAMA, please refer to this [page](https://github.com/dajiangliu/rareGWAMA).

------------------------------------------------------
## Citations
Liu DJ*†, Peloso GM*, Zhan X*, Holmen O*, Zawistowski M, Feng S, Nikpay M, Auer PL, Goel A, Zhang H, Peters U, Farrall M, Orho-Melander M, Kooperberg C, McPherson R, Watkins H, Willer CJ, Hveem, K, Melander O, Kathiresan S, Abecasis GR†    
**Meta-analysis of gene-level tests of rare variant association, Nature Genetics, 46, 200–204 (2014)**  
[doi: 10.1038/ng.2852.](https://www.nature.com/articles/ng.2852)

**Trans-ethnic Transcriptome-wide Association Study for Smoking Addiction in 1.3 Million Individuals Yields Insights into Tobacco Use Biology and Drug Repurposing**
(in preparation)

------------------------------------------------------
## Installing the rareGWAMA R package <a name="Installing-the-rareGWAMA-R-package"></a>

The package is hosted on github, which allows installation and update to be very easy. First, make sure you have the `mvtnorm` and `data.table` packages installed:

    install.packages("devtools")

And also, you need the latest version of [seqminer](https://github.com/zhanxw/seqminer):

    library(devtools)
    devtools::install_github("zhanxw/seqminer")

Then you could use:

    install_github("funfunchen/rareGWAMA")
    
With `library(rareGWAMA)`, your are ready to go!


------------------------------------------------------
## Quick tutorial <a name="quick-tutorial"></a>

1. The [wiki](https://github.com/funfunchen/rareGWAMA/wiki) has detailed tutorials on [input formats](https://github.com/funfunchen/rareGWAMA/wiki/2.-Input-files-and-arguments), [analysis](https://github.com/funfunchen/rareGWAMA/wiki/3.-Analysis) and [results interpretation](https://github.com/funfunchen/rareGWAMA/wiki/4.-Results-interpretation).
2. The methods are described in the papers [(citations above)](#citations)

### TESLA analysis <a name="Single-variant-tests"></a>

1.The very basic test is using:  

```{r}
res.gene.ii <- rareGWAMA.gene(score.stat.file,
                                imp.qual.file=imp.qual.file,
                                vcf.ref.file, refFileFormat="vcf.vbi",
                                anno=anno.ii, annoType=c('-'),
                                rvtest='TESLA',
                                ref.ancestry=ref.ancestry, trans.ethnic=TRUE, study.ancestry=study.ancestry,
                                maf.cutoff=1,
                                study.ref.panel=study.ref.panel,
                                chrVcfPrefix=c("chr",""), chrSumstatPrefix="chr",
                                pc.no=3,
                                af.pca=as.matrix(af.pca), af.pca.eqtl=gtex.pca,
                                gc=TRUE, maf.bin=maf.bin, gc.lambda=gc.merge,
                                regMat.lambda=0.0);
```  

please find more details in the wiki: [input formats](https://github.com/funfunchen/rareGWAMA/wiki/2.-Input-files-and-arguments) for the arguments:
> * score.stat.file: The file paths of all the score statistic files, which is a **vector object**;
> * imp.qual.file: Default is `NULL`. The file names of imputation quality, which is a **vector object**;
> * vcf.ref.file: The file paths of the reference panel files, which is a **list object**. Each list element is a **data.frame** which contains all the vcf files paths for this panel;
> * refFileFormat: `vcf` or `vcf.vbi`;
> * anno: eQTL weights and annotation file, which is a **data.frame** object;
> * annoType: The annotation types: could be `Nonsynonymous|Stop|Splice` or just `-`;
> * rvtest: The Rare-Variant Association Testing you want to use (i.e. 'VT', 'BURDEN', 'SKAT'), usging `TESLA` here;
> * ref.ancestry: Individuals' ancestry information, which is a **list object**;
> * trans.ethnic: `True` for multi-ethnic analysis;
> * study.ancestry: Ancestry information for each study, which is a **vector object**;
> * maf.cutoff: The minor allele frequency cut off, `1` as default;
> * study.ref.panel: `ref.panel` used for each study, which is a **vector object**;
> * chrVcfPrefix: The prefix of chromosome colomn used for each ref panel, which is a **vector object**;
> * chrSumstatPrefix: What the prefix is, usually is `chr`;
> * pc.no: The number of PCs are incorporated in the analysis;
> * af.pca: MDS information (or PCA) of all the studies, which is a **matrix object**. Please add `a column of 1` to represent the intercept;
> * af.pca.eqtl: The MDS information of the eQTL data, which is a **vector object**. Please add `1` to represent the intercept;
> * gc: Default is `FALSE`;
> * maf.bin: The minor allelle frequency used, which is a **matrix object**;
> * gc.lambda: GC values for each study, which is a **data.frame** object;
> * regMat.lambda: `0` as default;  

2.The out put should be as follows:  
`head(res$res.formatted))`  

```
GENE	RANGE	STAT_PC0	STAT_PC1	STAT_PC2	STAT_PC3	PVALUE_PC0	PVALUE_PC1	PVALUE_PC2	PVALUE_PC3	PVALUE_TETWAS	MAF_CUTOFF	NUM_VAR	TOTAL_MAF	POS_VAR	N	POS_SINGLE_MINP	BETA_SINGLE_MINP	SD_SINGLE_MINP	col_20
ENSG00000115364_MRPL19	2:74691309-76284257	0.53226	0.29468	0.02075	0.00972	0.466	0.587	0.885	0.921	0.768	1	5	16.5	2:74691309_A/G,2:74743266_C/T,2:74743454_G/A,2:74743589_G/T,2:75966100_C/T	367600	2:75966100_C/T	-0.00710230700790038	0.00264051650131425	NA
ENSG00000176204_LRRTM4	2:75785917-78591153	0.764	1.461	3.015	2.500	0.3820	0.2267	0.0825	0.1138	0.151	1	9	25	2:75785917_C/T,2:75811099_C/T,2:75924603_A/G,2:75929335_A/G,2:75934448_T/C,2:75947603_A/G,2:75964796_G/A,2:75968106_G/A,2:77706383_A/G	366578	2:77706383_A/G	-0.003566617901943	0.0014587152918681	NA
ENSG00000042445_RETSAT	2:85060154-85931869	5.78	3.38	3.84	3.08	0.0162	0.0662	0.0502	0.0792	0.0249	1	26	8.24	2:85060154_G/A,2:85183810_C/T,2:85186540_T/C,2:85195299_C/T,2:85195870_A/C,2:85333360_G/A,2:85338511_T/C,2:85340538_C/T,2:85345910_G/A,2:85346280_C/A,2:85346892_G/A,2:85347139_C/A,2:85441189_G/A,2:85441865_A/G,2:85649558_A/G,2:85922642_A/C,2:85925585_A/C,2:85925776_C/T,2:85926546_A/G,2:85926881_C/T,2:85927621_G/A,2:85928121_C/T,2:85929854_T/C,2:85930472_G/T,2:85930874_C/A,2:85931869_T/G	366411	2:85333360_G/A	0.00374031622935748	0.00146376973618988	NA
ENSG00000042493_CAPG	2:84485399-86235624	6.09	1.62	1.97	1.35	0.0136	0.2027	0.1605	0.2460	0.0289	1	30	4.81	2:84485399_A/G,2:84580261_T/C,2:84588891_A/G,2:84641177_A/G,2:84670377_A/G,2:84679046_A/G,2:84689275_C/T,2:84731535_G/A,2:84779586_C/T,2:84792782_G/T,2:84812439_C/T,2:84814512_A/G,2:84819794_G/T,2:84846573_G/A,2:85310189_G/A,2:85349849_G/A,2:85366983_A/G,2:85383920_T/G,2:85389635_T/C,2:85394936_T/C,2:85395194_T/C,2:85398099_T/G,2:85411200_A/G,2:85420466_C/T,2:85806164_G/A,2:86141681_G/A,2:86165645_A/G,2:86216085_C/T,2:86228369_C/T,2:86235624_G/A	366790	2:85411200_A/G	-0.00408101301034469	0.00146471774904106	NA
ENSG00000115486_GGCX	2:84890543-86411977	0.5747	0.8634	0.0126	0.0260	0.448	0.353	0.911	0.872	0.662	1	41	14.5	2:84890543_G/A,2:85041479_C/T,2:85199986_G/A,2:85210822_C/T,2:85253575_T/G,2:85254276_G/A,2:85540612_C/T,2:85578244_C/A,2:85581614_A/G,2:85581748_C/T,2:85582866_C/T,2:85584106_G/A,2:85587861_C/T,2:85591365_T/C,2:85725907_C/A,2:85734723_C/T,2:85858022_G/A,2:85860886_T/C,2:85861843_G/A,2:85862014_C/A,2:85865130_T/C,2:85867056_G/T,2:85868078_G/A,2:85868309_C/T,2:85869687_C/A,2:85873484_G/A,2:85873888_T/C,2:85876532_C/T,2:85877606_C/T,2:85878029_T/C,2:85878375_A/G,2:85883080_A/G,2:85938320_T/C,2:86285865_G/A,2:86288998_C/T,2:86291126_A/C,2:86314543_A/C,2:86325796_G/T,2:86359649_T/C,2:86368878_A/G,2:86411977_C/T	365220	2:85725907_C/A	0.0026933772332081	0.00179705864001883	NA
```

3.For demo data, please see `?rareGWAMA.gene`.


------------------------------------------------------
## Feedback/Contact <a name="Feedback/Contact"></a>

Questions and requests can be sent to
Github issue page ([link](https://github.com/dajiangliu/rareGWAMA/issues))
or
Dajiang Liu ([dajiang.liu@outlook.com](mailto:dajiang.liu@outlook.com "mailto:dajiang.liu@outlook.com")) and Fang Chen([fchen1@hmc.psu.edu](mailto:fchen1@hmc.psu.edu))



------------------------------------------------------
## References

<a name="myfootnote1">1</a>: Xiaowei Zhan, Youna Hu, Bingshan Li, Goncalo R. Abecasis, and Dajiang J. Liu       
**RVTESTS: An Efficient and Comprehensive Tool for Rare Variant Association Analysis Using Sequence Data**      
Bioinformatics 2016 32: 1423-1426. [doi:10.1093/bioinformatics/btw079](http://bioinformatics.oxfordjournals.org/content/32/9/1423.short)  ([PDF](http://bioinformatics.oxfordjournals.org/content/32/9/1423.full.pdf+html))

<a name="myfootnote2">2</a>: Yang L, Jiang S, Jiang B, Liu DJ, Zhan X  
**Seqminer2: an efficient tool to query and retrieve genotypes for statistical genetics analyses from biobank scale sequence dataset**
Bioinformatics. 2020 Oct 1;36(19):4951-4. [doi:10.1093/bioinformatics/btaa628](https://academic.oup.com/bioinformatics/article-abstract/36/19/4951/5881355?redirectedFrom=fulltext)