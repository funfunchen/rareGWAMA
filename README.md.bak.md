# A rareGWAMA R package

**Table of Contents**

- [Introduction](#introduction)
- [Citation](#citation)
- [Installation](#Installing-the-rareGWAMA-R-package)
- [Quick tutorial](#quick-tutorial)
    - [Single variant tests](#Single-variant-tests)
    - [Conditional single variant tests](#conditional-single-variant-tests)
    - [Gene based tests](#gene-based-tests)
- [Input files and arguments](#input-files)
    - [Score statistic files files (Summary statistics)](#Score-statistics-files)
    - [Imputation quality files](#Imputation-quality-files)
    - [VCF reference files](#vcf-reference-files)
    - [Annotation files](#annotation-files)
    - [Ancestry files](#ancestry-files)
- [Feedback/Contact](#Feedback/Contact)


## Introduction

**rareGWAMA** is a flexible software package for imputation based GWAS meta-analysis.   
It is developed and maintained by [Dajiang Liu's Group](https://dajiangliu.blog/).


------------------------------------------------------
## Citation
Liu DJ*†, Peloso GM*, Zhan X*, Holmen O*, Zawistowski M, Feng S, Nikpay M, Auer PL, Goel A, Zhang H, Peters U, Farrall M, Orho-Melander M, Kooperberg C, McPherson R, Watkins H, Willer CJ, Hveem, K, Melander O, Kathiresan S, Abecasis GR†    
**Meta-analysis of gene-level tests of rare variant association, Nature Genetics, 46, 200–204 (2014)**  
[doi: 10.1038/ng.2852.](https://www.nature.com/articles/ng.2852)



------------------------------------------------------
## Installing the rareGWAMA R package <a name="Installing-the-rareGWAMA-R-package"></a>

The package is hosted on github, which allows installation and update to be very easy. First, make sure you have the `mvtnorm` and `data.table` packages installed:

    install.packages("devtools")

And also, you need the latest version of [seqminer](https://github.com/zhanxw/seqminer):

    library(devtools)
    devtools::install_github("zhanxw/seqminer")

Then you could use:

    install_github("dajiangliu/rareGWAMA")
    
With `library(rareGWAMA)`, your are ready to go!


------------------------------------------------------
## Quick tutorial <a name="quick-tutorial"></a>

### Single variant tests <a name="Single-variant-tests"></a>

1.The very basic test is using:  

`res <- rareGWAMA.single(score.stat.file=study, imp.qual.file	
=imp.qual, tabix.range="1:11000-58000", alternative="two.sided", col.impqual=5, impQual.lb=0, impQualWeight=FALSE, weight="Npq+impQ",gc=FALSE, rmMultiAllelicSite=TRUE);`   

please find more details in the [input and arguments part](#input-files) for the arguments:
> * score.stat.file: The file names of score statistic files, which could be a **vector object**;
> * imp.qual.file: The file names of imputation quality, which could be a **vector object**;
> * tabix.range: The **tabix range**, which must be in quote and provided as a string like this: "1:11000-58000";
> * alternative: The alternative hypothesis. Default is two.sided;
> * col.impqual: The column number for the imputation quality score;
> * impQual.lb: The lower bound for the imputation quality. Variants with imputaiton quality less than impQual.lb will be labelled as missing;
> * impQualWeight: Using imputation quality as weight;
> * rmMultiAllelicSite: Default is TRUE. Multi-allelic sites will be removed from the analyses if set TRUE, and a variable posMulti will be output; The variant site with multiple alleles can be analyzed using rareGWAMA.single.multiAllele function;  

2.The out put should be as follows:  
`head(res$res.formatted))`
```
     POS           REF ALT  AF         STAT       PVALUE     BETA        SD         N        DIRECTION                          EFFECTIVE_N
[1,] "9:100000040" "G" "A"  "2.28e-05" "1.01e+00" "3.16e-01" " 7.29e-01" " 0.72624" "41634"  "XXXXXXXXXXXXXXXXXXXXXXXXXXXX++XX" "7925"
[2,] "9:100000056" "T" "TA" "7.07e-01" "1.62e+00" "2.04e-01" "-9.09e-03" " 0.00715" "47219"  "XXXXXXXXXXXXXXXX-XXXXXXXXXXX--XX" "46951"
[3,] "9:100000102" "G" "C"  "6.91e-04" "8.00e-01" "3.71e-01" "-7.78e-02" " 0.08691" "95869"  "+++---+-------++XX-++--X++----XX" "36378"
[4,] "9:100000172" "C" "T"  "2.26e-05" "2.25e-01" "6.35e-01" "-1.99e-01" " 0.41982" "125656" "XXXXXXXXXXXXXXXX-XXXXXXXXXXX--X+" "20434"
[5,] "9:100000177" "C" "G"  "2.16e-04" "8.38e-01" "3.60e-01" "-1.04e-01" " 0.11335" "179891" "-+---X+XX-+--+-+-X-+-+-X++-+--X+" "89179"
[6,] "9:100000187" "A" "T"  "1.13e-04" "6.78e-01" "4.10e-01" "-2.35e-01" " 0.28514" "54235"  "+-++-XX+--++--+-XX+++++X----XXXX" "21285"
```
3. For demo data, please see `?rareGWAMA.single`.

### Conditional single variant tests <a name="conditional-single-variant-tests"></a>
1.The command should be like:  

`res <- rareGWAMA.cond.single(score.stat.file=study, imp.qual.file=imp.qual, vcf.ref.file="{$your_path}/ALL.chr9.phase3.genotypes.vcf.gz", candidateVar="9:97018619", knownVar="9:100000172", alternative="two.sided", col.impqual=5, impQual.lb=0, impQualWeight=FALSE, weight="Npq+impQ", gc=FALSE, rmMultiAllelicSite=TRUE);`      
  
please find more details in the [input and arguments part](#input-files) for the arguments:
> * score.stat.file: The file names of score statistic files, which could be a **vector object**;
> * imp.qual.file: The file names of imputation quality, which could be a **vector object**;
> * vcf.ref.file: the file names of the reference panel file, e.g. could be downloaded from [1000 Genomes Project](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/). 
> * candidateVar: The tabix range;
> * knownVar: The known variant;
> * alternative: The alternative hypothesis. Default is two.sided;
> * col.impqual: The column number for the imputation quality score;
> * impQual.lb: The lower bound for the imputation quality. Variants with imputaiton quality less than impQual.lb will be labelled as missing;
> * impQualWeight: Using imputation quality as weight;
> * rmMultiAllelicSite: Default is TRUE. Multi-allelic sites will be removed from the analyses if set TRUE, and a variable posMulti will be output; The variant site with multiple alleles can be analyzed using rareGWAMA.single.multiAllele function;  

2.The out put should be as follows:  
`head(res$res.formatted))` 
```
    POS          REF ALT AF         STAT    PVALUE BETA    SD      N       numStability
[1,] "9:97018619" "T" "C" "8.82e-05" "0.093" "0.76" "0.131" "0.428" "54235" "0"
```
3. For demo data, please see `?rareGWAMA.cond.single`.

### Gene based tests <a name="gene-based-tests"></a>

**For more details, please see the [SHOWCASE](https://funfunchen.github.io/pages/rareGWAMA_page)**

1.The command should be like:  
```
res.gene <- rareGWAMA.gene(score.stat.file, imp.qual.file=imp.qual.file, vcf.ref.file, refFileFormat="vcf.vbi", anno=anno, annoType=c('Nonsynonymous','Stop_Gain',"Essential_Splice_site"), rvtest='VT', ref.ancestry=ref.ancestry, trans.ethnic=TRUE, study.ancestry=study.ancestry, maf.cutoff=0.01);
```

please find more details in the [input and arguments part](#input-files) for the arguments:
> * score.stat.file: The file names of score statistic files, which could be a **vector object**;
> * imp.qual.file: The file names of imputation quality, which could be a **vector object**;
> * vcf.ref.file: The file names of the reference panel file, e.g. could be downloaded from [1000 Genomes Project](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/). 
> * anno: Annotation file or list;
> * annoType: The annotation types you want to use;
> * rvtest: The method you want to use (i.e. 'VT', 'BURDEN', 'SKAT');
> * ref.ancestry: The ancestry information for each sample;
> * study.ancestry: The ancestry information for each study;
> * maf.cutoff: The cutoff for the MAF, could be 0.01, 0.05, 0.001, whatever you want;

2.The out put should be as follows:  
```
GENE    RANGE   STAT    P-VALUE MAF_CUTOFF      NUM_VAR TOTAL_MAF       POS_VAR N       POS_SINGLE_MINP BETA_SINGLE_MINP        SD_SINGLE_MINP
MATN2   8:97931370-98033638     2.19    0.353   0.009472856879857       3       0.0123  8:97961468_G/A,8:98018087_G/A,8:98021213_G/A    204783  8:98021213_G/A  0.01813449913
STK3    8:98767360-98767360     0.0894  0.765   0.00582326555223991     1       0.00582 8:98767360_C/T  235921  8:98767360_C/T  0.00505283659608632     0.0168969355386225
VPS13B  8:99121478-99853608     0.305   0.752   0.00630135839199519     2       0.0117  8:99778930_G/A,8:99820031_A/G   283684  8:99778930_G/A  0.0114882148633614      0.016
COX6C   8:99892015-99892015     0.349   0.555   0.00244114819976761     1       0.00244 8:99892015_G/A  231499  8:99892015_G/A  -0.0159448806371386     0.0269803627402012
```


------------------------------------------------------
## Input files and arguments <a name="input-files"></a> 

:point_right: **All the score statistics files and imputation quality files should be tabix indexed.**      
Say if you have such files: `study1.gz(score statistics files), study2.gz, study1.R2.gz(imputation quality files), study2.R2.gz`in your folder, and already added **tabix** to your `$PATH` environment variable, you could use:   

```for file in study*;do g=`zcat $file | head -200 | grep -n CHROM | cut -f1 -d":"`; tabix -f -S $g -s 1 -b 2 -e 2 $file; done```

to tabix each file.


### Score statistics files:  
If you use <sup>[1](#myfootnote1)</sup> [RVTESTS](https://github.com/zhanxw/rvtests), the output is ready to go:
```
CHROM   POS     REF     ALT     N_INFORMATIVE   AF      INFORMATIVE_ALT_AC      CALL_RATE       HWE_PVALUE      N_REF   N_HET   N_ALT   U_STAT  SQRT_V_STAT     ALT_EFFSIZE     PVALUE
1       10177   A       AC      2352    0.5     2352    1       0       0       2352    0       1.67496 2.51553 0.264695        0.505508
1       10235   T       TA      2352    0       0       1       1       2352    0       0       -0.108472       0.207841        -2.51104        0.601742
1       10352   T       TA      2352    0.5     2352    1       0       0       2352    0       0.665562        2.61389 0.0974122       0.799013
1       10539   C       A       2352    0       0       1       1       2352    0       0       -0.00020902     0.0626305       -0.0532862      0.997337
```

### Imputation quality files:
An example file has the following format:  
```
CHROM   POS     REF     ALT     Rsq
1       10177   A       AC      0.00581
1       10235   T       TA      0.00396
1       10352   T       TA      0.00608
1       10539   C       A       0.00154
1       10616   CCGCCGTTGCAAAGGCGCGCCG  C       0.02085
1       10642   G       A       0.00013reference allele
```
There are five columns, which are `Chromosome #`, `Position`, `Reference allele`, `Alternative allele` and `R square quality`. The higher of `R square quality`(the 5th column) value, the better of the genotype imputation quality.


### VCF reference files: <a name="vcf-reference-files"></a>
The file names of the reference panel file.  
You could download the files from **1000 Genomes Project**: <ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502>.   
Also, you could subset the ranges you want by using **tabix**, with command looks like:   

`tabix -h ALL.2of4intersection.20100804.genotypes.vcf.gz 2:39967768-39967768`   

For more details, please see: [How do I get a sub-section of a VCF file?](http://www.internationalgenome.org/faq/how-do-i-get-sub-section-vcf-file/);

### Annotation files: <a name="annotation-files"></a>
```
       chrom      pos ref alt         af          anno    gene
112207     1 32247432   G   A 0.03425920 Nonsynonymous FAM167B
112208     1 32247598   A   G 0.07993410          Exon FAM167B
112209     1 32247709   G   T 0.00645166        Intron FAM167B
112210     1 32248405   G   C 0.00290383 Nonsynonymous FAM167B
140157     1 41479019   C   T 0.09314330          Utr3    EDN2
140158     1 41479369   C   T 0.06422470          Utr3    EDN2
```

### Ancestry files: <a name="ancestry-files"></a>
1.ref.ancestry. 

You should have a original file (i.e. `ref.ancestry.ori`) as:
```
head(ref.ancestry.ori)
        fid ancestry study
1 samp1      HIM  HCHS
2 samp2      HIM  HCHS
3 samp3      HIM  HCHS
4 samp4      HIM  HCHS
5 samp5      HIM  HCHS
6 samp6      HIM  HCHS
```
Then, you use:
`ref.ancestry <- cbind(ref.ancestry.ori[,1], paste(ref.ancestry.ori[,2], ref.ancestry.ori[,3], sep=","))`

So, the final format should be a *matrix*:
```
     [,1]        [,2]
[1,] "samp1" "HIM,HCHS"
[2,] "samp2" "HIM,HCHS"
[3,] "samp3" "HIM,HCHS"
[4,] "samp4" "HIM,HCHS"
[5,] "samp5" "HIM,HCHS"
[6,] "samp6" "HIM,HCHS"
```

2.study.ancestry.  

Just a vector, as
```
[1] "AACAC" "AMISH" "ARIC"  "BAGS"  "CFS"   "CHS"
```


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
