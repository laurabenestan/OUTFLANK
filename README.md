# OUTFLANK

OUTFLANK is an R package that implements the method developed by Whitlock and Lotterhos (2015) to use likelihood on a trimmed distribution of FST values to infer the distribution of FST for neutral markers. See the [vignette](https://htmlpreview.github.io/?https://github.com/whitlock/OutFLANK/blob/master/inst/doc/OutFLANKAnalysis.html) for further explanations and escription of the package utilities.

## 1. Prepare your dataset in the right format

To be able to analyse your file with OUTFLANK, use the option **-012** [VCFTOOLS](http://vcftools.sourceforge.net) in the terminal. 

This option outputs three files:
1- a large matrix with suffix ".012", which contains the genotypes of each individual on a separate line. Genotypes are represented as 0, 1 and 2, where the number represent that number of non-reference alleles. Missing genotypes are represented by -1. 
2 - a file with suffix ".012.indv" details the individuals included in the main file. 
3 - a file, with suffix ".012.pos" details the site locations included in the main file.

```{r, engine = 'bash', eval = FALSE}
vcftools --vcf yourvcffile.vcf --012
```

## 2. Install OUTFLANK

In R download devtools and qvalue libraries required for OUTFLANK.
```{r}
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("qvalue")
```
Then, install OUTFLANK from Github and load the librairy.

```{r}
install.packages("devtools")
library(devtools)
install_github("whitlock/OutFLANK")
library(OutFLANK)
```

## 3. Import your dataset for OUTFLANK analysis
Import your .012 file to R workspace.
```{r}
SNPdata <- read.table("70indsbg.012", header=FALSE)
```

Remove the first column fo this file.
```{r}
SNPdata2 <- SNPdata[,-1]
```

Your dataset, here called SNPdata, should be a matrix with columns = loci, and rows = inds. 
We have to transpose the genotype matrix (SNPdata) to get it into OutFLANK format.
```{r}
SNPdata <- t(SNPdata)
head(SNPdata)
dim(SNPdata)
```

Load the individual population data. 
```{r}
inds <- read.table("70indsbg.012.indv", header=FALSE)
#inds2 <-inds[-1,]
head(inds)
dim(inds)
```

Add locus names.
```{r}
locusname <- as.character(read.table("70indsbg.012.pos", header=FALSE)
locusname
```

First, we calculate FST on all the loci in our dataset. Calculating FSTs, may take a few minutes...
```{r}
FstDataFrame <- MakeDiploidFSTMat(SNPdata2,locusname,inds)
head(FstDataFrame)
```

Check for missing data and get a vector indicating the amount of missing data for each locus.
```{r}
missing_data <- apply(SNPdata,2, function(i){sum(i==9)})
hist(missing_data, breaks=seq(0,100,1))
```

Visualize the Fst distribution across the SNPs
```{r}
plot(FstDataFrame$FST, FstDataFrame$FSTNoCorr, 
     xlim=c(-0.01,0.3), ylim=c(-0.01,0.3),
     pch=20)
abline(0,1)
```

