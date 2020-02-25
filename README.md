# OUTFLANK

OUTFLANK is a nice software to detect outliers.

## 1. Prepare your dataset
To use OUTFLANK, you can convert your vcf file to 012 file using vcftools. 
```{r, engine = 'bash', eval = FALSE}
vcftools --vcf yourvcffile.vcf --012
```

## 2. Install OUTFLANK

In R download devtools and qvalue libraries required for OUTFLANK.
```{r}
install.packages("devtools")
library(devtools)
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")
```
Then, install OUTFLANK from Github and load the librairy.

```{r}
install_github("whitlock/OutFLANK")
library(OutFLANK)
```

## 3. Import your dataset for OUTFLANK analysis
Import your 012 file to R workspace.
```{r}
SNPdata <- read.table("70indsbg.012", header=FALSE)
```

Remove the first column fo this file.
```{r}
SNPdata2 <- SNPdata[,-1]
```

Your dataset, here called SNPdata, should be a matrix with columns = loci, and rows = inds.
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
locusname <- as.character(1:2442)
locusname
```

Calculate Fst per locus.
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

