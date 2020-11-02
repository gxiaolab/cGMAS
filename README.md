# cGMAS

cGMAS (Concordance-based GMAS) is a method for predicting functional SNPs for GMAS (genetically-modulated alternative splicing) events. GMAS events have associated SNPs that can serve as tag SNPs. If a SNP is functional for GMAS, we expect to see a concordant SNP genotype and alternative splicing pattern accross a large number of individuals. We can quanitify the concordance between genotype and splicing pattern using a concordance score (Si). This method cannot distinguish true functional and neutal SNPs if they are in perfect LD. Also, the method requires a large number of individuals.

## Table of contents

A. [Preprocessing](#preprocessing)
- Identify GMAS events (GMAS SNPs and GMAS exons)
- Determine genotype of candidate functional SNPs

B. [Calculating concordance score (Si)](#calc-Si-score)

C. [Filtering for candidate functional SNPs for GMAS](#filtering-candidate-snps)
- Identify Si score peak
- Filter by mean Si scores and Si peak magnitude
- FDR corrections


## Preprocessing

- We can use the [ASARP pipeline](https://legacy.ibp.ucla.edu/research/xiao/Software_files/ASARP/doc/) to get high-quality tag SNPs and GMAS events as input for this pipeline.

## Calculating concordance score (Si)


## Filtering for candidate functional SNPs for GMAS
