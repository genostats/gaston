# gaston 1.1

## Corrections

* corrected C++ code (compilation errors on some distributions)

# gaston 1.2

## Improvements

* kinship matrix is now in single precision (faster) 

* various minor changes

# gaston 1.3

## Bug fixes

* corrected bug in HWE chi square (monomorphic SNP) 

## New features

* select.snps and select.inds (negative indices, presence of NAs)

# gaston 1.3.1

## Bug fixes

* fixed an (intermittent) bug (?!) in kinhsip matrix computations

# gaston 1.4

## New features

* introduced option gastion.auto.set.stats / updated vignette accordingly

* added 'hz' in ped stats

* new function reshape.GRM

* as.bed.matrix no longer is a method, but a function

* minor modifs in read.vcf

## Implementation details

* matrix4 now uses `uint8_t` array instead of char array (induced many small changes in the C++ code).

# gaston 1.4.5

## Bug fixes

* corrected a bug in method 'show'

## Implementation details

* better handling of constraints in AIREML methods 

* faster implementation of GRM

## New features

* new function genomic.sex

# gaston 1.4.6

## Implementation details

* Modified Imports/Depends fields to remove a warning about RcppParallel

* Modified VCF handling functions to follow the new interface of WhopGenome (patch by Tomas Kalibera and Ulrich Wittelsbuerger)

## New features

* New functions logistic.mm.aireml, score.variance.linear, score.variance.logistic, score.fixed.linear, score.fixed.logistic

* New argument get.P in lmm.aireml

* New features and arguments in association.test

* Minor modifications of the vignette (illustration of association.test)

# gaston 1.4.7

## Bug fixes

* Corrected (minor) bug in association.test (handling of monomorphic SNPs in some tests)

## New features

* New function LD.plot (LDheatmap dependence removed)

* Handling NAs in select.snp select.inds

* Added chr X, Y, MT support / modified set.stats, set.genomic.sex accordingly

* Added which.snps argument in LD.thin and GRM

* rbind and cbind now check individuals ids / snp ids ; rbind check reference alleles and perform reference inversions / strand flips if needed.

* Integrated code from CompQuadForm / removed dependence 

* Minor modifications of the vignette (update description of modified functions)

# gaston 1.4.8

## Bug fixes

* Corrected bug in set.stats (stats for chr Y)

* Finally removed CompQuadForm from dependence list !

## New features

* Minor modifs in read vcf, rbind, set.genomic.sex

* New function merge.inds (function finally not exported)

* New functions (short cuts) is.autosome, is.chr.x, is.chr.y, is.chr.mt

# gaston 1.4.9

## New features

* Modification of rbind, cbind to handle bed.matrices with different column names in @ped, @snps

* A new test in association.test

* Minor changes in GRM and LD.thin 

# gaston 1.5

## New features

* VCF files support

* Computation of the Dominance Matrix (DM)

* New functions test.snps test.inds which.snps which.inds

* Faster optimization with the diagonalization trick

* Auto set the #threads when loaded

# gaston 1.5.1 

## New features

* Improved functions association.test / random.pm

* cbind : now keep p, mu, sigma

* rbind : now handles alleles longer than 1

* New functions qqplot.pvalue, manhattan

* New functions SNP.duplicated, SNP.match

* New functions lmm.restricted.likelihood, lmm.profile.restricted.likelihood

## Implementation details

* Vignette in knitr

# gaston 1.5.2  

## Implementation details

* Remove some compilation warnings (-Wreorder)

* Temporarily disable TBB (UBSAN test) + function setThreadOptions

# gaston 1.5.3  

## Implementation details

* Cast in log() function (fix compilation error on solaris)

* Modified plot thinning in qqplot.pvalues

# gaston 1.5.4  

## New features

* New functions SNP.rm.duplicates and set.dist

* Argument 'by = ' in SNP.match

## Implementation details

* Re-enabled TBB

* Improved logistic regression code / small improvements in association.testa

# gaston 1.5.5  

## Bug fixes

* Corrected naugthy bug in logistic regression (in case of 20 iterations, still compute covariance matrix!!!) (bug introduced in 1.5.4)

## New features

* Started to write an FAQ

# gaston 1.5.6  

## Bug fixes

* Corrected bug in gwas lmm wald (PQL) (was sending back NaNs only :/ )

# gaston 1.5.7  

## Implementation details

* Corrected return; to return() as requested by the CRAN

* Improved code for snp_hash.h (from gaston.utils)

## New features

* Argument EM in logistic-mm-aireml 

# gaston 1.5.8  

## New features

* new function LD.clump

# gaston 1.5.9  

## New features

* Default argument max.dist = 500e3 in LD.thin

## Bug fixes

* Corrected bug in ld clump, thanks to valgrind and Pr Ripley <3

