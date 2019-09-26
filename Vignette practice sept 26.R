# bed file source: https://github.com/mebbert/Dark_and_Camouflaged_genes
library("GenomicDistributions", lib.loc="/Library/Frameworks/R.framework/Versions/3.6/Resources/library")

### Initialization
#BiocManager::install("rtracklayer")
library(rtracklayer)
bed_file <- import("illuminaRL100.hg38.camo.realign.sorted.bed", format="bed")
View(bed_file)
typeof(bed_file)
# This is type "S4" which is the same as the file that Nathan provided which is good
# using HG38 instead of HG19 as reference

### Chromosome distribution plots
x = calcChromBinsRef(bed_file, "hg38")
plotChromBins(x)

# showing two query sets at the same time
query2 = GenomicRanges::shift(bed_file, 1e6)
queryList = GRangesList(original=bed_file, shifted=query2)
x2 = calcChromBinsRef(queryList, "hg38")
plotChromBins(x2)

### Feature distance distribution plots
TSSdist = calcFeatureDistRefTSS(bed_file, "hg38")
plotFeatureDist(TSSdist, featureName="TSS")
# interesting- more reads on trans start site compared to Nathan's example

# plot two dictribution relative to features plots side by side
TSSdist2 = calcFeatureDistRefTSS(queryList, "hg38")
plotFeatureDist(TSSdist2)

# note- can check distances to features other than TSS as well
# to fabricate a feature to check dist from, shift our file by a normally distributed random number
featureExample = GenomicRanges::shift(bed_file, round(rnorm(length(bed_file), 0,1000)))
# use calcFeatureDist instead of calcFeatureDistRefTSS
fdd = calcFeatureDist(bed_file, featureExample)
plotFeatureDist(fdd)

### Partition distribution plots
gp = calcPartitionsRef(bed_file, "hg38")
plotPartitions(gp)

# 2 at a time
gp2 = calcPartitionsRef(queryList, "hg38")
plotPartitions(gp2)
