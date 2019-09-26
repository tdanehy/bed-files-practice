# source: https://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html#introduction

# package is foundation for representing genomic locations within the Bioconductor project
# builds upon IRanges, many other packages depend on it

library(GenomicRanges)
#GRanges class represents a collection of genomic ranges that each have a single start and end location on the genome.

# constructing a GRanges object
gr <- GRanges(
  seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
  ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
  strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
  score = 1:10,
  GC = seq(1, 0, length=10))
gr
# left and right hand are sep by |
# genomic coordinates (seqnames, ranges, and strand) on right
# metadata on right (score, GC, anything really)

seqnames(gr)
ranges(gr)
strand(gr)

# accessor function to extract without corresponsing metadata
granges(gr)

# Annotations for these coordinates
mcols(gr)
mcols(gr)$score

# info about the length of the sequences the ranges are aligned to
seqlengths(gr) <- c(249250621, 243199373, 198022430)
seqlengths(gr)

# accessing length and name
names(gr)
length(gr)

## Splitting and combining GRanges objects
# creates a GrangesList object
sp <- split(gr, rep(1:2, each=5))
sp
# can put back together with c and append
c(sp[[1]], sp[[2]])

## Subsetting GRanges objects
# act like subsetting a vector
gr[2:3]
# this specifies the two rows and one column to be shown
gr[2:3, "GC"]

# this replaces the second row of the gr with the values of the first row
singles <- split(gr, names(gr))
grMod <- gr
grMod[2] <- singles[[1]]
head(grMod, n=3)

#repeat, reverse, or select specific portions of GRanges objects
rep(singles[[2]], times = 3)
rev(gr)
head(gr,n=2)
tail(gr,n=2)

## Basic interval operations
g <- gr[1:3]
# appending one row onto row of 3 created above
g <- append(g, singles[[10]])
start(g)
end(g)
width(g)
range(g)
# methods for manipulating ranges: intra-range methods, inter-range methods, and between-range methods
# find flanking regions on either side of range (10 means include the 10 bases upstream)
flank(g, 10)
# to include the downstream bases:
flank(g, 10, start=FALSE)

# an intra-range method to move ranges by specified number of base pairs
shift(g, 5) # used this in GenomicDistributions Vignette!
# resize extends the ranges by a specific width
resize(g, 30)

#Inter-range methods- compares the ranges in a single GR obj
# reduce method will align the ranges and merge overlapping ranges
reduce(g)
## 102-112 and 103-113 are merged into 102-113
# find gaps with gaps method
gaps(g)
# disjoin is a collection of non-overlapping ranges
disjoin(g)
# coverage quantifies the degree of overlap
coverage(g)

#Between-range methods: relationships between different GRanges objects
# make new data frame g2 with top 2 lines of original
g2 <- head(gr, n=2)
union(g, g2)
intersect(g, g2)
# calc the asymmetric difference
setdiff(g, g2)

# when objects are parallel to eachother: element 1 of obj 1 is related to element 1 of obj 2
# functions all start with a p
# combine element-wise, req is that the number of elements in each G ranges is same
# also that they have the same seqnames and strand assignments

# simulate this:
g3 <- g[1:2]
ranges(g3[1]) <- IRanges(start=105, end=112)
punion(g2, g3)
pintersect(g2, g3)
psetdiff(g2, g3)

# below can be used to find a comprehensive list of methods
methods(class="GRanges")
