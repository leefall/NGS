#install.packages("ExomeDepth")
output.path<-"/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/All_wise_Test_Inhouse/OriBed/ExomeDepth/Result/"
library(ExomeDepth)
library("GenomicRanges")
data(exons.hg19)
pass = function(){
}
null.df<-data.frame()

#my.counts <- getBamCounts(bed.frame = exons.hg19,
                          #bam.files = my.bam,
                          #include.chr = FALSE,
                          #referenceFasta = fasta)

head(exons.hg19,n=10)
data(ExomeCount)
ExomeCount.dafr <- as(ExomeCount, 'data.frame')


ExomeCount.dafr<-data.frame(read.table("/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/All_wise_Test_Inhouse/OriBed/ExomeDepth/Target001.txt",
                                       header=TRUE,stringsAsFactors = FALSE))

ExomeCount.dafr$chromosome <- gsub(as.character(ExomeCount.dafr$seqnames),
                                   pattern = 'chr',
                                   replacement = '') ##remove the annoying chr letters

print(head(ExomeCount.dafr))
data(Conrad.hg19)
exons.hg19.GRanges <- GRanges(seqnames = exons.hg19$chromosome,
                              IRanges(start=exons.hg19$start,end=exons.hg19$end),
                              names = exons.hg19$name)

ExomeCount.mat <- as.matrix(ExomeCount.dafr[, grep(names(ExomeCount.dafr), pattern = 'SNU*')])
nsamples <- ncol(ExomeCount.mat)

#i<-1

for (i in 1:nsamples) {
my.choice <- select.reference.set (test.counts = ExomeCount.mat[,i],
                                   reference.counts = ExomeCount.mat[,-i],
                                   bin.length = (ExomeCount.dafr$end - ExomeCount.dafr$start)/1000,
                                   n.bins.reduced = 10000)

my.reference.selected <- apply(X = ExomeCount.mat[, my.choice$reference.choice, drop = FALSE],
                               MAR = 1,
                               FUN = sum)
message('Now creating the ExomeDepth object')
all.exons <- new('ExomeDepth',
                 test = ExomeCount.mat[,i],
                 reference = my.reference.selected,
                 formula = 'cbind(test, reference) ~ 1')
################ Now call the CNVs
all.exons <- CallCNVs(x = all.exons,
                      transition.probability = 10^-4,
                      chromosome = ExomeCount.dafr$chromosome,
                      start = ExomeCount.dafr$start,
                      end = ExomeCount.dafr$end,
                      name = ExomeCount.dafr$names)
########################### Now annotate the ExomeDepth object
if (identical(all.exons@CNV.calls,null.df)){
  pass 
} else {
all.exons <- AnnotateExtra(x = all.exons,
                           reference.annotation = exons.hg19.GRanges,
                           min.overlap = 0.0001,
                           column.name = 'exons.hg19')

#all.exons@CNV.calls

output.file <- paste(output.path,'/', colnames(ExomeCount.mat)[i], '.txt', sep = '')


write.table(file = output.file, x = all.exons@CNV.calls, row.names = FALSE,quote=FALSE,sep="\t")
}

}








