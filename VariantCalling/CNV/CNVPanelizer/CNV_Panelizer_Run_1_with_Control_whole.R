library(CNVPanelizer)
library(ggplot2)
bedFilepath<-"/storage/home/leefall2/mypro/Cancer_Panel_Package/CNVPanelizer.bed"
amplColumnNumber<-8

genomicRangesFromBed <- BedToGenomicRanges(bedFilepath,
                                           ampliconColumn = amplColumnNumber,
                                           split = "_")


genomicRangesFromBed
metadataFromGenomicRanges <- elementMetadata(genomicRangesFromBed)
geneNames = metadataFromGenomicRanges["geneNames"][, 1]
ampliconNames = metadataFromGenomicRanges["ampliconNames"][, 1]

sampleDirectory<-"/storage/home/leefall2/mypro/Cancer_Panel_Package/Bam/Target_001"
referenceDirectory<-"/mnt/QNAP/leefall2/Package/Cancer_core_Control/"

sampleFilenames <- list.files(path = sampleDirectory,
                              pattern = ".bam$",
                              full.names = TRUE)
# Vector with reference filenames
referenceFilenames <- list.files(path = referenceDirectory,
                                 pattern = ".bam$",
                                 full.names = TRUE)


removePcrDuplicates <- TRUE # TRUE is only recommended for Ion Torrent data



referenceReadCounts <- ReadCountsFromBam(referenceFilenames,
                                         genomicRangesFromBed,
                                         sampleNames = referenceFilenames,
                                         ampliconNames = ampliconNames,
                                         removeDup = removePcrDuplicates)
# Read the sample of interest data set
sampleReadCounts <- ReadCountsFromBam(sampleFilenames,
                                      genomicRangesFromBed,
                                      sampleNames = sampleFilenames,
                                      ampliconNames = ampliconNames,
                                      removeDup = removePcrDuplicates)


head(sampleReadCounts)

normalizedReadCounts <- CombinedNormalizedCounts(sampleReadCounts,
                                                 referenceReadCounts,
                                                 ampliconNames = ampliconNames)
# After normalization data sets need to be splitted again to perform bootstrap
samplesNormalizedReadCounts = normalizedReadCounts["samples"][[1]]
referenceNormalizedReadCounts = normalizedReadCounts["reference"][[1]]


replicates <- 10000
# Perform the bootstrap based analysis
bootList <- BootList(geneNames,
                     samplesNormalizedReadCounts,
                     referenceNormalizedReadCounts,
                     replicates = replicates)

# Estimate the background noise left after normalization
backgroundNoise <- Background(geneNames,
                              samplesNormalizedReadCounts,
                              referenceNormalizedReadCounts,
                              bootList,
                              replicates = replicates,
                              significanceLevel = 0.1,
                              robust = TRUE)


# Build report tables
reportTables <- ReportTables(geneNames,
                             samplesNormalizedReadCounts,
                             referenceNormalizedReadCounts,
                             bootList,
                             backgroundNoise)

#write.table(reportTables,"/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/CNVPanelizer/Batch1_with_Whole_Control_except_Case.txt",quote=FALSE,sep="\t")

aa<-PlotBootstrapDistributions(bootList, reportTables)
out.dir<-"/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/CNVPanelizer/batch1"
dir.create(out.dir, showWarnings = FALSE, recursive=TRUE)

for (i in 1:length(aa)){
pngNames<-strsplit(aa[[i]]$labels$title,"\\.")[[1]][1]
print(aa[[i]])
ggsave(paste(out.dir,"/",pngNames,".png",sep=""),width = 20, height = 20, unit="cm",type="cairo")
}


