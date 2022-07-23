library(CNVPanelizer)
bedFilepath<-"/storage/home/leefall2/mypro/Cancer_Panel_Package/CNVPanelizer.bed"
amplColumnNumber<-8

genomicRangesFromBed <- BedToGenomicRanges(bedFilepath,
                                           ampliconColumn = amplColumnNumber,
                                           split = "_")


genomicRangesFromBed
metadataFromGenomicRanges <- elementMetadata(genomicRangesFromBed)
geneNames = metadataFromGenomicRanges["geneNames"][, 1]
ampliconNames = metadataFromGenomicRanges["ampliconNames"][, 1]

sampleDirectory<-"/mnt/QNAP/leefall2/Package/Cancer_core_Target004"
referenceDirectory<-"/mnt/QNAP/leefall2/Package/Cancer_core_Control_1"

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

PlotBootstrapDistributions(bootList, reportTables)



reportTables


write.table(reportTables,"/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Coverage/Bed_Coverage/Merged_Bed/20Q/CNVPanelizer/Batch4_with_Control1.txt",quote=FALSE,sep="\t")

###


options(width=500)  # to show the entire table..
# to avoid have to print to other page..
numberOfGenesViewport = 20
sampleIndexToShow = 2

indexToHide <- which(colnames(reportTables[[1]]) %in% c("MeanRatio", "Passed"))

hidden.report.tables<-reportTables[[sampleIndexToShow]][1:numberOfGenesViewport, -c(indexToHide)]







reportTables$`IonXpress_005_R_2015_12_04_17_16_58_user_PROTON1-64-20151204_Cancer_Core_Target_004_20151204_Cancer_Core_Target_004_Re.bam`

###Discovery Function


BedToGenomicRanges


