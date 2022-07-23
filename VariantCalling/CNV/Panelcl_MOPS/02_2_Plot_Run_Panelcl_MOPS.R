#!/usr/bin/env R
library(panelcn.mops)

bed <- "/storage/home/leefall2/mypro/Cancer_Panel_Package/CNV/Panelcn_MOPS/Cancer_Panel.bed"

countWindows <- getWindows(bed,chr=TRUE) 

#bam <- "/storage/home/leefall2/mypro/Cancer_Panel_Package/Bam/Target_001/SNU005.bam" 
#test.data <- countBamListInGRanges(countWindows = countWindows, bam.files = bam, read.width = 150)

Control.list.file<-read.table("Control_File_list.txt",header=FALSE,stringsAsFactors=FALSE)

control.data <- countBamListInGRanges(countWindows = countWindows, bam.files = Control.list.file[,1], read.width = 150)





Case.list.file<-read.table("Case_File_list.txt",header=FALSE,stringsAsFactors=FALSE)



selectedGenes<-c("CDK4","CDK6","EGFR","ERBB2","FGFR1","FGFR2","FGFR3","MET","MYC")
for (i in Case.list.file[,1]){
bam <- i
test.data <- countBamListInGRanges(countWindows = countWindows, bam.files = bam, read.width = 150)
XandCB <- test.data
elementMetadata(XandCB) <- cbind(elementMetadata(XandCB), elementMetadata(control.data)) 

resultlist <- runPanelcnMops(XandCB, countWindows = countWindows, selectedGenes = selectedGenes)
sampleNames <- colnames(elementMetadata(test.data)) 
pngNames<-strsplit(sampleNames,"\\.")[[1]][1]

resulttable <- createResultTable(resultlist = resultlist, XandCB = XandCB, countWindows = countWindows, selectedGenes = selectedGenes, sampleNames = sampleNames) 

#resulttable[[1]][which((resulttable[[1]]$CN != "CN2"| is.na(resulttable[[1]]$CN))&resulttable[[1]]$lowQual!="lowQual"),]

#write.table(resulttable[[1]],paste("./ResultTable/",pngNames,".txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)
write.table(resulttable[[1]][which((resulttable[[1]]$CN != "CN2"| is.na(resulttable[[1]]$CN))&resulttable[[1]]$lowQual!="lowQual"),],paste("./ResultTable/",pngNames,".txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)

#png(paste("./All_Plot/",pngNames,".png",sep=""),width = 1920, height = 1280, unit='px')
#par(mfrow=c(3,3))

for (i in 1:length(selectedGenes)){
png(paste("./All_Plot/",pngNames,"_",selectedGenes[i],".png",sep=""),width = 1920, height = 1280, unit='px')
plotBoxplot(result = resultlist[[1]], sampleName = sampleNames[1], countWindows = countWindows, selectedGenes = selectedGenes, showGene = i)
dev.off()
}
}
#dev.off()

#}

#resulttable[[1]][which((resulttable[[1]]$CN != "CN2"| is.na(resulttable[[1]]$CN))&resulttable[[1]]$lowQual!="lowQual"),]
