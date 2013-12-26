######
### Given Barcodes files, generate paired barcode lookup file
######

################################
barcodes <- read.table(file.path('barcodes.txt'), header=T, as.is=T)
P7idx <- barcodes$Adapter == 'P7'
P5_barcode <- barcodes$Seq[!P7idx]
P7_barcode <- barcodes$Seq[P7idx]
## Code to generate lookupTable

barcodeIdentifier = c("Alpha", "Bravo","Charlie", "Delta", "Echo", "Foxtrot", "Golf", "Hotel", "India", "Juliett", "Kilo", "Lima","Mike","November", "Oscar", "Papa", "Quebec", "Romeo")

lookupTable <- data.frame(Barcode_name = paste(rep(barcodeIdentifier, each=48), 1:(length(barcodeIdentifier)*48), sep=""),
#                          P5=rep(paste("P5-Leela", 1:24, "CS1", sep="-"), length(barcodeIdentifier)*2),
                          P5_barcode=rep(P5_barcode, length(barcodeIdentifier)*2),
#                          P7=rep(paste("P7-Leela", length(barcodeIdentifier)*2, "CS2", sep="-"), each=24),
                          P7_barcode=rep(P7_barcode[1:(length(barcodeIdentifier)*2)], each=24))
rownames(lookupTable) <- lookupTable$Barcode
################################

write.table(lookupTable,"barcodeLookupTable.txt",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(lookupTable,"../python/barcodeLookupTable.txt",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
