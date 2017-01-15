######
### Given Barcodes files, generate paired barcode lookup file
######

################################
barcodes <- read.table(file.path('barcodes.txt'), header=T, as.is=T)
I1idx <- barcodes$Adapter == 'I1'
I2_barcode <- barcodes$Seq[!I1idx]
I1_barcode <- barcodes$Seq[I1idx]
## Code to generate lookupTable

barcodeIdentifier = c("Alpha", "Bravo","Charlie", "Delta", "Echo", "Foxtrot", "Golf", "Hotel", "India", "Juliett", "Kilo", "Lima","Mike","November", "Oscar", "Papa", "Quebec", "Romeo")

lookupTable <- data.frame(Barcode_name = paste(rep(barcodeIdentifier, each=48), 1:(length(barcodeIdentifier)*48), sep=""),
#                          I2=rep(paste("I2-Leela", 1:24, "CS1", sep="-"), length(barcodeIdentifier)*2),
                          I2_barcode=rep(I2_barcode, length(barcodeIdentifier)*2),
#                          I1=rep(paste("I1-Leela", length(barcodeIdentifier)*2, "CS2", sep="-"), each=24),
                          I1_barcode=rep(I1_barcode[1:(length(barcodeIdentifier)*2)], each=24))
rownames(lookupTable) <- lookupTable$Barcode
################################

write.table(lookupTable,"barcodeLookupTable.txt",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
write.table(lookupTable,"../python/barcodeLookupTable.txt",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
