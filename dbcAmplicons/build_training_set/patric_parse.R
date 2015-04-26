
# location of sync db
library(parallel)
ncores = 60
folder <- "/data/patric/genomes"

patricGenomes <- dir(path=folder)
patricGenomes <- patricGenomes[grep("_",patricGenomes)]

geneText <- "Small Subunit"


# ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
names <- read.table("../taxonomy_ncbi/names.dmp",sep="\t",comment.char="",quote="",as.is=TRUE)
names <- names[,seq(1,ncol(names),by=2)]
names <- names[which(names$V7 == "scientific name"),]

nodes <- read.table("../taxonomy_ncbi/nodes.dmp",sep="\t",comment.char="",quote="",as.is=TRUE)
nodes <- nodes[,seq(1,ncol(nodes),by=2)]


getFullPath <- function(genome){
  genus <- sapply(strsplit(genome,split="_"),"[[", 1L)
  species <- sapply(strsplit(genome,split="_"),"[[", 2L)
  getIsolate <- function(array){
    isolate <- sub(paste(array[2],array[3],sep="_"),"",array[1])
    isolate <- sub("^_","",isolate)
    isolate <- gsub("_"," ",isolate)
    isolate
  }
  isolate <- apply(cbind(genome,genus,species),1,getIsolate)
  if(species == "sp") species = "sp."
  
  accept_path <- c("superkingdom","phylum","class","order","family","genus","species")
  getLineage <- function(x) {
    #    cat(paste(names[match(x,names$V1),"V3"],nodes[match(x,nodes$V1),"V1"],nodes[match(x,nodes$V1),"V5"],"\n",sep=" "))
    if(names[match(x,names$V1),"V3"] != "Bacteria"){
      if (nodes[match(x,nodes$V1),"V5"] %in% accept_path){
        return(c(getLineage(nodes[match(x,nodes$V1),"V3"]),as.numeric(nodes[match(x,nodes$V1),"V1"])))
      }
      else{
        return(getLineage(nodes[match(x,nodes$V1),"V3"]))
      }
    }else{
      return(as.numeric(nodes[match(x,nodes$V1),"V1"]))
    }
  }
  taxid <- names[match(paste(genus,species,isolate),names$V3),"V1"]
  if(is.na(taxid)) {
    taxid <- names[match(paste(genus,species,sub(" ","_",isolate)),names$V3),"V1"]
    if (!is.na(taxid)) isolate=sub(" ","_",isolate)
  }
  if(is.na(taxid)) taxid <- names[match(paste(genus,species),names$V3),"V1"]
  if(is.na(taxid)) taxid <- names[match(paste(genus),names$V3),"V1"]
  if(is.na(taxid)) return(list("taxid"=NA,"lineage"=paste(rep(NA,8),collapse="_:_")))
  line <- getLineage(taxid)
  if (isolate == ""){
    isolate = NA
  }else{
    isolate = sub(" +$","",paste(genus,species,isolate))
  }
  return(list("taxid"=taxid,"lineage"=paste(c(names$V3[match(line[match(accept_path,nodes$V5[match(line,nodes$V1)])],names$V1)],isolate),collapse="_:_")))
}
## end taxonomy ##

extractPatricGene <- function(genome, gene="Small Subunit", pathToGenome="/data/patric/genomes"){
  # test
  # genome = patricGenomes[1]
  results = list()
  results$genome_name = genome
  results$taxonomy = getFullPath(genome)
  
  featureFile <- file.path(pathToGenome,genome,paste(genome,"PATRIC.features.tab",sep="."))
  sequenceFile <- file.path(pathToGenome,genome,paste(genome,"fna",sep="."))
  
  if(!(file.exists(featureFile)) | !(file.exists(sequenceFile))){
    results$status ="File Error"
    return(results)
  }
  
  if(file.info(sequenceFile)$size > 0){
    sequence <- try(readDNAStringSet(sequenceFile))
    if (inherits(sequence, "try-error")){ results$status = "Sequence Read Error"; return(results); }
  } else { results$status = "Sequence Size Error"; return(results); }
  results$genome_size <- sum(width(sequence))
  results$num_contigs <- length(sequence)
  
  if(file.info(featureFile)$size > 0){
    annot <- read.table(featureFile,sep="\t",header=TRUE,as.is=TRUE,quote="",comment.char="")
    if (inherits(annot, "try-error")) { results$status = "Feature Read Error"; return(results); }
  } else { results$status = "Feature Size Error"; return(results); }
  
  geneannot <- annot[intersect(grep(gene, annot$product),which(annot$feature_type=="rRNA")),]
  if (nrow(geneannot) == 0) { results$status = "No Gene Present"; return(results); }
  
  yankSequence <- function(genome,seqs,anno) {
    sizeOk <- apply(anno,1,function(x){
      wid = width(seqs[grep(paste0(x["accession"]," "), names(seqs))]) >= as.numeric(x["end_min"])
      if (length(wid) == 0) return(FALSE)
      return(wid)
    })
    if(sum(sizeOk) != nrow(anno))
      warning(paste(genome,anno$genome_name[!sizeOk],paste(anno$accession[!sizeOk],":",anno$start_max[!sizeOk],"-",anno$end_min[!sizeOk],sep=""),anno$product[!sizeOk],": does not exist or is past the end of the contig\n"))
    if (sum(sizeOk)){ ## If any are left
      tmp <- apply(anno[sizeOk,],1,function(x){
        seq <- subseq(x=seqs[grep(paste0(x["accession"]," "),names(seqs))],
                      start=as.numeric(x["start_max"]),
                      end=as.numeric(x["end_min"])); 
        if(as.vector(x["strand"] == "-")){
          seq <- reverseComplement(seq)};
        seq 
      })
      tmp <- DNAStringSet(sapply(tmp,as.character))
      names(tmp) <- paste(paste(anno$accession[sizeOk],":",anno$start_max[sizeOk],"-",anno$end_min[sizeOk],sep=""),genome,anno$product[sizeOk],sep="|")
      return(tmp)
    } else return(c())
  }
  results$seqs <- yankSequence(genome,sequence,geneannot)
  if (length(results$seqs) == 0) { results$status = "Gene Extraction Error"; return(results); }
  results$num_gene = length(results$seqs)
  results$status = "OK"
  
  return(results)
}

require("Biostrings")
geneResults <- mclapply(patricGenomes, function(x) {print(x);extractPatricGene(x,gene=geneText,pathToGenome=folder)},mc.cores = ncores)
# Failed
table(sapply(geneResults,inherits,"try-error"))
while(1){
  geneResults2 <- mclapply(patricGenomes[which(sapply(geneResults,inherits,"try-error"))], function(x) {print(x);extractPatricGene(x,gene=geneText,pathToGenome=folder)},mc.cores = ncores)
  geneResults[which(sapply(geneResults,inherits,"try-error"))] <- geneResults2
  table(sapply(geneResults,inherits,"try-error"))
  if(!any(sapply(geneResults2,inherits,"try-error")) | all(sapply(geneResults2,inherits,"try-error"))) break()
}

names(geneResults) <- patricGenomes

#final
failed <- which(sapply(geneResults,inherits,"try-error"))
geneResults <- geneResults[-failed]

geneStatus <- sapply(geneResults,"[[","status")
notOK <- which(geneStatus != "OK")

badGENOMES <- cbind(c(names(failed),names(notOK)),c(rep("Failed",length(failed)),geneStatus[notOK]))

write.table(badGENOMES,"failed_genomes.txt",sep="\t",row.names=F,col.names=F)

geneResults <- geneResults[-notOK]
seqs <- sapply(geneResults,function(x) x$seqs)
names(seqs) <- NULL
seqs <- do.call("c",seqs)

writeXStringSet(seqs,"TestSeqs.fasta")
seqs <- readDNAStringSet("TestSeqs.fasta")


cross_match_minmatch <- 8
cross_match_minscore <- 12
cross_match_screenfile <- "../primer2.fasta"

cross_match_call <- function(x) paste("/mnt/home/msettles/opt/bin/cross_match ",x," ",cross_match_screenfile, " -minmatch ",cross_match_minmatch," -minscore ", cross_match_minscore," -tags > ",x,".cmout",sep="")
res <- system(cross_match_call("TestSeqs.fasta"),intern = TRUE)

### read in and format cross_match output
"parse_cm" <- function(filenames)
{
  align <- sapply(filenames,function(x) {
    lines <- readLines(x)
    align <- grep("^ALIGNMENT",lines,value=T)
    align <- sub(" C ", " ",align)
    align
  })
  align<- unlist(align)
  cm_out <- matrix(unlist(strsplit(align,split=" +")),ncol=13,byrow=T)
  cm_out <- data.frame(cm_out,stringsAsFactors=FALSE)
  cm_out$FC <- "F"
  cm_out$FC[grep("(",cm_out$X11,fixed=TRUE)] <- "C"
  cm_out[cm_out$FC == "C",c("X11","X12","X13")] <- cm_out[cm_out$FC == "C",c("X13","X12","X11")]
  colnames(cm_out) <- c("Alignment","score","perc_sub","perc_del","perc_ins","read_id",
                        "read_start","read_end","read_remain","adapt","adapt_start","adapt_end","adapt_remain","FC")   
  cm_out$read_start <- as.numeric(cm_out$read_start)
  cm_out$read_end <- as.numeric(cm_out$read_end)
  cm_out$read_remain <- as.numeric(gsub("[()]","",cm_out$read_remain))
  cm_out$adapt_remain <- as.numeric(gsub("[()]","",cm_out$adapt_remain))
  cm_out$adapt_start <- as.numeric(cm_out$adapt_start)
  cm_out$adapt_end <- as.numeric(cm_out$adapt_end)
  cm_out$score <- as.numeric(cm_out$score)
  cm_out$perc_sub <- as.numeric(cm_out$perc_sub)
  cm_out$perc_del <- as.numeric(cm_out$perc_del)
  cm_out$perc_ins <- as.numeric(cm_out$perc_ins)
  
  cm_out$perc_sub <- round(cm_out$perc_sub *(cm_out$read_end - cm_out$read_start +1),digits=0)/100
  cm_out$perc_del <- round(cm_out$perc_del *(cm_out$read_end - cm_out$read_start +1),digits=0)/100
  cm_out$perc_ins <- round(cm_out$perc_ins *(cm_out$read_end - cm_out$read_start +1),digits=0)/100
  
  cm_out$read_len <- cm_out$read_end+cm_out$read_remain
  
  cm_out$err <- apply(cm_out[,c("perc_sub","perc_del","perc_ins","adapt_remain","adapt_start")],1,sum) -1
  cm_out
}

cmout <- parse_cm("TestSeqs.fasta.cmout")
cmout$primer <- sapply(strsplit(cmout$adapt,split="_"),"[[",1L)

tb_fa <- table(cmout$read_id,cmout$primer)

ord <- match(sapply(strsplit(names(seqs),split=" "),"[[",1L),rownames(tb_fa))
tb_fa <- as.data.frame(tb_fa[ord,])
seqs <- seqs[!is.na(ord)]
tb_fa <- tb_fa[!is.na(ord),]
tb_fa$width = width(seqs)
tb_fa$primers <- apply(tb_fa[,1:7],1,function(x) sum(x>0))

## only keep those with V1-V3 primer matches
v1v3 <- seqs[match( sapply(strsplit(rownames(tb_fa),split="\\|")[tb_fa$'8F' & tb_fa$'519F'],"[[",1L),sapply(strsplit(names(seqs),"\\|"),"[[",1L))]

writeXStringSet(unique(v1v3),file="AmpSeqs.V1.V3.fasta")