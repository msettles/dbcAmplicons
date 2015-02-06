#!/usr/bin/env Rscript

### Required packages
#source('http://bioconductor.org/biocLite.R')
#biocLite(c('optparse','parallel','ggplot2','RColorBrewer','reshape','plyr','Biostrings','ShortRead'))
###
suppressPackageStartupMessages(library("optparse"))

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to
# make_option(c("-h", "--help"), action="store_true", default=FALSE,
# help="Show this help message and exit")
option_list <- list(
    make_option(c("-p","--program"), type="character", default="consensus",
                help="comma separated list of the functions of consensus,ambiguities,occurrence [default %default]",
                dest="program"),
    make_option(c("-s", "--min-seq"), type="integer", default=5,
                help="minimum number of reads [default %default]",
                dest="min_seq"),
    make_option(c("-f", "--min-freq"), type="double", default=0.05,
                help="minimum frequecy of reads [default %default]",
                dest="min_freq"),
    make_option(c("-1", "--trim-1"), type="integer", default=0,
                help="number of bases to trim from the end of read 1 [default %default]",
                dest="trimOne"),
    make_option(c("-2", "--trim-2"), type="integer", default=0,
                help="number of bases to trim from the end of read 2 [default %default]",
                dest="trimTwo"),   
    make_option(c("-r", "--reuse"), action="store_true", default=FALSE,
                help="Reuse the reads within the folder [default %default]",
                dest="reuse"),
    make_option(c("-o", "--outputDir"), type="character", default=NULL,
                help="directory name to place output files [default [basename].reduced]",
                dest="output"),
    make_option(c("-c", "--cpus"), type="integer", default=0,
                help="number of processors to use, choose 0 to use all available cores [default %default]",
                dest="procs")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
parser <- OptionParser(usage = "%prog [options] basename",option_list=option_list)
arguments <- parse_args(parser, positional_arguments = 1)

#arguments <- list(options = list(program="consensus,ambiguities,occurrence",min_seq=5,min_freq=0.05,trimOne=0,trimTwo=0,reuse=FALSE,output="Diego_Chloroplast_R1-0_R2-0",procs=30),args="Diego_Chloroplast_primer")

opt <- arguments$options
basename <- arguments$args

avail_functions <- c("consensus","ambiguities","occurrence")
program_list = unlist(strsplit(opt$program,split=","))
if (!all(program_list %in% avail_functions)) {
    stop(paste("program list parameter must be 1) comma separated and 2) in the list",paste(avail_functions,collapse=","),sep=" "))
}

min_seq = opt$min_seq
min_freq = opt$min_freq
trimOne = opt$trimOne
trimTwo = opt$trimTwo

output = opt$output
if (is.null(output)) output <- paste(basename,"reduced",sep=".")

if(file.exists(output) & !opt$reuse){
    warning(paste("output directory '",output,"' already exists, deleting directory"),sep=" ")
    unlink(output,recursive=TRUE,force=TRUE)
}
if (!opt$reuse) dir.create(output,recursive=TRUE,mode="0777")

suppressPackageStartupMessages(library("parallel"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("reshape"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("Biostrings"))
suppressPackageStartupMessages(library("ShortRead"))

######################################################################
## prepareCore
##  Set up the numer of processors to use
## 
## Parameters
##  opt_procs: processors given on the option line
##  samples: number of samples
##  targets: number of targets
"prepareCore" <- function(opt_procs){
    # if opt_procs set to 0 then expand to samples by targets
    if( opt_procs == 0 ) opt_procs <- detectCores()
    if (opt_procs > detectCores())
        stop(paste("number of processors specified exceeds the number of available processors '",detectCores,"'",sep=" "))
    write(paste("Using",opt_procs,"processors",sep=" "),stdout())
    return(opt_procs)
}
procs <- prepareCore(opt$procs)

### first read in paired end files
R1 <- paste(basename,"_R1.fastq.gz",sep="")
R2 <- paste(basename,"_R2.fastq.gz",sep="")
if (!file.exists(R1) | !file.exists(R2)){
    R1 <- paste(basename,"_R1.fastq",sep="")
    R2 <- paste(basename,"_R2.fastq",sep="")
    if (!file.exists(R1) | !file.exists(R2)){
        stop(paste("cannot find input files:\n",R1,"\n",R2))
    }
}

if (trimOne != 0 | trimTwo != 0){
    if (!opt$reuse){
        write(paste("Trimming Sequences:\n",trimOne,"bases from Read1\n",trimTwo,"bases from Read2"),stdout())
        write(paste("Reading input files:\n",R1,"\n",R2),stdout())
        fq_r1 <- readFastq(R1)
        fq_r2 <- readFastq(R2)
        write(paste("Read in",length(fq_r1), " paired reads"),stdout())
        
        trimSRQ <- function(sr, trimL, trimR){
            trimR[trimR < 1] = 1
            trimL[trimR<trimL]  = trimR
            ShortReadQ(sread=subseq(sread(sr),start=trimL,end=trimR),
                             quality=new(Class=class(quality(sr)),quality=subseq(quality(quality(sr)),start=trimL,end=trimR)),
                             id=id(sr))
    }
    tryCatch({
        fq_r1 <- trimSRQ(fq_r1,1,width(fq_r1)-trimOne);
        fq_r2 <- trimSRQ(fq_r2,1,width(fq_r2)-trimTwo);
        R1 <- file.path(output,paste(basename,".trimmed_R1.fastq.gz",sep=""))
        R2 <- file.path(output,paste(basename,".trimmed_R2.fastq.gz",sep=""))
        write(paste("Writing trimmed input files:\n",R1,"\n",R2),stdout())
        writeFastq(fq_r1,R1)
        writeFastq(fq_r2,R2)
        },error=function(e) {print("cannot trim reads, make sure bases to trim does not exceed read lengths")})
        } else {
            R1 <- file.path(output,paste(basename,".trimmed_R1.fastq.gz",sep=""))
            R2 <- file.path(output,paste(basename,".trimmed_R2.fastq.gz",sep=""))            
    }
}

if (!opt$reuse){
    write(paste("Joining reads with flash"),stdout())
    flash_prefix = "flash"
    call <- paste("flash --max-overlap=600 --allow-outies -t",procs,"-x 0.25 -z -o", file.path(output,flash_prefix),R1,R2,sep=" " )
    flash_output <- system(call,intern = TRUE)
    
    flash_data = rep(NA,4)
    if(length(oflash <- which(flash_output=="[FLASH] Read combination statistics:")) != 0){
        flash_res <- flash_output[(oflash+1):(oflash+4)]
        flash_data <- as.numeric(sapply(strsplit(flash_res,split=" +|%"),"[[",4L))
        flash_data <- paste("\nReads_Combined:\t\t",flash_data[2],"\nReads_Uncombined:\t",flash_data[3],"\nCombined_Percentage:\t",flash_data[4],"%",sep="")
    }
    write(flash_data,stdout())
}
write(paste("Reading merged read files:\n",file.path(output,"flash.extendedFrags.fastq.gz"),"\n",file.path(output,"flash.notCombined_1.fastq.gz"),"\n",file.path(output,"flash.notCombined_2.fastq.gz")),stdout())

### read in input files      
fq <- readFastq(file.path(output,"flash.extendedFrags.fastq.gz"))
fq_r1 <- readFastq(file.path(output,"flash.notCombined_1.fastq.gz"))
fq_r2 <- readFastq(file.path(output,"flash.notCombined_2.fastq.gz"))

### process merged files
nms <- as.character(id(fq))
id <- sapply(strsplit(sapply(strsplit(nms,split=" "),"[[",2L),split=":"),"[[",4L)
primer <- sapply(strsplit(sapply(strsplit(nms,split=" "),"[[",2L),split=":"),"[[",5L)
occurrence = table(id,primer)

### process unmerged files
nms_p <- as.character(id(fq_r1))
id_p <- sapply(strsplit(sapply(strsplit(nms_p,split=" "),"[[",2L),split=":"),"[[",4L)
primer_p <- sapply(strsplit(sapply(strsplit(nms_p,split=" "),"[[",2L),split=":"),"[[",5L)
occurrence_p = table(id_p,primer_p)

uprimer <- sort(unique(c(primer,primer_p)))
uid <- sort(unique(c(id,id_p)))

## Plot read ratio 
uoccurrence <- matrix(0,nrow=length(uid),ncol=length(uprimer))
uoccurrence_p <- matrix(0,nrow=length(uid),ncol=length(uprimer))
rownames(uoccurrence) <- uid
rownames(uoccurrence_p) <- uid
colnames(uoccurrence) <- uprimer
colnames(uoccurrence_p) <- uprimer

uoccurrence[rownames(occurrence),colnames(occurrence)] <- occurrence
uoccurrence_p[rownames(occurrence_p),colnames(occurrence_p)] <- occurrence_p
ratio <- log10((uoccurrence+1)/(uoccurrence_p+1))
mratio <- melt(ratio)
colnames(mratio) <- c("SampleID","PrimerID","value")

jRdGyFun <- colorRampPalette(brewer.pal(n = 11, "RdGy"))
paletteSize <- 256
jRdGyPalette <- jRdGyFun(paletteSize)

axis.y.res = if(length(unique(mratio$SampleID))> 96) { element_blank() }else{element_text(size=8)}

p <- ggplot(mratio, aes(x = PrimerID, y = SampleID, fill = value)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=8),
              axis.text.y = if(length(unique(mratio$SampleID))> 96) { element_blank() }else{element_text(size=8)}) +
#     axis.text.y = element_text(size=8)) +
    geom_tile() +
    scale_fill_gradient2(low = jRdGyPalette[1],
                         mid = jRdGyPalette[paletteSize/2],
                         high = jRdGyPalette[paletteSize],
                         midpoint = 0,
                         name = "log10 Read Count") +
    labs(title=paste("Overlapping Read Bias\nPositive values indicate joined reads\nnegative values indecate unpaired reads\nCombined_Percentage:\t",flash_data[4],"%",sep=""))

png(file.path(output,"merged_read_results.png"),width=8,height=10.5,units="in",res=300)
print(p)
invisible(dev.off())

#####################################################################
### Analysis Functions BEGIN
#####################################################################

######################################################################
## haplotypes - NOTREADY
##  compute haplotypes for each sample/amplicon
## 
## Parameters
##  name: amplicon name
##  min_freq: minimum frequence to accept variant
##  min_seq: minimum number of sequences to accept variant
# "haplotypes" <-  function(name,min_freq,min_seq) {
# #    cat(name,"\n")      
#     seqs <- DNAStringSet()
#     
#     ssingle = c(); spaired=c();
#     if (name %in% names(splitfq)){
#         ssingle <- table(splitfq[[name]])[table(splitfq[[name]])/counts[[name]] >=min_freq & table(splitfq[[name]])>= min_seq]
#         if (length(ssingle) > 0){
#             ssingle_str <- DNAStringSet(names(ssingle))
#             names(ssingle_str) <- paste(name,"merged",ssingle, counts[[name]], signif(ssingle/counts[[name]],3),sep="|")
#             seqs <- c(seqs,ssingle_str)
#         }
#     }
#     if (name %in% names(splitfq_p)){
#         tmp <- splitfq_p[[name]]
#         lapply(split(tmp,width(tmp)), function(wtmp){    
#             abc<-alphabetByCycle(wtmp)
#             consensus <- rownames(abc)[apply(abc,2,which.max)]
#             sites <- apply(sweep(abc,MARGIN=2,STATS=colSums(abc),FUN="/"),2,max) < (1-0.25)
#             lapply(strsplit(as.character(wtmp),split=""),"[",sites)
#         }
#         sweep(abc,MARGIN=2,STATS=colSums(abc),FUN="/")
#         apply(sweep(abc,MARGIN=2,STATS=colSums(abc),FUN="/"),2,max)
#         summary(apply(sweep(abc,MARGIN=2,STATS=colSums(abc),FUN="/"),2,max))
# 
# 
# 
#         spaired <- table(splitfq_p[[name]])[table(splitfq_p[[name]])/counts[[name]] >=min_freq & table(splitfq_p[[name]])>= min_seq]
#         if (length(spaired) > 0){
#             spaired_str <- lapply(strsplit(names(spaired),split=".....",fixed = TRUE),DNAStringSet)
#             spaired_str <- do.call("c",spaired_str)
#             names(spaired_str) <- paste(name,c("read1","read2"),rep(spaired,each=2), counts[[name]], rep(signif(spaired/counts[[name]],3),each=2),sep="|")
#             seqs <- c(seqs,spaired_str)
#         }
#     }
#     seqs
# }
# 
######################################################################
## consensus
##  compute strict consensus sequence for each sample/amplicon
## 
## Parameters
##  name: amplicon name
##  min_freq: minimum frequence to accept variant
##  min_seq: minimum number of sequences to accept variant
"consensus" <- function(name,min_freq,min_seq) {
#    cat(name,", ")
    seqs <- DNAStringSet()
    maxsize = 0    
    if (name %in% names(splitfq)){
        tmp <- splitfq[[name]]
        wtmp = split(tmp,width(tmp))
        wtmp = wtmp[[which.max(sapply(wtmp,length))]]
        if(length(wtmp) > maxsize & length(wtmp) > min_seq & (length(wtmp)/counts[[name]]) > min_freq){
            abc<-alphabetByCycle(wtmp)
            consensus <- paste(rownames(abc)[apply(abc,2,which.max)],collapse="")
            error_rate=sum(1-apply(sweep(abc,MARGIN=2,STATS=colSums(abc),FUN="/"),2,max))/nchar(consensus)
            ssingle_str <- DNAStringSet(consensus)
            names(ssingle_str) <- paste(name,"merged", seq.int(1,length(ssingle_str)), length(wtmp), counts[[name]], round(length(wtmp)/counts[[name]],3),round(error_rate,3),sep="|")
            maxsize = length(wtmp)
            seqs <- ssingle_str
        }
    }
    if (name %in% names(splitfq_p)){
        tmp <- splitfq_p[[name]]
        wtmp = split(tmp,width(tmp))
        wtmp = wtmp[[which.max(sapply(wtmp,length))]]
        if(length(wtmp) > maxsize & length(wtmp) > min_seq & (length(wtmp)/counts[[name]]) > min_freq){
            abc<-alphabetByCycle(wtmp)
            consensus <- paste(rownames(abc)[apply(abc,2,which.max)],collapse="")
            error_rate=sum(1-apply(sweep(abc,MARGIN=2,STATS=colSums(abc),FUN="/"),2,max))/(nchar(consensus)-5) ## - 5 for Ns
            spaired_str <- lapply(strsplit(consensus,split=".....",fixed = TRUE),DNAStringSet)
            spaired_str <- do.call("c",spaired_str)
            names(spaired_str) <- paste(name,c("read1","read2"),rep(seq.int(1,length(spaired_str)/2),each=2),length(wtmp), counts[[name]], round(length(wtmp)/counts[[name]],3),round(error_rate,3),sep="|")            
            maxsize = length(wtmp)
            seqs <- spaired_str
        }
    }
    seqs
}

######################################################################
## ambiguities
##  compute ambiguity sequence for each sample/amplicon
## 
## Parameters
##  name: amplicon name
##  min_freq: minimum frequence to accept variant
##  min_seq: minimum number of sequences to accept variant
"ambiguities" <- function(name,min_freq,min_seq) {
#    cat(name,", ")
    seqs <- DNAStringSet()
    maxsize = 0    
    if (name %in% names(splitfq)){
        tmp <- splitfq[[name]]
        wtmp = split(tmp,width(tmp))
        wtmp = wtmp[[which.max(sapply(wtmp,length))]]
        if(length(wtmp) > maxsize & length(wtmp) > min_seq & (length(wtmp)/counts[[name]]) > min_freq){
            abc<-alphabetByCycle(wtmp)
            abcf<-sweep(abc,MARGIN=2,STATS=colSums(abc),FUN="/")
            abcA = abc
            abcA[abc < min_seq | abcf < min_freq] = 0
            consensus <- names(IUPAC_CODE_MAP)[match(apply(abcA,2,function(x) paste(rownames(abcA)[x>0],collapse="")),IUPAC_CODE_MAP)]
            consensus[is.na(consensus)] = "N"
            amb_bases <- sum(consensus %in% names(IUPAC_CODE_MAP[-c(1:4)]))
            consensus = paste(consensus,collapse="")
            error_rate=sum(1-apply(sweep(abc,MARGIN=2,STATS=colSums(abc),FUN="/"),2,max))/nchar(consensus)
            ssingle_str <- DNAStringSet(consensus)
            names(ssingle_str) <- paste(name,"merged", seq.int(1,length(ssingle_str)), length(wtmp), counts[[name]], round(length(wtmp)/counts[[name]],3), round(error_rate,3), amb_bases,sep="|")
            maxsize = length(wtmp)
            seqs <- ssingle_str
        }
    }
    if (name %in% names(splitfq_p)){
        tmp <- splitfq_p[[name]]
        wtmp = split(tmp,width(tmp))
        wtmp = wtmp[[which.max(sapply(wtmp,length))]]
        if(length(wtmp) > maxsize & length(wtmp) > min_seq & (length(wtmp)/counts[[name]]) > min_freq){
            abc<-alphabetByCycle(wtmp)
            abcf<-sweep(abc,MARGIN=2,STATS=colSums(abc),FUN="/")
            abcA = abc
            abcA[abc < min_seq | abcf < min_freq] = 0
            consensus <- names(IUPAC_CODE_MAP)[match(apply(abcA,2,function(x) paste(rownames(abcA)[x>0],collapse="")),IUPAC_CODE_MAP)]
            consensus[is.na(consensus)] = "N"
            consensus[which(abc['.',] > (length(wtmp)*0.9))] = "."  ## at least 90% of reads must have '.' at the join location
            amb_bases <- sum(consensus %in% names(IUPAC_CODE_MAP[-c(1:4)])) ## - 5 for Ns
            consensus = paste(consensus,collapse="")
            error_rate=sum(1-apply(sweep(abc,MARGIN=2,STATS=colSums(abc),FUN="/"),2,max))/(nchar(consensus)-5) ## - 5 for Ns
            spaired_str <- lapply(strsplit(consensus,split=".....",fixed = TRUE),DNAStringSet)
            spaired_str <- do.call("c",spaired_str)
            names(spaired_str) <- paste(name,c("read1","read2"),rep(seq.int(1,length(spaired_str)/2),each=2),length(wtmp), counts[[name]], round(length(wtmp)/counts[[name]],3),round(error_rate,3),amb_bases,sep="|")            
            maxsize = length(wtmp)
            seqs <- spaired_str
        }
    }
    seqs
}

######################################################################
## occurrence
##  output all sequences, for each sample/amplicon, which meet the min_freq and min_seq criteria
## 
## Parameters
##  name: amplicon name
##  min_freq: minimum frequence to accept variant
##  min_seq: minimum number of sequences to accept variant
"occurrence" <- function(name,min_freq,min_seq) {
#    cat(name,", ")
    seqs <- DNAStringSet()
    
    ssingle = c(); spaired=c();
    if (name %in% names(splitfq)){
        ssingle <- table(splitfq[[name]])[table(splitfq[[name]])/counts[[name]] >=min_freq & table(splitfq[[name]])>= min_seq]
        if (length(ssingle) > 0){
            ssingle_str <- DNAStringSet(names(ssingle))
            names(ssingle_str) <- paste(name,"merged", seq.int(1,length(ssingle_str)), ssingle, counts[[name]], signif(ssingle/counts[[name]],3),sep="|")
            seqs <- c(seqs,ssingle_str)
        }
    }
    if (name %in% names(splitfq_p)){
        spaired <- table(splitfq_p[[name]])[table(splitfq_p[[name]])/counts[[name]] >=min_freq & table(splitfq_p[[name]])>= min_seq]
        if (length(spaired) > 0){
            spaired_str <- lapply(strsplit(names(spaired),split=".....",fixed = TRUE),DNAStringSet)
            spaired_str <- do.call("c",spaired_str)
            names(spaired_str) <- paste(name,c("read1","read2"),rep(seq.int(1,length(spaired_str)/2),each=2),rep(spaired,each=2), counts[[name]], rep(signif(spaired/counts[[name]],3),each=2),sep="|")
            seqs <- c(seqs,spaired_str)
        }
    }
    seqs
}

#####################################################################
### Analysis Functions END
#####################################################################

## merge
merged <- DNAStringSet(paste(as.character(sread(fq_r1)),".....",as.character(reverseComplement(sread(fq_r2))),sep=""))
names(merged) <- nms_p

single <- sread(fq)
names(single) <- nms

joined = paste(id,primer,sep=":")
joined_p = paste(id_p,primer_p,sep=":")
counts <- table(c(joined,joined_p))

anames <- unique(c(joined,joined_p))
splitfq <- split(single,paste(id,primer,sep=":"))
splitfq_p <- split(merged,paste(id_p,primer_p,sep=":"))

sapply(program_list,function(program){
    write(paste("Running program:",program),stdout())
    freq_seqs <- mclapply(anames, function(nms) do.call(program,list(name=nms,min_freq=min_freq,min_seq=min_seq)), mc.cores = procs)
    
    redo <- sapply(freq_seqs,function(x) class(x) != "DNAStringSet")
    if (sum(redo) > 0){
        write(paste("There were",sum(redo),"failed amplicon analysis, trying to rerun"),stdout())
#        stop(paste("There were",sum(redo),"failed amplicon analysis, trying to rerun"))
    } 
    retry = 0
    while (any(redo) & retry < 5){
        freq_seqs[redo] <- mclapply(anames[redo], do.call(program,list(name=nms,freq_min=min_freq,min_seq=min_seq)), mc.cores = procs)
        redo <- sapply(freq_seqs,function(x) class(x) != "DNAStringSet")
        retry = retry+1
        write(paste("Retry number", retry,sep=" "),stdout())
    }
    if (any(redo)) stop("Failed to complete all amplicons, try reducing the number of processors.")
    write(paste("Finished analyzing amplicons"),stdout())
        
    names(freq_seqs) <- anames
    result_seqs <- sapply(freq_seqs,length)
    tt <- DNAStringSet(unlist(unname(sapply(freq_seqs,as.character))))
   
    onms <- names(tt)
    first <- sapply(strsplit(onms,split="|",fixed=TRUE),"[[",1L)
    oprimer <- sapply(strsplit(first,split=":"),function(x) tail(x,1))
    oid <- apply(cbind(paste(":",oprimer,sep=""),first),1,function(x) sub(x[1],"",x[2]))
    
    second <- sapply(strsplit(onms,split="|",fixed=TRUE),"[[",2L)

    write(paste("Writing Reads"),stdout())
    writeXStringSet(tt,file.path(output,paste(program,"reduced.fasta",sep=".")))

    dir.create(file.path(output,paste(program,"split_samples",sep=".")))
    split_tt <- split(tt,oid)
    mclapply(names(split_tt), function(x){
        writeXStringSet(split_tt[[x]],file.path(output,paste(program,"split_samples",sep="."),paste("Sample",x,"fasta",sep=".")))
    }, mc.cores = procs)
    dir.create(file.path(output,paste(program,"split_amplicon",sep=".")))
    split_tt <- split(tt,paste(oprimer,second,sep="."))
    mclapply(names(split_tt), function(x){
        writeXStringSet(split_tt[[x]],file.path(output,paste(program,"split_amplicon",sep="."),paste("Amplicon",x,"fasta",sep=".")))
    }, mc.cores = procs)
    
    write(paste("Producing final images"),stdout())
    #### PLOTTING RESULTS    
    png(file.path(output,paste(program,"freq_read_counts.png",sep=".")),height=8,width=10.5,units="in",res=300)
    hist(as.numeric(sapply(strsplit(onms,split="|",fixed=TRUE),"[[",5L)),breaks=200,main=paste("histogram of amplicon read counts\n",program," trimR1:",trimOne," trimR2:",trimTwo,"\n",sep=""),xlab="frequency")
    invisible(dev.off())
    png(file.path(output,paste(program,"freq_most_occuring.png",sep=".")),height=8,width=10.5,units="in",res=300)
    hist(as.numeric(sapply(strsplit(onms,split="|",fixed=TRUE),"[[",6L)),breaks=200,main=paste("histogram of chosen amplicon frequency\n",program," trimR1:",trimOne," trimR2:",trimTwo,"\n",sep=""),xlab="read counts")
    invisible(dev.off())
    if (program %in% c("consensus","ambiguities")){
        png(file.path(output,paste(program,"error_rate.png",sep=".")),height=8,width=10.5,units="in",res=300)
        hist(as.numeric(sapply(strsplit(onms,split="|",fixed=TRUE),"[[",7L)),breaks=200,main=paste("histogram of error rate in amplicons\n",program," trimR1:",trimOne," trimR2:",trimTwo,"\n",sep=""),xlab="error rate")
        invisible(dev.off())
        if (program %in% c("ambiguities")){
            png(file.path(output,paste(program,"ambiguities.png",sep=".")),height=8,width=10.5,units="in",res=300)
            hist(as.numeric(sapply(strsplit(onms,split="|",fixed=TRUE),"[[",8L)),breaks=200,main=paste("histogram of number of ambiguites in amplicons\n",program," trimR1:",trimOne," trimR2:",trimTwo,"\n",sep=""),xlab="number of ambiguity bases")
            invisible(dev.off())
        }
    }
    
    occurrence_a <- table(oid,oprimer)
    uoccurrence_a <- matrix(0,nrow=length(uid),ncol=length(uprimer))
    rownames(uoccurrence_a) <- uid
    colnames(uoccurrence_a) <- uprimer
    
    uoccurrence_a[rownames(occurrence_a),colnames(occurrence_a)] <- occurrence_a
    
    mratio <- melt(uoccurrence_a)
    colnames(mratio) <- c("SampleID","PrimerID","value")
    mratio$value = as.factor(mratio$value)
    
    jRdGyPalette <- brewer.pal(n = 4, "Set1")
    paletteSize <- 4
    
    p <- ggplot(mratio, aes(x = PrimerID, y = SampleID, fill = value)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=8), 
              axis.text.y = if(length(unique(mratio$SampleID))> 96) { element_blank() }else{element_text(size=8)}) +
        geom_tile() +
        scale_fill_brewer(palette="Set1") +
        #     scale_fill_gradient2(low = jRdGyPalette[1],
        #                          mid = jRdGyPalette[paletteSize/2],
        #                          high = jRdGyPalette[paletteSize],
        #                          midpoint = 0,
        #                          name = "Amplicon Count") +
        labs(title=paste("Resulting number of amplicons\n",program," trimR1:",trimOne," trimR2:",trimTwo,"\n",sep=""))
    
    png(file.path(output,paste(program,"Amplicons","Per","Sample","png",sep=".")),width=8,height=10.5,units="in",res=300)
    print(p)
    invisible(dev.off())
})

write("Finished",stdout())

