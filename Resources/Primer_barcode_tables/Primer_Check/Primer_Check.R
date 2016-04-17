## Primer check

F_adapt <- read.table("Forward_adapters.txt",sep="\t",header=T,skip=1,as.is=T)
R_adapt <- read.table("Reverse_adapters.txt",sep="\t",header=T,skip=1,as.is=T)

F_targ <- read.table("Forward_target.txt",sep="\t",header=T,skip=1,as.is=T)
R_targ <- read.table("Reverse_target.txt",sep="\t",header=T,skip=1,as.is=T)

library(ShortRead)
F_targ.seq <- DNAStringSet(F_targ$Oligonucleotide.Sequence)
names(F_targ.seq) <- F_targ$Oligo.Name

R_targ.seq <- DNAStringSet(R_targ$Oligonucleotide.Sequence)
names(R_targ.seq) <- R_targ$Oligo.Name

F_adapt.seq <- DNAStringSet(F_adapt$Oligonucleotide.Sequence)
names(F_adapt.seq) <- F_adapt$Oligo.Name

R_adapt.seq <- DNAStringSet(R_adapt$Oligonucleotide.Sequence)
names(R_adapt.seq) <- R_adapt$Oligo.Name

tagdiff[tagpairs[,1] != tagpairs[,2]] <- apply(tagpairs[tagpairs[,1] != tagpairs[,2],],1,stringDist, method="hamming")

