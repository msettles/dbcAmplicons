### Generate a QC report for results

libary(qrqc)

s1.fastq <- readSeqFile("Fish_Intestinal-Hagerman_R1.fastq.gz")
s2.fastq <- readSeqFile("Fish_Intestinal-Hagerman_R2.fastq.gz")
u.fastq <- readSeqFile("Fish_Intestinal-Hagerman.extendedFrags.fastq.gz")

layout(matrix