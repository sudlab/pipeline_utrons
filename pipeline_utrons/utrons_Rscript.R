#!/data/mb1cna/cgat/envs/sharc/bin/Rscript

#import necessary files
tx2gene<-read.delim(file='expression.dir/csvdb_files/tx2gene.txt', header = TRUE, sep = '\t', dec = '.')
novel_utrons_merged<-read.delim(file='expression.dir/csvdb_files/novel_utrons_ids.txt', header = TRUE, sep = '\t', dec = '.')
all_utrons_merged<-read.delim(file='expression.dir/csvdb_files/all_utrons_ids.txt', header=TRUE, sep='\t', dec='.')
db_artifacts<-read.delim(file='expression.dir/csvdb_files/60db_novel_utrons_ids.txt', header = TRUE, sep = '\t', dec = '.')

#rename match_gene_id to gene_id in tx2gene_db
names(tx2gene)<-c("transcript_id", "gene_id")

#define the quantification files
library(tximport)
dir<-"quantification.dir/"
IDs<-list.files(dir, pattern = "sf")
fileno<-length(IDs)

IDs<-IDs[1:fileno]
files<-file.path(dir, IDs) # in file.patch function you can use multiple arguments to structure your directory name. In here dir is you main bit of directory, IDs specify the second variable bit of directory and quant.sf specifies the last bit which is our count file
sample_names <- sapply(IDs, function(x) substr(x, 0, nchar(x)-3))
sample_names <- sapply(IDs, function(x) strsplit(x, "[.]")[[1]][1])
names(files) <- sample_names
all(file.exists(files)) #checks if you have all of the files it searched for

#tximport
library(GenomicFeatures)
library(readr)
library(tximport)

#Transcript-level estimates, txOut=TRUE avoids gene-level summarization
txi.tr<-tximport(files, type="salmon", txOut=TRUE, tx2gene = tx2gene, dropInfReps = TRUE, countsFromAbundance = "lengthScaledTPM") #dropInfReps gets rid of this error: "Error: package_version(minfo$salmon_version) >= "0.8.0" is not TRUE"

#Gene-level summarization with txOut=FALSE
txi.sum<-tximport(files, type="salmon", txOut=FALSE, tx2gene = tx2gene, dropInfReps = TRUE, countsFromAbundance = "lengthScaledTPM")
all.equal(txi.tr, txi.sum)

#save the $abundance (TPM) information
tpm.tr<-txi.tr$abundance
tpm.sum<-txi.sum$abundance

#remove the txi.tr and txi.sum
rm(txi.tr)
rm(txi.sum)

#convert the tpm files into data frames
tpm.tr<-data.frame(tpm.tr)
tpm.sum<-data.frame(tpm.sum)

#convert rownames to a first column and name it transcript_id/gene_id
library(tibble)
tpm.tr<-rownames_to_column(tpm.tr, var="transcript_id")
tpm.sum<-rownames_to_column(tpm.sum, var="gene_id")

#add gene_id matching the transcript_id and viceversa
tpm.tr<-merge(tpm.tr, tx2gene, by=c("transcript_id"), all.y = TRUE)
tpm.sum<-merge(tpm.sum, tx2gene, by=c("gene_id"), all.y = TRUE)

#bring the last column gene_id as first column, before transcript_id (if needed)
tpm.tr<-tpm.tr[,c("gene_id", setdiff(names(tpm.tr),"gene_id"))]
tpm.sum<-tpm.sum[,c("transcript_id", setdiff(names(tpm.tr),"transcript_id"))]
tpm.sum<-tpm.sum[,c("gene_id", setdiff(names(tpm.tr),"gene_id"))]

#filter the entries without any gene_id
tpm.sum<-with(tpm.sum, tpm.sum[!(gene_id == "" | is.na(gene_id)), ])
tpm.tr<-with(tpm.tr, tpm.tr[!(gene_id == "" | is.na(gene_id)), ])

#re-order the tpm.tr
library(plyr)
tpm.tr<-arrange(tpm.tr, gene_id)
tpm.sum<-arrange(tpm.sum, transcript_id)
tpm.sum<-arrange(tpm.sum, gene_id)

#filter the artefacts from the transcript txi
tpm.tr<-tpm.tr[!(tpm.tr$transcript_id %in% db_artifacts$match_transcript_id),]
tpm.sum<-tpm.sum[!(tpm.sum$transcript_id %in% db_artifacts$match_transcript_id),]

#built data frame with the fract. expression by dividing rows from transcripts to the rows from sum, by common transcript_id
fract.expr<-cbind(tpm.tr[1:2], tpm.tr[-(1:2)]/tpm.sum[match(tpm.tr$transcript_id, tpm.sum$transcript_id), -(1:2)])
fract.expr<-fract.expr[,c("gene_id", "transcript_id", setdiff(names(tpm.tr),c("gene_id", "transcript_id")))]

#Save RDAta for non-melted data frames
save.image(file = "expression.dir/nonmelted_data.RData")


#Melt transcript, sum and fract.expr tpm dataframes
library(reshape)
mtr<-melt(tpm.tr, id=c("gene_id", "transcript_id"))
msum<-melt(tpm.sum, id=c("gene_id", "transcript_id"))
mfract<-melt(fract.expr, id=c("gene_id", "transcript_id"))
#Re-order and rename colnames
mtr<-mtr[,c("variable", "transcript_id", "gene_id", setdiff(names(mtr),c("gene_id", "transcript_id", "variable")))]
names(mtr)<-c("Sample", "transcript_id", "gene_id", "tr.expr")
msum<-msum[,c("variable", "transcript_id", "gene_id", setdiff(names(msum),c("gene_id", "transcript_id", "variable")))]
names(msum)<-c("Sample", "transcript_id", "gene_id", "gene.expr")
mfract<-mfract[,c("variable", "transcript_id", "gene_id", setdiff(names(mfract),c("gene_id", "transcript_id", "variable")))]
names(mfract)<-c("Sample", "transcript_id", "gene_id", "fract.expr")

# Add columns from fract.expr, tpm.tr and tpm.sum together
mexpr<-merge(mtr, msum, by=c("Sample", "transcript_id", "gene_id"))
mdata<-merge(mexpr, mfract, by=c("Sample", "transcript_id", "gene_id"))
mdata<-arrange(mdata, Sample, gene_id)

#Save RDAta and the melted expression dataframe as a text file
save.image(file = "expression.dir/data.RData")
write.table(mdata, "expression.dir/utrons_expression.txt", sep = "\t", row.names = FALSE, quote = FALSE)

