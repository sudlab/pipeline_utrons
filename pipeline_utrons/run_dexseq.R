library(DEXSeq)
library(GenomicRanges)
library(optparse)
library(tidyr)
library(dplyr)
library(data.table)

opts = list(make_option(c("-p","--processes"),
                        type="integer",
                        dest="processes",
                        default=1,
                        help="Number of procsses to use. Using multiple processes requires the BiocParallel library"),
            make_option(c("-g", "--tx2gene"),
                      
                        dest="tx2gene", 
                        help="Two column file makking transcripts to gene. Expecting two columns: transcript_id and match_gene_id"),
          make_option(c("-O", "--outfile"),
                       
                        dest="outfile",
                        default="deseq_results.tsv",
                        help="output file for the DEXSeq results table"),
          make_option(c("-c", "--control"),
                      
                      dest="control", 
                      help="condition value for control condition, where condition is CANCERTYPE-CONDITION-REP"),
          make_option(c("-t", "--treatment"),
                      
                      dest="treatment",
                      help="condition value for treatment/test condition, where condition is CANCERTYPE-CONDITION-REP"))

args <- optparse::parse_args2(OptionParser(option_list=opts))
#args <- list(args=c("/shared/sudlab1/utrons/TCGA_GTEx/utrons/colon_utrons/salmon_quant.tsv"),
        #     options=list(processes=1, control="COAD-NO", treatment="COAD", outfile="dexseq_out.tsv",
        #                  tx2gene="/shared/sudlab1/utrons/TCGA_GTEx/utrons/colon_utrons/expression.dir/csvdb_files/tx2gene.txt"))
cat("Args are:")
print(args)
expression_data <- fread(args$args)
tx2gene <- read.delim(args$options$tx2gene)
expression_data_matrix = dcast(expression_data, Name ~ track, value.var="NumReads")
rn <- expression_data_matrix$Name
rm(expression_data)
expression_data_matrix <- as.matrix(expression_data_matrix[,-1])
rownames(expression_data_matrix) <- rn


col_data <- data.frame(sample=colnames(expression_data_matrix)) %>%
  extract(sample, into=c("cancer", "condition",  "rep"), regex="([^-]+)-(.+)-([^-]+\\.b)", remove = FALSE)

print(dim(expression_data_matrix))
expression_matrix <- expression_data_matrix[rownames(expression_data_matrix) %in% tx2gene$transcript_id,]
print(dim(expression_matrix))
print(head(rownames(expression_data_matrix)))
print(head(tx2gene$transcript_id))
tx2gene_sub <- tx2gene[tx2gene$transcript_id %in% rownames(expression_matrix),]
rownames(tx2gene) <- tx2gene$transcript_id
tx2gene_sub <- tx2gene[rownames(expression_matrix),]

col_data_sub <- col_data %>%
  filter(condition %in% c(args$options$control, args$options$treatment)) 

col_data_sub$condition <- factor(col_data_sub$condition, levels=c(args$options$control, args$options$treatment))

expression_matrix_sub <- expression_matrix[,col_data_sub$sample]
rm(expression_data_matrix)
rm(expression_matrix)
gc()

storage.mode(expression_matrix_sub) <- "integer"
print(dim(expression_matrix_sub))
print(dim(col_data_sub))
dxd <- DEXSeqDataSet(expression_matrix_sub,
                     sampleData = data.frame(col_data_sub),
                     design=~sample + condition + condition:exon, 
                     featureID=tx2gene_sub$transcript_id,
                     groupID=tx2gene_sub$match_gene_id)


if (args$options$processes > 1) {
  BPPARAM = MulticoreParam(workers=args$options$processes)
  dxd_results <- DEXSeq(dxd, BPPARAM=BPPARAM, quiet=F)
} else {
  dxd_results <- DEXSeq(dxd, quiet=F)
}

#save the results
save.image(paste0(args$options$outfile, ".RData"))
write.table(data.frame(dxd_results), paste0(args$options$outfile, ".tsv"),
            sep="\t",
            quote=FALSE,
            row.names = FALSE)
