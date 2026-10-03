library(rtracklayer)
library(GenomeInfoDb)
library(GenomicFeatures)
library(tximport)
gtf = import('/data/ChioriniCompCor/Pipeline/ATAC/CHM13/CHM13/chm13/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf')
chrominfo <- getChromInfoFromNCBI("T2T-CHM13v2.0")
seqlevels(gtf) <- setNames(chrominfo$SequenceName, chrominfo$RefSeqAccn)
seqinfo(gtf) <- Seqinfo(genome="T2T-CHM13v2.0")
txdb <- makeTxDbFromGRanges(gtf, taxonomyId=9606)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

files <- sort(Sys.glob(c('/data/ChioriniCompCor/Pipeline/RNA/CHM13/Output/BRA*/quant.sf','/data/ChioriniCompCor/Pipeline/RNA/CHM13/Output/MSG_*/quant.sf')))
filenames <- sort(unlist(lapply(Sys.glob(c('/data/ChioriniCompCor/Pipeline/RNA/CHM13/Output/BRA*','/data/ChioriniCompCor/Pipeline/RNA/CHM13/Output/MSG_*')), basename)))
names(files) <- filenames
metadata = read.table('meta.csv',sep=',',header=TRUE,row.names='X')
files <- files[rownames(metadata)]
sampleTable = data.frame(disease = factor(metadata$SjD.Diagnosis),
	ssa = factor(metadata$SSA), ssb = factor(metadata$SSB),
	ana = factor(metadata$ANA), rf = factor(metadata$RF),
	cohort = factor(metadata$Cohort))

txi <- tximport(files,type='salmon',tx2gene=tx2gene)
txi.tx <- tximport(files,type='salmon',txOut=TRUE)
txi_adj <- sva::ComBat_seq(txi$counts,batch=sampleTable$cohort,group=sampleTable$disease)
txi.tx_adj <- sva::ComBat_seq(txi.tx$counts,batch=sampleTable$cohort,group=sampleTable$disease)
write.csv(txi_adj,file='all_gene.csv')
write.csv(txi.tx_adj,file='all_transcript.csv')

library(DESeq2)
for(catvar in colnames(sampleTable)){
	files_nona = files[sampleTable[catvar] != '']
	sampleTable_nona = sampleTable[sampleTable[catvar] != '',]
	txi <- tximport(files_nona,type='salmon',tx2gene=tx2gene)
	dds <- DESeqDataSetFromTximport(txi,sampleTable_nona,as.formula(paste('~cohort+',catvar)))
	dds <- DESeq(dds)
	res <- results(dds,alpha=0.05)
	write.csv(as.data.frame(res[order(res$padj),]),file=paste('de/gene.',catvar,'.csv',sep=''))
	txi.tx <- tximport(files_nona,type='salmon',txOut=TRUE)
	dds <- DESeqDataSetFromTximport(txi.tx,sampleTable_nona,as.formula(paste('~cohort+',catvar)))
	dds <- DESeq(dds)
	res <- results(dds,alpha=0.05)
	write.csv(as.data.frame(res[order(res$padj),]),file=paste('de/tx.',catvar,'.csv',sep=''))
}

