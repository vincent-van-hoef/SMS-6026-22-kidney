process FASTP {
	
	tag "$id"
	publishDir "${params.outdir}/QC/fastp/", mode: 'copy', overwrite: true, pattern: '*.json'

	executor "slurm"
	cpus 3
	time 90.m
	module 'bioinfo-tools:fastp/0.23.2'
	clusterOptions '-A sens2022505'

	input:
	tuple val(id), path(reads)

	output:
	tuple val(id), path("*.json"), emit: json
	tuple val(id), path("${id}_1_clean.fastq.gz"), path("${id}_2_clean.fastq.gz"), emit: clean_fq

	script:
	"""
	ln -s ${reads[0]} ${id}_1.fastq.gz
	ln -s ${reads[1]} ${id}_2.fastq.gz

	fastp -i ${id}_1.fastq.gz \
		-I ${id}_2.fastq.gz \
		-o ${id}_1_clean.fastq.gz \
		-O ${id}_2_clean.fastq.gz \
		-p \
		-j ${id}_fastp.json \
		-h ${id}_fastp.html \
		-R "${id}_fastp_report" \
		-w 3
	"""
} 

process FASTQC {

	tag "$id"
	publishDir "${params.outdir}/QC/raw/", mode: 'copy', overwrite: true, pattern: '*.zip'

	executor "slurm"
	cpus 3
	time 90.m
	scratch false
	module 'bioinfo-tools:FastQC/0.11.8'
	clusterOptions '-A sens2022505'

	input:
	tuple val(id), path(reads)

	output:
	tuple val(id), path("*.html"), emit: html
	tuple val(id), path("*.zip"), emit: zip

	script:
	"""
	fastqc --threads $task.cpus $reads
	"""
}

process TRIMGALORE {

	tag "$id"

	executor "slurm"
	cpus 12
	time 120.m
	scratch false
	module 'bioinfo-tools:TrimGalore/0.6.1'
	clusterOptions '-A sens2022505'

	input:
	tuple val(id), path(reads)

	output:
	tuple val(id), path("*val*.fq.gz")		, emit: reads

	script:
	"""
	ln -s ${reads[0]} ${id}_1.fastq.gz
        ln -s ${reads[1]} ${id}_2.fastq.gz

	trim_galore -q 20 \
		--paired \
		--cores 4 \
		${id}_1.fastq.gz \
		${id}_2.fastq.gz
	"""
}

process FASTQC_TRIM {

	tag "$id"
	publishDir "${params.outdir}/QC/trimmed/", mode: 'copy', overwrite: true, pattern: '*.zip'

	executor "slurm"
	cpus 3
	time 90.m
	scratch false
	module 'bioinfo-tools:FastQC/0.11.8'
	clusterOptions '-A sens2022505'

	input:
	tuple val(id), path(reads)

	output:
	tuple val(id), path("*.html"), emit: html
	tuple val(id), path("*.zip"), emit: zip

	script:
	"""
	fastqc --threads $task.cpus $reads
	"""
}

process MULTIQC {

	publishDir "${params.outdir}/QC/trimmed/", mode: 'copy', overwrite: true, pattern: '*trimmed_multiqc*'
	publishDir "${params.outdir}/QC/raw/", mode: 'copy', overwrite: true, pattern: '*raw_multiqc*'
	publishDir "${params.outdir}/QC/fastp/", mode: 'copy', overwrite: true, pattern: '*fastp_multiqc*'

	executor "slurm"
	cpus 3
	time 90.m
	module 'bioinfo-tools:MultiQC'
	clusterOptions '-A sens2022505'

	input:
	path(fastp)
	path(fastqc)
	path(trimgalore)

	output:
	path("*")

	script:
	"""
	multiqc *_fastp.json -n fastp_multiqc
	multiqc *_R*_fastqc.zip -n raw_multiqc
	multiqc *_val_*_fastqc.zip -n trimmed_multiqc
	"""
}

process SALMON_INDEX {

	storeDir "/proj/sens2022505/analysis/db/"

	executor "slurm"
	cpus 6
	time 120.m
	module 'bioinfo-tools:Salmon/1.6.0'
	clusterOptions '-A sens2022505'

	input:
	path(transcripts)
	path(genome)

	output:
	path(salmon_index)

	script:
	"""
	grep '^>' $genome | cut -d ' ' -f 1 > decoys.txt
	sed -i.bak -e 's/>//g' decoys.txt
        cat $transcripts $genome > tx_decoy.fa
        salmon index \
        	--threads $task.cpus \
        	-t tx_decoy.fa \
        	-d decoys.txt \
        	-i salmon_index
	"""
}

process SALMON_QUANT {

	tag "$id"
	publishDir "${params.outdir}/quant/salmon/", mode: 'copy', overwrite: true

	executor "slurm"
	cpus 4
	time 120.m
	module 'bioinfo-tools:Salmon/1.6.0'
	clusterOptions '-A sens2022505'

	input:
	tuple val(id), path(freads), path(rreads)
	path(index)

	output:
	path("${id}")

	script:
	"""
	salmon quant \
		-i $index \
		-l A \
		-1 $freads \
		-2 $rreads \
		--validateMappings \
		--seqBias \
		--gcBias \
		-o ${id}
	"""
}

process CREATE_TX2GENE {

	input:
	path(txs)

	output:
	path("tx2gene.txt")

	script:
	"""
	cat $txs | grep '>' |cut -d ' ' -f1,4,7 > tmp
	paste <(cut -d '>' -f2 tmp | cut -d ' ' -f1) <(cut -d ' ' -f2 tmp | cut -d ':' -f2) >> tx2gene.txt
	"""
}

process CREATE_META {

	publishDir "${params.outdir}", mode: 'copy', overwrite: true

	input:
	path(ids)
	path(meta)

	output:
	path("meta.csv")

	script:
	"""
	awk -F ',' 'NR==FNR {id[\$1]; next} \$1 in id' $ids $meta > tmp.csv
	cat <(echo "sample,group") tmp.csv > meta.csv
	"""

}

process SALMON_TXIMPORT {

	publishDir "${params.outdir}/quant/salmon/merged", mode: 'copy', overwrite: true

	executor "slurm"
	cpus 2
	time 30.m
	module 'bioinfo-tools:R_packages/4.1.1'
	clusterOptions '-A sens2022505'

	input:
	path(txs2gene)
	path(quant)
	path(meta)

	output:
	path("txi.Rds")

	script:
	"""
	#!/usr/bin/env Rscript
	
	library(tximport)
	
	meta <- read.csv("${meta}", header = TRUE)
	tx2gene <- read.table("${txs2gene}")
		
	files <- file.path(meta\$sample, "quant.sf")
	names(files) <- gsub("/.*", "", files)
	txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

	stopifnot(identical(meta\$sample, colnames(txi\$counts)))
	
	saveRDS(txi, "txi.Rds")
	"""
}

process CREATE_DESEQ {

	executor "slurm"
	cpus 2
	time 30.m
	module 'bioinfo-tools:R_packages/4.1.1'
	clusterOptions '-A sens2022505'

	input:
	path(meta)
	path(txi)
	
	output:
	path("dds.Rds")

	script:
	"""
	#!/usr/bin/env Rscript

	library(DESeq2)

	meta <- read.csv("${meta}", header = TRUE)
	rownames(meta) <- meta\$sample
	txi <- readRDS("${txi}")

	stopifnot(identical(rownames(meta), colnames(txi\$counts)))

	dds <- DESeqDataSetFromTximport(txi, meta, ~ 1)	
	saveRDS(dds, "dds.Rds")
	"""

}

process QC_DESEQ {

	publishDir "${params.outdir}/DESeq2/QC", mode: 'copy', overwrite: true

	executor "slurm"
	cpus 12
	time 3.h
	module 'bioinfo-tools:R_packages/4.1.1'
	clusterOptions '-A sens2022505'

	input:
	path(dds)
	each gene_number
	
	output:
	path("PCA*.png")
	path("Heatmap*.png")

	script:
	"""
	#!/usr/bin/env Rscript

	library(DESeq2)
	library(ggplot2)
	library(pheatmap)

	dds <- readRDS("${dds}")

	# Diffferent Normalizations
	vsd <- vst(dds, blind = FALSE)
	ntd <- normTransform(dds)	
	
	# PCA plots
	plotPCA(vsd, intgroup = "group", ntop = $gene_number) + theme_bw()
	ggsave("PCA_ntop_${gene_number}_vsd.png")
	plotPCA(ntd, intgroup = "group", ntop = $gene_number) + theme_bw()
	ggsave("PCA_ntop_${gene_number}_ntd.png")

	# Heatmaps
	gene_no <- ${gene_number}/10
	dds <- estimateSizeFactors(dds)
	select <- order(rowVars(counts(dds, normalized = TRUE)), decreasing=TRUE)[1:gene_no]	
	df <- as.data.frame(colData(dds))[,"group", drop = FALSE]
	png(paste0("Heatmap_", gene_no, "_most_variable_genes_vsd.png"))
	pheatmap(assay(vsd)[select,], 
		cluster_rows = FALSE, 
		show_rownames = FALSE,
		cluster_cols = TRUE,
		annotation_col = df)
	dev.off()
	png(paste0("Heatmap_", gene_no, "_most_variable_genes_ntd.png"))
	pheatmap(assay(ntd)[select,], 
		cluster_rows = FALSE, 
		show_rownames = FALSE,
		cluster_cols = TRUE,
		annotation_col = df)
	dev.off()
	"""

}

process QC_DESEQ_HM {

	publishDir "${params.outdir}/DESeq2/QC", mode: 'copy', overwrite: true

	executor "slurm"
	cpus 2
	time 1.h
	module 'bioinfo-tools:R_packages/4.1.1'
	clusterOptions '-A sens2022505'

	input:
	path(dds)
	
	output:
	path("*.png")

	script:
	"""
	#!/usr/bin/env Rscript

	library(DESeq2)
	library(ggplot2)
	library(RColorBrewer)
	library(pheatmap)

	dds <- readRDS("${dds}")
	vsd <- vst(dds, blind = FALSE)

	sampleDists <- dist(t(assay(vsd)))
	sampleDistMatrix <- as.matrix(sampleDists)
	rownames(sampleDistMatrix) <- vsd\$group
	colnames(sampleDistMatrix) <- NULL
	colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
	png("sample_to_sample_distance.png")
	pheatmap(sampleDistMatrix,
		clustering_distance_rows=sampleDists,
		clustering_distance_cols=sampleDists,
		col=colors)
	dev.off()			
	"""
}

process RUN_DESEQ {

	publishDir "${params.outdir}/DESeq2/DE", mode: 'copy', overwrite: true

	executor "slurm"
	cpus 6
	time 1.h
	module 'bioinfo-tools:R_packages/4.1.1'
	clusterOptions '-A sens2022505'

	input:
	path(dds)
	
	output:
	path("dds_group.Rds")

	script:
	"""
	#!/usr/bin/env Rscript

	library(DESeq2)
	
	dds <- readRDS("${dds}")
	colData(dds)\$group <- as.factor(colData(dds)\$group)
	design(dds) <- ~ group
	dds <- DESeq(dds)	

	saveRDS(dds, "dds_group.Rds")
	"""
}


process EXTRACT_RES {

         tag "${contrast}"

         publishDir "${params.outdir}/DESeq2/DE/${contrast}", mode: 'copy'

         executor "slurm"
         cpus 4
         time 30.m
         module 'bioinfo-tools:R_packages/4.1.1'
         clusterOptions '-A sens2022505'

         input:
         tuple val(contrast), val(treat), val(contr), path(dds)

         output:
         path("*.xlsx")
         path("*.csv")
         path("*.png")

         script:
         """
         #!/usr/bin/env Rscript

         library(DESeq2)
         library(dplyr)
         library(tidyr)
         library(grid)
         library(ggplot2)
         library(ggfortify)
	 library(org.Hs.eg.db)

         dds <- readRDS("${dds}")
         res <- results(dds, contrast = c("group", "${treat}", "${contr}"))
	 res <- res %>%
        		as.data.frame() %>%
        		mutate(ensembl = gsub("\\\\.[0-9]*\$", "", rownames(.)),
               			symbol = mapIds(org.Hs.eg.db, keys = ensembl, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")) %>%
        		arrange(padj, pvalue)
         openxlsx::write.xlsx(res, file = "${contrast}.xlsx", overwrite = TRUE)
         write.csv(res, file = "${contrast}.csv")

         res_sig <- res %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)

         res_sig_up <- res_sig %>% filter(log2FoldChange > 1) %>% nrow
         res_sig_down <- res_sig %>% filter(log2FoldChange < -1) %>% nrow

         grob <- grobTree(textGrob(paste0(res_sig_up, " (log2FC > 1 & p.adj < 0.05)"), x = 0.57, y = 0.95, hjust = 0),
                             textGrob(paste0(res_sig_down, " (log2FC < -1 & p.adj < 0.05)"), x = 0.01, y = 0.95, hjust = 0))

         png("${contrast}_volcano.png")
         ggplot(res, aes(x=log2FoldChange, y=-log10(pvalue))) +
              geom_point(alpha = 0.2) +
              annotation_custom(grob) +
              geom_hline(yintercept = -log10(0.05)) +
              geom_vline(xintercept = 0) +
              geom_vline(xintercept = c(-1,1), col = "grey") +
              geom_point(data=res_sig, aes(log2FoldChange, y=-log10(pvalue)), col = "red", alpha = 0.2) +
              ggtitle("${contrast}") +
              theme_bw()
         dev.off()

	 plotGenes <- res %>% arrange(-log2FoldChange) %>% head()	
	 for(gene in rownames(plotGenes)){
		png(paste0(gene, "_plot_counts.png"))
		plotCounts(dds, gene = gene, intgroup = "group")
		dev.off()
	}
         """
 }
