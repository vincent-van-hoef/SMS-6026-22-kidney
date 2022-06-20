process FASTP {
	
	tag "$id"
	publishDir "${params.outdir}/QC/", mode: 'symlink', overwrite: true, pattern: '*.{html,json}'
	publishDir "${params.outdir}/fastq_clean/", mode: 'symlink', overwrite: true, pattern: '*.fastq.gz'

	executor "slurm"
	cpus 3
	time 90.m
	module 'bioinfo-tools:fastp/0.23.2'
	clusterOptions '-A sens2022505'

	input:
	tuple val(id), path(reads)

	output:
	tuple val(id), path("*.html"), emit: html
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
		-j ${id}.json \
		-h ${id}.html \
		-R "${id}_fastp_report" \
		-w 3
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
	publishDir "${params.outdir}/quant/salmon/", mode: 'symlink', overwrite: true

	executor "slurm"
	cpus 6
	time 60.m
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

process SALMON_TXIMPORT {

	publishDir "${params.outdir}/quant/salmon/merged", mode: 'symlink', overwrite: true

	executor "slurm"
	cpus 6
	time 30.m
	module 'bioinfo-tools:R_packages/4.1.1'
	clusterOptions '-A sens2022505'

	input:
	path(txs2gene)
	path(quant)

	output:
	path("txi.Rds")

	script:
	"""
	#!/usr/bin/env Rscript
	
	library(tximport)

	tx2gene <- read.table("${txs2gene}")
		
	files <- list.files(pattern = "quant.sf", recursive = TRUE)
	names(files) <- gsub("/.*", "", files)
	txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
	
	saveRDS(txi, "txi.Rds")
	"""
}
