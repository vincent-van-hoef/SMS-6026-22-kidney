#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
=================================
Analysis for #6026: RNAseq Kidney
=================================
*/

// Parameters
params.outdir	= "/proj/sens2022505/nobackup/results/"
params.meta	= "/proj/sens2022505/analysis/db/metadata.csv"
params.txs	= "/proj/sens2022505/analysis/db/Homo_sapiens.GRCh38.cdna.all.fa"
params.genome	= "/sw/data/iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa"

// Modules
include { 
	FASTP;
	SALMON_INDEX;
	SALMON_QUANT;
	CREATE_TX2GENE;
	SALMON_TXIMPORT;
	CREATE_DESEQ
} from "/proj/sens2022505/analysis/src/modules/functions" 

workflow{

	fq_ch = Channel
		.fromFilePairs(["/proj/sens2022505/raw_data/**/*R{1,2}*fastq.gz"])
		.map{ it -> [it[0].split("_S")[0], it[1]]}
		.take( 5 )

	FASTP(fq_ch)
	SALMON_INDEX(params.txs,
		params.genome)
	SALMON_QUANT(FASTP.out.clean_fq,
		SALMON_INDEX.out)
	CREATE_TX2GENE(params.txs)
	SALMON_TXIMPORT(CREATE_TX2GENE.out,
		SALMON_QUANT.out.collect(),
		params.meta)
	CREATE_DESEQ(params.meta,
		SALMON_TXIMPORT.out)	

}
