#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
=================================
Analysis for #6026: RNAseq Kidney
=================================
*/

// Parameters
params.outdir	= "/proj/sens2022505/nobackup/custom_analysis/results/"
params.meta	= "/proj/sens2022505/analysis/db/metadata.csv"
params.txs	= "/proj/sens2022505/analysis/db/Homo_sapiens.GRCh38.cdna.all.fa"
params.genome	= "/sw/data/iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa"
params.index	= "/sw/data/iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/STARIndex/"
params.junctions	= "/sw/data/iGenomes/Homo_sapiens/NCBI/GRCh38/Sequence/STARIndex//sjdbList.fromGTF.out.tab"
params.gtf	= "/sw/data/iGenomes/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf"

// Modules
include { 
	FASTQC;
	FASTQC_TRIM;
	TRIMGALORE;
	MULTIQC;
	STAR_ALIGN;
	SALMON_INDEX;
	SALMON_QUANT;
	CREATE_TX2GENE;
	CREATE_META;
	SALMON_TXIMPORT;
	QC_DESEQ;
	QC_DESEQ_HM;
	RUN_DESEQ;
	EXTRACT_RES;
	CREATE_DESEQ } from "/proj/sens2022505/analysis/src/modules/functions" 

workflow{

	Channel
		.fromFilePairs(["/proj/sens2022505/raw_data/**/*R{1,2}*fastq.gz"])
		.map{ it -> [it[0].split("_S")[0], it[1]]}
		.multiMap{ it ->
			reads: it
			ids: it[0]
			}
		.set{ fq_ch }

	FASTQC(fq_ch.reads)
	TRIMGALORE(fq_ch.reads)
	FASTQC_TRIM(TRIMGALORE.out.reads)
	MULTIQC(FASTQC.out.zip.collect({it[1]}),
		FASTQC_TRIM.out.zip.collect({it[1]}))
	STAR_ALIGN(TRIMGALORE.out.reads,
		params.index,
		params.gtf,
		params.junctions)
	SALMON_INDEX(params.txs,
		params.genome)
	SALMON_QUANT(TRIMGALORE.out.reads,
		SALMON_INDEX.out)
	CREATE_TX2GENE(params.txs)
	CREATE_META(fq_ch.ids.collectFile(name: "ids.csv", newLine: true),
		params.meta)
	SALMON_TXIMPORT(CREATE_TX2GENE.out,
		SALMON_QUANT.out.collect(),
		CREATE_META.out)
	CREATE_DESEQ(CREATE_META.out,
		SALMON_TXIMPORT.out)

	gene_nos = [250,500,1000,2000]
	QC_DESEQ(CREATE_DESEQ.out,
		gene_nos)

	QC_DESEQ_HM(CREATE_DESEQ.out)
	RUN_DESEQ(CREATE_DESEQ.out)

	treat = Channel.of([1, "Glomerular_IgAN"], [2, "RoK_IgAN"])
	control = Channel.of([1, "Glomerular_LD"], [2, "RoK_LD"])
	treat.join(control).map{a,b,c -> tuple(b + "_vs_" + c, b, c)}.combine(RUN_DESEQ.out).set{ contrasts }

	EXTRACT_RES(contrasts)
	
}

