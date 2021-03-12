#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MARK_DUPLICATES{	

	label 'picard'
	tag "$bam" // Adds name to job submission instead of (1), (2) etc.

	input:
		path(bam)
		val(outputdir)
		val(mark_duplicates_args)
		val(verbose)

	output:
		path "*bam", emit: bam
		path "*txt", emit: metrics
	 	publishDir "$outputdir", mode: "link", overwrite: true

	script:
		if (verbose){
			println ("[MODULE] MARK_DUPLICATES ARGS: " + mark_duplicates_args)
		}

		"""
		module load picard

		picard MarkDuplicates INPUT=${bam} OUTPUT=${bam}.dedup.bam ASSUME_SORTED=true METRICS_FILE=${bam}.MarkDuplicates.metrics.txt ${mark_duplicates_args}

		rename .bam.dedup .dedup *
		rename .bam.MarkDuplicates.metrics.txt .MarkDuplicates.metrics.txt *
		"""
}