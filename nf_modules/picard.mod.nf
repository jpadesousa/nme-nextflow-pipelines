#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MARK_DUPLICATES{	

	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	input:
		tuple val(name), path(reads)
		val(outputdir)
		val(mark_duplicates_args)
		val(verbose)

	output:
		tuple val(name), path("*fastqc*"), emit: all
	 	publishDir "$outputdir", mode: "link", overwrite: true

	script:
		if (verbose){
			println ("[MODULE] MARK_DUPLICATES ARGS: " + mark_duplicates_args)
		}

		"""
		module load picard
		picard MarkDuplicates ${reads}
		"""
}