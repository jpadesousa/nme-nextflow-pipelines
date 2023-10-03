#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.verbose  = true


/* ========================================================================================
    PROCESSES
======================================================================================== */
process SRADOWNLOADER {

	tag "$geoaccession" // Adds name to job submission
    label 'sradownloader'

	input:
	  	file(sra_metadata)
        val(geoaccession)
		val(outputdir)
		val(sradownloader_args)
		val(verbose)

	output:
	  	path "*gz", emit: fastq
		publishDir "$outputdir/fastq", mode: "copy", overwrite: true

	script:

		/* ==========
			Verbose
		========== */
		if (verbose){
			println ("[MODULE] SRADOWNLOADER ARGS: " + sradownloader_args)
		}

		"""
		module load sradownloader/3.8

		sradownloader_axel ${sradownloader_args} ${sra_metadata} 
		"""
}
