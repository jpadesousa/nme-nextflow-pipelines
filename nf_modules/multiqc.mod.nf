#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MULTIQC {

	label 'multiQC'

	input:
		path(file)
		val(outputdir)
		val(multiqc_args)
		val(verbose)

	output:
		path "*html", emit: html
		
		publishDir "$outputdir", mode: "link", overwrite: true

	script:

		if (verbose){
			println ("[MODULE] MULTIQC ARGS: " + multiqc_args)
		}

		"""
		module load multiqc
		multiqc $multiqc_args -x work .
		"""

}
