#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    PROCESSES
======================================================================================== */
process BEDTOOLS_GENOMECOV{	

	tag "$bam" // Adds name to job submission instead of (1), (2) etc.

	input:
		path(bam)
		val(outputdir)
		val(bedtools_genomecov_args)
		val(verbose)

	output:
		path "*bedgraph", emit: bedgraph

		publishDir "$outputdir/aligned/bedgraph", mode: "link", overwrite: true

    script:
	    // Verbose
		if (verbose){
			println ("[MODULE] BEDTOOLS GENOMECOV ARGS: " + bedtools_genomecov_args)
		}

		// bedtools genomecov parameters for the CUT&Tag pipeline
        if (params.cutntag) {
			bedtools_genomecov_args += " -bga " 
		}

		"""
		module load bedtools2

		bedtools genomecov ${bedtools_genomecov_args} -ibam ${bam} > ${bam}.bedgraph

		rename .bam.bedgraph .bedgraph *
    	"""
}
