#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process BISMARK2BEDGRAPH {

    label 'bismark2bedGraph'
	tag "$context_files" // Adds name to job submission instead of (1), (2) etc.

	input:
		path(context_files)
		val(outputdir)
		val(bismark2bedgraph_args)
		val(verbose)

	output:
		path "*{.bedgraph.gz,.bedgraph}",   emit: bedgraph
        path "*{.cov.gz,.cov}",             emit: coverage
		publishDir "$outputdir", mode: "link", overwrite: true

	script:
		outfile_basename = context_files.toString()  // Important to convert nextflow.processor.TaskPath object to String first
		outfile_basename = (outfile_basename - ~/.txt.gz$/)
        outfile_basename = (outfile_basename - ~/.txt$/)


		if (verbose){
			println ("[MODULE] BISMARK2BEDGRAPH ARGS: " + bismark2bedgraph_args)
		}

		"""
		module load bismark
		bismark2bedGraph --output ${outputdir}/bismark2bedGraph/${outfile_basename} ${context_files}
		"""
}