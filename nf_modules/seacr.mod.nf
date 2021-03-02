#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SEACR {	
    
	tag "$bedgraph" // Adds name to job submission instead of (1), (2) etc.

	input:
		path(bedgraph)
		val(outputdir)
		val(seacr_args)
		val(verbose)

	output:
		path "*bed", emit: bed
		publishDir "$outputdir", mode: "link", overwrite: true

    script:

		if (verbose){
			println ("[MODULE] SEACR ARGS: " + seacr_args)
		}

		if (seacr_args =~ /.*\.bedgraph.*|.*\.bg.*/){
			output_suffix = ".IgG"
		} else {
			if (seacr_args =~ /.*\d.*/){
            	output_suffix = ".FDR"
			} else {
				seacr_args = " 0.01 " + seacr_args
            	output_suffix = ".FDR_0.01"
			}
        }

        if(!(seacr_args =~ /.*norm.*|.*non.*/)){
            seacr_args += " norm "
        }

        if(!(seacr_args =~ /.*relaxed.*|.*stringent.*/)){
            seacr_args += " stringent "
        }

        bedgraph_name = bedgraph.toString() - ".bedgraph"

		"""
		module load seacr
        seacr ${bedgraph} ${seacr_args} "${bedgraph_name}${output_suffix}"
    	"""
}