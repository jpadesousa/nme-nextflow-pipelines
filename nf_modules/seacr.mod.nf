#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process SEACR {	

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


		if (bedgraph instanceof List) {
			files_command   = bedgraph[0] + " " + bedgraph[1]
			output_name     = bedgraph[0].toString() - ".bedgraph"
			output_suffix   = ""
			seacr_threshold = ""
		} else {
			files_command = bedgraph
			output_name   = bedgraph.toString() - ".bedgraph"

			if (seacr_args =~ /.*\d.*/){
				arr = seacr_args.split(" ", 2)
				seacr_threshold = arr[0]
            	output_suffix = ".FDR_" + arr[0]
			} else {
				seacr_threshold = "0.01"
           		output_suffix = ".FDR_0.01"
			}
		}



        if(!(seacr_args =~ /.*norm.*|.*non.*/)){
            seacr_normalization = "norm"
        } else {
			seacr_normalization = "non"
		}

        if(!(seacr_args =~ /.*relaxed.*|.*stringent.*/)){
            seacr_mode = "stringent"
        } else {
			seacr_mode = "relaxed"
		}



		"""
		module load seacr
        seacr ${files_command} ${seacr_threshold} ${seacr_normalization} ${seacr_mode} "${output_name}${output_suffix}"
    	"""
}