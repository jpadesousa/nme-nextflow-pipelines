#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MACS_CALLPEAK {	
    
	label 'macs_callpeak'
	tag "$bam" // Adds name to job submission instead of (1), (2) etc.

	input:
		val(bamList)
		val(outputdir)
		val(macs_callpeak_args)
		val(verbose)

	output:
		path "*bam", emit: bam
		publishDir "$outputdir", mode: "link", overwrite: true

    script:
		
		if (verbose){
			println ("[MODULE] MACS CALLPEAK ARGS: " + macs_callpeak_args)
		}

        // Effective genome size
        if (params.genome["name"] == 'GRCh37' || params.genome["name"] == 'GRCh38') {
			gsize = 'hs'
		}

        if (params.genome["name"] == 'GRCm38') {
			gsize = 'mm'
		}

        if (params.genome["name"] == 'WBcel235') {
			gsize = 'ce'
		}

        if (params.genome["name"] == 'BDGP6') {
			gsize = 'dm'
		}

        treatment = bamList.treatment
        control   = bamList.control

		"""
		module load macs
		macs3 callpeak -t ${treatment} -c ${control} -g ${gsize} $macs_callpeak_args
    	"""
}
