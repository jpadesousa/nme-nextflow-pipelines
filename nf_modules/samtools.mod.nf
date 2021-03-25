#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// setting if the bam file should be published
params.bam_output = true

process SAMTOOLS_SORT{	
    
	label 'samtools'
	tag "$bam" // Adds name to job submission instead of (1), (2) etc.

	input:
		path(bam)
		val(outputdir)
		val(samtools_sort_args)
		val(verbose)

	output:
		path "*bam", emit: bam

		publishDir "$outputdir/aligned/bam", mode: "link", overwrite: true, enabled: params.bam_output

    script:
		samtools_sort_options = samtools_sort_args
		
		if (verbose){
			println ("[MODULE] SAMTOOLS SORT ARGS: " + samtools_sort_args)
		}

		output_name = bam.toString() - ".bam"

		"""
		module load samtools
		samtools sort $samtools_sort_options $bam -o ${output_name}.sorted.bam
    	"""
}

process SAMTOOLS_INDEX{	
    
	label 'samtools'
	tag "$bam" // Adds name to job submission instead of (1), (2) etc.

	input:
		path(bam)
		val(outputdir)
		val(samtools_index_args)
		val(verbose)

	output:
		path "*.bai", emit: bai

		publishDir "$outputdir/aligned", mode: "link", overwrite: true

    script:
		samtools_index_options = samtools_index_args
		
		if (verbose){
			println ("[MODULE] SAMTOOLS INDEX ARGS: " + samtools_index_args)
		}
		
		"""
		module load samtools
		samtools index $samtools_index_options $bam
		"""
}