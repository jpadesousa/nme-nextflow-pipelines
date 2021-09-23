#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    PROCESSES
======================================================================================== */
process HISAT2 {

	label 'hisat2'
	tag "$name" // Adds name to job submission instead of (1), (2) etc.

	input:
		tuple val(name), path(reads)
		val(outputdir)
		val(hisat2_args)
		val(verbose)

	output:
		path "*bam",       emit: bam
		path "*stats.txt", emit: stats
		val(single_end)  , emit: single_end
		
		publishDir "$outputdir/aligned/bam",  mode: "link", overwrite: true, pattern: "*bam"
		publishDir "$outputdir/aligned/logs", mode: "link", overwrite: true, pattern: "*stats.txt"

	script:
		// Verbose
		if (verbose){
			println ("[MODULE] HISAT2 ARGS: " + hisat2_args)
		}

		// Default options
		hisat2_args = hisat2_args + " --no-unal --no-softclip "

		// File names
		readString = ""
		if (reads instanceof List) {
			readString  = "-1 " + reads[0] + " -2 " + reads[1]
			hisat2_args = hisat2_args + " --no-mixed --no-discordant"
			single_end  = false
		}
		else {
			readString  = "-U " + reads
			single_end  = true
		}

		// Index
		index = params.genome["hisat2"]

		// Splices
			// TODO: need to add a check if the splice-site infile is present or not, and leave out this parameter otherwise
		splices = " --known-splicesite-infile " + params.genome["hisat2_splices"]

		// Basename
		hisat_name = name + "_" + params.genome["name"]

		"""
		module load hisat2 samtools

		hisat2 -p ${task.cpus} ${hisat2_args} -x ${index} ${splices} ${readString}  2>${hisat_name}_hisat2_stats.txt | samtools view -bS -F 4 -F 8 -F 256 -> ${hisat_name}_hisat2.bam
		"""
}
