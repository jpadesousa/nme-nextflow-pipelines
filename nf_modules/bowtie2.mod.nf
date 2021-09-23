#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.cutntag    = false
params.bam_output = true // setting if the bam file should be published


/* ========================================================================================
    PROCESSES
======================================================================================== */
process BOWTIE2 {

		label 'bowtie2'
		tag "$name" // Adds name to job submission instead of (1), (2) etc.

	input:
		tuple val(name), path(reads)
		val(outputdir)
		val(bowtie2_args)
		val(verbose)

	output:
		path "*bam",  	   emit: bam
		path "*stats.txt", emit: stats

		publishDir "$outputdir/aligned/bam",  mode: "link", overwrite: true, pattern: "*bam", enabled: params.bam_output
		publishDir "$outputdir/aligned/logs", mode: "link", overwrite: true, pattern: "*stats.txt"

	script:
		// Verbose
		if (verbose){
			println ("[MODULE] BOWTIE2 ARGS: " + bowtie2_args)
		}

		// Default options
		bowtie2_args +=  " --no-unal " // We don't need unaligned reads in the BAM file

		// File names
		readString = ""
		if (reads instanceof List) {
			readString = "-1 " + reads[0] + " -2 " + reads[1]
			bowtie2_args += " --no-discordant --no-mixed " // just output properly paired reads
		}
		else {
			readString = "-U " + reads
		}

		// CUT&Tag
		if (params.cutntag) {
			bowtie2_args += " --local --very-sensitive --phred33 -I 10 -X 700 "
		}

		// Index
		index = params.genome["bowtie2"]

		// Basename
		bowtie_name = name + "_" + params.genome["name"]

		"""
		module load bowtie2 samtools

		bowtie2 -x ${index} -p ${task.cpus} ${bowtie2_args} ${readString}  2>${bowtie_name}_bowtie2_stats.txt | samtools view -bS -F 4 -F 8 -F 256 -> ${bowtie_name}_bowtie2.bam
		"""
}
