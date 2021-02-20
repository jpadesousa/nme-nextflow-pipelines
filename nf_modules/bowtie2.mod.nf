#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// parameters passed in by specialised pipelines
params.cutntag = false

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
		publishDir "$outputdir", mode: "link", overwrite: true

	script:
		readString = ""
		
		if (verbose){
			println ("[MODULE] BOWTIE2 ARGS: " + bowtie2_args)
		}

		// Options we add are
		bowtie_options = bowtie2_args
		bowtie_options +=  " --no-unal " // We don't need unaligned reads in the BAM file

		if (reads instanceof List) {
			readString = "-1 " + reads[0] + " -2 " + reads[1]
			bowtie_options += " --no-discordant --no-mixed " // just output properly paired reads
		}
		else {
			readString = "-U " + reads
		}

		if (params.cutntag) {
			bowtie_options += " --local --very-sensitive --phred33 -I 10 -X 700 " // Bowtie parameters for the CUT&Tag pipeline
		}

		index = params.genome["bowtie2"]
		bowtie_name = name + "_" + params.genome["name"]

		"""
		module load bowtie2 samtools
		bowtie2 -x ${index} -p ${task.cpus} ${bowtie_options} ${readString}  2>${bowtie_name}_bowtie2_stats.txt | samtools view -bS -F 4 -F 8 -F 256 -> ${bowtie_name}_bowtie2.bam
		"""

}
