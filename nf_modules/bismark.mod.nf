#!/usr/bin/env nextflow
nextflow.enable.dsl=2


/* ========================================================================================
    DEFAULT PARAMETERS
======================================================================================== */
params.pbat 	     = false
params.unmapped   	 = false
params.singlecell    = ""
params.read_identity = ""


/* ========================================================================================
    PROCESSES
======================================================================================== */
process BISMARK {
	
	label 'bismark'
	tag "$name" // Adds name to job submission instead of (1), (2) etc.
		
    input:
	    tuple val(name), path(reads)
		val (outputdir)
		val (bismark_args)
		val (verbose)

	output:
	    tuple val(name), path ("*bam"), emit: bam
		path "*report.txt", 			emit: report
		
		// we always pass back the original name so we can use .join() later on, e.g. for bismark2bedGraph
		tuple val(name), path ("*unmapped_reads_1.fq.gz"), optional: true, emit: unmapped1
		tuple val(name), path ("*unmapped_reads_2.fq.gz"), optional: true, emit: unmapped2

		publishDir "$outputdir/aligned/bam", 	 mode: "link", overwrite: true, pattern: "*bam"
		publishDir "$outputdir/aligned/logs",    mode: "link", overwrite: true, pattern: "*report.txt"
		publishDir "$outputdir/unaligned/fastq", mode: "link", overwrite: true, pattern: "*.fq.gz"

    script:
		// Verbose
		if (verbose){
			println ("[MODULE] BISMARK ARGS: " + bismark_args)
		}

		// Single-cell
		if (params.singlecell){
			bismark_args += " --non_directional "
		}
		else {
		
		}

		// PBAT
		if (params.pbat){
			bismark_args += " --pbat "
		}

		// File names
		readString = ""
		if (reads instanceof List) {
			readString = "-1 " + reads[0] + " -2 " + reads[1]
		}
		else {
			readString = reads
		}

		// Index
		index = "--genome " + params.genome["bismark"]

		// Basename
			/// adds Genome build and aligner to output name	
		unmapped_name = ''	
		if (params.read_identity == "1" || params.read_identity == "2"){

			if (params.read_identity == "1"){
				unmapped_name = name + "_unmapped_R1"
			}
			else {
				unmapped_name = name + "_unmapped_R2"
			}

			if (bismark_args =~ /-hisat/){ // if HISAT2 was given on the command line
				bismark_name = unmapped_name + "_" + params.genome["name"] + "_bismark_ht2"
			}
			else { // default is Bowtie 2
				bismark_name = unmapped_name + "_" + params.genome["name"] + "_bismark_bt2"
			}
		}
		else {

			if (bismark_args =~ /-hisat/){ // if HISAT2 was given on the command line
				bismark_name = name + "_" + params.genome["name"] + "_bismark_ht2"
			}
			else { // default is Bowtie 2
				bismark_name = name + "_" + params.genome["name"] + "_bismark_bt2"
			}
		}	
		
		"""
		module load bismark

		bismark --parallel 1 --basename $bismark_name $index $bismark_args $readString
		"""
}