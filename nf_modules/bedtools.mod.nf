nextflow.enable.dsl=2

process BEDTOOLS_GENOMECOV{	
    
	label 'bedtools'
	
	tag "$bam" // Adds name to job submission instead of (1), (2) etc.

	input:
		path(bam)
		val (outputdir)
		val (bedtools_genomecov_args)
		val (verbose)

	output:
		// path "*report.txt", emit: report
		path "*bam",        emit: bam

		publishDir "$outputdir",
		mode: "link", overwrite: true

    script:
		bedtools_genomecov_options = bedtools_genomecov_args
		
		if (verbose){
			println ("[MODULE] BEDTOOLS GENOMECOV ARGS: " + bedtools_genomecov_args)
		}

        if (params.cutntag) {
			bedtools_genomecov_options += " -ibam -bga " // Bedtools genomecov parameters for the CUT&Tag pipeline
		}

        //fasta = params.genome["fasta"] + "/*.fa"

		"""
		module load bedtools2

		bedtools genomecov $bedtools_genomecov_options $bam
    	"""
}