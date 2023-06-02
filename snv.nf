/*
 * Somatic Variant Calling Routine
 * 
 */

params.correctedBam = "SRR23350360.aligned.merged.marked.sorted.corrected.bsqr.bam"
params.referencefasta = "/home/jon/ref/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"
params.outdir = "/home/jon/samples/dev"
params.referenceVcf = "/home/jon/ref/gnomad.genomes.v3.1.2.sites.chr17.vcf.bgz"
params.referenceVcfIndex = "/home/jon/ref/gnomad.genomes.v3.1.2.sites.chr17.vcf.bgz.tbi"

process INDEX {
	input:
	path input
	
	output:
	path "${input}.bai"
	
	script:
	"""
	samtools index ${input} -o ${input}.bai
	"""
}

process MUTECT {
	publishDir "$params.outdir/mutect_out", mode: 'copy'
	
    input:
    path input
	path index
	
    output:
    path "${input.simpleName}.vcf.gz"

    script:
    """
	gatk Mutect2 \
	-R ${params.referencefasta} \
	-I ${input} \
	-O "${input.simpleName}.vcf.gz" \
	--tmp-dir ~/samples/dev/temp \
	--java-options "-Xmx16G" 
    """
}

process GETPILEUPSUMMARIES {
    input:
    path input
	path refVcf
	path refVcfIndex
	path inputIndex
	
    output:
    path "${input.simpleName}.pileups.table"

    script:
    """
	gatk GetPileupSummaries -I ${input} -V ${refVcf} -L ${refVcf} -O "${input.simpleName}.pileups.table"
    """
}

process GETCONTAMINATIONTABLE {
    input:
    path pileupTable
	
    output:
    path "${pileupTable.simpleName}.contamination.table"

    script:
    """
	gatk CalculateContamination \
		-I ${pileupTable} \
		-O "${pileupTable.simpleName}.contamination.table"
    """
}

workflow CALCULATECONTAMINATION {
	take:
		inputBam
		input_index_ch
	main:
		pileups_table_ch = GETPILEUPSUMMARIES(inputBam, params.referenceVcf, params.referenceVcfIndex, input_index_ch)
		contamination_table_ch = GETCONTAMINATIONTABLE(pileups_table_ch)
	emit:
		contamination_table_ch
}

workflow {
	main:
		input_index_ch = INDEX(params.correctedBam)
		MUTECT(params.correctedBam, input_index_ch)
		CALCULATECONTAMINATION(params.correctedBam, input_index_ch)
}