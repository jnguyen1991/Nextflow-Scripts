/*
 *
 * Preprocesssing routine for fastq data
 * Converts to bam and goes through quality preprocessing routine
 *
 */
params.reads = "/home/jon/samples/dev/SRR23350360_{1,2}.fastq.gz"
params.reference = "/home/jon/ref/"
params.refvcftable = "/home/jon/ref/hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf"
params.outdir = "/home/jon/samples/dev"
params.ref_folder = "/home/jon/ref/"

process BWAMEM {
    input:
    path ref_folder
	tuple val(sample_id), path(reads)
	
    output:
    path "${sample_id}.aligned.bam"

    script:
    """
    bwa mem "$ref_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"  ${reads[0]} ${reads[1]} > "${sample_id}.aligned.bam"
    """
}

process FASTQTOSAM {
	input:
	tuple val(sample_id), path(reads)
	
	output:
	path "${sample_id}.unaligned.bam"
	
	script:
	"""
	java -jar ~/bin/picard.jar FastqToSam F1=${reads[0]} F2=${reads[1]} O="${sample_id}.unaligned.bam" SM=${sample_id}
	"""
}

process MERGEBAM {	
	input:
	path aligned_bam
	path unaligned_bam
	path ref_folder
	
	output:
	path "${aligned_bam.baseName}.merged.bam"
	
	script:
	"""
	java -jar ~/bin/picard.jar MergeBamAlignment ALIGNED=${aligned_bam} UNMAPPED=${unaligned_bam} O="${aligned_bam.baseName}.merged.bam" R="$ref_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"
	"""
}

process MARKDUPLICATES {
	input:
	path input
	
	output:
	path "${input.baseName}.marked.bam"
	
	script:
	"""
	java -jar ~/bin/picard.jar MarkDuplicates \
      I=${input} \
      O="${input.baseName}.marked.bam" \
      M=marked_dup_metrics.txt
	"""
}

process SORTSAM {
	input:
	path input
	
	output:
	path "${input.baseName}.sorted.bam"
	
	script:
	"""
	java -jar ~/bin/picard.jar SortSam \
		  I=${input} \
		  O="${input.baseName}.sorted.bam" \
		  SORT_ORDER=coordinate
	"""
}

process BASERECALIBRATOR {
	input:
	path input
	path ref_folder
	
	output:
	path "${input.baseName}.recal_data.table"
	
	script:
	"""
	gatk BaseRecalibrator \
		-I ${input} \
		-R "$ref_folder/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta" \
		--known-sites ${params.refvcftable} \
		-O "${input.baseName}.recal_data.table"
	"""
}

process ADDREADGROUPINFO {
	input:
	path input

	output:
	path "${input.baseName}.corrected.bam"
	
	script:
	"""
	 java -jar ~/bin/picard.jar AddOrReplaceReadGroups \
		   I=${input} \
		   O="${input.baseName}.corrected.bam" \
		   RGLB=lib1 \
		   RGPL=ILLUMINA \
		   RGPU=unit1 \
		   RGSM=${input.simpleName}
	"""
}

process APPLYBSQR {
	publishDir "$params.outdir/bam_files", mode: 'copy'
	
	input:
	path input
	path recal_table

	output:
	path "${input.baseName}.bsqr.bam"
	
	script:
	"""
	gatk ApplyBQSR \
		-R "${params.ref_folder}/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta" \
		-I ${input} \
		--bqsr-recal-file ${recal_table} \
		-O "${input.baseName}.bsqr.bam"
	"""
}



workflow MAPTOREFERENCE {
	main:
		Channel
			.fromFilePairs(params.reads, checkIfExists: true)
			.set { read_pairs_ch }	
		aligned_bam_ch = BWAMEM(params.reference, read_pairs_ch)
		unaligned_bam_ch = FASTQTOSAM(read_pairs_ch)
		merged_bam_ch = MERGEBAM(aligned_bam_ch,unaligned_bam_ch,params.reference)
	emit:
		merged_bam_ch
}



workflow MARKDUPLICATESANDSORT {
	take: alignedbam
	main:
		marked_bam_ch = MARKDUPLICATES(alignedbam)
		sort_bam_ch = SORTSAM(marked_bam_ch)
	emit:
		sort_bam_ch
}

workflow BQSR {
	take: 
		sortedbam
	main:
		corrected_bam_ch = ADDREADGROUPINFO(sortedbam)
		recalibration_table_ch = BASERECALIBRATOR(corrected_bam_ch, params.reference)
		bsqr_bam_ch = APPLYBSQR(corrected_bam_ch, recalibration_table_ch)
}

workflow {
	main:
		MAPTOREFERENCE()
		MARKDUPLICATESANDSORT(MAPTOREFERENCE.out)
		BQSR(MARKDUPLICATESANDSORT.out)
}