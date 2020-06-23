/*
 * Reading the input files for RNASeq pipeline
 */

params.reads = ""
params.genome = ""
params.gtf = ""
params.outputDir = ""

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { read1_ch; raw_reads_trimgalore; raw_reads_for_hisat2 }

/*
 * Generating QC reports for Raw reads
 */

process RawFASTQC {
	tag "FASTQC on $sample_id"
	publishDir "${params.outputDir}/RawFASTQC", mode: 'copy',
		saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename" }
    	input:
    	set sample_id, file(reads) from read1_ch

    	output:
	file "fastqc_${sample_id}_output" into fastqc_results_ch

    	script:
    	"""
		mkdir fastqc_${sample_id}_output
	        fastqc -f fastq -q ${params.reads} -o fastqc_${sample_id}_output
    	"""
}

/*
 * Generating hisat2 genome index 
 */
process BuildHISat2Index {
	tag "$genome.baseName"
	memory '8 GB'
	publishDir "${params.outputDir}/ReferenceIndex/",mode: 'copy', pattern: ""

	input:
	path genome from params.genome
	
	output:
	path 'genome.index*' into index_ch
	
	"""
	hisat2-build -f ${genome} genome.index
	"""
}

/*
 * Trimming raw reads
 */
process Trimming {
    publishDir "${params.outputDir}/TrimmedReads/", mode: 'copy', pattern: "*trimmed.fq.gz"
    publishDir "${params.outputDir}/TrimmedReads/", mode: 'copy', pattern: "*_trimming_report.txt"

    input:
	set val(name), file(reads) from raw_reads_trimgalore
	
	output:
    	set val(name), file("*.fq.gz") into trimgalore_reads
	set val(name), file("*.fq.gz") into trimgalore_reads_for_alignment
    	file "*trimming_report.txt" into trimgalore_results

	script:
	"""
        trim_galore --phred33 --gzip --stringency 3 --length 25 --trim-n ${params.reads}
        """
}

/*
 * Generating QC reports for trimmed reads
 */
process Trimmed_Reads_FASTQC {
	tag "FASTQC on $sample_id"
        publishDir "${params.outputDir}/TrimmedFASTQC", mode: 'copy',
		saveAs: { filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename" }
	input:
        set sample_id, path(reads) from trimgalore_reads
	println "trimgalore_reads"
        output:
        file "trimmed_fastqc_${sample_id}_output" into fastqc_trimmed_results_ch

        script:
        """
                mkdir trimmed_fastqc_${sample_id}_output
                fastqc -f fastq -q $reads -o trimmed_fastqc_${sample_id}_output
        """
}

/*
 * Generating alignment file using HISat2
 */
process Hisat2_align {
	tag "$name"
	cpus 2
	memory '8 GB'
	publishDir "${params.outputDir}/HISat2_Output/hisat2_mapstats", mode: 'copy', pattern: "*hisat2_mapstats.txt"
	publishDir "${params.outputDir}/HISat2_Output/", mode: 'copy', pattern: "*sam"	
	input:
	path genome from params.genome
	path index from index_ch
	set val(name), file(reads) from raw_reads_for_hisat2

	output:
	set val(name), file("*.sam") into hisat2_sam
	file("*hisat2_mapstats.txt") into hisat2_mapstats

	"""
	hisat2 -p 4 --very-sensitive -x genome.index -1 ${reads[0]} -2 ${reads[1]} --new-summary > ${name}.sam 2> ${name}.hisat2_mapstats.txt 
	"""
}

/*
 * Generating BAM file from SAM file
 * Sorting BAM file
 * Generating the index for BAM file
 */
process SAM2BAM {
        tag "$name"
        publishDir "${params.outputDir}/HISat2_Output/", mode: 'copy', pattern: "*sorted.bam"
        publishDir "${params.outputDir}/HISat2_Output/", mode: 'copy', pattern: "*sorted.bai"
        publishDir "${params.outputDir}/HISat2_Output/", mode: 'copy', pattern: "*flagstat"

        input:
        set val(name), file(mapped_sam) from hisat2_sam

        output:
        set val(name), file("${name}.sorted.bam") into sorted_bam_for_feature_count
        set val(name), file("${name}.sorted.bam.bai") into sorted_bam_indices_for_feature_count
        set val(name), file("${name}.flagstat") into bam_flagstat
	set val(name), file("${name}.sorted.bam") into sorted_bam_for_qualimap
        set val(name), file("${name}.sorted.bam.bai") into sorted_bam_indices_for_qualimap

        script:
        """
        samtools view -@ 8 -bS -o ${name}.bam ${mapped_sam}
        samtools sort ${name}.bam ${name}.sorted
        samtools index ${name}.sorted.bam ${name}.sorted.bam.bai
        samtools flagstat ${name}.sorted.bam > ${name}.flagstat
        """
}

/*
 * Generating gene count file from sorted BAM file
 */
process Feature_count{
        tag "$name"
        publishDir "${params.outputDir}/Feature_counts/" , mode: 'copy', pattern: "*.t*"

        input:
	path gtf from params.gtf
        set val(name), file(bam_file) from sorted_bam_for_feature_count
        file(bam_indices) from sorted_bam_indices_for_feature_count

        output:
	file "*.gene.featurecount.txt" into geneCounts
        file "*.gene.featurecounts.txt.summary" into featureCounts_logs

	script:
	"""
	featureCounts -F GTF -a ${gtf} -g gene_id -o ${name}.gene.featurecount.txt ${bam_file}
        """
 }

/*
 * Generating QC reports for sorted BAM file
 */
process QC_Using_Qualimap{
        tag "$name"
        publishDir "${params.outputDir}/Qualimap_QC/" , mode: 'copy', pattern: "*.pdf"

        input:
	path gtf from params.gtf
	set name, file(bam_file2) from sorted_bam_for_qualimap
        file(bam_indices2) from sorted_bam_indices_for_qualimap
	
        output:
        set val(name), file("*.pdf") into qualimap_output

        script:
        """
        qualimap rnaseq -s -a proportional -bam ${bam_file2} -p non-strand-specific -gtf ${gtf} -outformat PDF --java-mem-size=6G
        """
}
workflow.onComplete { 
	log.info ( workflow.success ? "Successfully Completed ..Done!" : "Something in the pipeline went wrong" )
}
