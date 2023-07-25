nextflow.enable.dsl=2


// Initialise variables


def help_message() {
  log.info """
    This pipeline was developed to generate consensus sequences of Mpox virus from Oxford Nanopore amplicon sequencing
    files. Here we combine tools such as chopper, minimap2, samtools consensus, and medaka to generate a polished
    final consensus file.

    Usage:
    . nextflow mpox.nf --reference reference_genome.fa --fastq ONT_reads.fastq --primer_scheme primer.scheme.bed --threads [NUM]

    Required Parameters:
    --reference --> Specify the reference_genome fasta file
    --fastq --> Specify the ONT read fastq file(s)
    --sra --> Specify SRA Project / Sample Accession numbers e.g 'PRJAA000001'
    --primer_scheme --> Specify the amplicon primer scheme bed file
    --threads --> Specify the max threads that the pipeline can use

    Optional Parameters:
    --intermediate --> Specify whether intermediate files are saved, such as corrected read and alignment files [FALSE]
    --help --> Prints out this help message [FALSE]
    --ncbi_key --> Provide NCBI API key if you want to access SRA files
    --[qc, draft, stats, align, med]_dir ––> Provide other directory paths for quality controlled reads, alignment,
      pipeline statistics and Medaka outputs

    """
}

// Show help message
if (params.help) {
    help_message()
    exit 0
}


// This process takes the raw read files and applies size and quality filtering
process read_qc {
    // thread allocation
    cpus params.threads

    // Output Dir
    publishDir "${params.qc_dir}/${file(sample).getSimpleName()}", pattern: "${sample_id}.fasta*"

    input:
      tuple val(sample_id), path(sample), path(scheme), path(reference)

    output:
      tuple val(sample_id), path('*fasta*'), path(scheme), path(reference), path('*.txt')

    script:
    if(sample_id == 'blank'){ sample_id = file(sample).getSimpleName()}
    """
    # Running NanoStat to extract the median read quality and sets chopper q threshold to 90% of the median quality value
    q_threshold=\$(nanoq -i ${sample} -sv | grep 'Median read quality:' | awk '{ threshold = \$4 * 0.9; printf "%.1f\\n", threshold }')

    # Setting size filtering limits
    amplicon_len=\$(awk 'NR % 2 == 1 {start = \$2} NR % 2 == 0 {end = \$3; len = end - start; sum += len; count++} END {if (count > 0) printf "%.0f\\n", sum / count}' ${scheme})
    min_chop=\$((amplicon_len - 500))
    max_chop=\$((amplicon_len + 500))

    if [[ ${sample} == *.fastq.gz ]]; then
        nanoq -l "\$min_chop" -m "\$max_chop" -q "\$q_threshold" -i ${sample} -svv -r "${sample_id}.txt" -O g -o "${sample_id}.fasta.gz"
    else
        nanoq -l "\$min_chop" -m "\$max_chop" -q "\$q_threshold" -i ${sample} -svv -r "${sample_id}.txt" -O u -o "${sample_id}.fasta"
    fi
    """
}

process mini_sam {
    // thread allocation
    cpus params.threads

    // Output Dir
    publishDir "${params.draft_dir}/${file(sample).getSimpleName()}", pattern: "mini_sam.fasta", mode: "copy"
    publishDir "${params.draft_dir}/${file(sample).getSimpleName()}", pattern: "mini_sam.bam"
    publishDir "${params.stats_dir}/${file(sample).getSimpleName()}", pattern: "mini_sam_stats.txt"

    input:
      tuple val(sample_id), path(sample), path(scheme), path(reference), path(stats)
      
    output:
      path '*.fasta', emit: fasta
      path '*.bam'
      path '*.txt', optional: true

    script:
    if(sample_id == null){ sample_id = file(sample).getSimpleName()}
    """
    # Alignment step

    minimap2 -a -t 3 -x map-ont ${reference} ${sample} | \

    # Samtools sorting, clipping primers, and converting to bam format
    samtools ampliconclip -u --both-ends -f clip_stats.txt -b ${scheme} - | \
    samtools sort -u -O bam -@ ${task.cpus} -o mini_sam.bam
    
    # Appending amplicon clip and alignment stats
    cp ${stats} mini_sam_stats.txt
    echo >> mini_sam_stats.txt
    cat clip_stats.txt >> mini_sam_stats.txt
    samtools stats temp/aln.bam -c 1,1000,25 | awk '
    /^SN/ { if (!sn_flag) { print "Summary Numbers"; sn_flag = 1 } print substr(\$0, index(\$0, \$2)); next }
    /^ID/ { if (!id_flag) { print "\\nIndel Distribution"; id_flag = 1 } print substr(\$0, index(\$0, \$2)); next }
    /^COV/ { if (!cov_flag) { print "\\nCoverage Distribution (increments of 25)"; cov_flag = 1 } print substr(\$0, index(\$0, \$2)); next }' >> mini_sam_stats.txt

    # Draft consensus to kick off Medaka
    samtools consensus -o mini_sam.fasta mini_sam.bam
    """
}

process mini_vc {
    // thread allocation
    cpus params.threads

    // Output Dir
    publishDir "${params.draft_dir}/${file(sample).getSimpleName()}", pattern: "mini_vc.fasta", mode: "copy"
    publishDir "${params.draft_dir}/${file(sample).getSimpleName()}", pattern: "mini_vc.bam"
    publishDir "${params.stats_dir}/${file(sample).getSimpleName()}", pattern: "mini_vc_stats.txt"
    
    input:
      tuple val(sample_id), path(sample), path(scheme), path(reference), path(stats)

    output:
      path '*.fasta', emit: fasta
      path '*.bam'
      path '*.txt', optional: true

    script:
    if(sample_id == null){ sample_id = file(sample).getSimpleName()}
    """
    # Alignment step

    minimap2 -a -t 3 -x map-ont ${reference} ${sample} | \

    # Samtools sorting, clipping primers, and converting to bam format
    samtools ampliconclip -u --both-ends -f clip_stats.txt -b ${scheme} - | \
    samtools sort -u -O bam -@ ${task.cpus} -o mini_vc.bam
    
    # Appending amplicon clip and alignment stats
    cp ${stats}  mini_vc_stats.txt
    echo >> mini_vc_stats.txt
    cat clip_stats.txt >> mini_vc_stats.txt
    samtools stats temp/aln.bam -c 1,1000,25 | awk '
    /^SN/ { if (!sn_flag) { print "Summary Numbers"; sn_flag = 1 } print substr(\$0, index(\$0, \$2)); next }
    /^ID/ { if (!id_flag) { print "\\nIndel Distribution"; id_flag = 1 } print substr(\$0, index(\$0, \$2)); next }
    /^COV/ { if (!cov_flag) { print "\\nCoverage Distribution (increments of 25)"; cov_flag = 1 } print substr(\$0, index(\$0, \$2)); next }' >> mini_vc_stats.txt

    # Draft consensus to kick off Medaka
    viral_consensus -i mini_vc.bam -r ${reference} -o mini_vc.fasta
    """
}


process ngmlr_sam {
    // thread allocation
    cpus params.threads

    // Output Dir
    publishDir "${params.draft_dir}/${file(sample).getSimpleName()}", pattern: "ngmlr_sam.fasta", mode: "copy"
    publishDir "${params.draft_dir}/${file(sample).getSimpleName()}", pattern: "ngmlr_sam.bam"
    publishDir "${params.stats_dir}/${file(sample).getSimpleName()}", pattern: "ngmlr_sam_stats.txt"

    input:
      tuple val(sample_id), path(sample), path(scheme), path(reference), path(stats)

    output:
      path '*.fasta', emit: fasta
      path '*.bam'
      path '*.txt', optional: true

    script:
    if(sample_id == null){ sample_id = file(sample).getSimpleName()}
    """
    # Alignment step

    ngmlr -x ont -t ${task.cpus} --bam-fix -r ${reference} -q ${sample} | \

    samtools ampliconclip -u --both-ends -f clip_stats.txt -b ${scheme} - | \

    samtools sort -u -O bam -@ ${task.cpus} -o ngmlr_sam.bam
    
    # Appending amplicon clip and alignment stats
    cp ${stats}  ngmlr_sam_stats.txt
    echo >> ngmlr_sam_stats.txt
    cat clip_stats.txt >> ngmlr_sam_stats.txt

    samtools stats temp/aln.bam -c 1,1000,25 | awk '
    /^SN/ { if (!sn_flag) { print "Summary Numbers"; sn_flag = 1 } print substr(\$0, index(\$0, \$2)); next }
    /^ID/ { if (!id_flag) { print "\\nIndel Distribution"; id_flag = 1 } print substr(\$0, index(\$0, \$2)); next }
    /^COV/ { if (!cov_flag) { print "\\nCoverage Distribution (increments of 25)"; cov_flag = 1 } print substr(\$0, index(\$0, \$2)); next }' >> ngmlr_sam_stats.txt

    samtools consensus -o ngmlr_sam.fasta ngmlr_sam.bam
    """
}

process ngmlr_vc {
    // thread allocation
    cpus params.threads

    // Output Dir
    publishDir "${params.draft_dir}/${file(sample).getSimpleName()}", pattern: "ngmlr_vc.fasta", mode: "copy"
    publishDir "${params.draft_dir}/${file(sample).getSimpleName()}", pattern: "ngmlr_vc.bam"
    publishDir "${params.stats_dir}/${file(sample).getSimpleName()}", pattern: "ngmlr_vc_stats.txt"

    input:
      tuple val(sample_id), path(sample), path(scheme), path(reference), path(stats)

    output:
      path '*.fasta', emit: fasta
      path '*.bam'
      path '*.txt', optional: true

    script:
    if(sample_id == null){ sample_id = file(sample).getSimpleName()}
    """
    # Alignment step

    ngmlr -x ont -t ${task.cpus} --bam-fix -r ${reference} -q ${sample} | \

    samtools ampliconclip -u --both-ends -f clip_stats.txt -b ${scheme} - | \

    samtools sort -u -O bam -@ ${task.cpus} -o ngmlr_vc.bam

    # Appending amplicon clip and alignment stats
    cp ${stats}  ngmlr_vc_stats.txt
    echo >> ngmlr_vc_stats.txt
    cat clip_stats.txt >> ngmlr_vc_stats.txt
    samtools stats temp/aln.bam -c 1,1000,25 | awk '
    /^SN/ { if (!sn_flag) { print "Summary Numbers"; sn_flag = 1 } print substr(\$0, index(\$0, \$2)); next }
    /^ID/ { if (!id_flag) { print "\\nIndel Distribution"; id_flag = 1 } print substr(\$0, index(\$0, \$2)); next }
    /^COV/ { if (!cov_flag) { print "\\nCoverage Distribution (increments of 25)"; cov_flag = 1 } print substr(\$0, index(\$0, \$2)); next }' >> ngmlr_vc_stats.txt

    viral_consensus -i ngmlr_vc.bam -r ${reference} -o ngmlr_vc.fasta
    """
}


// This process will generate a statistics text file using Nanostats on amplicons, raw reads, filtered
// reads, and the alignment.


process medaka_1 {

    input:
      path raw_reads
      path draft

    output:
      tuple path(raw_reads), path("medaka/*.fasta")

    script:
    def sample_id = file(raw_reads).getSimpleName()
    """
    medaka_consensus -i ${raw_reads} -d ${draft} -t 2
    """
}

process medaka_2 {

    input:
      tuple path(raw_reads), path(med_one)

    output:
      tuple path(raw_reads), path("medaka/*.fasta")

    script:
    def sample_id = file(raw_reads).getSimpleName()
    """
    medaka_consensus -i ${raw_reads} -d ${med_one} -t 2
    """
}

process medaka_3 {
    // Output Dir
    publishDir "${params.med_dir}", pattern: "${file(sample).getSimpleName()}/*.fasta", mode: "move"

    input:
      tuple path(raw_reads), path(med_two)

    output:
      path "${sample_id}/*.fasta"

    script:
    def sample_id = file(raw_reads).getSimpleName()
    """
    medaka_consensus -i ${raw_reads} -d ${med_two} -o ${sample_id}/ -t 2
    """
}



workflow align {
  // Setting channels
  scheme=Channel.fromPath(params.primer_scheme, checkIfExists: true)
  reference=Channel.fromPath(params.reference, checkIfExists: true)

  // these two options creates a channel depending on if locally stored fastq files
  // are supplied or if files are downloaded from the NCBI SRA database.

  // the fastq condition creates a tuple that has a blank id value and the fastq file path,
  // the reason for the blank is that the SRA channel automatically creates a tuple with
  // an ID and SRA file path. The blank id for fastq maintains the input cardinality for
  // the processes.

  if(params.fastq){
    sample=Channel.fromPath(params.fastq, checkIfExists: true)
    sample_tuple=Channel.value('blank').combine(sample)}

  if(params.sra){
      sample_tuple=Channel.fromSRA(params.sra, apiKey: params.ncbi_key)}

  // raw read tuple so that the scheme and reference file gets consumed with every sample id and file tuple
  raw_samples=sample_tuple.combine(scheme).combine(reference)


  // This set of pipped processes takes the raw_samples tuple and pipes it through the fastq
  // correcting process and then output results to each alignment process.
  raw_samples | read_qc | (mini_sam & mini_vc & ngmlr_sam & ngmlr_vc)

  emit:
    mini_sam.out.fasta
}

workflow medaka {
  // takes the output from the align workflow
  take: draft

  main:
    // This workflow runs the draft alignment through medaka polishing three times and
    // then output is saved in the working directory in a directory called Medaka results

    medaka_1(raw_reads=Channel.fromPath(params.fastq), draft) | medaka_2 | medaka_3
}

// Overall workflow that runs through both sub-workflows
workflow {
  align()
  medaka(align.out)
}
