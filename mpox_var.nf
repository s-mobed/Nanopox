nextflow.enable.dsl=2


// Initialise variables


def help_message() {
  log.info """
    This pipeline was developed to generate consensus sequences of Mpox virus from Oxford Nanopore amplicon sequencing
    files. Here we combine tools such as chopper, minimap2, samtools consensus, and medaka to generate a polished
    final consensus file.

    Usage:
    nextflow run s-mobed/Nanopox --reference reference_genome.fa --fastq ONT_reads.fastq --primer_scheme primer.scheme.bed --threads [NUM]

    Required Parameters:
    --reference --> Specify the reference_genome fasta file
    --fastq --> Specify the ONT read fastq file(s)
    --sra --> Specify SRA Project / Sample Accession numbers e.g 'PRJAA000001'
    --primer_scheme --> Specify the amplicon primer scheme bed file

    Optional Parameters:
    --threads --> Specify the max threads that the pipeline can use [1]
    --method --> Selects the alignment method to use, mini_sam is the default [mini_sam, mini_vc, ngmlr_sam, ngmlr_vc, all]
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

    // I/O
    input:
      tuple val(sample_id), path(sample), path(scheme), path(reference)

    output:
      tuple val(sample_id), path('*fasta*'), path(scheme), path(reference), path("${sample_id}.txt")

    script:
    if(sample_id == 'blank'){ sample_id = file(sample).getSimpleName()}
    """
    # Running Nanoq to get raw read statistics
    nanoq -i ${sample} -r "raw.txt" -svv

    # Running Nanoq to extract the median read quality and sets chopper q threshold to 50% of the median quality value
    q_threshold=\$(nanoq -i ${sample} -sv | grep 'Median read quality:' | awk '{ threshold = \$4 * 0.5; printf "%.1f\\n", threshold }')

    # Setting size filtering limits
    amplicon_len=\$(awk 'NR % 2 == 1 {start = \$2} NR % 2 == 0 {end = \$3; len = end - start; sum += len; count++} END {if (count > 0) printf "%.0f\\n", sum / count}' ${scheme})
    min_chop=\$((amplicon_len - 500))
    max_chop=\$((amplicon_len + 500))

    if [[ ${sample} == *.fastq.gz ]]; then
        nanoq -l "\$min_chop" -m "\$max_chop" -q "\$q_threshold" -i ${sample} -svv -r "filtered.txt" -O g -o "${sample_id}.fasta.gz"
    else
        nanoq -l "\$min_chop" -m "\$max_chop" -q "\$q_threshold" -i ${sample} -svv -r "filtered.txt" -O u -o "${sample_id}.fasta"
    fi

    # stitching stats for output
    echo "Raw read statistics" > ${sample_id}.txt

    cat raw.txt >> ${sample_id}.txt

    echo "Filtered read statistics" >> ${sample_id}.txt

    cat filtered.txt >> ${sample_id}.txt

    echo "Amplicon clipping statistics" >> ${sample_id}.txt
    """
}

// Alignment Processes
process minimap {
    // thread allocation
    cpus params.threads

    // I/O
    input:
      tuple val(sample_id), path(sample), path(scheme), path(reference), path(stats)
      
    output:
      tuple val(sample_id), path(sample), path('mini.sam'), path(scheme), path(reference), path(stats)

    script:
    if(sample_id == null){ sample_id = file(sample).getSimpleName()}
    """
    minimap2 -a -t 3 -x map-ont ${reference} ${sample} -o mini.sam
    """
}


// Amplicon primer clipping
process sam_clip {
    // thread allocation
    cpus params.threads

    // Output Dir
    publishDir "${params.draft_dir}/${sample_id}", pattern: "*.bam"
    publishDir "${params.stats_dir}/${sample_id}", pattern: "*_stats.txt"

    // I/O
    input:
      tuple val(sample_id), path(sample), path(sam_file), path(scheme), path(reference), path(stats)
      
    output:
      tuple val(sample_id), path(sample), path('*.bam'), path(reference)

    script:
    """
    # Samtools sorting, clipping primers, and converting to bam format
    samtools ampliconclip -u --both-ends -f clip_stats.txt -b ${scheme} ${sam_file} | \
    samtools sort -u -O bam -@ ${task.cpus} -o mini_sam.bam
    
    # Appending amplicon clip and alignment stats
    cp ${stats} mini_sam_stats.txt
    echo >> mini_sam_stats.txt
    cat clip_stats.txt >> mini_sam_stats.txt
    samtools stats *.bam -c 1,1000,25 | awk '
    /^SN/ { if (!sn_flag) { print "Summary Numbers"; sn_flag = 1 } print substr(\$0, index(\$0, \$2)); next }
    /^ID/ { if (!id_flag) { print "\\nIndel Distribution"; id_flag = 1 } print substr(\$0, index(\$0, \$2)); next }
    /^COV/ { if (!cov_flag) { print "\\nCoverage Distribution (increments of 25)"; cov_flag = 1 } print substr(\$0, index(\$0, \$2)); next }' >> mini_sam_stats.txt
    """
}

// Consensus Calling Processes

process sam_con {
    // thread allocation
    cpus params.threads

    // Output Dir
    publishDir "${params.draft_dir}/${sample_id}", pattern: "*.fasta", mode: "copy"

    // I/O
    input:
      tuple val(sample_id), path(sample), path(bam_file), path(reference)
      
    output:
      tuple val(sample_id), path(sample), path(bam_file), path(reference), path('*.fasta')

    script:
    """
    # Draft consensus to kick off Medaka
    samtools consensus -o ${sample_id}.fasta ${bam_file}
    """
}

// Polishing Processes

process medaka_1 {
    // thread allocation
    cpus params.threads

    // I/O
    input:
      tuple val(sample_id), path(sample), path(bam_file), path(reference), path(draft)

    output:
      tuple val(sample_id), path(sample), path(bam_file), path(reference), path("medaka/*.fasta")

    script:
    """
    medaka_consensus -i ${sample} -d ${draft} -t ${task.cpus}
    """
}

process medaka_2 {
    // thread allocation
    cpus params.threads

    // I/O
    input:
      tuple val(sample_id), path(sample), path(bam_file), path(reference), path(med_one)

    output:
      tuple val(sample_id), path(sample), path(bam_file), path(reference), path("medaka/*.fasta")

    script:
    """
    medaka_consensus -i ${sample} -d ${med_one} -t ${task.cpus}
    """
}

process medaka_3 {
    // thread allocation
    cpus params.threads

    // Output Dir
    publishDir "${params.med_dir}", pattern: "${file(sample).getSimpleName()}/*.fasta", mode: "copy"

    // I/O
    input:
      tuple val(sample_id), path(sample), path(bam_file), path(reference), path(med_two)

    output:
      tuple val(sample_id), path(bam_file), path(reference), path("${file(sample).getSimpleName()}/*.fasta")

    script:
    """
    medaka_consensus -i ${sample} -d ${med_two} -o ${sample_id}/ -t ${task.cpus}
    """
}

// Coverage Plot process
process plot {
    // I/O
    input:
      tuple val(sample_id), path(mini_bam), path(ngmlr_bam), path(primer_scheme)

    script:
    """
    coverage_amplicon_plot.sh -p ${primer_scheme} -m ${mini_bam} -n ${ngmlr_bam} -s ${sample_id}
    """
}

// Variant calling processes
process bcfcall {
    // thread allocation
    cpus params.threads

    // I/O
    input:
      tuple val(sample_id), path(sample), path(bam_file), path(reference), path(fasta)

    output:
      tuple val(sample_id), path('*.vcf'), path(reference)

    script:
    if(sample_id == null){ sample_id = file(sample).getSimpleName()}
    """
    bcftools mpileup -Ou -f ${reference} ${bam_file} | \
    bcftools call -mv --ploidy 1 -Ov -o ${sample_id}.calls.vcf
    """
}

process snpann {
    // thread allocation
    cpus params.threads

    // Output Dir
    publishDir "${params.draft_dir}/${sample_id}", pattern: "*.ann.vcf", mode: "copy"

    // I/O
    input:
      tuple val(sample_id), path(calls), path(reference)

    output:
      tuple val(sample_id), path('*.ann.vcf')

    script:
    if(sample_id == null){ sample_id = file(sample).getSimpleName()}
    """
    snpEff ann mpox -v ${calls} > ${sample_id}.calls.ann.vcf
    """
}

process vcf_parser {
    // thread allocation
    cpus params.threads

    // Output Dir
    publishDir "${params.draft_dir}/${sample_id}", pattern: "*.txt", mode: "copy"

    // I/O
    input:
      tuple val(sample_id), path(ann_calls)

    output:
      path '*.txt'

    script:
    template 'annotated_vcf_parser.py'
}


// Read Loading and Filtering Workflow

workflow qc {
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
  raw_samples | read_qc

  // Named outputs to send to alignment workflows
  emit:
    read_qc.out
}

// Alignment Workflow
workflow m2 {
  take:
    qc_output
    
  main:
    minimap(qc_output)

  // Named outputs to send to medaka workflows
  emit:
    minimap.out
}


// Consensus Workflow
workflow consensus {
  take:
    aln_output

  main:
    sam_clip(aln_output) | sam_con

  // Named outputs to send to medaka workflows
  emit:
    sam_con.out
}


// Polishing Workflow
workflow medaka {
  // takes the output from the align workflow
  take:
    aln_output

  main:
    // This workflow runs the draft alignment through medaka polishing three times and
    // then output is saved in the working directory in a directory called Medaka results
    medaka_1(aln_output) | medaka_2 | medaka_3

  emit:
    medaka_3.out

}


// Overall Workflow
workflow {

  qc | m2 | consensus

  if (params.method == 'all' || params.method == 'polish_only' || params.method == 'variant_only'){
    if (params.method == 'polish_only'){
    medaka(consensus.out)
}
    else {
      println "Variant calling only method was chosen"
}
    if (params.method == 'variant_only'){
    bcfcall(consensus.out) | snpann | vcf_parser
}
    else {
    println "Polishing only method was chosen"
}
}
  else {
    println "Alignment only method was chosen"
}
}