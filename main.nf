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
      tuple val(sample_id), path(sample), path(scheme), path(reference), path(stats), val(prefix)
      
    output:
      tuple val(sample_id), path(sample), path('mini.sam'), path(scheme), path(reference), path(stats), val(prefix)

    script:
    if(sample_id == null){ sample_id = file(sample).getSimpleName()}
    """
    minimap2 -a -t 3 -x map-ont ${reference} ${sample} -o mini.sam
    """
}

process ngmlr {
    // thread allocation
    cpus params.threads

    // I/O
    input:
      tuple val(sample_id), path(sample), path(scheme), path(reference), path(stats), val(prefix)

    output:
      tuple val(sample_id), path(sample), path('ngmlr.sam'), path(scheme), path(reference), path(stats), val(prefix)

    script:
    if(sample_id == null){ sample_id = file(sample).getSimpleName()}
    """
    ngmlr -x ont -t ${task.cpus} --bam-fix -r ${reference} -q ${sample} | \

    # Fixes NGMLR issue where it sets MAPQ value to lowest int32 value and stops ampliconclip from working since MAPQ cannot be negative
    # Due to the way it calculates MAPQ if a value exceeds the int32 range, integer overflow occurs and results in
    # it wrapping to the lowest value -2,147,483,648. Therefore I set the MAPQ value to 60, the highest value for NGMLR

    awk -F '\t' -v OFS='\t' '!/^@/ && \$5 < 0 { \$5 = 60 } { print }' > ngmlr.sam
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
      tuple val(sample_id), path(sample), path(sam_file), path(scheme), path(reference), path(stats), val(prefix), val(suffix)
      
    output:
      tuple val(sample_id), path(sample), path('*.bam'), path(reference), val(prefix), val(suffix)

    script:
    """
    # Samtools sorting, clipping primers, and converting to bam format
    samtools ampliconclip -u --both-ends -f clip_stats.txt -b ${scheme} ${sam_file} | \
    samtools sort -u -O bam -@ ${task.cpus} -o ${prefix}${suffix}.bam
    
    # Appending amplicon clip and alignment stats
    cp ${stats} ${prefix}${suffix}_stats.txt
    echo >> ${prefix}${suffix}_stats.txt
    cat clip_stats.txt >> ${prefix}${suffix}_stats.txt
    samtools stats *.bam -c 1,1000,25 | awk '
    /^SN/ { if (!sn_flag) { print "Summary Numbers"; sn_flag = 1 } print substr(\$0, index(\$0, \$2)); next }
    /^ID/ { if (!id_flag) { print "\\nIndel Distribution"; id_flag = 1 } print substr(\$0, index(\$0, \$2)); next }
    /^COV/ { if (!cov_flag) { print "\\nCoverage Distribution (increments of 25)"; cov_flag = 1 } print substr(\$0, index(\$0, \$2)); next }' >> ${prefix}${suffix}_stats.txt
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
      tuple val(sample_id), path(sample), path(bam_file), path(reference), val(prefix), val(suffix)
      
    output:
      tuple val(sample_id), path(sample), path('*.fasta'), val(prefix), val(suffix)

    script:
    """
    # Draft consensus to kick off Medaka
    samtools consensus -o ${prefix}${suffix}.fasta ${bam_file}
    """
}

process vir_con {
    // thread allocation
    cpus params.threads

    // Output Dir
    publishDir "${params.draft_dir}/${sample_id}", pattern: "*.fasta", mode: "copy"

    // I/O
    input:
      tuple val(sample_id), path(sample), path(bam_file), path(reference), val(prefix), val(suffix)
      
    output:
      tuple val(sample_id), path(sample), path('*.fasta'), val(prefix), val(suffix)

    script:
    """
    # Draft consensus to kick off Medaka
    viral_consensus -q 1 -d 1 -i ${bam_file} -r ${reference} -o ${prefix}${suffix}.fasta
    """
}


// Polishing Processes

process medaka_1 {
    // thread allocation
    cpus params.threads

    // I/O
    input:
      tuple val(sample_id), path(sample), path(draft), val(prefix), val(suffix)

    output:
      tuple val(sample_id), path(sample), path("medaka/*.fasta"), val(prefix), val(suffix)

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
      tuple val(sample_id), path(sample), path(med_one), val(prefix), val(suffix)

    output:
      tuple val(sample_id), path(sample), path("medaka/*.fasta"), val(prefix), val(suffix)

    script:
    """
    medaka_consensus -i ${sample} -d ${med_one} -t ${task.cpus}
    """
}

process medaka_3 {
    // thread allocation
    cpus params.threads

    // Output Dir
    publishDir "${params.med_dir}", pattern: "${file(sample).getSimpleName()}/${process}/*.fasta", mode: "copy"

    // I/O
    input:
      tuple val(sample_id), path(sample), path(med_two), val(prefix), val(suffix)

    output:
      path "${file(sample).getSimpleName()}/${prefix}${suffix}/*.fasta"

    script:
    """
    medaka_consensus -i ${sample} -d ${med_two} -o ${sample_id}/${prefix}${suffix} -t ${task.cpus}
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

// Alignment Workflows
workflow mini_wf {
  take:
    qc_output
    
  main:
    input_tuple=qc_output.combine(channel.value('mini_'))
    minimap(input_tuple)

  // Named outputs to send to medaka workflows
  emit:
    minimap.out
}

workflow ngmlr_wf {
  take:
    qc_output

  main:
    input_tuple=qc_output.combine(channel.value('ngmlr_'))
    ngmlr(input_tuple)

  // Named outputs to send to medaka workflows
  emit:
    ngmlr.out
}

// Consensus Workflows
workflow clip_1 {
  take:
    aln_output

  main:
    input_tuple=aln_output.combine(channel.value('sam'))
    sam_clip(input_tuple)

  // Named outputs to send to medaka workflows
  emit:
    sam_clip.out
}

workflow clip_2 {
  take:
    aln_output

  main:
    input_tuple=aln_output.combine(channel.value('vc'))
    sam_clip(input_tuple)

  // Named outputs to send to medaka workflows
  emit:
    sam_clip.out
}

workflow samtools_con_1 {
  take:
    aln_output

  main:
    input_tuple=aln_output.combine(channel.value('sam'))
    sam_clip(input_tuple) | sam_con

  // Named outputs to send to medaka workflows
  emit:
    sam_con.out
}

workflow samtools_con_2 {
  take:
    aln_output

  main:
    input_tuple=aln_output.combine(channel.value('sam'))
    sam_clip(input_tuple) | sam_con

  // Named outputs to send to medaka workflows
  emit:
    sam_con.out
}

workflow viral_con_1 {
  take:
    aln_output

  main:
    input_tuple=aln_output.combine(channel.value('vc'))
    sam_clip(input_tuple) | vir_con

  // Named outputs to send to medaka workflows
  emit:
    vir_con.out
}

workflow viral_con_2 {
  take:
    aln_output

  main:
    input_tuple=aln_output.combine(channel.value('vc'))
    sam_clip(input_tuple) | vir_con

  // Named outputs to send to medaka workflows
  emit:
    vir_con.out
}

// Polishing Workflows

workflow medaka_mini_sam {
  // takes the output from the align workflow
  take:
    aln_output

  main:
    // This workflow runs the draft alignment through medaka polishing three times and
    // then output is saved in the working directory in a directory called Medaka results
    medaka_1(aln_output) | medaka_2 | medaka_3
}

workflow medaka_mini_vc {
  // takes the output from the align workflow
  take:
    aln_output

  main:
    // This workflow runs the draft alignment through medaka polishing three times and
    // then output is saved in the working directory in a directory called Medaka results
    medaka_1(aln_output) | medaka_2 | medaka_3
}

workflow medaka_ngmlr_sam {
  // takes the output from the align workflow
  take:
    aln_output

  main:
    // This workflow runs the draft alignment through medaka polishing three times and
    // then output is saved in the working directory in a directory called Medaka results
    medaka_1(aln_output) | medaka_2 | medaka_3
}

workflow medaka_ngmlr_vc {
  // takes the output from the align workflow
  take:
    aln_output

  main:
    // This workflow runs the draft alignment through medaka polishing three times and
    // then output is saved in the working directory in a directory called Medaka results
    medaka_1(aln_output) | medaka_2 | medaka_3
}


workflow {
  qc()
  
  if ( params.method =~ 'mini_*' || params.method == 'all'){
    mini_wf(qc.out)

    if ( params.method == 'mini_sam' || params.method == 'all'){
      samtools_con_1(mini_wf.out) | medaka_mini_sam
    }

    if ( params.method == 'mini_vc' || params.method == 'all'){
      viral_con_1(mini_wf.out) | medaka_mini_vc
    }
  }
  
  if ( params.method =~ 'ngmlr_*' || params.method == 'all'){
    ngmlr_wf(qc.out)

    if ( params.method == 'ngmlr_sam' || params.method == 'all'){
      samtools_con_2(ngmlr_wf.out) | medaka_ngmlr_sam
    }

    if ( params.method == 'ngmlr_vc' || params.method == 'all'){
      viral_con_2(ngmlr_wf.out) | medaka_ngmlr_vc
    }
  }
}