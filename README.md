![Nanopox_logo](https://github.com/s-mobed/Nanopox/assets/91598812/f98e6d73-324a-43ef-8075-cf901def55b5)


# Nanopox
This is the repository of my Msc BioInformatics summer project. The goal of this project is to generate consensus genomes of Mpox from Oxford Nanopore long read data from a tiled amplicon sequencing strategy. 

# Pipeline Installation

To install the pipeline, the user would first need to install Nextflow by downloading it directly using this command: 

    curl -fsSL https://get.nextflow.io | bash

It also can be installed through the Bioconda channel:

    conda install -c bioconda nextflow

The user would also need to download Docker on their system to fetch the containers needed for the pipeline to execute. As Nextflow has integrated Docker support, there are no further steps needed. This is the list of Docker images used in the pipeline, which can be found in the Nextflow config file:
 
-	smobed/nanoq:latest
-	staphb/minimap2:latest
-	staphb/samtools:1.12
-	staphb/samtools:1.17
-	ontresearch/medaka:v1.8.0
-	cyverse/ngmlr:0.2.7
-	niemasd/viral_consensus

However, if having the Docker images pre-installed before running the pipeline is preferred by the user. Using the Docker pull command with names of the images from the list above, wil pre-install them.

Lastly, the pipeline can be pulled from the Github repository it sits in, using this command:

    nextflow pull s-mobed/Nanopox/

This command will download the Nextflow script and the configuration file necessary to run the pipeline.



# Pipeline execution
To execute the complete pipeline, employ the following command:

    nextflow run s-mobed/Nanopox --reference reference_genome.fa --fastq ONT_reads.fastq  --primer_scheme primer.scheme.bed –threads [NUM] 

These flags are essential for proper pipeline execution. Running the command without these flags or solely with the --help flag will display the help message, offering usage details and optional flags. If multiple FASTQ files need processing in a single run (e.g., when files are in a folder named Data), the input flag would be --fastq ‘Data/*.fastq’. The inclusion of quotes is vital for Nextflow to interpret it as a global path matcher.

The pipeline also integrates an alternative input format within Nextflow. Utilising the --sra flag instead of the --fastq flag, along with the optional --ncbi_key flag, users can extract files from the NCBI SRA database by providing project or accession numbers.
nextflow run s-mobed/Nanopox --reference reference_genome.fa --sra ‘PRJNA925815’  --primer_scheme primer.scheme.bed --threads [NUM] --ncbi_key ‘0123456789abcdef’. When dealing with multiple SRA IDs, supply them as a list: “['ERR000001', 'ERR000002', 'ERR000003']”.

There are four combinations of tools used in this pipeline using alignment tools Minimap2 (Li, 2018; https://github.com/lh3/minimap2) and NGMLR (Sedlazeck et al., 2018; https://github.com/philres/ngmlr), and the consensus callers are Samtools consensus (http://www.htslib.org/doc/samtools-consensus.html), and Viral consensus (Moshiri et al., 2023; https://github.com/niemasd/ViralConsensus). 

The combinations are:
-	Minimap2 + Samtools consensus  (mini_sam)
-	NGMLR + Samtools consensus  (ngmlr_sam)
-	Minimap2 + Viral consensus (mini_vc)
-	NGMLR + Viral consensus (ngmlr_vc)

The pipeline will only execute mini_sam by default, if users prefer to execute other or all combinations, they can use the --method flag with the following options: 
- mini_sam 
- mini_vc 
- ngmlr_sam
- ngmlr_vc 
- all

# Variant Calling Alternative Method

Variant calling wasn't implemented in the final version of project code, but was implemented shortly afterwards. To use this alternative  pipeline. If you pull the separate branch with this command:

    nextflow pull -r variant_call s-mobed/Nanopox

and then you just add the -r flag like this, when running the alternative:

    nextflow run s-mobed/Nanopox -r variant_call --reference reference_genome.fa --fastq 'Data/*.fastq.gz' --primer_scheme primer.scheme.bed --threads [NUM] --method variant_only

To simplify the pipeline, I removed the NGMLR and Viral consensus tools and now the pipeline only runs the mini_sam combination from before for the alignment and draft consensus calling steps. Therfore, I’ve changed the method flag so there are these methods:

- alignment_only
- polish_only (alignment + polish) 
- variant_only (alignment + variant calling)
- all

For the variant calling step, I use these docker images, and have been included in the alternative nextflow.config file.

- staphb/bcftools:1.12
- smobed/snpeff_mpox:latest
- smobed/python_mpox:latest

# Acknowledgements

I would like to thank my supervisor Dr. Sreenu Vattipally (https://github.com/vbsreenu) for all the guidance and help throughout my Masters Project. I would like to thank my brother, Dean Mobed (https://www.linkedin.com/in/dean-mobed-75231b198/), for the logo of the pipeline.
