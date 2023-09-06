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
- staphb/bcftools:1.12
- smobed/snpeff_mpox:latest
- smobed/python_mpox:latest
 

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

![image](https://github.com/s-mobed/Nanopox/assets/91598812/e247d209-32ae-4d10-8b4a-cff9bc5550ab)

