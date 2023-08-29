#!/usr/bin/env python3
import sys

info_len = []
high_var_type = []
mod_var_type = []
low_var_type = []

file = "${ann_calls}"
output_file1 = "high_impact_mutations.txt"
output_file2 = "moderate_impact_mutations.txt"
output_file3 = "low_impact_mutations.txt"

# SNPEff Effect Seq. Ontology acronyms dict
name_codes = {
    "coding_sequence_variant": "CDSV", "chromosome": "CHR", "duplication": "DUP",
    "inversion": "INV", "inframe_insertion": "INF_INS", "disruptive_inframe_insertion": "DII",
    "inframe_deletion": "INF_D", "disruptive_inframe_deletion": "DID", "downstream_gene_variant": "DGV",
    "exon_variant": "EXV", "exon_loss_variant": "EXLV", "frameshift_variant": "FSV",
    "gene_variant": "GNV", "eature_ablation": "FAB", "gene_fusion": "GFS", "bidirectional_gene_fusion": "BGF",
    "rearranged_at_DNA_level": "RAD", "intergenic_region": "IRV",
    "conserved_intergenic_variant": "CIV", "conservative_inframe_deletion": "CID",
    "intragenic_variant": "ITV", "intron_variant": "INTV",
    "conserved_intron_variant": "CINTV", "miRNA": "MIR", "missense_variant": "MSV", "initiator_codon_variant": "ICV",
    "stop_retained_variant": "STOPRV", "protein_protein_contact": "PPC", "structural_interaction_variant": "SIV",
    "rare_amino_acid_variant": "RAV", "splice_acceptor_variant": "SAV", "splice_donor_variant": "SDV",
    "splice_region_variant": "SRV", "stop_lost": "STOP-", "5_prime_UTR_premature_": "5UTR_PRE",
    "start_codon_gain_variant": "START+", "start_lost": "START-", "stop_gained": "STOP+",
    "synonymous_variant": "SYV", "start_retained": "STRT_R", "transcript_variant": "TRANSV",
    "regulatory_region_variant": "RRV", "upstream_gene_variant": "UGV", "3_prime_UTR_variant": "3UTRV",
    "3_prime_UTR_truncation + exon_loss": "3UTR_TRUNC_EXL", "5_prime_UTR_variant": "5UTRV",
    "5_prime_UTR_truncation + exon_loss_variant": "5UTR_TRUNC_EXL", "sequence_feature + exon_loss_variant": "SFV_EXL"
}

# VCF INFO codes to description dict
info_dict = {
    "IDV": "Maximum number of raw reads supporting an indel",
    "IMF": "Maximum fraction of raw reads supporting an indel",
    "DP": "Raw read depth",
    "VDB": "Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",
    "RPBZ": "Mann-Whitney U-z test of Read Position Bias (closer to 0 is better)",
    "MQBZ": "Mann-Whitney U-z test of Mapping Quality Bias (closer to 0 is better)",
    "BQBZ": "Mann-Whitney U-z test of Base Quality Bias (closer to 0 is better)",
    "MQSBZ": "Mann-Whitney U-z test of Mapping Quality vs Strand Bias (closer to 0 is better)",
    "SCBZ": "Mann-Whitney U-z test of Soft-Clip Length Bias (closer to 0 is better)",
    "SGB": "Segregation based metric",
    "MQ0F": "Fraction of MQ0 reads (smaller is better)",
    "AC": "Allele count in genotypes for each ALT allele, in the same order as listed",
    "AN": "Total number of alleles in called genotypes",
    "DP4": "Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases",
    "MQ": "Average mapping quality"
}


def variant_extractor(info_list):
    # Stripping allele line
    info = info_list[1:]

    # Splitting the info string into chunks of 15 sections
    info = [info[i:i + 15] for i in range(0, len(info), 15)]


    if info[0][1] == 'HIGH':
        high_var_type.append(list(set(chunk[0] for chunk in info)))

    elif info[0][1] == 'MODERATE':
        mod_var_type.append(list(set(chunk[0] for chunk in info)))

    else:
        low_var_type.append(list(set(chunk[0] for chunk in info)))


with open(file) as vcf:
    # This for loop aggregates all the SNPEff variant calls in the file
    for line in vcf:
        # Skipping file header lines
        if not line.startswith('#'):

            # Splitting line to better process each entry
            split_line = line.split('\t')

            # Extracting INFO data
            INFO = split_line[7].split(';')

            for _ in range(len(INFO)):
                # Variant extractor
                if INFO[_].startswith("ANN"):
                    variant_extractor(INFO[_].split('|'))

with open(file) as vcf:
    with open(output_file1, "w") as high_file, open(output_file2, "w") as moderate_file, open(output_file3, "w") as low_file:
        # This finds all the unique SNPEff calls and makes a list
        high_unique_vars = list(set(string for sublist in high_var_type for string in sublist))
        mod_unique_vars = list(set(string for sublist in mod_var_type for string in sublist))
        low_unique_vars = list(set(string for sublist in low_var_type for string in sublist))


        # Creating a dictionary of output file names
        unique_vars = {
            tuple(high_unique_vars): high_file,
            tuple(mod_unique_vars): moderate_file,
            tuple(low_unique_vars): low_file
        }

        # Iterate over the dictionary items and write the results to the corresponding output files
        # Here, it'll print all the SNPEff calls present in the vcf file and change them to their acronyms
        # This legend will be printing at the beginning of the output file for reference
        for var_list, output_file in unique_vars.items():
            output_file.write('SNPEff variant calls acronyms\n')
            output_file.write('-----------------------------\n')
            output_file.write('Below are all the types of variant calls present in the vcf file provided\n\n')

            for call in var_list:
                if '&' in call:
                    double_code = call.split('&')
                    spaced_call = ' & '.join(code for code in double_code)
                    double_acronym = ' & '.join(name_codes[code] for code in double_code)
                    output_file.write(spaced_call + ':' + double_acronym + '\n')
                else:
                    output_file.write(call + ':' + name_codes[call] + '\n')

        for line in vcf:
            # Skipping file header lines
            if not line.startswith('#'):
                # Splitting line to better process each entry
                split_line = line.split('\t')

                # Extracting INFO data
                INFO = split_line[7].split(';')

                # Checking impact of top variant per entry
                for _ in range(len(INFO)):
                    # Variant extractor
                    if INFO[_].startswith("ANN"):
                        impact = INFO[_].split('|')[2]

                # This section simply prints the original bcftools output for every entry
                fmt = "CHROM\tPOS\tREF\tALT\tQUAL".split('\t')
                line_list = [elem for i, elem in enumerate(split_line) if i not in (2, 6, 7, 8, 9)]
                concat = [': '.join(f_i) for f_i in zip(fmt, line_list)]
                final = '|'.join(concat)

                # Setting outputfile to the respective one based on the impact variable
                if impact == 'HIGH':
                    output_file = high_file

                elif impact == 'MODERATE':
                    output_file = moderate_file

                else:
                    output_file = low_file

                # Final printing
                output_file.write('\n')
                dash = '-' * len(final) + '\n'
                output_file.write(dash)
                output_file.write(final + '\n')
                info_header = ' SnpEff Variant Information '
                half_dash = '-' * (int((len(final) - len(info_header)) / 2))
                output_file.write(half_dash + info_header + half_dash + '\n' + '\n')

                # There is different spacing depending on the variant recording, so this handles the difference
                for _ in range(len(INFO)):
                    for key in info_dict.keys():
                        if INFO[_].split('=')[0] == key:
                            var_stats = INFO[_].split('=')
                            output_file.write(info_dict[key] + '= ' + var_stats[1] + '\n')

                    output_file.write('\n')

                    # INDEL variant handler
                    if INFO[_].startswith("INDEL"):
                        output_file.write('The variant is an INDEL\n\n')

                    # Annotation processor
                    if INFO[_].startswith("ANN"):
                        # Format sections as list
                        fmt = "Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | " \
                              "Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / " \
                              "CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO".split('| ')

                        # Filling in empty cells with N/A
                        info = [elem if elem else 'N/A' for elem in INFO[_].split('|')]

                        # Extracting the unique allele info and then trimming it
                        allele = info[0]
                        info = info[1:]

                        # Splitting the info string into chunks of 15 sections
                        info = [info[i:i + 15] for i in range(0, len(info), 15)]

                        for chunk in info:
                            call = chunk[0]
                            if '&' in call:
                                double_code = call.split('&')
                                chunk[0] = ' & '.join(name_codes[code] for code in double_code)
                            else:
                                chunk[0] = name_codes[call]

                        # Finding the longest string length and then applying a ljust to make the output evenly spaced
                        max_len_fmt = max(len(element) for element in fmt)
                        max_len_info = max(len(element) for chunk in info for element in chunk)

                        fmt = [element.ljust(max_len_fmt) for element in fmt]
                        info = [[element.ljust(max_len_info) for element in chunk] for chunk in info]

                        # chunked_info is the 3-level nested list that contains info but in chunks of length 7 to
                        # correct for vcf mutation entries with more than 7 mutation variant calls that made the final
                        # printout illegible

                        chunked_info = []  # Initialsing the 3rd nested list

                        for seven_wide in range(0, len(info), 7):
                            chunked_info.append(info[seven_wide:seven_wide + 7])

                        # Final printing
                        output_file.write('Allele : ' + allele + '\n')  # prints the single allele line

                        # Takes each length 7 chunk and transposes the nested info list to align it with
                        # the format column when printing
                        for _ in range(round(len(info) / 7)):
                            transposed_a = zip(*chunked_info[_])
                            concat = [' | '.join([elem] + list(sublist)) + ' |' for elem, sublist in
                                      zip(fmt, transposed_a)]

                            for i in concat:
                                output_file.write(i)
                                output_file.write('\n')

                            output_file.write('\n')

                        output_file.write('\n')

                    # Loss of Function variant handler
                    if INFO[_].startswith('LOF'):
                        # Statement
                        output_file.write('Predicted loss of function effects for this variant:\n')
                        lof_fmt = 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | ' \
                                  'Percent_of_transcripts_affected'.split(' | ')

                        # Lists all the associated stats seen above in lof_fmt
                        for position in range(4):
                            output_file.write(lof_fmt[position] + ':' + INFO[_][5:-1].split('|')[position] + '\n')

                    # Nonsense mediated decay variant handler
                    if INFO[_].startswith('NMD'):
                        # Statement
                        output_file.write('Predicted nonsense mediated decay effects for this variant:\n')
                        lof_fmt = 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | ' \
                                  'Percent_of_transcripts_affected'.split(' | ')

                        # Lists all the associated stats seen above in lof_fmt
                        for position in range(4):
                            output_file.write(lof_fmt[position] + ':' + INFO[_][5:-1].split('|')[position] + '\n')
