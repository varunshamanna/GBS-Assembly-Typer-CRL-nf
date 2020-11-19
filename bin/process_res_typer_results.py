#!/usr/bin/env python
import argparse
import sys
import os
import re
import glob
import subprocess
from Bio import SeqIO

# TODO Add column headings to the BIN file
# TODO it does the variant resistance typing (v2) but no obvious output.
#      The one that looked vaguely similar was the file genes_ResFinder_results and the same for ARGANNOT.
#      we should think of a way to output in a table if not done already

drugRes_Col = {
    'TET': 'neg',
    'EC': 'neg',
    'FQ': 'neg',
    'OTHER': 'neg',
}

Res_Targets = {
    'ERM': 'neg',
    'LNUB': 'neg',
    'LSA': 'neg',
    'MEF': 'neg',
    'TET': 'neg',
    'CAT': 'neg',
    'GYRA': 'neg',
    'PARC': 'neg',
    '23S1': 'neg',
    '23S3': 'neg',
    'RPOB1': 'neg',
    'RPOB2': 'neg',
    'RPOB3': 'neg',
    'RPOB4': 'neg',
}

EOL_SEP = "\n"


def extract_seq_by_id(lookup, fasta_file):
    """ Extract a sequence from the given fasta file by feature id """

    rf = None
    try:
        rf = open(fasta_file, "r")
        for record in SeqIO.parse(rf, 'fasta'):
            if record.id == lookup:
                return lookup + EOL_SEP + record.seq.strip() + EOL_SEP
    finally:
        rf.close()

    return None


def freebayes_prior_fix(bam_file, ref_file, target):
    sam_file = bam_file.replace('.bam', '.sam')

    os.system("samtools view -h " + bam_file + " > " + samfile)
    os.system("cat " + sam_file + " | grep -E \"^\@HD|^\@SQ.*" + target + "|^\@PG\" > CHECK_target_seq.sam")
    os.system("awk -F'\t' '\$3 == \"" + target + "\" {print \$0}' " + sam_file + " >> CHECK_target_seq.sam")
    os.system("samtools view -bS CHECK_target_seq.sam > CHECK_target_seq.bam")
    os.system("samtools index CHECK_target_seq.bam CHECK_target_seq.bai")

    ref_seq = extract_seq_by_id(target, ref_file)

    rf = None
    try:
        rf = open("CHECK_target_ref.fna", "w")
        rf.write(ref_seq + EOL_SEP)
    finally:
        rf.close()

    os.system("freebayes -q 20 -p 1 -f CHECK_target_ref.fna CHECK_target_seq.bam -v CHECK_target_seq.vcf")
    os.system("bgzip CHECK_target_seq.vcf")

    os.system("tabix -p vcf CHECK_target_seq.vcf.gz")
    extract_seq = subprocess.check_output("echo \"" + ref_seq + "\" | vcf-consensus CHECK_target_seq.vcf.gz", shell=True)

    os.system("rm CHECK_target*")

    return extract_seq.rstrip(EOL_SEP)


def codon2aa(codon):

    codon = codon.upper()

    codon_dict = {
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*', 'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G'
    }

    try:
        result = codon_dict[codon]
    except KeyError:
        if re.search(r"GC.", codon):
            result = 'A'
        elif re.search(r"GG.", codon):
            result = 'G'
        elif re.search(r"CC.", codon):
            result = 'P'
        elif re.search(r"AC.", codon):
            result = 'T'
        elif re.search(r"GT.", codon):
            result = 'V'
        elif re.search(r"CG.", codon):
            result = 'R'
        elif re.search(r"TC.", codon):
            result = 'S'
        else:
            result = 'x'
            print("Bad codon " + codon + "!!\n")

    return result


def six_frame_translate(seq_input, opt_f):
    pass


def derive_presence_absence_targets(input_file):
    with open(input_file, 'r') as fd:
        # Skip header row
        next(fd)
        # Process file lines
        for line in fd:
            fields = line.split('\t')
            gene = fields[2]
            allele = fields[3]
            depth = float(fields[5])

            if depth >= 10:
                if re.search(r"ERM|LNUB|LSA|MEF", allele):
                    if drugRes_Col["EC"] == "neg":
                        drugRes_Col["EC"] = gene
                    else:
                        new_val = drugRes_Col["EC"] + ":" + gene
                        drugRes_Col["EC"] = new_val

                if re.search(r"TET", allele):
                    if drugRes_Col["TET"] == "neg":
                        drugRes_Col["TET"] = gene
                    else:
                        new_val = drugRes_Col["TET"] + ":" + gene
                        drugRes_Col["TET"] = new_val

                if re.search(r"OTHER", allele):
                    if drugRes_Col["OTHER"] == "neg":
                        drugRes_Col["OTHER"] = gene
                    else:
                        new_val = drugRes_Col["OTHER"] + ":" + gene
                        drugRes_Col["OTHER"] = new_val

                if re.search(r"ERM", allele):
                    Res_Targets["ERM"] = "pos"
                elif re.search(r"LNUB", allele):
                    Res_Targets["LNUB"] = "pos"
                elif re.search(r"LSA", allele):
                    Res_Targets["LSA"] = "pos"
                elif re.search(r"MEF", allele):
                    Res_Targets["MEF"] = "pos"
                elif re.search(r"TET", allele):
                    Res_Targets["TET"] = "pos"
                elif re.search(r"CAT", allele):
                    Res_Targets["CAT"] = "pos"
                elif re.search(r"PARC", allele):
                    Res_Targets["PARC"] = "pos"
                elif re.search(r"GYRA", allele):
                    Res_Targets["GYRA"] = "pos"
                elif re.search(r"23S1", allele):
                    Res_Targets["23S1"] = "pos"
                elif re.search(r"23S3", allele):
                    Res_Targets["23S3"] = "pos"
                elif re.search(r"RPOB1", allele):
                    Res_Targets["RPOB1"] = "pos"
                elif re.search(r"RPOBN", allele):
                    Res_Targets["RPOB2"] = "pos"

            elif re.search(r"RPOBN", allele):
                Res_Targets["RPOB3"] = "pos"
            elif re.search(r"RPOBN", allele):
                Res_Targets["RPOB4"] = "pos"


def derive_presence_absence_targets_for_arg_res(input_file):
    """ For ARG-ANNOT / ResFinder """

    with open(input_file, 'r') as fd:
        # Skip header row
        next(fd)
        # Process file lines
        for line in fd:
            fields = line.split('\t')
            gene = fields[2]
            allele = fields[3]
            depth = float(fields[5])

            if depth >= 10:
                if re.search(r"ERM", allele):
                    if Res_Targets["ERM"] == "neg":
                        if drugRes_Col["EC"] == "neg":
                            drugRes_Col["EC"] = "ERM"
                        else:
                            drugRes_Col["EC"] = drugRes_Col["EC"] + ":ERM"

                        Res_Targets["ERM"] = "pos";

                elif re.search(r"LNU", allele):
                    if Res_Targets["LNUB"] == "neg":
                        if drugRes_Col["EC"] == "neg":
                            drugRes_Col["EC"] = "LNU"
                        else:
                            drugRes_Col["EC"] = drugRes_Col["EC"] + ":LNU"

                        Res_Targets["LNUB"] = "pos"

                elif re.search(r"LSA", allele):
                    if Res_Targets["LSA"] == "neg":
                        if drugRes_Col["EC"] == "neg":
                            drugRes_Col["EC"] = "LSA"
                        else:
                            drugRes_Col["EC"] = drugRes_Col["EC"] + ":LSA"

                        Res_Targets["LSA"] = "pos"

                elif re.search(r"MEF", allele):
                    if Res_Targets["MEF"] == "neg":
                        if drugRes_Col["EC"] == "neg":
                            drugRes_Col["EC"] = "MEF"
                        else:
                            drugRes_Col["EC"] = drugRes_Col["EC"] + ":MEF"

                        Res_Targets["MEF"] = "pos"

                elif re.search(r"TET", allele):
                    if Res_Targets["TET"] == "neg":
                        if drugRes_Col["TET"] == "neg":
                            drugRes_Col["TET"] = "TET"
                        else:
                            drugRes_Col["TET"] = drugRes_Col["TET"] + ":TET"

                        Res_Targets["TET"] = "pos"

                elif re.search(r"CAT", allele):
                    if Res_Targets["CAT"] == "neg":
                        if drugRes_Col["OTHER"] == "neg":
                            drugRes_Col["OTHER"] = "CAT"
                        else:
                            drugRes_Col["OTHER"] = drugRes_Col["OTHER"] + ":CAT";

                        Res_Targets["CAT"] = "pos"

                else:
                    if drugRes_Col["OTHER"] == "neg":
                        drugRes_Col["OTHER"] = gene
                    else:
                        drugRes_Col["OTHER"] = drugRes_Col["OTHER"] + ":" + gene


def run(srst2_gbs_output, srst2_argannot_output, srst2_resfinder_output):
    TEMP_RES_bam = glob.glob("RES_*.sorted.bam")
    TEMP_RES_fullgene = glob.glob("RES_*__fullgenes__*__results.txt")
    RES_bam = TEMP_RES_bam[0]
    RES_full_name = TEMP_RES_fullgene[0]
    print("res bam is: " + RES_bam + ", res full gene is " + RES_full_name)
    RES_vcf = RES_bam.replace(".bam", ".vcf")
    RES_bai = RES_bam.replace(".bam", ".bai")

    TEMP_ARG_fullgene = glob.glob("ARG_*__fullgenes__*__results.txt")
    ARG_full_name = TEMP_ARG_fullgene[0]
    TEMP_RESFI_fullgene = glob.glob("RESFI_*__fullgenes__*__results.txt")
    RESFI_full_name = TEMP_RESFI_fullgene[0]
    merged_net = "ARG-RESFI_fullgenes_results.txt"

    # Merge the arg-annot and resfinder results
    os.system("tail -n+2 {} > {}".format(ARG_full_name, merged_net))
    os.system("tail -n+2 {} >> {}".format(RESFI_full_name, merged_net))

    derive_presence_absence_targets(RES_full_name)
    print(Res_Targets)


def get_arguments():
    parser = argparse.ArgumentParser(description='Modify SRST2 sequence typing output files.')
    parser.add_argument('--srst2_gbs', '-g', dest='srst2_gbs_output', required=True,
                        help='Input SRST2 output for the GBS reference database.')
    parser.add_argument('--srst2_argannot', '-a', dest='srst2_argannot_output', required=True,
                        help='Input SRST2 output for the ARG-ANNOT reference database.')
    parser.add_argument('--srst2_resfinder', '-r', dest='srst2_resfinder_output', required=True,
                        help='Input SRST2 output for the ResFinder reference database.')
    parser.add_argument('--output', '-o', dest='output', required=True,
                        help='Output results filename.')
    parser.add_argument('--output_bin', '-b', dest='output_bin', required=True,
                        help='Output BIN results filename.')

    return parser


def main():
    args = get_arguments().parse_args()
    run(args.srst2_gbs_output, args.srst2_argannot_output, args.srst2_resfinder_output)


if __name__ == "__main__":
    sys.exit(main())
