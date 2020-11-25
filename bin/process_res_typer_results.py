#!/usr/bin/env python
import argparse
import sys
import os
import re
import glob
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

class nSeq(str): # Nucleotide sequence
    pass

class aSeq(str): # Amino acid sequence
    pass

drugRes_Col = {
    'TET': 'neg',
    'EC': 'neg',
    'FQ': 'neg',
    'OTHER': 'neg',
}

drugToClass = {
    'ERM': 'EC',
    'LNU': 'EC',
    'LSA': 'EC',
    'MEF': 'EC',
    'TET': 'TET',
    'CAT': 'OTHER',
    'PARC': 'FQ',
    'GYRA': 'FQ',
    '23S1': 'EC',
    '23S3': 'EC',
    'RPOBGBS-1': 'OTHER',
    'RPOBGBS-2': 'OTHER',
    'RPOBGBS-3': 'OTHER',
    'RPOBGBS-4': 'OTHER',
}

Res_Targets = {
    'ERM': 'neg',
    'LNU': 'neg',
    'LSA': 'neg',
    'MEF': 'neg',
    'TET': 'neg',
    'CAT': 'neg',
}

GBS_Res_Targets = {
    'GYRA': 'neg',
    'PARC': 'neg',
    '23S1': 'neg',
    '23S3': 'neg',
    'RPOBGBS-1': 'neg',
    'RPOBGBS-2': 'neg',
    'RPOBGBS-3': 'neg',
    'RPOBGBS-4': 'neg',
}

GBS_Res_var = {
    'GYRA': '', #10
    'PARC': '', #14
    '23S1': '', #0
    '23S3': '', #1
    'RPOBGBS-1': '', #6
    'RPOBGBS-2': '', #7
    'RPOBGBS-3': '', #8
    'RPOBGBS-4': '', #9,
}

geneToRef = defaultdict(lambda: '')
geneToRef.update({
    'PARC': aSeq('HPHGDSSIYDAMVRMSQ'),
    'GYRA': aSeq('VMGKYHPHGDSSIYEAMVRMAQWW'),
    '23S1': nSeq('GTTACCCGCGACAGGACGGAAAGACCCCATGGAG'),
    '23S3': nSeq('CGGCACGCGAGCTGGGTTCAGAACGTCGTGAGACAGTTCGGTCCCTATCCGTCGCGGGCG'),
    'RPOBGBS-1': aSeq('FGSSQLSQFMDQHNPLSELSHKRRLSALGPGGL'),
    'RPOBGBS-2': aSeq('VSQLVRSPGV'),
    'RPOBGBS-3': aSeq('FTVAQANSKLNEDGTFAEEIVMGRHQGNNQEFPSSI'),
    'RPOBGBS-4': aSeq('LIDPKAPYVGT')
})

geneToTargetSeq = defaultdict(lambda: '')
geneToTargetSeq.update({
    'PARC': '7__PARCGBS__PARCGBS-1__7',
    'GYRA': '5__GYRAGBS__GYRAGBS-1__5',
    '23S1': '11__23S1__23S1-1__11',
    '23S3': '12__23S3__23S3-3__12',
    'RPOBGBS-1': '16__RPOBgbs__RPOBgbs-1__16',
    'RPOBGBS-2': '17__RPOBgbs__RPOBgbs-2__17',
    'RPOBGBS-3': '18__RPOBgbs__RPOBgbs-3__18',
    'RPOBGBS-4': '19__RPOBgbs__RPOBgbs-4__19'
})

EOL_SEP = "\n"
MIN_DEPTH = 10


def extract_seq_by_id(lookup: str, fasta_file: str) -> str:
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


def freebayes_prior_fix(bam_file: str, ref_file: str, target: str) -> str:
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


def codon2aa(codon: str) -> str:

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
            print("Bad codon " + codon + "!!")

    return result


def extract_frame_aa(sequence: str, frame: int) -> str:
    """
    :param sequence: dna bases
    :param frame: frame number to extract
    :return: Amino acid translates frame
    """

    bases = sequence
    if frame > 3:
        seq = Seq(sequence)
        bases = str(seq.reverse_complement())
        frame -= 3

    i = frame - 1
    protein = ""
    while i < len(bases)-2:
        protein += codon2aa(bases[i:i+3])
        i += 3

    return protein


def six_frame_translate(seq_input: str, frame: int) -> str:
    """
    :param seq_input: Fasta feature including id and sequence lines
    :param frame: Codon number:
    :return: protein translation
    """

    if frame < 1 or frame > 6:
        raise IndexError("Frame number argument is out of bounds: " + str(frame))

    lines = seq_input.splitlines()
    dna = ""
    for line in lines:
        if line.startswith('>'):
            continue
        else:
            dna += line.strip()

    return extract_frame_aa(dna, frame)


def update_presence_absence_target(gene, allele, depth, drug_res_col_dict, res_target_dict, gbs_res_target_dict):
    if depth >= MIN_DEPTH:

        for gene_name in res_target_dict.keys():
            if re.search(gene_name, allele.upper()):
                drugCat = drugToClass[gene_name]

                if drug_res_col_dict[drugCat] == "neg":
                    drug_res_col_dict[drugCat] = gene + '(' + allele + ')'
                else:
                    new_val = drug_res_col_dict[drugCat] + ":" + gene + '(' + allele + ')'
                    drug_res_col_dict[drugCat] = new_val

                res_target_dict[gene_name] = "pos"

        for gene_name in gbs_res_target_dict.keys():
            if re.search(gene_name, allele):

                gbs_res_target_dict[gene_name] = "pos"


def derive_presence_absence_targets(input_file: str):
    try:
        with open(input_file, 'r') as fd:
            # Skip header row
            next(fd)

            # Process file lines
            for line in fd:
                fields = line.split('\t')
                gene = fields[2]
                allele = fields[3]
                depth = float(fields[5])
                update_presence_absence_target(gene, allele, depth, drugRes_Col, Res_Targets, GBS_Res_Targets)
    except IOError:
        print('Cannot open {}.'.format(filename))


def update_presence_absence_target_for_arg_res(gene, allele, depth, drug_res_col_dict, res_target_dict):
    if depth >= MIN_DEPTH:

        other = 1
        for gene_name in res_target_dict.keys():
            if re.search(gene_name, allele.upper()):
                other = 0
                drugCat = drugToClass[gene_name]
                if res_target_dict[gene_name] == "neg":
                    if drug_res_col_dict[drugCat] == "neg":
                        drug_res_col_dict[drugCat] = gene + '(' + allele + ')'
                    else:
                        drug_res_col_dict[drugCat] = drug_res_col_dict[drugCat] + ':' + gene + '(' + allele + ')'

                    res_target_dict[gene_name] = "pos"

        if other:
            if drug_res_col_dict["OTHER"] == "neg":
                drug_res_col_dict["OTHER"] = gene + '(' + allele + ')'
            else:
                drug_res_col_dict["OTHER"] = drug_res_col_dict["OTHER"] + ":" + gene + '(' + allele + ')'


def derive_presence_absence_targets_for_arg_res(input_files: list):
    """ For ARG-ANNOT / ResFinder """
    for input_file in input_files:
        try:
            with open(input_file, 'r') as fd:
                # Skip header row
                next(fd)
                # Process file lines
                for line in fd:
                    fields = line.split('\t')
                    gene = fields[2]
                    allele = fields[3]
                    depth = float(fields[5])
                    update_presence_absence_target_for_arg_res(gene, allele, depth, drugRes_Col, Res_Targets)
        except IOError:
            print('Cannot open {}.'.format(filename))


def find_mismatches(seq_diffs, query_Seq, ref_Seq):
    for resi in range(len(query_Seq)):
        if query_Seq[resi] != ref_Seq[resi]:
            seq_diffs.append(ref_Seq[resi] + str(resi+1) + query_Seq[resi])
    return seq_diffs


def get_seq_diffs(query_Seq, ref_Seq):
    if type(ref_Seq) == aSeq:
        query_Seq = six_frame_translate(query_Seq, 1)
    seq_diffs = []
    if query_Seq != ref_Seq:
        seq_diffs = find_mismatches(seq_diffs, query_Seq, ref_Seq)
    return seq_diffs


def update_GBS_Res_var(gene_name, seq_diffs, bin_res_arr):
    if seq_diffs:
        bin_res_arr[gene_name] = gene_name + '-' + ','.join(seq_diffs)
    else:
        bin_res_arr[gene_name] = gene_name


def update_drug_res_col_dict(gene_name, seq_diffs, drugRes_Col, drugToClass):
    drugClass = drugToClass[gene_name]
    if seq_diffs:
        gene_var = gene_name + '-' + ','.join(seq_diffs)
    else:
        gene_var = gene_name
    if drugRes_Col[drugClass] == 'neg':
        drugRes_Col[drugClass] = gene_var
    else:
        new_value = drugRes_Col[drugClass] + ':' + gene_var
        drugRes_Col[drugClass] = new_value


def get_consensus_seqs(consensus_seqs_file):
    consensus_seq_dict = defaultdict(lambda: '')
    try:
        with open(consensus_seqs_file, 'r') as fd:
            for line in fd:
                if line[0] == '>':
                    seq_name = line.split('>')[1].split('\n')[0]
                else:
                    consensus_seq_dict[seq_name] = line.split('\n')[0]
    except IOError:
        print('Cannot open {}.'.format(filename))

    return consensus_seq_dict


def get_gene_names_from_consensus(consensus_seq_dict):
    gene_names = []
    for gene_name in geneToTargetSeq.keys():
        if geneToTargetSeq[gene_name] in consensus_seq_dict.keys():
            gene_names.append(gene_name)
    return gene_names


def get_variants(consensus_seqs: str):

    consensus_seq_dict = get_consensus_seqs(consensus_seqs)
    gene_names = get_gene_names_from_consensus(consensus_seq_dict)
    for gene_name in gene_names:
        if GBS_Res_Targets[gene_name] == "pos" and geneToTargetSeq[gene_name] and geneToRef[gene_name]:
            seq_diffs = get_seq_diffs(consensus_seq_dict[geneToTargetSeq[gene_name]], geneToRef[gene_name])
            update_GBS_Res_var(gene_name, seq_diffs, GBS_Res_var)
            update_drug_res_col_dict(gene_name, seq_diffs, drugRes_Col, drugToClass)


def write_output(output_pref):
    try:
        with open(output_pref + '_res_incidence.txt', 'w') as out_inc:
            out_inc.write('\t'.join(Res_Targets.keys()) + '\t' + '\t'.join(GBS_Res_Targets.keys()) + '\n')
            out_inc.write('\t'.join(Res_Targets.values()) + '\t' + '\t'.join(GBS_Res_Targets.values()) + '\n')

        with open(output_pref + '_res_gbs_variants.txt', 'w') as out_var:
            out_var.write('\t'.join(GBS_Res_var.keys()) + '\n')
            out_var.write('\t'.join(GBS_Res_var.values()) + '\n')

        with open(output_pref + '_res_alleles.txt', 'w') as out_allele:
            out_allele.write('\t'.join(drugRes_Col.keys()) + '\n')
            out_allele.write('\t'.join(drugRes_Col.values()) + '\n')
    except IOError:
        print('Cannot open filename starting "{}"'.format(output_pref))


def run(args):
    # Get presence/absence of genes
    if args.srst2_gbs_fg_output is not None:
        print("Running GBS presence/absence...")
        derive_presence_absence_targets(args.srst2_gbs_fg_output)

    if args.srst2_argannot_fg_output is not None and args.srst2_resfinder_fg_output is not None:
        print("Running ARGANNOT/ResFinder presence/absence...")
        derive_presence_absence_targets_for_arg_res([args.srst2_argannot_fg_output, args.srst2_resfinder_fg_output])
    elif args.srst2_argannot_fg_output is not None:
        derive_presence_absence_targets_for_arg_res([args.srst2_argannot_fg_output])
    elif args.srst2_resfinder_fg_output is not None:
        derive_presence_absence_targets_for_arg_res([args.srst2_resfinder_fg_output])

    # Get alleles
    if args.srst2_gbs_cs_output is not None:
        print("Running GBS variants...")
        get_variants(args.srst2_gbs_cs_output)

    write_output(args.output)


def get_arguments():
    parser = argparse.ArgumentParser(description='Modify SRST2 sequence typing output files.')
    parser.add_argument('--srst2_gbs_fullgenes', dest='srst2_gbs_fg_output', required=False,
                        help='Input SRST2 fullgenes output for the GBS reference database.')
    parser.add_argument('--srst2_gbs_consensus', dest='srst2_gbs_cs_output', required=False,
                        help='Input freebayes consensus sequence output for the GBS reference database.')
    parser.add_argument('--srst2_argannot_fullgenes', dest='srst2_argannot_fg_output', required=False,
                        help='Input SRST2 fullgenes output for the ARG-ANNOT reference database.')
    parser.add_argument('--srst2_resfinder_fullgenes', dest='srst2_resfinder_fg_output', required=False,
                        help='Input SRST2 fullgenes output for the ResFinder reference database.')
    parser.add_argument('--output_prefix', dest='output', required=True,
                        help='Output prefix of filename.')

    return parser


def main():
    args = get_arguments().parse_args()
    run(args)


if __name__ == "__main__":
    sys.exit(main())
