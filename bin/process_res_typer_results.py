#!/usr/bin/env python3
import argparse
import sys
import os
import re
import glob
import subprocess
from Bio.Seq import Seq
from collections import defaultdict

class nSeq(str): # Nucleotide sequence
    pass

class aSeq(str): # Amino acid sequence
    pass

# Drug Class Resistance dictionary
drugRes_Col = {
    'TET': 'neg',
    'EC': 'neg',
    'FQ': 'neg',
    'OTHER': 'neg',
}

# Gene to Drug Class Resistance lookup dictionary
geneToClass = {
    'CAT': 'OTHER',
    'ERMB': 'EC',
    'ERMT': 'EC',
    'FOSA': 'OTHER',
    'LNUB': 'EC',
    'LSAC': 'EC',
    'MEFA': 'EC',
    'MPHC': 'EC',
    'MSRA': 'OTHER',
    'MSRD': 'OTHER',
    'TETB': 'TET',
    'TETL': 'TET',
    'TETM': 'TET',
    'TETO': 'TET',
    'TETS': 'TET',
    'SUL2': 'OTHER',
    'PARC': 'FQ',
    'GYRA': 'FQ',
    '23S1': 'EC',
    '23S3': 'EC',
    'RPOBGBS-1': 'OTHER',
    'RPOBGBS-2': 'OTHER',
    'RPOBGBS-3': 'OTHER',
    'RPOBGBS-4': 'OTHER',
}

# Other Resistance Targets dictionary
Res_Targets = {
    'CAT': 'neg',
    'ERMB': 'neg',
    'ERMT': 'neg',
    'FOSA': 'neg',
    'LNUB': 'neg',
    'LSAC': 'neg',
    'MEFA': 'neg',
    'MPHC': 'neg',
    'MSRA': 'neg',
    'MSRD': 'neg',
    'TETB': 'neg',
    'TETL': 'neg',
    'TETM': 'neg',
    'TETO': 'neg',
    'TETS': 'neg',
    'SUL2': 'neg',
}

# GBS Resistance Targets dictionary
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

# GBS Gene Resistance Variants dictionary
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

# Reference sequence dictionary
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

# Gene to target name lookup dictionary
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


def codon2aa(codon):
    """Translate codons to amino acids"""

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


def extract_frame_aa(sequence, frame):
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


def six_frame_translate(seq_input, frame):
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


def update_presence_absence_target(gene, allele, depth, gbs_res_target_dict):
    """Update presence/absence for GBS Targets dictionary"""
    if depth >= MIN_DEPTH:

        for gene_name in gbs_res_target_dict.keys():
            if re.search(gene_name, allele):

                gbs_res_target_dict[gene_name] = "pos"


def derive_presence_absence_targets(input_file):
    """Find gene presence/absence for the GBS resistance database"""
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
                update_presence_absence_target(gene, allele, depth, GBS_Res_Targets)
    except IOError:
        print('Cannot open {}.'.format(input_file))


def update_presence_absence_target_for_arg_res(gene, allele, depth, drug_res_col_dict, res_target_dict):
    """Update presence/absence for Other Resistance Targets dictionary"""
    if depth >= MIN_DEPTH:

        other = 1
        for gene_name in res_target_dict.keys():
            if re.search(gene_name, "".join(re.split("[^a-zA-Z0-9]*", allele)).upper()):
                other = 0

                if res_target_dict[gene_name] == "neg":
                    res_target_dict[gene_name] = "pos"

                drugCat = geneToClass[gene_name]
                if drug_res_col_dict[drugCat] == "neg":
                    drug_res_col_dict[drugCat] = gene + '(' + allele + ')'
                else:
                    drug_res_col_dict[drugCat] = drug_res_col_dict[drugCat] + ':' + gene + '(' + allele + ')'

        if other:
            if drug_res_col_dict["OTHER"] == "neg":
                drug_res_col_dict["OTHER"] = gene + '(' + allele + ')'
            else:
                drug_res_col_dict["OTHER"] = drug_res_col_dict["OTHER"] + ":" + gene + '(' + allele + ')'


def derive_presence_absence_targets_for_arg_res(input_files):
    """Find gene presence/absence for other resistance databases"""
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
            print('Cannot open {}.'.format(input_file))


def find_mismatches(seq_diffs, query_Seq, ref_Seq):
    """Find mismatches between query and reference sequences"""
    for resi in range(len(query_Seq)):
        if query_Seq[resi] != ref_Seq[resi]:
            seq_diffs.append(ref_Seq[resi] + str(resi+1) + query_Seq[resi])
    return seq_diffs


def get_seq_diffs(query_Seq, ref_Seq):
    """Get SNP variants"""
    if type(ref_Seq) == aSeq:
        query_Seq = six_frame_translate(query_Seq, 1)
    seq_diffs = []
    if query_Seq != ref_Seq:
        seq_diffs = find_mismatches(seq_diffs, query_Seq, ref_Seq)
    return seq_diffs


def update_GBS_Res_var(gene_name, seq_diffs, bin_res_arr):
    """Update GBS Gene Resistance Variants dictionary with gene variants"""
    if seq_diffs:
        bin_res_arr[gene_name] = gene_name + '-' + ','.join(seq_diffs)
    else:
        bin_res_arr[gene_name] = gene_name


def update_drug_res_col_dict(gene_name, seq_diffs, drugRes_Col, geneToClass):
    """Update Drug Resistance Class dictionary with GBS gene variants"""
    drugClass = geneToClass[gene_name]
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
    """Get consensus sequences from FASTA file into dictionary"""
    consensus_seq_dict = defaultdict(lambda: '')
    try:
        with open(consensus_seqs_file, 'r') as fd:
            for line in fd:
                if line[0] == '>':
                    seq_name = line.split('>')[1].split('\n')[0]
                else:
                    if consensus_seq_dict[seq_name] == '':
                        consensus_seq_dict[seq_name] = line.split('\n')[0]
                    else:
                        tmp_seq = consensus_seq_dict[seq_name] + line.split('\n')[0]
                        consensus_seq_dict[seq_name] = tmp_seq
    except IOError:
        print('Cannot open {}.'.format(consensus_seqs_file))

    return consensus_seq_dict


def get_gene_names_from_consensus(consensus_seq_dict):
    """Get the gene name identifiers from the consensus sequence IDs"""
    gene_names = []
    for gene_name in geneToTargetSeq.keys():
        if geneToTargetSeq[gene_name] in consensus_seq_dict.keys():
            gene_names.append(gene_name)
    return gene_names


def get_variants(consensus_seqs):
    """Get resistance gene variants from freebayes consensus GBS sequences"""
    consensus_seq_dict = get_consensus_seqs(consensus_seqs)
    gene_names = get_gene_names_from_consensus(consensus_seq_dict)
    for gene_name in gene_names:
        if GBS_Res_Targets[gene_name] == "pos" and geneToTargetSeq[gene_name] and geneToRef[gene_name]:
            seq_diffs = get_seq_diffs(consensus_seq_dict[geneToTargetSeq[gene_name]], geneToRef[gene_name])
            update_GBS_Res_var(gene_name, seq_diffs, GBS_Res_var)
            update_drug_res_col_dict(gene_name, seq_diffs, drugRes_Col, geneToClass)


def write_output(content, output_filename):
    """Write table content to output file"""
    try:
        with open(output_filename, 'w') as out:
            out.write(content)
    except IOError:
        print('Cannot open filename starting "{}"'.format(output_filename))


def create_output_contents(final_dict):
    """Create tab-delimited table from dictionary"""
    final = sorted(final_dict.items(), key=lambda item: item[0], reverse=False)
    content = ''
    for n, item in enumerate(final):
        if n == len(final)-1:
            content += item[0] + '\n'
        else:
            content += item[0] + '\t'
    for n, item in enumerate(final):
        if n == len(final)-1:
            content += item[1] + '\n'
        else:
            content += item[1] + '\t'
    return content


def run(args):

    # Set minimum read depth
    global MIN_DEPTH
    MIN_DEPTH = args.min_depth

    # Get presence/absence of genes
    derive_presence_absence_targets(args.srst2_gbs_fg_output)

    if args.srst2_other_fg_output is not None:
        derive_presence_absence_targets_for_arg_res(args.srst2_other_fg_output)
        Res_Targets.update(GBS_Res_Targets)
        inc_out = create_output_contents(Res_Targets)
    else:
        inc_out = create_output_contents(GBS_Res_Targets)

    # Get variants
    get_variants(args.srst2_gbs_cs_output)
    var_out = create_output_contents(GBS_Res_var)

    # Get alleles for all drug classes
    allele_out = create_output_contents(drugRes_Col)

    # Write incidence output
    write_output(inc_out, args.output + '_res_incidence.txt')
    # Write gbs variant output
    write_output(var_out, args.output + "_res_gbs_variants.txt")
    # Write allele output
    write_output(allele_out, args.output + "_res_alleles.txt")


def get_arguments():
    parser = argparse.ArgumentParser(description='Modify SRST2 sequence typing output files.')
    parser.add_argument('--srst2_gbs_fullgenes', dest='srst2_gbs_fg_output', required=True,
                        help='Input SRST2 fullgenes output for the GBS reference database.')
    parser.add_argument('--srst2_gbs_consensus', dest='srst2_gbs_cs_output', required=True,
                        help='Input freebayes consensus sequence output for the GBS reference database.')
    parser.add_argument('--srst2_other_fullgenes', dest='srst2_other_fg_output', required=False,
                        help='Input SRST2 fullgenes outputs for other references databases.',
                        nargs = '*')
    parser.add_argument('--min_read_depth', dest='min_depth', required=True, type=int, default=30,
                        help = 'Minimum read depth where mappings with fewer reads are excluded.')
    parser.add_argument('--output_prefix', dest='output', required=True,
                        help='Output prefix of filename.')

    return parser


def main():
    args = get_arguments().parse_args()
    run(args)


if __name__ == "__main__":
    sys.exit(main())
