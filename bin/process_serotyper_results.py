#!/usr/bin/env python
import argparse, sys


def write_line(gene, gene_dict, out):
    serotype = gene_dict[gene]
    status = 'identical'
    if serotype != '':
        status = 'imperfect'
        out.write(serotype[1]+'\t'+serotype[0]+'='+status+'\t'+serotype[0]+'\t'+serotype[3]+'\n')
    del gene_dict[gene]

def make_gene_dict(input_file, depth_threshold):
    gene_dict = dict()
    with open(input_file, 'r') as fg_file:
        next(fg_file) # Skip header row
        for line in fg_file:
            feature = line.split('\t')
            if float(feature[5]) > depth_threshold:
                gene_dict[feature[2]] = feature[2:-1]
    return gene_dict

def write_outfile(gene_dict, out_file):
    with open(out_file, 'w') as out:
        out.write('Matched_Allele'+'\t'+'Match_Type'+'\t'+'Serotype'+'\t'+'AvgDepth'+'\n')
        if 'III' in gene_dict.keys():
            write_line('III', gene_dict, out)
        if 'II' in gene_dict.keys():
            write_line('II', gene_dict, out)
        for key in gene_dict:
            write_line(key, gene_dict, out)

def get_arguments():
    parser = argparse.ArgumentParser(description='Modify fullgenes output of SRST2.')
    parser.add_argument('--input', '-i', dest='input', required=True,
                        help='Input fullgenes results tab file.')
    parser.add_argument('--output', '-o', dest='output', required=True,
                        help='Output filename.')
    parser.add_argument('--depth-threshold', '-d', dest='depth', default = 10, type=int,
                        help='Include SNPs greater than depth threshold. Default: 10.')

    return parser

def main():
    args = get_arguments().parse_args()
    gene_dict = make_gene_dict(args.input, args.depth)
    write_outfile(gene_dict, args.output)


if __name__ == "__main__":
    sys.exit(main())
