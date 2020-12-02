#!/usr/bin/env python3
import argparse, sys

def write_line(gene, gene_dict, out):
    """Write serotype to out stream of file"""
    serotype = gene_dict[gene]
    status = 'identical'
    if serotype[4] != '':
        status = 'imperfect'
    out.write(serotype[1]+'\t'+serotype[0]+'='+status+'\t'+serotype[0]+'\t'+serotype[3]+'\n')

def make_gene_dict(input_file, depth_threshold):
    """Get features from SRST2 input file into dictionary depending on read depth threshold"""
    gene_dict = dict()
    with open(input_file, 'r') as fg_file:
        next(fg_file) # Skip header row
        for line in fg_file:
            feature = line.split('\t')
            if float(feature[5]) > depth_threshold:
                gene_dict[feature[2]] = feature[2:-1]
    return gene_dict

def write_outfile(gene_dict, out_file):
    """Write serotype, match type status and average read depth to output file"""
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
    parser.add_argument('--srst2_output', '-s', dest='id', required=True,
                        help='Input fullgenes results tab file.')
    parser.add_argument('--sero_db', '-b', dest='db', required=True,
                        help='Input fullgenes results tab file.')
    parser.add_argument('--output', '-o', dest='output', required=True,
                        help='Output filename.')
    parser.add_argument('--depth-threshold', '-d', dest='depth', default = 10, type=int,
                        help='Include SNPs greater than depth threshold. Default: 10.')

    return parser

def main():
    args = get_arguments().parse_args()
    db_name = ' '.join(args.db.split('.')[:-1])

    # Specift fullgenes file path from ID and database name
    fullgenes_file = args.id + '__fullgenes__' + db_name + '__results.txt'

    # Get feature dictionary
    gene_dict = make_gene_dict(fullgenes_file, args.depth)

    # Write tab-delimited output file with serotype features
    write_outfile(gene_dict, args.output)


if __name__ == "__main__":
    sys.exit(main())
