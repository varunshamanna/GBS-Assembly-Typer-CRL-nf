#!/usr/bin/env python
import argparse
import sys
import re


# Protein to family lookup
proteinFamilies = {
    'ALP': 'ALPH',
    'RIB': 'ALPH',
    'SRR': 'SRR',
    'PI': 'PILI',
    'HVGA': 'HVGA',
}

featureCol = {
    "ALPH": "neg",
    "SRR": "neg",
    "PILI": "neg",
    "HVGA": "neg",
}

binFeatureCol = {
    "HVGA": "-",
    "PI1": "-",
    "PI2A1": "-",
    "PI2A2": "-",
    "PI2B": "-",
    "SRR1": "-",
    "SRR2": "-",
    "ALP1": "-",
    "ALP23": "-",
    "ALPHA": "-",
    "RIB": "-",
}


def update_protein_presence_absence(
        gene, allele, min_depth, depth, feature_col_dict, bin_feature_col_dict, protein_family_dict):
    """Update presence/absence"""

    if depth >= min_depth:

        for protein in protein_family_dict.keys():
            if re.search(protein, "".join(re.split("[^a-zA-Z0-9]*", allele)).upper()):

                family = protein_family_dict[protein]
                if feature_col_dict[family] == "neg":
                    feature_col_dict[family] = gene
                else:
                    feature_col_dict[family] = feature_col_dict[family] + ':' + gene

        for protein in bin_feature_col_dict.keys():
            if re.search(protein, "".join(re.split("[^a-zA-Z0-9]*", allele)).upper()):
                bin_feature_col_dict[protein] = "+"


def derive_presence_absence(input_file, min_depth, processor):
    """Find surface protein gene presence/absence for GBS surface database"""

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
                processor(gene, allele, min_depth, depth, featureCol, binFeatureCol, proteinFamilies)
    except IOError:
        print('Cannot open {}.'.format(input_file))


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
    """Top level run method"""
    db_name = ' '.join(args.db.split('.')[:-1])

    # Specify fullgenes file path from ID and database name
    fullgenes_file = args.fullgenes_file_id + '__fullgenes__' + db_name + '__results.txt'

    # Get presence/absence of genes
    derive_presence_absence(fullgenes_file, args.min_depth, update_protein_presence_absence)

    feature_out = create_output_contents(featureCol)
    bin_feature_out = create_output_contents(binFeatureCol)

    # Write incidence output
    write_output(feature_out, args.output + '_surface_protein_incidence_sample.txt')
    # Write gbs variant output
    write_output(bin_feature_out, args.output + "_surface_protein_variants_sample.txt")


def get_arguments():
    """Parse surface typer command line arguments"""
    parser = argparse.ArgumentParser(description='Modify fullgenes output of SRST2.')
    parser.add_argument('--srst2_gbs_fullgenes', '-s', dest='fullgenes_file_id', required=True,
                        help='Input fullgenes results tab file.')
    parser.add_argument('--surface_db', '-b', dest='db', required=True,
                        help='Surface proteins database file.')
    parser.add_argument('--output_prefix', '-o', dest='output', required=True,
                        help='Output prefix for filename.')
    parser.add_argument('--min_read_depth', '-d', dest='min_depth', default=30, type=float,
                        help='Minimum read depth where mappings with fewer reads are excluded. Default: 30.')
    return parser


def main():
    run(get_arguments().parse_args())


if __name__ == "__main__":
    sys.exit(main())
