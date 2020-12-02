#!/usr/bin/env python3
from __future__ import print_function
import argparse, sys


def write_output(lines, output):
    with open(output, 'w') as out:
        out.write(lines)


def get_content_with_id(id, file):
    content = 'ID' + '\t'
    count = 0
    with open(file, 'r') as f:
        for line in f:
            if count:
                content += id + '\t' + line
            else:
                content += line
            count += 1
    return content


def get_sero_res_contents(id, sero_file, res_file):
    header = 'ID' + '\t' + 'Serotype' + '\t'
    type = ''
    res_incidence = ''
    with open(sero_file, 'r') as sero:
        next(sero) # Skip header row
        count = 0
        for line in sero:
            if count:
                type = type + ";" + line.split('\t')[2]
            else:
                type = id + '\t' + line.split('\t')[2]
            count += 1

    with open(res_file, 'r') as res:
        count = 0
        for line in res:
            if count:
                results = line.split('\n')[0].split('\t')
                results_out = []
                for result in results:
                    if result == 'pos':
                        results_out.append('+')
                    elif result == 'neg':
                        results_out.append('-')
                res_incidence = '\t'.join(results_out)
            else:
                header += line
            count += 1

    return(header + type + '\t' + res_incidence + '\n')


def get_arguments():
    parser = argparse.ArgumentParser(description='Combine resistance gene typer and serotyper results.')
    parser.add_argument('--id', '-i', dest='id', required=True,
                        help='Sample ID')
    parser.add_argument('--serotyper_results', '-s', dest='sero', required=True,
                        help='Input SeroType results tab file.')
    parser.add_argument('--res_incidence_results', '-r', dest='inc', required=True,
                        help='Input resistance typing incidence results file.')
    parser.add_argument('--res_alleles_results', '-a', dest='alleles', required=True,
                        help='Input resistance typing alleles results file.')
    parser.add_argument('--res_variants_results', '-v', dest='variants', required=True,
                        help='Input resistance typing variants results file.')
    parser.add_argument('--output', '-o', dest='output', required=True,
                        help='Output prefix.')
    return parser


def main():
    args = get_arguments().parse_args()
    sero_res_output_lines = get_sero_res_contents(args.id, args.sero, args.inc)
    write_output(sero_res_output_lines, args.output + "_sero_res_incidence.txt")

    res_alleles_output_lines = get_content_with_id(args.id, args.alleles)
    write_output(res_alleles_output_lines, args.output + "_id_alleles.txt")

    res_variants_output_lines = get_content_with_id(args.id, args.variants)
    write_output(res_variants_output_lines, args.output + "_id_variants.txt")


if __name__ == "__main__":
    sys.exit(main())