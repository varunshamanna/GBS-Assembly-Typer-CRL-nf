#!/usr/bin/env python
from __future__ import print_function
import argparse, sys


def print_sample_line(id, sero_file, res_file):
    with open(sero_file, 'r') as sero:
        next(sero) # Skip header row
        count = 0
        for line in sero:
            if count:
                type = type + ";" + line.split('\t')[2]
            else:
                type = line.split('\t')[2]
            count += 1

    with open(res_file, 'r') as res:
        for line in res:
            res_incidence = '\t'.join(line.split(','))

    print(type + '\t' + res_incidence, end = '')


def get_arguments():
    parser = argparse.ArgumentParser(description='Combine resistance gene typer and serotyper results.')
    parser.add_argument('--id', '-i', dest='id', required=True,
                        help='Sample ID')
    parser.add_argument('--serotyper_results', '-s', dest='sero', required=True,
                        help='Input SeroType results tab file.')
    parser.add_argument('--res_typer_results', '-r', dest='res', required=True,
                        help='Input BIN results file.')
    return parser

def main():
    args = get_arguments().parse_args()
    print_sample_line(args.id, args.sero, args.res)


if __name__ == "__main__":
    sys.exit(main())
