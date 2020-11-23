#!/usr/bin/env python
import argparse, sys, re

def in_line(line, target):
    if (bool(re.search('^@HD|^\@SQ.*' + target + '|^@PG', line))==True):
        return True
    elif (line.split('\t')[2] == target):
        return True
    else:
        return False


def write_sam_file(sam_file, target, out_sam):
    with open(out_sam, 'w') as out:
        with open(sam_file, 'r') as sam:
            for line in sam:
                if in_line(line, target):
                    out.write(line)


def get_arguments():
    parser = argparse.ArgumentParser(description='Get targets from sam file.')
    parser.add_argument('--sam_file', '-s', dest='sam', required=True,
                        help='Input sam file.')
    parser.add_argument('--target', '-t', dest='target', required=True,
                        help='Input target sequence name.')
    parser.add_argument('--output', '-o', dest='output', required=True,
                        help='Output sam file.')
    return parser


def main():
    args = get_arguments().parse_args()
    write_sam_file(args.sam, args.target, args.output)


if __name__ == "__main__":
    sys.exit(main())
