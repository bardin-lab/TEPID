#! /usr/bin/python

from argparse import ArgumentParser
import sys

parser = ArgumentParser(description='Invert the TE deletion calls to give a consistent data format between TE insertions and deletions')
parser.add_argument('-s', '--samples', help='list of all sample names', nargs="+", required=True)
parser.add_argument('-d', '--deletions', help='merged TEPID deletions', required=True)
parser.add_argument('-r', '--reference', help='reference sample name, eg Col-0', required=True)
parser.add_argument('-o', '--output', help='output file name', required=True)
options = parser.parse_args()


def filter_del(options):
    with open(options.deletions, 'r') as dels, open(options.output, 'w+') as outfile:
        sample_names = options.samples
        for line in dels:
            line = line.strip().split('\t')
            accessions = line[5]
            sys.stderr.write(accessions)
            sys.stderr.write(",".join(sample_names))
            accessions = accessions.split(',')
            coords = line[:4]
            temp = [options.reference]
            te = line[4]
            for sample in sample_names:
                if sample not in accessions:
                    temp.append(sample)
                else:
                    pass
            coords.pop(3)  # remove strand
            info = '\t'.join(coords) + '\t' + te + '\t' + ','.join(temp) + '\n'
            outfile.write(info)


if __name__ == "__main__":
    filter_del(options)
