#! /usr/bin/env python

import os
from argparse import ArgumentParser

def create_master_dict(sample, fname):
    with open(fname, 'r') as masterfile:
        x = 0
        master_dict = {}
        for line in masterfile:
            field = line.rsplit()
            if not line[0] == 'ins_chr':
                coords = '\t'.join(field[:5])
                master_dict[x] = {'coords': coords, 'accessions': [sample]}
                x += 1
        return master_dict


def merge_deletions(master, fname, sample):
    with open(fname, 'r') as infile:
        for line in infile:
            field = line.rsplit()
            coords = '\t'.join(field[:5])
            i = len(master)-1
            x = 0
            while x <= i:
                if master[x]['coords'] == coords:
                    master[x]['accessions'].append(sample)
                    break
                elif x == i:
                    master[x+1] = {'coords': coords, 'accessions': [sample]}
                    break
                else:
                    x += 1


def save_deletions(master, outf):
    with open(outf, 'w+') as outfile:
        for key, value in master.iteritems():
            accessions = set(value['accessions'])
            outfile.write('{c}\t{a}\n'.format(c=value['coords'], a=','.join(accessions)))

def get_name_from_filename(filename):
    return os.path.basename(filename).rsplit('.', 1)[0]

if __name__ == "__main__":

    parser = ArgumentParser(description='Merge TE deletions calls')
    parser.add_argument('-o', '--output', help="File to write merged deletions to.", required=True)
    parser.add_argument('-i', '--input', help='all files that should be merged', nargs="+", required=True)
    options = parser.parse_args()

    first_file = options.input[0]
    first_samplename = get_name_from_filename(first_file)
    master_dictionary = create_master_dict(first_samplename, first_file)
    for filename in options.input[1:]:
        samplename = get_name_from_filename(filename)
        merge_deletions(master_dictionary, filename, samplename)
    save_deletions(master_dictionary, options.output)
