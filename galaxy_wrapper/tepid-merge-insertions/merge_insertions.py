from argparse import ArgumentParser
import os
import tempfile
import pandas as pd
import pybedtools

COLUMNS = ['ins_chrom',
           'ins_start',
           'ins_end',
           'ref_chrom',
           'ref_start',
           'ref_end',
           'agi',
           'accession',
           'cluster']

def overlap(start1, stop1, start2, stop2, d=50):
    """returns True if sets of coordinates overlap. Assumes coordinates are on same chromosome"""
    return start1 <= stop2+d and stop1 >= start2-d


def merge(i, insertion, result):
    if len(result) == 0:
        result[i] = insertion
    else:
        if not can_merge(insertion, result):
            result[i] = insertion


def can_merge(insertion, result):
    """
    Merges insertions and returns True if all requirements are met
    """
    for j, master_insertion in result.items():
        if insertion['agi'] & master_insertion['agi']:
            if overlap(master_insertion['ins_start'], master_insertion['ins_end'], insertion['ins_start'],insertion['ins_end']):
                # Adjusting the insertion start (doesn't really do anything?!)
                if len(insertion['agi']) < len(master_insertion['agi']):
                    ref_start = master_insertion['ref_start']
                else:
                    ref_start = insertion['ref_start']
                if master_insertion['ins_chrom'] == insertion['ins_chrom'] and insertion['ref_chrom'] == master_insertion['ref_chrom'] and ref_start == master_insertion['ref_start']:
                   result[j]['accession'] = result[j]['accession'] | (insertion['accession'])
                   result[j]['agi'] = result[j]['agi'] | (insertion['agi'])
                return True
    return False


def inner_merge(s):
    result = {}
    for i, insertion in s.items():
        merge(i, insertion, result)
    return result.values()


def reduce_and_cluster(inputfiles):
    """
    Read in inputfiles using pandas, write additional column with sample identifier,
    sort and cluster using pybedtools and return dataframe.
    """
    usecols = [0,1,2,3,4,5,6]  # skip col 7, which contains the read support id
    tables = [pd.read_table(f, header=None) for f in inputfiles]
    sample_ids = [os.path.basename(f).rsplit('.')[0] for f in inputfiles]
    for sample_id, df in zip(sample_ids, tables):
        df[7] = sample_id
    merged_table = pd.concat(tables)
    tfile = tempfile.NamedTemporaryFile()
    merged_table.to_csv(tfile, sep="\t", header=None, index=False)
    tfile.flush()
    bedfile = pybedtools.BedTool(tfile.name).sort().cluster(d=50)
    df = bedfile.to_dataframe()
    df.columns = COLUMNS
    # Split comma separated agi values and make set
    df['agi'] = [set(v.split(',')) for v in df['agi'].values]
    df['accession'] = [set(str(v).split(',')) for v in df['accession'].values]
    return df


def split_clusters(df):
    """
    clusters as defined by bedtools allow for 50 nt distance. This means that
    clusters can be many kb large, so we check each individual insertion in
    the cluster against the other insertions. We split the clusters based on
    whether the overlap and TE identity criteria are fulfilled (so a
    different TE would lead to a split in the clusters)
    """
    groups = df.groupby('cluster')
    nested_list = [inner_merge(group.transpose().to_dict()) for _, group in groups]
    return pd.DataFrame([i for n in nested_list for i in n])[COLUMNS]


def write_output(df, output):
    # Turn sets back to comma-separated values
    df['agi'] = [",".join(agi) for agi in df['agi']]
    df['accession'] = [",".join(acc) for acc in df['accession']]
    df.to_csv(output, sep="\t",header=None, index=None)


def main(inputfiles, output):
    df = reduce_and_cluster(inputfiles)
    df = split_clusters(df)
    write_output(df, output)


if __name__ == "__main__":
    parser = ArgumentParser(description='Merge TE insertions calls')
    parser.add_argument('-o', '--output', help='output file', required=True)
    parser.add_argument('-i', '--input', help='Insertion files to merge', nargs="+", required=True)
    options = parser.parse_args()

    main(inputfiles=options.input, output=options.output)
