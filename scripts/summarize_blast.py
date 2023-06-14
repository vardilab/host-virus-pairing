#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from argparse import ArgumentParser



def load_data(path,file,db):
    list_columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore".split(" ")
    print("reading file " + file + "...")
    data = pd.read_table(path+file, sep = "\t", names = list_columns)
    #data['cell'] = file.replace(".metapr2.tsv","")
    data['cell'] =  data['qseqid'].str.split("i\d+_", expand = True)[1]
    # in case of pr2 open:
    if db == "PR2":
        data[["sseqid","rrna_type","rrna_loc","clone","Domain","Supergroup","Division","Class","Order","Family","Genus","Species"]] = data['sseqid'].str.split("|",expand = True)
    # in case of metapr2 open:
    elif db == "metaPR2":
        data[["sseqid","Supergroup","Division","Class","Order","Family","Genus","Species"]] = data['sseqid'].str.split("|",expand = True)
    return data


def subset_results(sample_data):

    print("filtering by length and identity...")
    sample_data = sample_data[(sample_data['pident'] >= 99)]
    sample_data = sample_data[(sample_data['length'] >= 100)]
    print("sorting data...")
    sample_data = sample_data.sort_values(by = ['qseqid','bitscore'],ascending = False)

    print("Choosing only best hit...")
    sample_data['rank'] = sample_data.groupby(['qseqid'])['bitscore'].rank("min",ascending = False)
    final_run = sample_data[sample_data['rank'] == 1]
    return final_run



def replace(data,level, name, new_name):
    group = set(data[data[level] == name]['Family'])
    value_to_replace = {x : new_name for x in group}
    return value_to_replace


def get_args():
    parser = ArgumentParser(description="Summarize blast results: pick only the best hit with identity of >= 99 percent and an alignment length of >100 bp")
    parser.add_argument("--data_dir", help = "Path to directory containing the blast results", required=True)
    parser.add_argument("--database",help = "The name of 18s rRNA database (either 'PR2' or 'metaPR2')", required = True)
    
    return parser.parse_args()

def main():
    args = get_args()
    path = args.data_dir
    db = args.database
    file = '/all_cells.transcripts.edit.' + db + '.tsv'
    df = load_data(path,file,db)
    df_dedup = subset_results(df)


    # annotate cells based on their best contig
    idx = df_dedup.groupby(['cell'])['bitscore'].transform(max) == df_dedup['bitscore']
    df_top = df_dedup[idx]
    df_top['Annotation'] = df_top['Family'] 
    df_dedup.to_csv(path+"/transcripts.summary.99." + db + ".tsv", sep = "\t")

if __name__ == "__main__":
    main()
  