#!/usr/bin/env python
# coding: utf-8


import numpy as np
import pandas as pd
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser(description="Choose a subset of cells that pass a criteria")
    parser.add_argument("--data_dir", help = "Path to directory containing the combined raw data", required=True)
    parser.add_argument("--sample_table",help = "A tab-delimited table of sample names and its raw fastq file", required = True)
    parser.add_argument("--output_folder",help = "Name of folder for output file" , required=True)
    # Threshold values
    parser.add_argument("--sum",type=int,help = "Threshold value: The sum of UMIs in a cell should be equal or higher than this number, default:1", default = 1)
    parser.add_argument("--count",type=int,help = "Threshold value: The total number of genes expressed in a cell should be equal or higher than this number, default:1", default = 1)
    parser.add_argument("--exp",type=int,help = "Threshold value: The most highly expressed gene should have a number of UMIs equal or higher than this number, default:1", default = 1)
  
    return parser.parse_args()

def load_data(path):
    
    data = pd.read_pickle(path+"data_raw.pickle.gz")
    #data_magic = pd.read_pickle(path+"/data_magic_nd.gz")
  #  metadata = pd.read_pickle(path+"metadata.pickle.gz")
    #wells_cells = pd.read_table(path+"/wells_cells.txt")
    return data



def load_genes(path,file):
    filename = path+file
    with open (filename,"r") as f:
        content_list = f.read().splitlines()
    return content_list
    


# load data
def main():
    args = get_args()
    base_path = args.data_dir

    data = load_data(base_path)

    data_new = data.fillna(0)

    sum, count, exp = args.sum, args.count, args.exp
    new_df = subset(data_new,sum, count, exp)
    table_samples = pd.read_table(args.sample_table,sep = "\t", header = None,index_col = 0, squeeze = True)
    dict_samples = table_samples.to_dict()
    save_cells(new_df,dict_samples,args.output_folder)

def subset(data, sum, count, exp):

    df = pd.DataFrame(index = data.index)

    df['count'] = data.select_dtypes(np.number).gt(0).sum(axis=1)
    df['sum'] = data.sum(axis = 1)
    df['max'] = data.max(axis = 1)


    df['sample'] = df.index
    df['sample'] = df['sample'].str.split(".",1, expand = True)[1]
    df = df.set_index(pd.Series(df.index).str.split("-", expand = True)[0])


    df_sum = df[(df['max'] >= exp) & (df['count'] >= count) & (df['sum'] >= sum)]
    return df_sum


def save_cells(df_sum,dict_samples,output_folder):
    df_cells_coexpressed = (df_sum['sample']).replace(dict_samples)
    df_cells_coexpressed.to_csv(output_folder+"/infected_cells.tsv", sep = "\t",header=None)


if __name__ == "__main__":
    main()
  
