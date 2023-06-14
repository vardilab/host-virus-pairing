#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import os
from argparse import ArgumentParser

def get_args():
    parser = ArgumentParser(description="Combine pickle files into one file")
    parser.add_argument("--data_dir", help="Path to directory containing the folders of the 10X outputs", required=True)
    parser.add_argument("--list_folders",help = "List of folders containing the 10X files (in pickle format); write the namer of the parent folder, not the '/outs/filtered_feature_bc_matrix' path", required = True)
    parser.add_argument("--output_folder",help = "Name of folder for combined output file" ,default= "Combined_R")
    
    parser.add_argument('--raw', action='store_true')
    parser.add_argument('--not-raw', dest='raw', action='store_false')
    parser.set_defaults(raw=True)
    return parser.parse_args()

def load_data(path,file,raw):

    data_dir = path+ file + "/outs/filtered_feature_bc_matrix"
    if raw:
        data = pd.read_pickle(os.path.join(data_dir,"data_raw.pickle.gz"))
        return(data)
    else:
        data = pd.read_pickle(os.path.join(data_dir,"data.pickle.gz"))
        metadata = pd.read_pickle(os.path.join(data_dir,"metadata.pickle.gz"))
        return(data,metadata)

def load_list(list_samples,path,raw):
    l_data = []
    l_metadata = []
    for file in list_samples:
        # load data
        print("Loading",file+"...")
        if raw:
            data = load_data(path,file,raw)
        else:
            data,metadata = load_data(path,file,raw = False)
            metadata.index =  metadata.index + '.' + file
        # add sample column to metadata
            metadata['sample'] = file
            l_metadata.append(metadata)
        # set new index as a combination of index and sample name
        data.index =  data.index + '.' + file

        l_data.append(data)
    if raw:
        return(l_data)
    else:
        return(l_data,l_metadata)

def main():

    args = get_args()
    raw = args.raw
    list_samples = [str(item) for item in args.list_folders.split(',')]

    path = args.data_dir
    

    out_folder = args.output_folder
    if raw:
        print("Processing raw data...")
        list_data = load_list(list_samples,path,raw)
        # combine 
        combined_data = pd.concat(list_data,axis=0)
        combined_data = combined_data.fillna(0)
        save_files_raw(path,out_folder,combined_data)
    else:
        print("Processing data...")
        list_data,list_metadata = load_list(list_samples,path,raw = False)
        # combine 
        combined_data = pd.concat(list_data,axis=0)
        combined_data = combined_data.fillna(0)
        combined_metadata = pd.concat(list_metadata)
        save_files(path,out_folder,combined_data,combined_metadata)



def save_files_raw(path,out_folder,combined_data):
    os.mkdir(path + "/" + out_folder)
    print(combined_data.shape)
    combined_data.to_pickle(path + "/" + out_folder + "/data_raw.pickle.gz")
  

def save_files(path,out_folder,combined_data,combined_metadata):
    os.mkdir(path + "/" + out_folder)
    print(combined_data.shape)
    combined_data.to_pickle(path + "/" + out_folder + "/data.pickle.gz")
    combined_metadata.to_pickle(path + "/" + out_folder + "/metadata.pickle.gz")


if __name__ == "__main__":
    main()
  
    

