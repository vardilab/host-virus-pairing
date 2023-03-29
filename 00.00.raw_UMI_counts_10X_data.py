#!/bin/python3

import os
from argparse import ArgumentParser
import pandas as pd
import scprep

def main():
    args = get_args()
    print('Loading data')
    data = load_data(args.data_dir, args.file_type)
    # remove clutter from gene names
    new_columns = pd.Series(data.columns).str.split('\t', expand=True)[0]
    data.columns = new_columns.to_list()
    print('Removing empty cells/features')
    data = data[data.sum(axis=1) > 0]
    data = data.T[data.sum() > 0].T
    print('Data shape (cells, features)', data.shape)
    print('Saving raw data file')
    data.to_pickle(os.path.join(args.data_dir, "data_raw.pickle.gz"))

def get_args():
    parser = ArgumentParser(description="Generate raw UMI counts for 10X single-cell data")
    parser.add_argument("--data_dir", help="Path to direcory with 10X outputs: matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)", required=True)
    parser.add_argument("--file_type", help="Note if the input count matrix is in format mtx, csv, or tsv (default: mtx)", default='mtx')
    #parser.add_argument("--out_file", help="path to file for saving output", required=True)
    return parser.parse_args()

def load_data(datadir, filetype):
    """
    :type datadir: object
    """
    filetypes = ['10X', '10X_HDF5', '10X_zip', 'csv', 'fcs', 'mtx', 'tsv']
    assert (filetype in filetypes), "Please enter file format from the following {}".format(filetypes)

    if filetype == 'mtx':
        data = scprep.io.load_mtx(os.path.join(datadir, 'matrix.mtx.gz'),
                                  gene_names=os.path.join(datadir, 'features.tsv.gz'),
                                  cell_names=os.path.join(datadir, 'barcodes.tsv.gz'),
                                  cell_axis="column")

    elif filetype == 'csv':
        data = scprep.io.load_csv(os.path.join(datadir, 'data.csv'))

    elif filetype == 'tsv':
        data = scprep.io.load_csv(os.path.join(datadir, 'data.tsv'))

    return data

if __name__ == '__main__':
    main()
