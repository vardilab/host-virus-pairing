#!/bin/python3

from argparse import ArgumentParser
import pandas as pd
import phate
import umap
import scprep
import os
import tasklogger
from sklearn.manifold import TSNE

def main():
    args = get_args()
    print('Loading metadata')
    metadata = pd.read_pickle(os.path.join(args.data_dir, args.metadata_file))
    if args.pca_data == -1:
        print('Loading data')
        data = pd.read_pickle(os.path.join(args.data_dir, args.data_file))
        print("Calculating PCA")
        data_pca, metadata = calc_pca(data, metadata, int(args.pca_components))
        data_pca.to_pickle(os.path.join(args.data_dir, args.data_pca_output))
    else:
        print("Loading PCA data")
        data_pca = pd.read_pickle(os.path.join(args.data_dir, args.pca_data))
        try:
            metadata = pd.concat([metadata, data_pca[['PC1','PC2','PC3']]], axis=1)
        except:
            metadata = pd.concat([metadata, data_pca[['PC1','PC2','PC3']]], axis=1,suffixes='.scprep')
        
    metadata = calc_tsne(data_pca, metadata, int(args.tsne_perplexity), int(args.pca_components))
    metadata = calc_umap(data_pca, metadata, int(args.pca_components),args.min_dist,args.spread)
    metadata = calc_phate(data_pca, metadata, int(args.pca_components))
    print('Saving dimentionality reduction data')
    metadata.to_pickle(os.path.join(args.data_dir, args.metadata_output))
            

def get_args():
    parser = ArgumentParser(description="Dimentionality reduction for single-cell data.")
    parser.add_argument("--data_dir", help="path to direcory with preprocessed single-cell data, including file names: data.pickle.gz and metadata.pickle.gz", required=True)
    parser.add_argument("--data_file", help="name of data file as appear in data_dir (default: data.pickle.gz)", default='data.pickle.gz')
    parser.add_argument("--pca_data", help="name of the pca data to use for dimensionality reduction. If not mentioned, it will produce a new dataset", default = -1)
    parser.add_argument("--metadata_file", help="name of metadata file as appear in data_dir (default: metadata.pickle.gz)", default='metadata.pickle.gz')
    parser.add_argument("--pca_components", help="Number of components to calculate for the PCA (default: 50)", default=50)
    parser.add_argument("--tsne_perplexity", help="float. The perplexity is related to the number of nearest neighbors that is used in other manifold learning algorithms. Larger datasets usually require a larger perplexity. Consider selecting a value between 5 and 50. Different values can result in significantly different results (default: 30)", default=30)
    parser.add_argument("--metadata_output", help="name of output metadata file (default: metadata_dimentionality_reduction.pickle.gz)", default='metadata_dimentionality_reduction.pickle.gz')
    parser.add_argument("--data_pca_output", help="name of output pca data file (default: data_pca.pickle.gz)", default='data_pca.pickle.gz')
    parser.add_argument("--min_dist", type=float, help="value of the min_dist variable in UMAP (default: 0.1)", default=0.1)
    parser.add_argument("--spread", type=float, help="value of the spread variable in UMAP  (default: 1)", default=1)
    return parser.parse_args()

def calc_pca(data, metadata, pca_components):
    with tasklogger.log_task('PCA components on {} cells'.format(metadata.shape[0])):
        data_pca = scprep.reduce.pca(data, n_components=pca_components, method='svd',eps=0.1, seed = 2023)
        try:
            metadata = pd.concat([metadata, data_pca[['PC1','PC2','PC3']]], axis=1)
        except:
            metadata = pd.concat([metadata, data_pca[['PC1','PC2','PC3']]], axis=1,suffixes='.scprep')
        #try:
        #    metadata = metadata.join(data_pca[['PC1','PC2','PC3']])
        #except:
        #    metadata = metadata.join(data_pca[['PC1','PC2','PC3']], rsuffix='.scprep')
        return data_pca, metadata

def calc_tsne(data_pca, metadata, tsne_perplexity, pca_components):
    with tasklogger.log_task('t-SNE on {} cells'.format(data_pca.shape[0])):
        # Fitting tSNE. Change the perplexity here.
        tsne_op = TSNE(n_components=3, perplexity=tsne_perplexity)
        data_tsne = tsne_op.fit_transform(data_pca.iloc[:,:pca_components])
	# Put output into a dataframe
        data_tsne = pd.DataFrame(data_tsne, index=data_pca.index)
        print(data_tsne.shape)
        data_tsne.columns = ['TSNE1', 'TSNE2', 'TSNE3']
        try:
            metadata = pd.concat([metadata, data_tsne], axis=1)
         #   metadata = metadata.merge(data_tsne)
        except:
            metadata = pd.concat([metadata, data_tsne], axis=1, suffixes='.scprep')
          #  metadata = metadata.join(data_tsne, rsuffix='.scprep')
        return metadata

def calc_umap(data_pca, metadata, pca_components,min_dist, spread):
    with tasklogger.log_task('UMAP on {} cells'.format(data_pca.shape[0])):
        ## Calculate UMAP and add results to metadata file
    #    min_dist, spread = min_dist.astype(float), spread.astype(float)
        umap_op = umap.UMAP(min_dist = min_dist, spread = spread,n_components=3, n_neighbors = 7)
        data_umap = umap_op.fit_transform(data_pca.iloc[:,:pca_components])
        print(data_umap.shape)
        data_umap = pd.DataFrame(data_umap, index=data_pca.index)
        print(data_umap.shape,metadata.shape)
        data_umap.columns = ['UMAP1', 'UMAP2', 'UMAP3']
        try:
            metadata = pd.concat([metadata, data_umap], axis=1)
         #   metadata = metadata.join(data_umap)
        except:
            metadata = pd.concat([metadata, data_umap], axis=1, suffixes='.scprep')
         #   metadata = metadata.join(data_umap, rsuffix='.scprep')
        return metadata

def calc_phate(data_pca, metadata, pca_components):
    with tasklogger.log_task('PHATE on {} cells'.format(data_pca.shape[0])):
        ## Calculate PHATE and add results to metadata file
        phate_op = phate.PHATE(n_components=3)
        data_phate = phate_op.fit_transform(data_pca.iloc[:,:pca_components])
        print(data_phate.shape)
        data_phate = pd.DataFrame(data_phate, index=data_pca.index)
        print(data_phate.shape,metadata.shape)
        data_phate.columns = ['PHATE1', 'PHATE2', 'PHATE3']
        try:
            metadata = pd.concat([metadata, data_phate], axis=1)
        #    metadata = metadata.join(data_phate)
        except:
            metadata = pd.concat([metadata, data_phate], axis=1,suffixes='.scprep')
           # metadata = metadata.join(data_phate, rsuffix='.scprep')
        print(metadata.shape)
        return metadata

if __name__ == '__main__':
    main()

