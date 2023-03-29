#!/bin/python3

from argparse import ArgumentParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import sklearn
import sklearn.cluster
import sklearn.manifold

import scprep
import os
import tasklogger

import graphtools as gt
import phenograph
import louvain

def main():
    args = get_args()
    print('Loading data')
    data_pca = pd.read_pickle(os.path.join(args.data_dir, args.data_pca))
    metadata = pd.read_pickle(os.path.join(args.data_dir, args.metadata_file))
    ## Cluster cells with different algorithms
    phenograph_clusters, k = clust_phenogrph(data_pca)
    G = gt.Graph(data_pca)
    G_igraph = G.to_igraph()
    louvain_clusters = clust_louvain(data_pca, G_igraph)
    kmeans_clusters = clust_kmeans(data_pca, k)
    spectral_clusters = clust_spectral(data_pca, k, G)
    dim = args.dimensionality
    data_phate = metadata[['{}1'.format(dim),'{}2'.format(dim)]]
    ## Reordering clusters by PHATE coordinates
    clusterings = {'Phenograph':phenograph_clusters,
            'Louvain':louvain_clusters,
            'KMeans':kmeans_clusters,
            'Spectral':spectral_clusters}
    for alg in clusterings:
        cl_nu = scprep.utils.sort_clusters_by_values(clusterings[alg], data_phate.iloc[:,0])
        clusterings[alg] = cl_nu
    metadata = pd.concat([metadata,pd.DataFrame(clusterings, index=data_pca.index)],axis = 1)
    print('Saving clustering data')
    metadata.to_pickle(os.path.join(args.data_dir, str(args.metadata_file).replace(".pickle.gz","_clusters.pickle.gz")))

def get_args():
    parser = ArgumentParser(description="Clustering single-cell data.")
    parser.add_argument("--data_dir", help="path to direcory with single-cell PCA data, including file names: data_pca.pickle.gz and metadata_dim_reduction.pickle.gz", required=True)
    parser.add_argument("--data_pca", help="name of data file as appear in data_dir (default: data_pca.pickle.gz)", default='data_pca.pickle.gz')
    parser.add_argument("--metadata_file", help="name of metadata file as appear in data_dir (default: metadata_dim_reduction.pickle.gz)", default='metadata_dim_reduction.pickle.gz')
    parser.add_argument("--dimensionality", help="name of dimentionality reduction algorithm to cluster by (default: PHATE)", default='PHATE')
  
    return parser.parse_args()

def clust_phenogrph(data_pca):
    phenograph_clusters, _, _ = phenograph.cluster(data_pca)
    return phenograph_clusters, len(set(phenograph_clusters))

def clust_louvain(data_pca, G_igraph):
    with tasklogger.log_task("Louvain"):
        partition = louvain.find_partition(G_igraph, louvain.RBConfigurationVertexPartition, weights="weight", resolution_parameter=1)
        louvain_clusters = np.array(partition.membership)
    return louvain_clusters

def clust_kmeans(data_pca, k):
    with tasklogger.log_task("KMeans"):
        kmeans_clusters = sklearn.cluster.KMeans(n_clusters=k).fit_predict(data_pca)
    return kmeans_clusters

def clust_spectral(data_pca, k, G):
    with tasklogger.log_task("Spectral clustering"):
        spec_op = sklearn.cluster.SpectralClustering(n_clusters=k, affinity='precomputed')
        spectral_clusters = spec_op.fit_predict(G.K)
    return spectral_clusters

if __name__ == '__main__':
    main()
