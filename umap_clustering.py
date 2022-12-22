import pandas as pd
import umap

from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, silhouette_score

from time import time


def umap_clustering_best(sample, df_data, n_clusters=20):
    """
    Best clustering method
    :param sample: MPFs
    :param df_data: NSC, SMILES
    :param n_clusters: Default 7 clusters
    :return: [NSC, SMILES, CLuster_ID]
    """
    t0 = time()
    x_red = umap.UMAP(n_neighbors=100, min_dist=0.0,
                      n_components=2, metric='jaccard',
                      random_state=42).fit_transform(sample)
    clustering = AgglomerativeClustering(linkage='ward', n_clusters=n_clusters)
    clustering.fit(x_red)
    tf = time() - t0
    # Assign cluster ID
    df_clusters = assign_cluster_id(df_data, clustering)
    # Metrics
    s1 = silhouette_score(x_red, clustering.labels_, metric='euclidean')
    c1 = calinski_harabasz_score(x_red, clustering.labels_)
    d1 = davies_bouldin_score(x_red, clustering.labels_)
    df_metrics = pd.DataFrame(data=[tf, s1, c1, d1],
                              columns=['Time', 'Silhouette', 'CH score', 'DB score'])
    return df_metrics, df_clusters


def assign_cluster_id(df_data, cluster_id):
    print('Cluster ID')
    df_data['Cluster_ID'] = cluster_id.labels_
    return df_data




