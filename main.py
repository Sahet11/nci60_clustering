from load_dataset import load_data
from umap_clustering import umap_clustering_best


def main(name_file, n_clusters=7):
    print('Best clustering method')
    # load NCI-60 data
    x_sample, df_data = load_data(name_file=name_file)
    df_metrics, df_clustering = umap_clustering_best(sample=x_sample, df_data=df_data, n_clusters=n_clusters)
    # Save results
    df_metrics.to_csv('clustering_quality_metrics.csv')
    df_clustering.to_csv('clustering_id.csv')
    return None







