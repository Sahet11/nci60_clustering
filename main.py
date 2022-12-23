from load_dataset import get_data_cluster, check_output_folder
from umap_clustering import umap_clustering_best


def main(name_file_nci60, name_file_smi, dir_file, n_clusters=7, outliers=True):
    print('Best clustering method')
    # load NCI-60 data
    df_data, x_sample = get_data_cluster(name_file_nci60, name_file_smi, dir_file, outliers)
    print('\nClustering...')
    df_metrics, df_clustering = umap_clustering_best(sample=x_sample, df_data=df_data, n_clusters=n_clusters)
    # Save results
    print('Saving output')
    if outliers:
        dir_out = dir_file + '/Results/non_outliers_molecules/'
    else:
        dir_out = dir_file + '/Results/all_molecules/'
    check_output_folder(dir_out)
    file_metrics = dir_out + 'clustering_quality_metrics_k' + str(n_clusters) + '.csv'
    df_metrics.to_csv(file_metrics)
    file_cluster = dir_out + 'clustering_id_k' + str(n_clusters) + '.csv'
    df_clustering.to_csv(file_cluster)
    return None


dir_main = '/home/hernadez/Documents/data_test_biomolecules'
file_nci60 = 'CANCER60GI50.LST'   # _Oct2020
file_smiles = 'Chem2D_Jun2016.smi'

main(name_file_nci60=file_nci60, name_file_smi=file_smiles, dir_file=dir_main, n_clusters=7)







