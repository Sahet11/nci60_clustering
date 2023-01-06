import os.path
import copy
import pandas as pd
import numpy as np

from rdkit import RDLogger
from rdkit import Chem
from rdkit.Chem import SaltRemover
from molvs import Standardizer
from rdkit.Chem import AllChem
from rdkit import DataStructs


def load_data(name_file, dir_file):
    file_in = dir_file + '/' + name_file
    df = pd.read_csv(file_in, delimiter=',')
    print('\nFile: ', name_file, '\n - ORIGINAL DATA - \n\t Shape: ', df.shape)
    print('\t Columns: ', df.columns.tolist())
    print('\t Number of cell lines: ', len(df.CELL.unique()))
    print('\t Unique NSC: ', len(df.NSC.unique()))
    print('\t Number of unique NSC - cell line pairs: ', df.groupby(['NSC', 'CELL']).size().shape)
    return df


def load_smiles_file(name_file, dir_file):
    file_in = dir_file + '/' + name_file
    df_smile = pd.read_csv(file_in, sep='\t', header=None)
    df_smile.columns = ['SMILE', 'NSC']
    df_smile = df_smile.iloc[:, [1, 0]]
    return df_smile


def get_mean(df):
    key_info = df["NSC"].map(str) + '_' + df['CELL']
    df.insert(0, 'ID', key_info)
    df_aux = df.loc[:, ['NLOGGI50', 'ID']]
    df_mean = df_aux.groupby('ID', as_index=False).mean()
    df_mean.columns = ['ID', 'NLOGGI50_N']
    return df_mean


def get_mfp(df_smiles):
    mol_list = df_smiles.SMILE.unique()
    fingerprints, fingerprints1 = [], []
    smiles_list = []
    std_smiles = []
    error_smiles = []
    RDLogger.DisableLog('rdApp.*')
    for m in mol_list:
        try:
            m1 = copy.copy(m)
            m = Chem.MolFromSmiles(m)
            remover = SaltRemover.SaltRemover()  # remove salt
            m = remover.StripMol(m)
            s = Standardizer()  # standardize molecule
            m = s.standardize(m)
            fp1 = AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024)
            fingerprints1.append(fp1)
            fingerprints.append(list(fp1))
            smiles_list.append(m1)
            std_smiles.append(Chem.MolToSmiles(m))
        except:
            error_smiles.append(m1)
    col_name = ['MFP_' + str(i) for i in range(1024)]
    df_mfp_aux = pd.DataFrame(data=fingerprints, columns=col_name)
    df_std_smile = pd.DataFrame(data=std_smiles, columns=['STD_SMILE'])
    df_smile_valid = pd.DataFrame(data=smiles_list, columns=['SMILE'])
    df_mfp = pd.concat([df_smile_valid, df_std_smile, df_mfp_aux], axis=1)
    return df_mfp, fingerprints1


def remove_few_data(df):
    to_remove = ['M19-MEL ', 'DMS 114 ', 'DLD-1 ', 'KM20L2 ', 'SNB-78 ', 'LXFL 529 ', 'DMS 273 ', 'XF 498 ', 'HOP-18 ',
                 'RXF-631 ', 'MDA-MB-468 ', 'SN12K1 ', 'P388 ', 'P388/ADR ', 'U-251/H.Fine ', 'A-172/H.Fine ',
                 'U-87/H.Fine ', 'T98G/H.Fine ', 'LN-229/H.Fine ', 'U118/H.Fine ', 'U-373/H.Fine ', 'SK-BR-3 ',
                 'MAXF 401 ', 'DB ', 'SW-156 ', 'COLO 741 ', 'SW-1573 ', 'UISO-BCA-1 ', 'COLO 746 ', 'CXF 264L ',
                 'MEXF 514L ', 'TK-164 ', 'HT ', 'RL ', 'H1299p53RE29 ', 'HT29p53RE22 ', 'MCF7-E6 ', 'RKO Waf1 ',
                 'RKOp53RE1 ', 'T47D FOS1 ', 'T47D NFkB15 ', 'HCT-116/CMV-2 ', 'HCT-116/E6-1 ', 'HCT-116/P ',
                 'HCT-116/P21/A ', 'HCT-116/P21/B ', 'T47D ERE4 ', 'HCT-116/CMV-1 ', 'HCT-116/PV ', 'HCT-116/P21/C ',
                 'TSU-PRI ', 'HCT-116/E6-2 ', 'JCA-1 ', 'ND-1 ', 'CACO-2 ', 'ES-2 ', 'MCF7/ATCC ', 'MLI-059 ',
                 'RPMI-7951 ', 'RXF 486L ', 'SF-767 ', 'SMS-KCNR ', 'SW 1088 ', 'UABLG22 ', 'UABMEL3 ', 'UOK-57 ',
                 'ZR-75-1 ', 'A-204 ', 'CALU-1 ', 'MLI-045 ', 'MLI-076 ', 'NB4 ', 'CHA-59 ', 'LXFS 650L ', 'MLI-019 ',
                 'A-C/EBP 3 ', 'A-CREB 1 ', 'A-CREB 2 ', 'A-FOS 2 ', 'A-FOS 3 ', 'A-JUN 1 ', 'A-JUN 3 ', 'A431 ',
                 'SW 1783 ', 'ZR-75-30 ', 'Mar-Bel ', 'WI-38 ', 'CCD-19LU ', 'MDA-MB-435S ', 'TE85 ', 'VDSO/CMV-9 ',
                 'VDSO/P ', 'NYH/ICRF-187-1 ', 'VDSO/CMV-8 ', 'VDSO/E6-18 ', 'VDSO/E6-19 ', 'CHO ', 'CHO/159-1 ', 'NYH ']
    df_60_cell_lines = df[~df['CELL'].isin(to_remove)]
    x_col = df_60_cell_lines["CELL"].replace({'A549/ATCC ': "A549_ATCC", 'OVCAR-8 ': 'OVCAR-8', 'SW-620 ': 'SW-620',
                                              'U251 ': 'U251', 'SF-295 ': 'SF-295', 'NCI-H23 ': 'NCI-H23',
                                              'KM12 ': 'KM12', 'SN12C ': 'SN12C', 'HCT-116 ': 'HCT-116',
                                              'HCT-15 ': 'HCT-15', 'HT29 ': 'HT29', 'SF-268 ': 'SF-268',
                                              'COLO 205 ': 'COLO_205', 'SNB-19 ': 'SNB-19',
                                              'ACHN ': 'ACHN', 'UACC-257 ': 'UACC-257', 'IGROV1 ': 'IGROV1',
                                              'UO-31 ': 'UO-31', 'SK-MEL-28 ': 'SK-MEL-28', '786-0 ': '786-0',
                                              'M14 ': 'M14', 'MOLT-4 ': 'MOLT-4', 'SK-MEL-5 ': 'SK-MEL-5',
                                              'UACC-62 ': 'UACC-62', 'OVCAR-5 ': 'OVCAR-5', 'OVCAR-3 ': 'OVCAR-3',
                                              'K-562 ': 'K-562', 'NCI-H460 ': 'NCI-H460', 'HOP-62 ': 'HOP-62',
                                              'TK-10 ': 'TK-10', 'NCI-H322M ': 'NCI-H322M', 'OVCAR-4 ': 'OVCAR-4',
                                              'SK-OV-3 ': 'SK-OV-3', 'EKVX ': 'EKVX', 'CCRF-CEM ': 'CCRF-CEM',
                                              'LOX IMVI ': 'LOX_IMVI', 'CAKI-1 ': 'CAKI-1', 'SF-539 ': 'SF-539',
                                              'SK-MEL-2 ': 'SK-MEL-2', 'RPMI-8226 ': 'RPMI-8226',
                                              'NCI-H226 ': 'NCI-H226', 'SNB-75 ': 'SNB-75', 'MALME-3M ': 'MALME-3M',
                                              'NCI-H522 ': 'NCI-H522', 'HL-60(TB) ': 'HL-60(TB)', 'RXF 393 ': 'RXF_393',
                                              'HCC-2998 ': 'HCC-2998', 'HOP-92 ': 'HOP-92', 'A498 ': 'A498',
                                              'SR ': 'SR', 'NCI/ADR-RES ': 'NCI_ADR-RES', 'MCF7 ': 'MCF7',
                                              'MDA-MB-435 ': 'MDA-MB-435', 'PC-3 ': 'PC-3', 'DU-145 ': 'DU-145',
                                              'MDA-MB-231/ATCC ': 'MDA-MB-231_ATCC',  'HS 578T ': 'HS_578T',
                                              'T-47D ': 'T-47D', 'BT-549 ': 'BT-549', 'MDA-N ': 'MDA-N'})
    df_60_cell_lines.insert(0, 'Cell line', x_col)
    f, c = df_60_cell_lines.shape
    print('\n Final data: \n\t Cell lines: ', len(df_60_cell_lines.CELL.unique()),
          '\n\t Unique smiles used: ', len(df_60_cell_lines.SMILE.unique()),
          '\n\t Unique NSC used: ', len(df_60_cell_lines.NSC.unique()),
          '\n\t Number of NSC-CELL pairs: ', f)
    return df_60_cell_lines


def closest_distance_similarity(fps):
    simis = []
    nfps = len(fps)
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        max_sims = np.max(np.array(sims))
        simis.append(max_sims)
    return simis


def check_output_folder(output_directory):
    """
    :param output_directory: folder path
    :return: None
    """
    if not(os.path.isdir(output_directory)):
        os.makedirs(output_directory, exist_ok=True)
    return None


def get_outlier_molecules(df_aux):
    df_tanimoto = pd.DataFrame(data=df_aux['NSC'])
    df_fps, fingerprints1 = get_mfp(df_smiles=df_aux)
    tanimoto_sims_vect = closest_distance_similarity(fingerprints1)
    tanimoto_sims_vect.insert(0, 1)
    df_tanimoto['Tanimoto similarity'] = tanimoto_sims_vect
    df_tanimoto['outlier_05'] = np.where(df_tanimoto['Tanimoto similarity'] <= 0.5, True, False)
    df_new = df_aux.merge(df_tanimoto, on='NSC')
    df_nonoutliers = df_new[df_new.loc[:, 'outlier_05'] == False]
    print('\tNon-outlier molecules: ', len(df_nonoutliers))
    df_fps_non, fingerprints_non = get_mfp(df_smiles=df_nonoutliers)
    df_fps_non.insert(0, 'NSC', df_nonoutliers['NSC'].values)
    return df_fps_non, fingerprints_non


def get_sample(df, outliers):
    df_aux = df.drop_duplicates(subset=['SMILE'], ignore_index=True)
    col_name = ['MFP_' + str(i) for i in range(1024)]
    if outliers:
        print('\nRemoving outliers molecules...')
        df_nonoutliers, fingerprints_non = get_outlier_molecules(df_aux)
        nsc_col = df_nonoutliers['NSC']
        df_fps_non, fingerprints_non = get_mfp(df_smiles=df_nonoutliers)
        df_fps_non.insert(0, 'NSC', nsc_col)
        df_sample = df_fps_non.loc[:, col_name]
        sample = df_sample.to_numpy()
        df_nsc_smile = df_fps_non.loc[:, ['NSC', 'SMILE']]
    else:
        df_sample = df_aux.loc[:, col_name]
        sample = df_sample.to_numpy()
        df_nsc_smile = df_aux.loc[:, ['NSC', 'SMILE']]
    return df_nsc_smile, sample


def get_data_cluster(name_file_nci60, name_file_smi, dir_file, outliers):
    df = load_data(name_file=name_file_nci60, dir_file=dir_file)
    print('\nPreprocessing data...')
    df_mean = get_mean(df)

    # Merge both dataframes and drop duplicates
    df_with_mean = df.merge(df_mean, how='left', left_on=['ID'], right_on=['ID'])
    df_with_mean.drop_duplicates(['NSC', 'CELL'], inplace=True)
    df_with_mean.drop(columns=['NLOGGI50'], inplace=True)  # remove the previous pGI50 value

    # Read smiles file
    df_smi = load_smiles_file(name_file=name_file_smi, dir_file=dir_file)
    df_mean_smi = df_with_mean.merge(df_smi, how='left', left_on=['NSC'], right_on=['NSC'])
    # Remove instances where there is no SMILES information
    df_final = df_mean_smi[df_mean_smi.loc[:, 'SMILE'].notnull()]

    # MFP
    df_mfp_bits, fp1 = get_mfp(df_smiles=df_final)
    df_to_process_aux = df_final.merge(df_mfp_bits, how='left', left_on=['SMILE'], right_on=['SMILE'])
    df_to_process = df_to_process_aux[df_to_process_aux.loc[:, 'SMILE'].notnull()]
    df_to_process.dropna(inplace=True)

    # Drop the cell lines with few data
    df_60c = remove_few_data(df_to_process)

    # get sample
    df_info, sample = get_sample(df_60c, outliers)
    return df_info, sample






