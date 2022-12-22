import os.path
import pandas as pd


def load_data(name_file):
    df_data = pd.read_csv()
    x_sample = df_data.to_numpy()
    return x_sample, df_data


def check_output_folder(output_directory):
    """

    :param output_directory: folder path
    :return: None
    """
    if not(os.path.isdir(output_directory)):
        os.makedirs(output_directory, exist_ok=True)
    return None

