# nci60_clustering

## Getting Started
Python implementation of the main clustering results for the paper entitled:
_On the best way to cluster NCI-60 molecules_


## Prerequisites
- Python 3.9.7
- [venv](https://docs.python.org/3/tutorial/venv.html) 

## Installation

1. Create the virtual environment
```bash
python -m venv .venv
```

2. Activate the virtual environment
   1. On Windows, run:
   ```bash
   .venv\Scripts\activate.bat
   ```
   2. On Linux or MacOs, run:
   ```bash
   .source .venv/bin/activate
   ```
3. Install all the necessary packages 
```bash
python -m pip install -r requirements.txt
```

## Usage

1. Download and store the from the [NCI-60 Growth Inhibition Data](https://wiki.nci.nih.gov/display/NCIDTPdata/NCI-60+Data+Download+-+Previous+Releases). 
   1. The required files contain the endpoints calculated from concentration curves ("CANCER60GI50_Oct2020.LST", for instance) and SMILES ("Chem2D_Jun2016.smi", for instance). 
      Other releases of both files are also available for download.
2. In the file _main.py_, replace the following lines.
   1. Replace with the directory  where the downloaded data is stored
   ```python
   dir_working = '/home/hernadez/Documents/NCI60_data/' 
   ```
   2. Replace only if the files are different from those suggested in 1.
   ```python
   file_nci60 = 'CANCER60GI50_Oct2020.LST'  
   file_smiles = 'Chem2D_Jun2016.smi' 
   ```
   
3. By default, the number of clusters (k) has been set to 7, and the removal of outliers has been requested. Both values can be modified in the following lines.
```python
k = 7  # number of clusters
outliers = True  # to remove outliers
```


To run:
```python
python main.py
```

Output:
As a result, a folder containing the clustering assignment ([NSC, SMILES, Cluster ID]) and the corresponding clustering quality metrics will be created.

### License


### Contact

