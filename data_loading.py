from typing import Tuple 

import pandas as pd

def load_data_and_counts(path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads any of the count tables and returns the full dataframe,
    as well as a dataframe just containing counts.
    """
    
    df = pd.read_csv(path, sep='\t')
    df = df.dropna(axis=1)
    counts = df.select_dtypes(exclude='object') # Drop all columns with strings

    return df, counts

def load_mibig_metadata(path: str='data/metadata/mibig_metadata.tsv') -> pd.DataFrame:
    """
    Loads mibig metadata and incorporates extra annotations.
    """
    m = pd.read_csv(path, sep='\t')

    m.loc[m['Type'] == 'lantipeptide', 'Type'] = 'Lanthipeptide'
    m.loc[m['Type'] == 'Lantipeptide', 'Type'] = 'Lanthipeptide'
    m.Name = [x.lower() for x in m.Name]

    with open('./data/labels/extra_antibiotics') as f:
        extra_antibiotics = [x.strip() for x in f.readlines()]
    with open('./data/labels/extra_antifungal') as f:
        extra_antifungals = [x.strip() for x in f.readlines()]
    with open('./data/labels/extra_siderophore') as f:
        extra_siderophores = [x.strip() for x in f.readlines()]

    m.loc[m['Name'].isin(extra_antibiotics), 'Antibacterial'] = 1
    m.loc[m['Name'].isin(extra_antifungals), 'Antifungal'] = 1
    m.loc[m['Name'].isin(extra_siderophores), 'Siderophore'] = 1

    return m
