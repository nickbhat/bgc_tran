from typing import Tuple, Dict

import numpy as np
import pandas as pd


def load_data_and_counts(path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Loads any of the count tables and returns the full dataframe,
    as well as a dataframe just containing counts.
    """

    df = pd.read_csv(path, sep='\t')
    df = df.dropna(axis=1)
    # Drop all columns with strings
    counts = df.select_dtypes(exclude='object')

    return df, counts


def load_mibig_metadata(path: str = 'data/metadata/mibig2_metadata.tsv') -> pd.DataFrame:
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


def load_base_data():
    mibig_metadata = load_mibig_metadata()
    cath_metadata = pd.read_csv(
        './data/metadata/final_fun_families.tsv', sep='\t')

    mibig_pfam = pd.read_csv(
        './data/pfam/mibig_pfam20.tsv', sep='\t').dropna(how='all', axis='columns')
    mibig_pfam['ABC2_membrane_3'] = mibig_pfam['ABC2_membrane_3'] + \
        mibig_pfam['ABC2_membrane_3.1']
    del mibig_pfam['ABC2_membrane_3.1']

    mibig_cath = pd.read_csv(
        './data/cathdb/mibig_cathdb.tsv', sep='\t').dropna(how='all', axis='columns')
    mibig_cath = mibig_cath[['BGC'] + list(cath_metadata['FunFam'])]
    return mibig_metadata, mibig_pfam, mibig_cath


def get_subset(fam, idx, metadata):
    # Grabs subset of bgcs matching some metadata criteria
    return fam[fam.BGC.isin(metadata[idx].BGC)]


def load_siderophore(metadata, pfam, cath: pd.DataFrame) -> Dict[str, Dict[str, pd.DataFrame]]:
    """Loads siderophore data split by gram and protein family"""
    siderophore = {'pfam': {}, 'cath': {}}
    idx = metadata.Gram == 1
    siderophore['pfam']['gram+'] = get_subset(pfam, idx, metadata).drop(
        ['BGC', 'Order'], axis=1), metadata.Siderophore[idx]
    siderophore['cath']['gram+'] = get_subset(cath, idx, metadata).drop(
        ['BGC'], axis=1), metadata.Siderophore[idx]
    idx = metadata.Gram == 0
    siderophore['pfam']['gram-'] = get_subset(pfam, idx, metadata).drop(
        ['BGC', 'Order'], axis=1), metadata.Siderophore[idx]
    siderophore['cath']['gram-'] = get_subset(cath, idx, metadata).drop(
        ['BGC'], axis=1), metadata.Siderophore[idx]

    return siderophore


def load_weight(metadata, pfam, cath, bins=[500., 1500.]) -> Dict[str, Dict[str, pd.DataFrame]]:
    """Loads molecular weight data split by gram and protein family."""
    meta_with_struct = metadata[pd.notnull(metadata.Structure)]
    bins = [0] + bins + [10000.]
    weight = {'pfam': {}, 'cath': {}}

    idx = metadata.Gram == 1
    weight['pfam']['gram+'] = get_subset(pfam, idx, meta_with_struct).drop(
        ['BGC', 'Order'], axis=1), pd.cut(metadata.MW[idx], bins, labels=False)
    weight['cath']['gram+'] = get_subset(cath, idx, meta_with_struct).drop(
        ['BGC'], axis=1), pd.cut(metadata.MW[idx], bins, labels=False)
    idx = metadata.Gram == 0
    weight['pfam']['gram-'] = get_subset(pfam, idx, meta_with_struct).drop(
        ['BGC', 'Order'], axis=1), pd.cut(metadata.MW[idx], bins, labels=False)
    weight['cath']['gram-'] = get_subset(cath, idx, meta_with_struct).drop(
        ['BGC'], axis=1), pd.cut(metadata.MW[idx], bins, labels=False)
    
    return weight
