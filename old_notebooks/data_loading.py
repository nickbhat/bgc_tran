from typing import Dict

import pandas as pd
from pathlib import Path


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

    m = m[m.Orf_length <= 75]  # Get rid of long BGC due to inflation

    return m


def _load_pfam(path: Path):
    d = pd.read_csv(path, sep='\t').dropna(how='all', axis='columns')
    # Data bug, should consider fixing in file.
    d['ABC2_membrane_3'] += d['ABC2_membrane_3.1']
    del d['ABC2_membrane_3.1']
    d = d.drop(['Order'], axis=1)

    return d


def _load_cath(path: Path, meta: pd.DataFrame):
    d = pd.read_csv(path, sep='\t').dropna(how='all', axis='columns')
    d = d[['BGC'] + list(meta['FunFam'])]

    return d


def _load_sbp(path: Path):
    d = pd.read_csv(path, sep='\t').dropna(how='all', axis='columns')
    d = d.drop(['Order'], axis=1)

    return d


def _load_non_transporter_pfam(path: Path):
    d = pd.read_csv(path, sep='\t').dropna(how='all', axis='columns')
    d = d.drop(['Unnamed: 0'], axis=1)

    return d


def load_data(meta: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """Central function for loading all datasets.

    Takes the following steps in loading data.
        1. Loads all datasets with necessary column modifications but no row modifications.
        2. Creates an index of examples with both pfam and cath transporters.
        3. Creates an index of examples with Gram status 0 or 1.
        4. Keep only examples with both transporters that are also gram +/-.
        5. For each dataset d:
            a. drop examples that are not in 4.
            b. append Gram status to d.
            c. Create a row of all 0's for any examples in 4, but not in d. 
    """
    # Step 1.
    data = {}
    path = Path('./data/pfam/mibig_pfam20.tsv')
    if not path.exists():
        raise FileNotFoundError(path)
    data['pfam'] = _load_pfam(path)

    cath_meta = pd.read_csv('./data/metadata/final_fun_families.tsv', sep='\t')
    path = Path('./data/cathdb/mibig_cathdb.tsv')
    if not path.exists():
        raise FileNotFoundError(path)
    data['cath'] = _load_cath(path, cath_meta)

    path = Path('./data/pfam/mibig2_sbp_hits.tsv')
    if not path.exists():
        raise FileNotFoundError(path)
    data['sbp'] = _load_sbp(path)

    path = Path('./data/metadata/biosynthetic_mibig2.tsv')
    if not path.exists():
        raise FileNotFoundError(path)
    data['biosynthetic'] = _load_non_transporter_pfam(path)

    # Step 2.
    shared_bgcs = set(data['cath'].BGC).intersection(set(data['pfam'].BGC))

    # Step 3.
    # This index implicitly contains other filtering, such as orf length filtering.
    # Any filtering done on metadata will apply to datasets here.
    # We could consider making this more explicit.
    gram_bgcs = set(meta[meta.Gram <= 1].BGC)

    # Step 4.
    filter_bgcs = shared_bgcs.intersection(gram_bgcs)
    print(f'Using {len(filter_bgcs)} BGCs with both pfam and cath transporters, as well as Gram +/-.')

    # Step 5.
    for key, d in data.items():
        # Step 5.a.
        tmp = d[d.BGC.isin(filter_bgcs)]

        # Step 5.b.
        all_zero_bgcs = filter_bgcs.difference(set(tmp.BGC))
        num_zero_examples = len(all_zero_bgcs)
        print(f'Dataset {key} has {num_zero_examples} examples with all zeros.')

        # Add examples if some are missing
        if num_zero_examples > 0:
            # Create zero df
            features = tmp.columns.drop('BGC')
            zeros = pd.DataFrame(0, index=range(num_zero_examples), columns=features)
            # Add BGC names and merge 
            zeros['BGC'] = all_zero_bgcs
            tmp = pd.concat([tmp, zeros], ignore_index=True, sort=False)
        
        tmp = tmp.sort_values('BGC').reset_index(drop=True)

        # Step 5.c.
        tmp = tmp.join(meta[['BGC', 'Gram']].set_index('BGC'), on='BGC')
        data[key] = tmp

    return data
