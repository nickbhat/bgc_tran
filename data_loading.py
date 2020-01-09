import pandas as pd


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
