import pandas as pd
from sklearn.metrics import precision_recall_curve, average_precision_score


def siderophore_preprocess(data: pd.DataFrame,
                           meta: pd.DataFrame,
                           verbose: bool = False) -> pd.DataFrame:
    """Filters dataframe and appends siderophore label.

    Keeps examples with either antibacterial, antifungal, or siderophore label OR
    a non-null activity label.
    """
    major_classes_idx = (meta.Antibacterial.isin([1]) |
                         meta.Antifungal.isin([1]) |
                         meta.Siderophore.isin([1]))
    activity_idx = ~meta.Activities.isnull()
    siderophore_idx = major_classes_idx | activity_idx
    sid_meta = meta[siderophore_idx]

    # Filter to valid examples
    data = data[data.BGC.isin(sid_meta.BGC)]

    # Append label
    data = data.join(sid_meta[['BGC', 'Siderophore']
                              ].set_index('BGC'), on='BGC')
    data = data.rename(columns={'Siderophore': 'label'})

    if verbose:
        print(f'There are {len(data)} examples for siderophore prediction.')
        print(f'Num positives: {(data.label == 1).sum()}')
        print(f'Num negatives: {(data.label == 0).sum()}')

    return data


def mw_preprocess(data: pd.DataFrame,
                  meta: pd.DataFrame,
                  cutoff: int = 1000,
                  verbose: bool = False) -> pd.DataFrame:
    """Filters dataframe and appends molecular weight binned label."""
    valid_struct_idx = ~meta.Structure.isnull()
    struct_meta = meta[valid_struct_idx]
    bins = [0] + [cutoff] + [10000]
    x = struct_meta.MW.copy()
    struct_meta = struct_meta.assign(label=pd.cut(x, bins, labels=False))

    # Filter to valid examples
    data = data[data.BGC.isin(struct_meta.BGC)]

    # Append label
    data = data.join(struct_meta[['BGC', 'label']].set_index('BGC'), on='BGC')

    if verbose:
        print(f'There are {len(data)} examples for binned weight prediction.')
        print(f'Num > {cutoff}: {(data.label > 0).sum()}')
        print(f'Num <= {cutoff}: {(data.label == 0).sum()}')

    return data


def fit_classifier(X: pd.DataFrame,
                            y: pd.Series,
                            clf):
    clf = clf.fit(X, y)
    preds = clf.predict_proba(X)
    aupr = average_precision_score(y, preds[:, 1])
    precision, recall, _ = precision_recall_curve(y, preds[:, 1])

    return clf, precision, recall, aupr
