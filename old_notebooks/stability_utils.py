from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import precision_recall_curve

from plotting_utils import plot_pr_curve


def fit_cross_val_and_plot_pr(X, y, clf_dict, num_splits=5, num_repeats=5, title='Cross validation PR curve for all repeats'):
    """Does repeated cross validation, plots pr curves for each cv and each classifier.

    Note that clf_dict is a Dict[str, sklearn classifier]
    """
    kf = StratifiedKFold(n_splits=num_splits)
    predictions = {k: [] for k in clf_dict}
    labels = {k: [] for k in clf_dict}

    cm = plt.get_cmap('tab20')
    colors = [cm(1. * i/len(clf_dict)) for i in range(len(clf_dict))]

    plt.figure(figsize=(8, 8))
    for i in range(num_repeats):
        for train, test in kf.split(X, y):
            for k, clf in clf_dict.items():
                X_train, y_train = X.iloc[train], y.iloc[train]
                X_test, y_test = X.iloc[test], y.iloc[test]

                clf.fit(X_train, y_train)
                preds = clf.predict_proba(X_test)
                predictions[k].extend(preds[:, 1])
                labels[k].extend(y_test)

        for j, k in enumerate(clf_dict):
            preds = predictions[k]
            trues = labels[k]
            precision, recall, _ = precision_recall_curve(trues, preds)
            # Hilarious hack to prevent legend duplication
            if i > 0:
                k = None
            plot_pr_curve(precision, recall, color=colors[j], legend=k)

    plt.title(title)
    plt.legend()
    plt.show()


def create_rule(feature, threshold, feature_names):
    return f'{feature_names[feature]} <= {threshold}'


def visualize_tree_rules(X, y, clf, feature_names, num_splits=5, num_repeats=5, title='next'):
    """Does repeated cross validation, keeps list of rules, then plots rules with counts."""
    kf = StratifiedKFold(n_splits=num_splits)
    rules = []
    num_split_examples = []
    depth = clf.tree_.max_depth

    for i in range(num_repeats):
        for train, test in kf.split(X, y):
            X_train, y_train = X.iloc[train], y.iloc[train]

            clf.fit(X_train, y_train)
            features = clf.tree_.feature[:depth]
            thresholds = clf.tree_.threshold[:depth]

            rules.append([create_rule(f, t, feature_names)
                          for f, t in zip(features, thresholds)])
            # This requires some fickle ordering of tree_ properties to work. Likely place for bug if
            # we deviate from original case (e.g. increasing tree depth further).
            num_split_examples.append(
                clf.tree_.weighted_n_node_samples[::-1][:depth])

    def flatten(l): return [item for sublist in l for item in sublist]

    reordered_rules = [[r[i] for r in rules] for i in range(depth)]
    reordered_nums = [[n[i] for n in num_split_examples] for i in range(depth)]

    mapped_nums = {r: n for r, n in list(
        zip(flatten(reordered_rules), flatten(reordered_nums)))}
    mean_nums = {r: np.mean(n) for r, n in mapped_nums.items()}

    print(title)
    for level in range(depth):
        c = Counter(reordered_rules[level])
        labels, counts = zip(*c.items())

        idx = np.arange(len(labels))
        fig, ax = plt.subplots()
        ax_count = ax.twinx()
        ax.barh(idx, counts, 0.8)
        ax.set_title(f'Rule occurrence in level {level+1} of tree.')
        ax.set_yticks(idx)
        ax.set_yticklabels(labels)
        ax.set_ylabel('Rule')
        ax.set_xlabel('Number of times rule appeared.')
        ax_count.set_yticks(idx*(1/len(labels)) + 0.1)
        ax_count.set_yticklabels([mean_nums[l] for l in labels])
        ax_count.set_ylabel('Mean number of splits per rule')
        plt.show()
    print('='*100)


def visualize_lasso_coeffs(X, y, clf, feature_names, num_repeats=5, num_splits=5, title='Effect sizes merged across all repeats and all folds'):
    kf = StratifiedKFold(n_splits=num_splits)
    dfs = []

    coefs = []
    for i in range(num_repeats):
        for j, (train, test) in enumerate(kf.split(X, y)):
            X_train, y_train = X.iloc[train], y.iloc[train]
            clf.fit(X_train, y_train)
            coefs.append(np.squeeze(clf.coef_))

    coefs = np.array(coefs)
    results = pd.DataFrame(coefs)
    results = results.rename(
        columns={i: f for i, f in enumerate(feature_names)})
    boring_idx = (results.mean() < 0.01) & (results.mean() > -0.01) & (results.std() < 0.1)
    results = results.loc[:, ~boring_idx]
    results = results.melt(var_name='feature', value_name='vals')
    dfs.append(results)

    plt.figure(figsize=(12, 22))
    sns.violinplot(x='vals', y='feature', data=pd.concat(dfs))
    plt.ylabel('Protein Families')
    plt.xlabel('Effect Size')
    plt.title(title)
    plt.show()
