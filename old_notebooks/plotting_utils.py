import matplotlib.pyplot as plt


def plot_pr_curve(precision, recall, color='black', legend=''):
    plt.step(recall, precision, alpha=1., color=color,
             where='post', label=legend)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.ylim([0.0, 1.05])
    plt.xlim([0.0, 1.0])
