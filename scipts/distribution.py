import pandas as pd
import numpy as np
import json
import os
import helper
from matplotlib import pyplot as plt


def calc_percentiles(df):
    """Calculate percentiles of pRRx/pRRx% parameters in AF and SR, return dict.

    Args:
        df (pd.DataFrame): DF with values of pRRx/pRRx% parameters

    Returns:
        (dict): Dictionary with values of percentiles in AF and SR
    """
    params = helper.get_params(df)
    # Parameter types (perc, ms)
    groups = [group for group in params if len(params[group]) > 0]
    perc_dict = {g: dict() for g in groups}
    for group in groups:
        features = [x[0] for x in params[group]]
        perc_dict[group] = {f: {'AF': {}, 'SR': {}} for f in features}
        # Calculate percentiles (0-100)
        for i, feature in enumerate(features):
            for is_af, rhythm, in ((0, 'SR'), (1, 'AF')):
                val = df[df['is_af'] == is_af][feature].values
                perc = list(np.arange(101)/100)
                q = list(np.quantile(val, perc))
                perc_dict[group][feature][rhythm] = {
                    'centiles': perc,
                    'values': q}
    return perc_dict


def plot_distr(perc_dict, db, x_sec, fig_dir):
    """Plot medians, and 10-90 and 25-75 percentile ranges
    of pRRx/pRRx% parameters in AF and SR

    Args:
        perc_dict (dict): Dictionary with parameter names and percentiles in AF and SR
        db (str): Acronym of the database
        x_sec (int): Length of RR sequence [s]
        fig_dir: Write directory for images
    """
    for group in perc_dict.keys():
        ax = helper.default_fig()
        title = f'Quartiles, 10-th and 90-th percentile of {group}'
        if group == 'pRRx':
            xlabel = 'Threshold x [ms]'
        else:
            xlabel = 'Threshold x [%]'
        features = list(perc_dict[group].keys())
        x_thr = helper.get_x_thr(features)
        for rhythm, color in (('SR', 'tab:blue'), ('AF', 'tab:orange')):
            # Values of percentiles
            val = np.array([perc_dict[group][f][rhythm]['values'] for f in features])
            # Median line
            p50 = val[:, 50]
            helper.plot_line(
                ax, xval=x_thr, yval=p50, label=rhythm, xlim=(0, x_thr[-1]),
                ylim=(0, 100), xlabel=xlabel, ylabel='[%]',
                title=title, color=color, marker='o', markersize=2,
                linewidth=2)
            # Bands: 25-75 and 10-90 percentile ranges
            for perc, lw, alpha in [(25, 1, 0.6), (10, 1, 0.3)]:
                p_low = val[:, perc]
                p_high = val[:, 100-perc]
                ax.plot(x_thr, p_low, '-o', markersize=1, color=color,
                        alpha=alpha, linewidth=lw)
                ax.plot(x_thr, p_high, '-o', markersize=1, color=color,
                        alpha=alpha, linewidth=lw)
                ax.fill_between(x_thr, p_low, p_high, alpha=alpha, color=color)
        ax.legend()
        fname = f"distr_{group}_{db}_vs_threshold_{x_sec}s.png"
        helper.save_fig(fig_dir, fname)
        plt.show()


if __name__ == '__main__':
    prrx_dir = '../data/processed'
    fig_dir = '../reports/images/distr'
    x_sec = 60
    for db in ['ltafdb', 'afdb']:
        df = pd.read_csv(os.path.join(prrx_dir, db, f"prrx_{db}_{x_sec}s.csv"))
        perc_dict = calc_percentiles(df)
        if not os.path.exists(fig_dir):
            os.makedirs(os.path.join(fig_dir))
        plot_distr(perc_dict, db, x_sec, fig_dir)
    