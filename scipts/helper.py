import re
import os
from matplotlib import pyplot as plt
import numpy as np


def get_params(df, col_name=None):
    """Get grouped parameter names from pRRx/pRRx% dataframe

    Args:
        df (pd.DataFrame): DF with pRRx/pRRx% parameters as column names or a single column
        col_name (str, optional): column with pRRx/pRRx% names. Defaults to None.

    Returns:
        (list): List of parameters
    """
    regex = {
        'pRRx':    r'(?<=pRR)\d+.\d+(?![%\d])',
        'pRRx%':   r'(?<=pRR)\d+\.*\d*(?=%)',
    }
    params = {}
    if col_name is None:
        lst = df.columns
    else:
        lst = df[col_name]
    for item in lst:
        for key in regex:
            x = re.search(regex[key], item)
            if x is not None:
                x = x.group()
                if key in params:
                    params[key].append(item)
                else:
                    params[key] = [item]
    return params


def default_fig():
    """Make a default figure

    Returns:
        plt.axis: Axis of a single-axis figure
    """
    fig, ax = plt.subplots(figsize=(8, 5))
    return ax


def get_x_thr(lst):
    """Get thresholds from a list of pRRx/pRRx% parameters

    Args:
        lst (list): List of pRRx/pRRx% parameters

    Returns:
        (np.array): Thresholds (ms or %)
    """
    return np.array([float(re.search(r'\d+(\.\d+)*', f).group()) for f in lst])


def plot_line(ax, xval, yval, label=None, xlim=None, ylim=None,
              xlabel=None, ylabel=None, title=None, color=None,
              marker=None, markersize=1, linewidth=1):
    ax.plot(xval, yval, color=color, marker=marker, label=label,
            markersize=markersize, linewidth=linewidth)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.minorticks_on()
    ax.grid(b=True, which='major', linestyle='-')
    ax.grid(b=True, which='minor', linestyle='--', linewidth=0.5)


def save_fig(fig_dir, fname):
    """Save figure to a file with specified name and directory

    Args:
        fig_path (str): Write directory
        fname (str): Filename
    """
    plt.savefig(os.path.join(fig_dir, fname),
                bbox_inches="tight", dpi=300)