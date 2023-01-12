import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt
from roc_analysis import read_opt_cutoff
import helper


def conf_mat(y_true, y_pred):
    """GregoryMorse implementation of confusion matrix
    from https://github.com/scikit-learn/scikit-learn/issues/15388"""
    return np.array([[np.sum(~y_true & ~y_pred), np.sum(~y_true & y_pred)],  # TN, FP
                    [np.sum(y_true & ~y_pred), np.sum(y_true & y_pred)]])    # FN, TP


def bootstrap(prrx_dir, cutoff_dir, db, x_sec, cutoff_method, N, boot_dir):
    df_prrx = helper.read_prrx(prrx_dir, db, x_sec)
    dfs_cutoff = read_opt_cutoff(cutoff_dir, db, x_sec, cutoff_method)
    params = helper.get_params(df_prrx)
    for group, features in params.items():  # pRRx, pRRx%
        if not os.path.exists(boot_dir):
            os.makedirs(boot_dir)
        features = params[group]
        tn = {f: [] for f in features}
        fp = {f: [] for f in features}
        fn = {f: [] for f in features}
        tp = {f: [] for f in features}
        cutoffs = dfs_cutoff[group][f'{x_sec} s'].values
        n_sampl = len(df_prrx.index)
        for iter in range(N):
            print(f'{db.upper()} ({group}), {iter = } / {N}')
            x_df = df_prrx[features].sample(n=n_sampl, replace=True)
            X = x_df.values
            y_true = df_prrx['is_af'].loc[x_df.index].values
            y_pred = X >= cutoffs
            
            for i, feature in enumerate(features):
                tn_, fp_, fn_, tp_ = conf_mat(y_true, y_pred[:, i]).ravel()
                tn[feature].append(tn_)
                fp[feature].append(fp_)
                fn[feature].append(fn_)
                tp[feature].append(tp_)
        fname = f"bootstrap_{group}_{x_sec}s_N={N}_cutoff_{cutoff_method}.xlsx"
        writer = pd.ExcelWriter(
            os.path.join(boot_dir, db, fname),
            engine='openpyxl'
            )
        for i, (metric, metric_name) in enumerate([(tn, "TN"), (fp, "FP"),
                                                   (fn, "FN"), (tp, "TP")]):
            pd.DataFrame.from_dict(metric).to_excel(
                writer, sheet_name=metric_name, index=False)
        writer.save()


def calculate_95ci(boot_dir, db, x_sec, N, group):
    # Calculate 95% CI from bootstrap results
    fname = f"bootstrap_{group}_{x_sec}s_N={N}.xlsx"
    xls_path = os.path.join(boot_dir, db, fname)
    sheets = pd.read_excel(xls_path, engine='openpyxl', sheet_name=None)
    features = sheets['TP'].columns.values
    tp = sheets['TP'].values
    fp = sheets['FP'].values
    tn = sheets['TN'].values
    fn = sheets['FN'].values
    scores = {
        "Accuracy": 100 * (tp+tn)/(tp+tn+fp+fn+1e-9),
        "Sensitivity": 100 * (tp)/(tp+fn+1e-9),
        "Specificity": 100 * (tn)/(tn+fp+1e-9),
        "PPV": 100 * (tp)/(tp+fp+1e-9),
        "NPV": 100 * (tn)/(tn+fn+1e-9),
        "F1-score": 100 * 2 * tp / (2*tp + fp + fn),
        "DOR": tp*tn/(fp*fn+1e-9)
    }
    medians = {}
    p2_5 = {}
    p97_5 = {}
    for metric_name in scores.keys():
        metric_arr = scores[metric_name]
        medians[metric_name] = np.median(metric_arr, axis=0)
        p2_5[metric_name] = np.percentile(metric_arr, 2.5, axis=0)
        p97_5[metric_name] = np.percentile(metric_arr, 97.5, axis=0)
    # Save percentiles (2.5, 50, 97.5) in Excel
    perc_dir = os.path.join(boot_dir, db, 'percentiles')
    if not os.path.exists(perc_dir):
        os.makedirs(perc_dir)
    writers = {percentile: pd.ExcelWriter(
        os.path.join(perc_dir, f"{percentile}_bootstrap_{group}_{x_sec}s_N={N}.xlsx"),
        engine='openpyxl'
    ) for percentile in ['medians', 'p2_5', 'p97_5']}
    for metric in medians:
        for percentile in ['medians', 'p2_5', 'p97_5']:
            df = pd.DataFrame.from_dict(
                {"Length [s]": [x_sec]})
            df[features] = np.array(eval(percentile)[metric])
            df.to_excel(
                writers[percentile],
                sheet_name=metric,
                index=False
            )
    for percentile in writers.keys():
        writers[percentile].save()    


def plot_boot(db, group, N, x_sec, boot_dir, fig_dir):
    color_tab = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple']
    if '%' in group:
        thr_unit = '%'
    else:
        thr_unit = 'ms'
    # Read percentiles (2.5, 50, 97.5) from Excel
    dfs_p2_5 = pd.read_excel(os.path.join(
        boot_dir, db, 'percentiles',
        f'p2_5_bootstrap_{group}_{x_sec}s_N={N}.xlsx'), sheet_name=None)
    dfs_p50 = pd.read_excel(os.path.join(
        boot_dir, db, 'percentiles',
        f'medians_bootstrap_{group}_{x_sec}s_N={N}.xlsx'), sheet_name=None)
    dfs_p97_5 = pd.read_excel(os.path.join(
        boot_dir, db, 'percentiles',
        f'p97_5_bootstrap_{group}_{x_sec}s_N={N}.xlsx'), sheet_name=None)
    # Plot medians and 95% CI
    ax = helper.default_fig()
    ylim_min = 100
    for i, metric in enumerate(['Accuracy', 'Sensitivity', 'Specificity',
                                'PPV', 'NPV', ]):
        features = [col for col in dfs_p50[metric].columns if 'Length' not in col]
        x_thr = helper.get_x_thr(features)
        [p2_5] = dfs_p2_5[metric][features].values
        [p50] = dfs_p50[metric][features].values
        [p97_5] = dfs_p97_5[metric][features].values
        ylim_step = 5
        if np.floor(np.min(p2_5)/ylim_step) * ylim_step < ylim_min:
            ylim_min = np.floor(np.min(p2_5)/ylim_step) * ylim_step
        ylim_max = 100
        helper.plot_line(
            ax, xval=x_thr, yval=p50, label=metric,
            xlim=(0, np.ceil(x_thr[-1]/5)*5), ylim=(ylim_min, ylim_max),
            #xlim=(0, 25), ylim=(80, 100),
            xlabel=f'Threshold x [{thr_unit}]',
            ylabel='Score [%]', title='Classification metrics',
            color=color_tab[i], marker='.', markersize=1.5,
            linewidth=1
            )
        ax.fill_between(x_thr, p2_5, p97_5, alpha=0.3, color=color_tab[i])
    ax.legend()
    fname = f'metrics_boot_N={N}_{group.replace("%", "_perc")}_{x_sec}s_{db}.png'
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    helper.save_fig(fig_dir, fname)
    plt.show()
    # DOR
    metric = 'DOR'
    features = [col for col in dfs_p50[metric].columns if 'Length' not in col]
    x_thr = helper.get_x_thr(features)
    [p2_5] = dfs_p2_5[metric][features].values
    [p50] = dfs_p50[metric][features].values
    [p97_5] = dfs_p97_5[metric][features].values
    ylim_step = 50
    ylim_min = np.floor(np.min(p2_5)/ylim_step) * ylim_step
    ylim_max = np.ceil(np.max(p97_5)/ylim_step) * ylim_step
    ax = helper.default_fig()
    helper.plot_line(
        ax, xval=x_thr, yval=p50, label=metric,
        xlim=(0, np.ceil(x_thr[-1]/5)*5), ylim=(ylim_min, ylim_max),
        xlabel=f'Threshold x [{thr_unit}]',
        ylabel='Metric [%]', title='Diagnostic odds ratio',
        color=color_tab[0], marker='.', markersize=1,
        linewidth=0.5
        )
    ax.fill_between(x_thr, p2_5, p97_5, alpha=0.3, color=color_tab[0])
    fname = f'DOR_boot_N={N}_{group.replace("%", "_perc")}_{x_sec}s_{db}.png'
    helper.save_fig(fig_dir, fname)
    plt.show()


def get_scores_from_confusion_matrix_sheets(sheets, param):
    tp = sheets['TP'][param].values
    fp = sheets['FP'][param].values
    tn = sheets['TN'][param].values
    fn = sheets['FN'][param].values
    scores = {
        "Accuracy": 100 * (tp+tn)/(tp+tn+fp+fn+1e-9),
        "Sensitivity": 100 * (tp)/(tp+fn+1e-9),
        "Specificity": 100 * (tn)/(tn+fp+1e-9),
        "PPV": 100 * (tp)/(tp+fp+1e-9),
        "NPV": 100 * (tn)/(tn+fn+1e-9),
        "F1-score": 100 * 2 * tp / (2*tp + fp + fn),
        "DOR": tp*tn/(fp*fn+1e-9)
    }
    return scores


def plot_compare_scores_distr(db, group_param_label,
                              N, x_sec, boot_dir, fig_dir, suptitle=None):
    scores = {}
    for (group, param, label) in group_param_label:
        fname = os.path.join(
            boot_dir, db, f"bootstrap_{group}_{x_sec}s_N={N}.xlsx")
        sheets = pd.read_excel(fname, engine='openpyxl', sheet_name=None)
        scores[param] = get_scores_from_confusion_matrix_sheets(sheets, param)
    color_tab = ['tab:orange', 'tab:blue', 'tab:green', 'tab:red', 'tab:purple']
    fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharey=True)
    for i, metric in enumerate(['Accuracy', 'Sensitivity', 'Specificity', 'DOR']):
        ax = axs[i//2, i % 2]
        for j, (group, param, label) in enumerate(group_param_label):
            score = scores[param][metric]
            ax.hist(score, bins=50, label=label, alpha=0.7, color=color_tab[j])
            ax.set_title(metric)
            ax.minorticks_on()
            ax.grid(b=True, which='major', linestyle='-')
            ax.grid(b=True, which='minor', linestyle='--', linewidth=0.5)
            ax.legend()
            if metric != 'DOR':
                ax.set_xlabel('[%]')
    plt.suptitle(suptitle, size=16)
    param_str = '_vs_'.join([l for g, p, l in group_param_label])
    param_str = param_str.replace('.', '_').replace('%', 'perc')
    fname = f'metrics_{param_str}_N={N}_{x_sec}s_{db}.png'
    plt.tight_layout()
    helper.save_fig(fig_dir, fname)
    plt.show()


if __name__ == '__main__':
    N = 5000
    prrx_dir = '../data/processed'
    cutoff_dir = '../reports/excel/roc'
    x_sec = 60
    cutoff_method = 'youden'
    boot_dir = '../reports/excel/boot'
    fig_dir = '../reports/images/boot'
    for db in ['ltafdb']:  # 'afdb'
        bootstrap(prrx_dir, cutoff_dir, db, x_sec, cutoff_method, N,
                  boot_dir)
        for group in ['pRRx', 'pRRx%']:
            calculate_95ci(boot_dir, db, x_sec, N, group)
            plot_boot(db, group, N, x_sec, boot_dir, fig_dir)
        # 5. Compare distributions of metrics for pRR5%, pRR31 and pRR50
        plot_compare_scores_distr(
            db, group_param_label=(('pRRx', 'pRR31.25', 'pRR31'),
                                   ('pRRx%', 'pRR3.25%', 'pRR3.25%')),
            N=N, x_sec=x_sec, boot_dir=boot_dir,
            fig_dir=fig_dir)
        plot_compare_scores_distr(
            db, group_param_label=(('pRRx', 'pRR54.6875', 'pRR50'),
                                   ('pRRx', 'pRR31.25', 'pRR31')),
            N=N, x_sec=x_sec, boot_dir=boot_dir,
            fig_dir=fig_dir)
