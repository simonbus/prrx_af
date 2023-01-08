import pandas as pd
import numpy as np
import os
from roc_analysis import read_opt_cutoff
import helper


def conf_mat(y_true, y_pred):
    """GregoryMorse implementation of confusion matrix
    from https://github.com/scikit-learn/scikit-learn/issues/15388"""
    return np.array([[np.sum(~y_true & ~y_pred), np.sum(~y_true & y_pred)],  # TN, FP
                    [np.sum(y_true & ~y_pred), np.sum(y_true & y_pred)]])    # FN, TP


def bootstrap(prrx_dir, cutoff_dir, db, x_sec, cutoff_method, num_iter, boot_dir):
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
        for iter in range(num_iter):
            print(f'{db.upper()} ({group}), {iter = } / {num_iter}')
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
        writer = pd.ExcelWriter(
                os.path.join(boot_dir, f"bootstrap_{group.replace('%', '_perc')}_{x_sec}s_N={num_iter}_cutoff_{cutoff_method}.xlsx"),
                engine='openpyxl'
            )
        for i, (metric, metric_name) in enumerate([(tn, "TN"), (fp, "FP"),
                                                   (fn, "FN"), (tp, "TP")]):
            pd.DataFrame.from_dict(metric).to_excel(
                writer, sheet_name=metric_name, index=False)
        writer.save()


if __name__ == '__main__':
    num_iter = 5000
    prrx_dir = '../data/processed'
    cutoff_dir = '../reports/excel/roc'
    x_sec = 60
    cutoff_method = 'youden'
    boot_dir = '../reports/excel/boot'
    for db in ['ltafdb']:  # 'afdb
        bootstrap(prrx_dir, cutoff_dir, db, x_sec, cutoff_method, num_iter,
                  boot_dir)
