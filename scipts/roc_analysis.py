import os
import pandas as pd
import numpy as np
from sklearn import metrics
import helper


def calc_auc(df, x_sec):
    """Calculate area under ROC curve (AUC)fo pRRx/pRRx% parameters

    Args:
        df (pd.DataFrame): DF with pRRx/pRRx% parameters
        x_sec (int): RR series length [s]

    Returns:
        dict: Dictionary of DFs with AUC of pRRx/pRRx% parameters
    """
    params = helper.get_params(df)
    dfs_auc = {}
    # Parameter types (pRRx, pRRx%)
    groups = [group for group in params if len(params[group]) > 0]
    for group in groups:
        features = params[group]
        df_auc = pd.DataFrame()

        df_auc['feature'] = features
        aucs = []
        y = df['is_af'].values

        # Calculate AUC
        for i, feature in enumerate(features):
            scores = df[feature].values
            auc = metrics.roc_auc_score(y, scores)
            if auc < 0.5:
                auc = 1 - auc
            aucs.append(auc)
        aucs = np.array(aucs)
        df_auc[f"{x_sec} s"] = aucs
        dfs_auc[group] = df_auc
    return dfs_auc


def auc_prrx_to_excel(prrx_dir, db, x_sec, auc_dir):
    """Calculate AUC of pRRx/pRRx% and save in a single Excel file

    Args:
        prrx_dir (str): Directory with pRRx/pRRx% in CSV format
        db (str): Acronym of the database
        x_sec (int): Length of RR sequence [s]
        auc_dir (str): Write directory for AUC
    """
    # 1. Read data (pRRx_ms or pRRx_%)
    df = pd.read_csv(os.path.join(prrx_dir, db, f"prrx_{db}_{x_sec}s.csv"))
    # 2. Calculate AUC
    dfs_auc = calc_auc(df, x_sec)
    # 3. Save results to Excel (single file)
    fname = f"auc_prrx_{db}_{x_sec}s.xlsx"
    writer = pd.ExcelWriter(os.path.join(auc_dir, fname), engine='openpyxl')
    for group in dfs_auc:
        dfs_auc[group].to_excel(writer, sheet_name=f'{group}', index=False)
    writer.save()


def find_optimal_cutoff(df, method):
    """Find optimal cutoff of pRRx/pRRx% using specified method

    Args:
        df (pd.DataFrame): DF with values of pRRx/pRRx% parameters
        method (str, optional): Method of optimal threshold calculation.
            Can be 'youden', 'dor_max' or '01_criterion'.

    Returns:
        (dict): Dict of DFs with cutoffs
    """
    dfs_cutoff = {}
    y_true = df['is_af'].values
    params = helper.get_params(df)
    groups = params.keys()  # Parameter types (pRRx, pRRx%)
    for group in groups:
        opt_cutoff = []
        features = params[group]
        for feature in features:
            negative = 0
            y_score = df[feature].values
            auc = metrics.roc_auc_score(y_true, y_score)
            if auc < 0.5:  # Classify as 1 if value < cutoff
                negative = 1
                y_score = -y_score
            # Get points on ROC curve
            fpr, tpr, thresholds = metrics.roc_curve(y_true, y_score)
            tnr = 1 - fpr
            fnr = 1 - tpr
            # 1. Optimal threshold maximizes informedness (TPR-FPR)
            # Also known as Youden's statistic
            informedness = tpr - fpr
            # 2. Optimal threshold closest to (0, 1) on ROC curve
            dist_01 = np.sqrt((1-tpr)**2 + fpr**2)
            # 3. Diagnostic Odds Ratio (DOR)
            dor = (tpr*tnr) / (fpr*fnr+1e-9)
            if method == 'youden':
                idx_opt = informedness.argmax()
            elif method == 'dor_max':
                idx_opt = np.nanargmax(dor)
                print(dor[idx_opt])
            elif method == '01_criterion':
                idx_opt = dist_01.argmin()
            else:
                print('Wrong method')
                return
            cutoff = thresholds[idx_opt]
            if negative:
                cutoff = -cutoff
                y_score = -y_score
            opt_cutoff.append(cutoff)
        dfs_cutoff[group] = pd.DataFrame({
            'features': features,
            f'{x_sec} s': opt_cutoff,
            })        
    return dfs_cutoff


def cutoff_prrx_to_excel(prrx_dir, db, method, x_sec, cutoff_dir):
    """Calculate optimal cutoff of pRRx/pRRx% and save in Excel.

    Args:
        prrx_dir (str): Directory with pRRx/pRRx% in CSV format.
        db (str): Acronym of the database.
        method (str): Method of optimal threshold calculation.
            Can be 'youden', 'dor_max' or '01_criterion'.
        x_sec (int): Length of RR sequence [s].
        cutoff_dir (str): Write directory for cutoffs.

    """
    # 1. Read data (pRRx, pRRx_%)
    df = pd.read_csv(os.path.join(prrx_dir, db, f"prrx_{db}_{x_sec}s.csv"))
    # 2. Calculate cutoffs
    dfs_cutoff = find_optimal_cutoff(df, method)
    # 3. Save in a single Excel file
    fname = f"cutoffs_prrx_{method}_{db}_{x_sec}s.xlsx"
    writer = pd.ExcelWriter(os.path.join(cutoff_dir, fname), engine='openpyxl')
    for group in dfs_cutoff:
        dfs_cutoff[group].to_excel(
            writer, sheet_name=f'{group}', index=False)
    writer.save()


if __name__ == '__main__':
    prrx_dir = '../data/processed'
    roc_dir = '../reports/excel/roc'
    x_sec = 60
    for db in ['ltafdb', 'afdb']:
        if not os.path.exists(roc_dir):
            os.makedirs(os.path.join(roc_dir))
        # auc_prrx_to_excel(prrx_dir, db, x_sec, roc_dir)
        cutoff_prrx_to_excel(
            prrx_dir, db, method='youden', x_sec=60, cutoff_dir=roc_dir)
