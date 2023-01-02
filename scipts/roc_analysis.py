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
        features = [x[0] for x in params[group]]
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


if __name__ == '__main__':
    prrx_dir = '../data/processed'
    auc_dir = '../reports/excel/roc'
    x_sec = 60
    for db in ['ltafdb', 'afdb']:
        if not os.path.exists(auc_dir):
            os.makedirs(os.path.join(auc_dir))
        auc_prrx_to_excel(prrx_dir, db, x_sec, auc_dir)
