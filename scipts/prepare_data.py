import os
import pandas as pd
import numpy as np
from pathlib import Path
import wfdb


def ann_from_file(full_path):
    """Read annotations from a single WFDB file and return them in a dictionary.
    
    Args:
        full_path (str): Full path to file (w/o extension)
    
    Returns:
        [dict]: QRS and rhythm annotations
    """
    # QRS annotation file
    ann_qrs = wfdb.io.rdann(
        record_name=full_path,  # without extension
        extension="qrs")
    # Corresponding rhythm annotation file
    ann_atr = wfdb.io.rdann(
        record_name=full_path,  # without extension
        extension="atr")
    # QRS locations
    qrs_loc = ann_qrs.sample / ann_qrs.fs
    # Rhythm annotation location
    atr_locs = ann_atr.sample / ann_qrs.fs  # [s]
    atr_rhythms = ann_atr.aux_note
    atr_loc = [loc for loc, rhythm in zip(atr_locs, atr_rhythms)
               if len(rhythm) > 0]
    atr_rhythm = [rhythm for loc, rhythm in zip(atr_locs, atr_rhythms)
                  if len(rhythm) > 0]
    return {
        "qrs_loc": qrs_loc,
        "atr_loc": atr_loc,
        "atr_rhythm": atr_rhythm
    }


def get_qrs_w_rhythm(ann):
    """Get a dictionary of QRS locations with corresponding rhythms
    
    Args:
        ann (dict): each element has a filename as key and is a dict returned by ann_from_file()
    
    Returns:
        dict: keys are filenames, each element is a dict
            * 'sr' and 'af' keys  
            * both elements are lists of 1D NumPy arrays with QRS locations 
    """
    qrs_dict = {}
    rec_names = ann.keys()
    for rec_num, name in enumerate(rec_names):
        sr = []
        af = []
        print(f"Associating rhythms with QRS's in record ({rec_num+1} / {len(rec_names)})", end='\r')

        num_rhythms = len(ann[name]["atr_loc"])
        # Iterate through rhythms
        for i, t in enumerate(ann[name]["atr_loc"]):  # [:-1]):
            # Rhythm type
            qrs = ann[name]["qrs_loc"]
            r = ann[name]["atr_rhythm"][i]
            # Next rhythm start
            if i < num_rhythms - 1:
                t_next = ann[name]["atr_loc"][i+1]
            else:
                t_next = qrs[-1]
            if r == "(N":
                sr.append(np.copy(qrs[
                    np.where((qrs >= t) & (qrs < t_next))]
                ))
            elif r == "(AFIB":
                af.append(np.copy(
                    qrs[np.where((qrs >= t) & (qrs < t_next))]))
            qrs_dict[name] = {"sr": sr, "af": af}
    print('')
    return qrs_dict


def save_qrs_to_files(qrs_dict, out_dir):
    """Save QRS of continuous AF and SR segments to CSV files
    
    Args:
        qrs_dict (dict): annotations dict returned by get_qrs_w_rhythm()
        out_dir (str): Directory to save CSV
    """
    for rec_name in qrs_dict:
        # Make directory for each record
        folder = os.path.join(out_dir, rec_name, 'qrs')
        if not os.path.exists(folder):
            os.makedirs(folder)
        # SR segments
        for i, sr in enumerate(qrs_dict[rec_name]["sr"]):
            df = pd.DataFrame(sr)
            if len(df.index > 0):
                df.to_csv(
                    folder + f"/{rec_name}_sr_{i}.csv",
                    index=False, header=False)
        # AF segments
        for i, af in enumerate(qrs_dict[rec_name]["af"]):
            df = pd.DataFrame(af)
            if len(df.index > 0):
                df.to_csv(
                    folder + f"/{rec_name}_af_{i}.csv",
                    index=False, header=False)


def prepare_qrs(rec_dir, db):
    """From all WFDB records in a directory, get QRS locations and rhythm info.
    Save QRS locations of continuous AF and SR segments to CSV files.
    
    Args:
        rec_dir (str): Directory with WFDB records
        db (str): abbreviation that will later be used to identify database
    """

    # #### 1. READ DATA ####

    # Names of all QRS annotation files in directory
    print(f'{db = }')
    file_lst = os.listdir(rec_dir)
    file_lst = [fn for fn in file_lst if fn.lower().endswith('.qrs')]

    # #### 2. GET QRS AND RHYTHMS FROM ANNOTATIONS ####

    # Read annotations (rhythm and QRS) to 'ann' dictionary
    # with record names as keys
    ann = {}
    for i, name in enumerate(file_lst):
        print(f'Reading annotations from record {i+1}/{len(file_lst)}', end='\r')
        name = Path(name).stem  # filename without extension
        ann[name] = ann_from_file(os.path.join(rec_dir, name))

    # For all records, divide QRS's to AF and SR
    # based on rhythm annotations and save in CSVs
    print('')
    qrs_dict = get_qrs_w_rhythm(ann)
    save_qrs_to_files(qrs_dict, f"../data/interim/{db}/qrs")


def prrx_dict_from_rr(rr, prr_r, prr_perc_r):
    """Get dict of pRRx / pRRx% parameters from RR sequence.

    Args:
        rr (np.array): 1D array - continuous sequence of RR intervals
        prr_r (np.array): 1D array - range of pRRx thresholds [ms]
        prr_perc_r (np.array): 1D array - range of pRRx% thresholds [%]

    Returns:
        (dict): parameters/metainfo as keys, single value for each item
    """
    rr = rr[(rr >= 0.24) & (rr <= 3.0)]  # Filter RR's [s]
    sd = rr[1:] - rr[:-1]
    prrx_dict = {}
    prrx_dict['mean RR'] = np.mean(rr)
    if len(rr) <= 3:
        return None
    # pRRx
    if prr_r is not None:
        for x in np.arange(prr_r[0], prr_r[1], prr_r[2]):
            prrx = 100 * np.count_nonzero((np.abs(sd) >= x/1000)) / len(sd)
            prrx_dict[f'pRR{x}'] = prrx
    # pRRx%
    if prr_perc_r is not None:
        for x_perc in np.arange(prr_perc_r[0], prr_perc_r[1], prr_perc_r[2]):
            prrx_perc = 100 * np.count_nonzero((np.abs(sd) >= rr[:-1] * x_perc/100)) / len(sd)
            prrx_dict[f'pRR{x_perc}%'] = prrx_perc
    return prrx_dict


def qrs_dir_to_prrx_df(qrs_dir, x_sec, prr_r, prr_perc_r):
    """Prepare Dataframe with assorted pRRx/pRRx% parameters from all records in a directory

    Args:
        qrs_dir (str): root directory of QRS location info in CSV format
        x_sec (double): length of RR sequence [s]
        prr_r (np.array): 1D array - range of pRRx thresholds [ms]
        prr_perc_r (np.array): 1D array - range of pRRx% thresholds [%]

    Returns:
        (pd.DataFrame): Dataframe with assorted pRRx/pRRx% parameters from all records
    """
    # List of records
    rec_names = sorted(os.listdir(qrs_dir))
    print(f"Record names: {rec_names}")
    prrx_dict = {}
    meta = ['rec_name', 't_start', 'len_sec', 'num_rr', 'is_af']
    for col in meta:
        prrx_dict[col] = []
    prrx_dict['mean RR'] = []
    for rn in rec_names:
        # List of segments (CSV files with QRS)
        qrs_segm = sorted(os.listdir(os.path.join(qrs_dir, rn, 'qrs')))
        # Calculate RR for each QRS segment
        for qrs_name in qrs_segm:
            # Get AF/SR label from QRS file name
            if 'af' in qrs_name:
                is_af = 1
            else:
                is_af = 0
            df = pd.read_csv(os.path.join(qrs_dir, rn, 'qrs', qrs_name))
            qrs = df.values.flatten()
            start_idx = 0
            # Find x_sec long series of QRS
            for idx, r in enumerate(qrs):
                if qrs[idx] - qrs[start_idx] > x_sec:
                    # HRV from x_sec long segment
                    rr = qrs[start_idx+1:idx] - qrs[start_idx:idx-1]
                    prrx_single = prrx_dict_from_rr(rr, prr_r, prr_perc_r)
                    if prrx_single is not None:
                        prrx_dict['rec_name'].append(rn)
                        prrx_dict['t_start'].append(qrs[start_idx])
                        prrx_dict['len_sec'].append(x_sec)
                        prrx_dict['num_rr'].append(idx - 1 - start_idx)
                        prrx_dict['is_af'].append(is_af)
                        for feature in prrx_single:
                            if feature in prrx_dict:
                                prrx_dict[feature].append(prrx_single[feature])
                            else:
                                prrx_dict[feature] = [prrx_single[feature]]
                    start_idx = idx - 1
    prrx_df = pd.DataFrame.from_dict(prrx_dict)
    return prrx_df


def prepare_prrx(db, fs, x_sec, qrs_dir, prrx_dir):
    """Calculate pRRx/pRRx% parameters and save in CSV files

    Args:
        db (str): Acronym of the database
        fs (double): Sampling frequency [fs]
        x_sec (int): Length of RR sequence [s]
        qrs_dir (str): Root directory of QRS location info in CSV format
        prrx_dir (str): Directory for pRRx/pRRx% data files in CSV format
    """
    # Make directory for results if it doesn't exist
    if not os.path.exists(os.path.join(prrx_dir, db)):
        os.makedirs(os.path.join(prrx_dir, db))
    # pRRx_ms, pRRx_%
    prr_step_ms = 1000 / fs  # pRRx threshold x increment step
    prrx_step_perc = 0.25  # pRRx% threshold x% increment step
    prrx_df = qrs_dir_to_prrx_df(
        qrs_dir=os.path.join(qrs_dir, f'{db}/qrs'),
        x_sec=x_sec,
        prr_r=(prr_step_ms, 201, prr_step_ms),
        prr_perc_r=(prrx_step_perc, 25.1, prrx_step_perc),
    )
    print(prrx_df.head())
    prrx_df.to_csv(
        os.path.join(prrx_dir, f"{db}/prrx_{db}_{x_sec}s.csv")
        )


if __name__ == '__main__':
    # 1. Save QRS locations of continuous AF and SR segments to CSV files
    for db in ['ltafdb', 'afdb']:
        rec_dir = f'../data/raw/{db}/1.0.0'
        prepare_qrs(rec_dir, db)
    # 2. Calculate pRRx and pRRx%
    x_sec = 60  # Length of RR segments [s]
    for db, fs in [['ltafdb', 128], ['afdb', 128]]:
        prepare_prrx(
            db, fs, x_sec,
            qrs_dir='../data/interim',
            prrx_dir='../data/processed')
