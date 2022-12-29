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
    for name in rec_names:
        sr = []
        af = []
        print(f"\n++++++++++++\nRecord name: {name}:")

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
                print("LAST RHYTHM")
            print(f"*** Rhythm no. {i} at {t} [s]: {r}. Ends at {t_next} [s]")
            if r == "(N":
                sr.append(np.copy(qrs[
                    np.where((qrs >= t) & (qrs < t_next))]
                ))
                print(f"Appended {len(sr[-1])} QRS's to SR")
            elif r == "(AFIB":
                af.append(np.copy(
                    qrs[np.where((qrs >= t) & (qrs < t_next))]))
                print(f"Appended {len(af[-1])} QRS's to AF")
            qrs_dict[name] = {"sr": sr, "af": af}
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
            df.to_csv(
                folder + f"/{rec_name}_sr_{i}.csv",
                index=False, header=False)
        # AF segments
        for i, af in enumerate(qrs_dict[rec_name]["af"]):
            df = pd.DataFrame(af)
            df.to_csv(
                folder + f"/{rec_name}_af_{i}.csv",
                index=False, header=False)


def prepare_qrs(rec_dir, db_name):
    """From all WFDB records in a directory, get QRS locations and rhythm info.
    Save QRS locations of continuous AF and SR segments to CSV files.
    Args:
        rec_dir (str): Directory with WFDB records
        db_name (str): abbreviation that will later be used to identify database
    """

    # #### 1. READ DATA ####

    # Names of all QRS annotation files in directory
    print(f'{db = }')
    file_lst = os.listdir(rec_dir)
    file_lst = [fn for fn in file_lst if fn.lower().endswith('.qrs')]
    print(file_lst)

    # #### 2. GET QRS AND RHYTHMS FROM ANNOTATIONS ####

    # Read annotations (rhythm and QRS) to 'ann' dictionary
    # with record names as keys
    ann = {}
    for name in file_lst[:2]:
        print(name)
        name = Path(name).stem  # filename without extension
        ann[name] = ann_from_file(os.path.join(rec_dir, name))

    # For all records, divide QRS's to AF and SR
    # based on rhythm annotations and save in CSVs
    qrs_dict = get_qrs_w_rhythm(ann)
    save_qrs_to_files(qrs_dict, f"../data/interim/{db_name}/qrs")


if __name__ == '__main__':
    # 1. Save QRS locations of continuous AF and SR segments to CSV files
    for db in ['ltafdb', 'afdb']:
        rec_dir = f'D:/Matlab_data/physionet/databases/{db}/1.0.0'
        print(rec_dir)
        prepare_qrs(rec_dir, db)
    # 2. Calculate pRRx and pRRx%
