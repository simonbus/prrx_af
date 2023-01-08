import os
import pandas as pd
import numpy as np


def rr_stats(qrs_dir, x_sec):
    """Info about very long and very short RR intervals

    Args:
        qrs_dir (str): root directory of QRS location info in CSV format
        x_sec (double): length of RR sequence [s]
    """
    # List of records
    rec_names = sorted(os.listdir(qrs_dir))
    print(f"{db = }")
    rr_dict = {rhythm: {rr_type: 0 for rr_type in
               ['total', '<240 ms', '<400 ms', '>1500 ms', '>3000 ms']}
               for rhythm in ['af', 'sr']}
    for rn in rec_names:
        # List of segments (CSV files with QRS)
        qrs_segm = sorted(os.listdir(os.path.join(qrs_dir, rn, 'qrs')))
        # Calculate RR for each QRS segment
        for qrs_name in qrs_segm:
            if 'af' in qrs_name:
                rhythm = 'af'
            else:
                rhythm = 'sr'
            df = pd.read_csv(os.path.join(qrs_dir, rn, 'qrs', qrs_name))
            qrs = df.values.flatten()
            start_idx = 0
            # Find x_sec long series of QRS
            for idx, r in enumerate(qrs):
                if qrs[idx] - qrs[start_idx] > x_sec:
                    # Update numbers of long and short RR intervals
                    rr = qrs[start_idx+1:idx] - qrs[start_idx:idx-1]
                    rr_dict[rhythm]['total'] += len(rr)
                    rr_dict[rhythm]['<240 ms'] += np.sum(rr < 0.24)
                    rr_dict[rhythm]['<400 ms'] += np.sum(rr < 0.4)
                    rr_dict[rhythm]['>1500 ms'] += np.sum(rr > 1.5)
                    rr_dict[rhythm]['>3000 ms'] += np.sum(rr > 3.0)
                    start_idx = idx - 1
    for rhythm in ['af', 'sr']:
        print(rhythm.upper())
        d = rr_dict[rhythm]
        for rr_type in d.keys():
            print(f'{rr_type}: {d[rr_type]} ({100*d[rr_type] / d["total"]:.2f}%)')


if __name__ == '__main__':
    x_sec = 60
    for db in ['ltafdb', 'afdb']:
        rr_stats(f'../data/interim/{db}/qrs', x_sec)