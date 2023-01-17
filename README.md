# prrx_af

This repository can be used for analysis of the diagnostic properties of [pRRx](https://www.mdpi.com/2077-0383/11/19/5702) and [pRRx%](https://www.mdpi.com/2077-0383/12/2/687) parameters in atrial fibrillation (AF) detection.
pRRx and pRRx% are families of heart rate variability (HRV) parameters.

## How to use this repository?

### 1. Download code and databases

1. Clone the repository.
    ```
    git clone git@github.com:simonbus/prrx_af.git
    ```
2. Download databases from Physionet ([Long Term AF Database (LTAFDB)](https://physionet.org/content/ltafdb/1.0.0/) and [MIT-BIH Atrial Fibrillation Database (AFDB)](https://physionet.org/content/afdb/1.0.0/))
    * Unzip the files.
    * Copy the content of the folders to `data/raw/ltafdb/1.0.0` and `data/raw/afdb/1.0.0`, respectively.

### 3. Use the demo Jupyter notebook to analyze the data

* [scipts/prrx_analysis.ipynb](scipts/prrx_analysis.ipynb)

### 4. Modify the analysis parameters

Different parameters of analysis can be changes, such as:
* `x_sec` - length of the analyzed RR series
* `N` - number of repetitions in nonparametric bootstrap
* `cutoff_method` - method of choosing the optimal cutoff value