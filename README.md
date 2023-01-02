# prrx_af

This is a repository for the analysis of diagnostic properties of pRRx and pRRx% parameters in atrial fibrillation detection. pRRx and pRRx% are families of heart rate variability (HRV) parameters.

## How to use this repository?

### 1. Download code and databases

1. Clone the repository.
    ```
    git clone git@github.com:simonbus/prrx_af.git
    ```
2. Download and unzip databases from Physionet
    * [MIT-BIH Atrial Fibrillation Database](https://physionet.org/content/afdb/1.0.0/)
    * [Long Term AF Database](https://physionet.org/content/ltafdb/1.0.0/)

### 2. Calculate pRRx and pRRx% parameters

* [scipts/prepare_data.py](scipts/prepare_data.py)

### 3. Analyze the data

1. [scipts/distribution.py](scipts/distribution.py)