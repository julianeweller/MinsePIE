# Setting up an environment to run MinsePIE
The instructions below enable you to create a minimal conda-based virtual environment to run MinsePIE. 
If you have trouble running MinsePIE, please check that you have the correct versions of the packages installed (see below).

## Install conda
On MacOS, run:
```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
sh Miniconda3-latest-MacOSX-x86_64.sh
```

On Linux, run:
```
curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
```

## Create conda environment
```
conda create --name pie38 python=3.8 --no-default-packages
conda activate pie38
```

## Install packages
```
conda install sys regex pandas
conda install -c conda-forge biopython xgboost
conda install -c bioconda viennarna
pip install pandarallel
```
## Run MinsePIE
See main [README](/README.md).

Don't forget to deactivate your conda environment when you are finished running MinsePIE
`conda deactivate`

## Environment packages
If you are having issues running MinsePIE, check the output of `conda list` matches the minimal MinsePIE environment below.
```
#Name                    Version                   Build  Channel
_py-xgboost-mutex         2.0                       cpu_0    conda-forge
biopython                 1.79             py38h96a0964_1    conda-forge
ca-certificates           2021.10.8            h033912b_0    conda-forge
colored                   1.4.3                    pypi_0    pypi
dill                      0.3.4                    pypi_0    pypi
gooey                     1.0.8.1                  pypi_0    pypi
joblib                    1.1.0              pyhd8ed1ab_0    conda-forge
libblas                   3.9.0           12_osx64_openblas    conda-forge
libcblas                  3.9.0           12_osx64_openblas    conda-forge
libcxx                    12.0.1               habf9029_0    conda-forge
libffi                    3.4.2                h0d85af4_5    conda-forge
libgfortran               5.0.0           9_3_0_h6c81a4c_23    conda-forge
libgfortran5              9.3.0               h6c81a4c_23    conda-forge
liblapack                 3.9.0           12_osx64_openblas    conda-forge
libopenblas               0.3.18          openmp_h3351f45_0    conda-forge
libxgboost                1.5.0                h4a89273_1    conda-forge
libzlib                   1.2.11            h9173be1_1013    conda-forge
llvm-openmp               12.0.1               hda6cdc1_1    conda-forge
ncurses                   6.2                  h2e338ed_4    conda-forge
numpy                     1.21.4           py38h49b9922_0    conda-forge
openssl                   3.0.0                h0d85af4_2    conda-forge
pandarallel               1.5.4                    pypi_0    pypi
pandas                    1.3.5            py38ha53d530_0    conda-forge
pillow                    8.4.0                    pypi_0    pypi
pip                       21.3.1             pyhd8ed1ab_0    conda-forge
psutil                    5.8.0                    pypi_0    pypi
py-xgboost                1.5.0            py38h50d1736_1    conda-forge
pygtrie                   2.4.2                    pypi_0    pypi
python                    3.8.12          h43ca1e7_2_cpython    conda-forge
python-dateutil           2.8.2              pyhd8ed1ab_0    conda-forge
python_abi                3.8                      2_cp38    conda-forge
pytz                      2021.3             pyhd8ed1ab_0    conda-forge
readline                  8.1                  h05e3726_0    conda-forge
regex                     2021.11.10       py38h96a0964_0    conda-forge
scikit-learn              1.0.1            py38h37f3bb3_3    conda-forge
scipy                     1.7.3            py38hd329d04_0    conda-forge
setuptools                59.6.0           py38h50d1736_0    conda-forge
six                       1.16.0             pyh6c4a22f_0    conda-forge
sqlite                    3.37.0               h23a322b_0    conda-forge
threadpoolctl             3.0.0              pyh8a188c0_0    conda-forge
tk                        8.6.11               h5dbffcc_1    conda-forge
viennarna                 2.4.18           py38he0002f0_0    bioconda
wheel                     0.37.0             pyhd8ed1ab_1    conda-forge
wxpython                  4.1.1                    pypi_0    pypi
xgboost                   1.5.0            py38hbb4f172_1    conda-forge
xz                        5.2.5                haf1e3a3_1    conda-forge
zlib                      1.2.11            h9173be1_1013    conda-forge
```
