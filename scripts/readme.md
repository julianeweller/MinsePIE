# Setting up your environment to run MinsePIE

## Create environment
conda create --name pie2 python=3.8 --no-default-packages
conda activate pie2

## Install packages
conda install sys
conda install regex
conda install pandas
conda install -c conda-forge xgboost
pip install pandarallel
conda install -c bioconda viennarna
conda install -c conda-forge biopython
