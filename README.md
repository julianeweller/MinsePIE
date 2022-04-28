# MinsePIE  &nbsp; :pie:
Find the code for the manuscript resubmission in [scripts](https://github.com/julianeweller/MinsePIE/blob/main/scripts/modelling.py).

### Modelling insertion efficiency for Prime Insertion Experiments

![Alt Text](img/minsepie_animation.gif)
</br></br>
Writing short sequences into the genome with prime eiditng  faciliates protein tagging, correction of pathogenic deletions and many more exciting applications. We studied the features that influence insertion efficiency and built a model to predict insertion rates based on the insert sequence. This helps users to choose optimal contructs for DNA insertion with prime editing. 

The provided model "MinsePIE.sav" was trained on 22974 events: a libary of 2,666 insert sequences up to 69 nt in length in four genomic sites (CLYBL, EMX1, FANCF, HEK3) in three human cell lines, using the PE2 prime editing system.

**System requirements**

- Python 3.8
- Python packages: argparse (1.4.0), more_itertools (8.12.0),[biopython (1.79)](https://biopython.org/wiki/Download), scikit-learn (0.24.2), scipy (1.5.3), [XGBoost (1.5.0)](https://xgboost.readthedocs.io/en/latest/install.html), pandas (1.3.4), [pandarallel (1.5.4)](https://github.com/nalepae/pandarallel), regex (2021.8.3), [RNAlib-2.4.18](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/examples_python.html)


If you encounter problems setting up the environment or packages, please check out the detailed description for installing the packages in the [scripts folder](https://github.com/julianeweller/MinsePIE/tree/main/scripts).

## Usage guide

The MinsePIE tools are constantly improving. Therefore, it is recommended to run clone the github repository and update it frequently:

```
# clone
git clone https://github.com/julianeweller/MinsePIE.git
# update
git pull
# install MinsePIE: go into the folder with setup.py
pip install .

```

Here is an example on how to use minsepie in python:
```
import minsepie
minsepie.predict(['TGTCA'], pbs = 'CAGACTGAGCACG', ha = 'TGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCA', spacer = 'GGCCCAGACTGAGCACGTGA', mmr = 0, outdir = "./")

```


## Reference

Predicting efficiency of writing short sequences into the genome using prime editing </br>
Jonas Koeppel, Elin Madli Peets, Juliane Weller, Ananth Pallaseni, Fabio Liberante, Leopold Parts </br>
bioRxiv https://www.biorxiv.org/content/10.1101/2021.11.10.468024v1 </br>
doi: https://doi.org/10.1101/2021.11.10.468024
