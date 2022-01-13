# MinsePIE  &nbsp; :pie:

### Modelling insertion efficiency for Prime Insertion Experiments

![Alt Text](img/minsepie_animation.gif)
</br></br>
Writing short sequences into the genome with prime eiditng  faciliates protein tagging, correction of pathogenic deletions and many more exciting applications. We studied the features that influence insertion efficiency and built a model to predict insertion rates based on the insert sequence. This helps users to choose optimal contructs for DNA insertion with prime editing. 

The provided model "MinsePIE.sav" was trained on 22974 events: a libary of 2,666 insert sequences up to 69 nt in length in four genomic sites (CLYBL, EMX1, FANCF, HEK3) in three human cell lines, using the PE2 prime editing system.

**System requirements**

- Python 3.8
- Python packages: argparse (1.4.0), [biopython (1.79)](https://biopython.org/wiki/Download), scikit-learn (0.24.2), scipy (1.5.3), [XGBoost (1.4.0)](https://xgboost.readthedocs.io/en/latest/install.html), pandas (1.3.4), [pandarallel (1.5.4)](https://github.com/nalepae/pandarallel), regex (2021.8.3), [RNAlib-2.4.18](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/examples_python.html)

If you encounter problems setting up the environment or packages, please check out the detailed description for installing the packages in the [scripts folder](https://github.com/julianeweller/MinsePIE/tree/main/scripts).

## Usage guide

The MinsePIE tools are constantly improving. Therefore, it is recommended to run clone the github repository and update it frequently:

```
# clone
git clone https://github.com/julianeweller/MinsePIE.git
# update
git pull
```


### Single insert mode
Input:
- sequence to be inserted (5' to 3' on target site)
- primer binding site and reverse transcriptase template (if target site on forward strand: respective 5' to 3' sequence on target site, otherwise the reverse complement in 5' to 3'; without insert sequence)
- MMR status of the cell line (default: MMR deficient)
- optional: expected mean and standard deviation for editing events in the experimental setup

Output:
- Predicted insertion efficiency score (z-score)
- optional: Conversion of z-score into insertion efficiency [%] based on expected mean and standard deviation for editing events in the experimental setup

Command line:
```
python minsepie.py -i [insert sequence] -p [PBS sequence] -r [RTT sequence]  -m [mmr status of cellline, default: 0] -a [optional: expected mean editing efficiency] -s [optional: expected standard deviation]
```
If you encounter problems setting up the environment or packages, please check out the detailed description for installing the packages in the [scripts folder](https://github.com/julianeweller/MinsePIE/tree/main/scripts).

**Example**

We want to insert ATAACTTCGATAATGTGATGCTATACGAAGTTAT into the HEK3 locus in HEK293T cells. The chosen primer binding site (PBS) is CAGACTGAGCACG and the reverse transcriptase template (RTT) is TGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCA. Please note that this corresponds to the target site, the pegRNA is in reverse complement to these sequences. 

Optionally, since we are not only interested in the predicted Z-score, but also want to get a feeling for what insertion rate this corresponds to, we input our expected average insertion rate and standard deviation based on previous experiences. During screens, we measured a mean insertion rate of 4.86 and a standard deviation of 4.28.

<img src="img/input_explanation.png" />

Command line:
```
python minsepie.py -i ATAACTTCGATAATGTGATGCTATACGAAGTTAT -p CAGACTGAGCACG -r TGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCA -a 4.86 -s 4.28
```
Output:
```
MMR status of cell line is considered as 0 (MMR deficient)
Insertion of ATAACTTCGATAATGTGATGCTATACGAAGTTAT 
Z-score: -0.4582706391811371 
Scaled score based on provided mean and standard deviation 2.898601664304733
```
### Batch mode
The batch mode is useful for evaluating different insert sequences into the same locus and where the pegRNA therefore has the same primer binding site and reverse transcription template for all insert sequences.

Input:
- File with list of sequence to be inserted (5' to 3' on target site), ".csv" or ".tsv"
- Output directory for table with prediction results
- primer binding site and reverse transcriptase template (if target site on forward strand: respective 5' to 3' sequence on target site, otherwise the reverse complement in 5' to 3'; without insert sequence)
- MMR status of the cell line (default: MMR deficient)
- optional: expected mean and standard deviation for editing events in the experimental setup

Output:
File with table that includes predicted insertion efficiency score (z-score) and optionally the conferion of the z-score into insertion efficiency [%] based on expected mean and standard deviation for editing events in the experimental setup

Command line:
```
python minsepie-batchinsert.py -i [path to input file] -p [PBS sequence] -r [RTT sequence] -o [path to output directory] -m [mmr status of cellline, default: 0] -a [optional: expected mean editing efficiency] -s [optional: expected standard deviation]
```

**Example**
In this example we want to add a His6-tag into the HEK3 locus in HEK293T cells. We compare 64 codon variants for the His6 tag with a combination of CAC and CAT codons (see example/his6_codonvariants.csv which serves as input). 

<p align="center">
  <img src="img/his6_codonvariants_table.png" />
</p>

The chosen primer binding site (PBS) is CAGACTGAGCACG and the reverse transcriptase template (RTT) is TGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCA. Please note that this corresponds to the target site, the pegRNA is in reverse complement to these sequences (see single insert mode). 

Optionally, since we are not only interested in the predicted Z-score, but also want to get a feeling for what insertion rate this corresponds to, we input our expected average insertion rate and standard deviation based on previous experiences. During screens, we measured a mean insertion rate of 4.86 and a standard deviation of 4.28.

Command line:
```
python minsepie-batchinsert.py -i his6_codonvariants.csv -p cagactgagcacg -r TGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCA -o . -a 4.86 -s 4.28
```
Output:
The three inserts with the highest predicted insertion rate will be displayed.
```
MMR status of cell line is considered as 0 (MMR deficient)
                insert    zscore  percIns_predicted
0   CACCACCACCACCACCAC  2.288961          14.656754
1   CATCACCACCACCACCAC  1.784901          12.499376
32  CACCACCACCACCACCAT  1.784901          12.499376
Prediction results are saved as his6_codonvariants_minsepie.csv.
```
The full prediction is saved as table in the defined output directory (see example/his6_codonvariants_predictions.csv).
<p align="center">
  <img src="img/his6_codonvariants_predictions.png" />
</p>

### Running MinsePIE with a graphical user interface
To facilitate the usage, we provide a simple user interface that allows the same functionalities as the command line tools. To use it, please install [Gooey v1.0.8] (https://pypi.org/project/Gooey/) and [Pythonw v.3.0.3] (https://pypi.org/project/pythonw/). To open the user interface run the following command:

```
pythonw minsepie_gui.py
```
An application will open. First choose between single and batch mode, then fill in the required information for running MinsePIE. 

<p align="center">
  <img src="img/gui_singlemode.png" width="300"/> 
  <img src="img/gui_batchmode.png" width="300"/>
</p>

If you encounter problems setting up the environment or packages, please check out the detailed description for installing the packages in the [scripts folder](https://github.com/julianeweller/MinsePIE/tree/main/scripts).

## Reference

Predicting efficiency of writing short sequences into the genome using prime editing </br>
Jonas Koeppel, Elin Madli Peets, Juliane Weller, Ananth Pallaseni, Fabio Liberante, Leopold Parts </br>
bioRxiv https://www.biorxiv.org/content/10.1101/2021.11.10.468024v1 </br>
doi: https://doi.org/10.1101/2021.11.10.468024
