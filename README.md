# MinsePIE  &nbsp; :pie:
**Modelling insertion efficiency for Prime Insertion Experiments**

Writing short sequences into the genome with prime eiditng  faciliates protein tagging, correction of pathogenic deletions and many more exciting applications. We studied the features that influence insertion efficiency and built a model to predict insertion rates based on the insert sequence. This helps users to choose optimal contructs for DNA insertion with prime editing. 

The provided model "MinsePIE.sav" was trained on 22974 events: a libary of 2,666 insert sequences up to 69 nt in length in four genomic sites (CLYBL, EMX1, FANCF, HEK3) in three human cell lines, using the PE2 prime editing system.

esigned a library of 2,666 sequences up to 69 nt in length and measured the frequency of their insertion into four genomic sites in three human cell lines, using different prime editor systems. 

Also, this is maybe in supplemental or the paper somewhere, but is the .sav file the model that was trained on ALL the data in the paper? It might be worth putting that in the readme, something like "Model included here was trained on 2,600+ pegRNA events targeting 293T & HAP1 at CLYBL etc.". Can one use the scripts to train one's own model in future? Looking at the code it appears so? 

**System requirements**

- Python 3.8
- Python packages: biopython (1.79), scikit-learn (0.24.2), scipy (1.5.3), XGBoost (1.4.0), pickle (3.10.0), RNAlib-2.4.18: https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/examples_python.html

**Usage guide**

Input:
- sequence to be inserted
- primer binding site and reverse transcriptase templates
- MMR status of the cell line (default: MMR deficient)
- optional: expected mean and standard deviation for editing events in the experimental setup

Output:
- Predicted insertion efficiency score (z-score)
- optional: Conversion of z-score into insertion efficiency [%] based on expected mean and standard deviation for editing events in the experimental setup

Command line:

python minsepie.py -i [insert sequence] -p [PBS sequence] -r [RTT sequence]  -m [mmr status of cellline, default: 0] -a [optional: expected mean editing efficiency] -s [optional: expected standard deviation]

**Reference**

Predicting efficiency of writing short sequences into the genome using prime editing </br>
Jonas Koeppel, Elin Madli Peets, Juliane Weller, Ananth Pallaseni, Fabio Liberante, Leopold Parts </br>
bioRxiv https://www.biorxiv.org/content/10.1101/2021.11.10.468024v1; doi: https://doi.org/10.1101/2021.11.10.468024
