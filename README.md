# MinsePIE  &nbsp; :pie:
**Modelling insertion efficiency for Prime Insertion Experiments**

Writing short sequences into the genome with prime eiditng  faciliates protein tagging, correction of pathogenic deletions and many more exciting applications. We studied the features that influence insertion efficiency and built a model to predict insertion rates based on the insert sequence. This helps users to choose optimal contructs for DNA insertion with prime editing. 

**System requirements**

- Python 3.8
- Python packages: biopython (1.79), scikit-learn (0.24.2), scipy (1.5.3), XGBoost (1.4.0)

**Usage guide**

Input:
- sequence to be inserted
- primer binding site and reverse transcriptase templates

Output:
- Prediction insertion efficiency (z-score)


