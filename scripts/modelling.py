import os
from scipy import stats
import sys
import logging
import numpy as np
import pandas as pd
from pandarallel import pandarallel

import xgboost
import sklearn
import pickle
import matplotlib.pyplot as plt

import pandas as pd
import logging
import sys
from scipy import stats
import numpy as np
import sklearn
import xgboost
import RNA
import random
from joblib import delayed, Parallel

from sklearn.model_selection import train_test_split
from scipy.cluster import hierarchy as sch
from sklearn import linear_model
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import RandomForestRegressor

import shap
import regex as re
from Bio.SeqUtils import MeltingTemp as mt
import RNA
import pandas as pd
import numpy as np
from Bio.Data.IUPACData import ambiguous_dna_values
from pandarallel import pandarallel
import math
from Bio import Align
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

#############################################################
###########             INITIALIZE                ########### 
#############################################################

pandarallel.initialize()
verbosity = 'info'
workdir = '../'
workdir_input = os.path.join(workdir, 'files/input')
workdir_output = os.path.join(workdir, 'files/output') # Change to your prefered output directory if you don't want to overwrite the files here
setseed = 3

# Input
screendata = os.path.join(workdir_input, 'data_screen.csv')
canadian_path = os.path.join(workdir_input,'data_canadian.csv')
barnacle_path = os.path.join(workdir_input,'data_barnacle.csv')
padding_path = os.path.join(workdir_input,'data_padding.csv')
pegRNA_mapping = os.path.join(workdir_input, 'Table_S2_pegRNAs.csv')

VF_baseline_data_insert = os.path.join(workdir_input, 'VF_baseline_meanstd.tsv')
VF_baseline_data_RTT = os.path.join(workdir_input, 'VF_baseline_RTT_meanstd.tsv')
Tm_baseline_data_insert = os.path.join(workdir_input,'Tm_baseline_meanstd.tsv')
Tm_baseline_data_RTT = os.path.join(workdir_input, 'Tm_baseline_RTT_meanstd.tsv')
VF_baseline_canadian_data_RTT = os.path.join(workdir_input, 'VF_baseline_canadian_RTT_meanstd.tsv')
VF_baseline_barnacle_data_RTT = os.path.join(workdir_input, 'VF_baseline_barnacle_RTT_meanstd.tsv')
VF_baseline_padding_data_RTT = os.path.join(workdir_input, 'VF_baseline_padding_RTT_meanstd.tsv')

# Output
data_featurized = os.path.join(workdir_output, 'data_screen_onehot.csv')
data_canadian_featurized = os.path.join(workdir_output, 'data_canadian_screen_onehot.csv')
data_barnacle_featurized = os.path.join(workdir_output, 'data_barnacle_screen_onehot.csv')
data_padding_featurized = os.path.join(workdir_output, 'data_padding_screen_onehot.csv')
datapath_train = os.path.join(workdir_output, 'data_train.csv')
datapath_test = os.path.join(workdir_output, 'data_test.csv')

model_path = os.path.join(workdir, 'models/MinsePIE_v3.sav')

path_corrmatrix_all = os.path.join(workdir_output, 'model_features_corrmatrix_all.csv')
path_corrmatrix_core = os.path.join(workdir_output, 'model_features_corrmatrix_core.csv')
model_crossvalidation_path = os.path.join(workdir_output, 'model_architecture_crossvalidation.csv')
model_addfeatures_path = os.path.join(workdir_output, 'model_features_sequentialadd.csv')
model_rmvfeatures_path = os.path.join(workdir_output, 'model_features_rmv.csv')
model_sites_crossvalidation_path = os.path.join(workdir_output, 'model_sites_crossvalidation.csv')

model_performance_path = os.path.join(workdir_output,'model_performance.csv')
model_performance_ca_path = os.path.join(workdir_output,'model_performance_ca.csv')
model_performance_ba_path = os.path.join(workdir_output,'model_performance_ba.csv')
model_performance_pa_path = os.path.join(workdir_output,'model_performance_pa.csv')
model_shap_impact_path = os.path.join(workdir_output,'model_shap_impact.pdf')
model_shap_importance_path = os.path.join(workdir_output,'model_shap_importance.pdf')

model_vffeatures_path = os.path.join(workdir_output,'model_performance_VF.csv')
model_mmrfeatures_path = os.path.join(workdir_output,'model_performance_mmr.csv')


# Set up logging based on the verbosity level set by the command line arguments:
logging.basicConfig(format='%(levelname)s: %(message)s', level=verbosity.upper())

#############################################################
###########             Functions                 ########### 
#############################################################
def error(msg, exit_code=1):
    """A function to quit with an error."""
    logging.error(msg)
    sys.exit(exit_code)

def writeDataFile(x, filename, datatype):
    """A function to write data to a given file, with correct logging"""
    try:
        with open(filename, 'w') as output_handle: output_handle.write(x)
    except Exception as err: error('failed to write {} data to {} ({})'.format(datatype, filename, err))

def enhance_feature_exploration(df, seq='insert', rtt='RTT', pbs='PBS', guide = 'spacer'):
    # sequences
    df['ext'] = df[pbs] + df[seq] + df[rtt]
    df['PBSinsert'] = df[pbs] + df[seq]
    # length features
    df['div3'] = df.length.apply(lambda x: True if x%3 == 0 else False)
    df['smaller4'] = df[seq].apply(lambda x: get_smallerI(x, 3))
    df['smaller7'] = df[seq].apply(lambda x: get_smallerI(x, 6))
    df['smaller13'] = df[seq].apply(lambda x: get_smallerI(x, 12))
    df['longer40'] = df[seq].apply(lambda x: get_higherI(x, 40))
    # nucleotide features
    df['N1'] = df[seq].apply(lambda x: x[0])
    df['N2'] = df[seq].apply(lambda x: x[1] if (len(x) >= 2) else np.NaN)
    df['N3'] = df[seq].apply(lambda x: x[2] if (len(x) >= 3) else np.NaN)
    df['N4'] = df[seq].apply(lambda x: x[3] if (len(x) >= 4) else np.NaN)
    df['NN'] = df[seq].apply(lambda x: x[-1])
    df['N123'] = df[seq].apply(lambda x: x[0:3])
      # Find runs
    df['Arun'] = df[seq].apply(lambda x: get_Nrun(x, 'A'))
    df['Crun'] = df[seq].apply(lambda x: get_Nrun(x, 'C'))
    df['Trun'] = df[seq].apply(lambda x: get_Nrun(x, 'T'))
    df['Grun'] = df[seq].apply(lambda x: get_Nrun(x, 'G'))
    df['Crun_maxlen'] = df[seq].apply(lambda x: get_Nrun_max(x, 'C'))
    df['Trun_maxlen'] = df[seq].apply(lambda x: get_Nrun_max(x, 'T'))
    df['Grun_maxlen'] = df[seq].apply(lambda x: get_Nrun_max(x, 'G'))
    # Alignment
    df['nick2match'] = df.apply(lambda x: pair_bases_pos(x[guide], x['insertRTT'], 1), axis=1)
    df['nick3match'] = df.apply(lambda x: pair_bases_pos(x[guide], x['insertRTT'], 2), axis=1)
    df['align_InsSpc'] = df.parallel_apply(lambda x: score_alignment(x[seq],x[guide]), axis=1)
    df['maxmicrohom_InsHA'] = df.apply(lambda x: length_maxmicrohomology(x[seq], x[rtt]),axis=1)
    df['loops_intact'] = df[seq].parallel_apply(scaffold_intact)
    df['loop1_intact'] = df['loops_intact'].apply(lambda x: x[0])
    df['loop2_intact'] = df['loops_intact'].apply(lambda x: x[1])
    # Structure
    df['VF_insert'] = df[seq].parallel_apply(get_vf)
    df['Tm_insert'] = df[seq].apply(get_tm)
    df['VF_ext'] = df['ext'].parallel_apply(get_vf)
    df['VF_PBSinsert'] = df['PBSinsert'].parallel_apply(get_vf)
    df['Tm_RTT'] = df['insertRTT'].apply(get_tm)
    return df

def splitdata(data, mainexperiment = 'HEK3 HEK293T', fx = 'sequence_original', fy1 = 'percIns_z', fy2 = 'insbins', seed = 12, test_size = 0.1):
    # I want to have each insert sequence only either in train and test data set, but they are used for different target sites and therefore repeat in the dataset --> just use the HEK3 HEK293T site to assign sequences initially
    seqQ = data.loc[data.experiment == mainexperiment]
    # Split into length bin and take a representative subset for each bin into test and train
    seqQ_dict = {}
    bins = [5, 10, 15, 20, 25, 30, 40, 50, 60]
    for i in bins:
        if int(i) == 60:
            seqQ_dict["seqQ"+str(int(i))] = seqQ.loc[seqQ.bin.astype(int) >= 60].copy()
        else:
            seqQ_dict["seqQ"+str(int(i))] = seqQ.loc[seqQ.bin == int(i)].copy()
            
    # For each length bin, take a representative subset and assign this to the test or train dataset
    train_seq = []
    test_seq = []
    for i in seqQ_dict:
        seqQ_dict[i][fy2] = pd.qcut(seqQ_dict[i][fy1], q=5, labels=range(0, 5))
        X = seqQ_dict[i][fx]
        y1 = seqQ_dict[i][fy1]
        y2 = seqQ_dict[i][fy2]
        X_train_temp, X_test_temp, y_train_temp, y_test_temp = train_test_split(X, y1, test_size= test_size, random_state=seed, stratify = y2)
        # add them to test and training set
        train_seq = train_seq + list(X_train_temp)
        test_seq = test_seq + list(X_test_temp)
        
    # Add data not included in the HEK3 HEK293T experimeriment to the dataset based on their percIns
    seqQ = data[~data[fx].isin(train_seq+test_seq)].drop_duplicates(subset=[fx])
    seqQ[fy1] = pd.qcut(seqQ[fy1], q=5, labels=range(0, 5))
    X = seqQ[fx]
    y1 = seqQ[fy1]
    y2 = seqQ[fy1]
    X_train_temp, X_test_temp, y_train_temp, y_test_temp = train_test_split(X, y1, test_size= test_size, random_state=seed, stratify = y2)

    # add them to test and training set
    train_seq = train_seq + (list(X_train_temp))
    test_seq = test_seq + (list(X_test_temp))
    logging.info(f"Final Training set sequences: {len(train_seq)}, Test set sequences: {len(test_seq)}")
    
    # Merge data from all target sites depending on sequences
    X_train = data[data[fx].isin(train_seq)]
    X_test = data[data[fx].isin(test_seq)]
    logging.info(f"Original dataset length: {data.shape[0]}, Final Training set:  {X_train.shape[0]}, Test set:  {X_test.shape[0]}")
    
    return(X_train, X_test)

def cluster_corr(corr_array, inplace=False):
    """
    Rearranges the correlation matrix, corr_array, so that groups of highly 
    correlated variables are next to eachother 
    
    Parameters
    ----------
    corr_array : pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix 
        
    Returns
    -------
    pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix with the columns and rows rearranged

    Adapted from https://wil.yegelwel.com/cluster-correlation-matrix/. 
    """
    pairwise_distances = sch.distance.pdist(corr_array)
    linkage = sch.linkage(pairwise_distances, method='complete')
    cluster_distance_threshold = pairwise_distances.max()/2
    idx_to_cluster_array = sch.fcluster(linkage, cluster_distance_threshold, 
                                        criterion='distance')
    idx = np.argsort(idx_to_cluster_array)
    
    if not inplace:
        corr_array = corr_array.copy()
    
    if isinstance(corr_array, pd.DataFrame):
        return corr_array.iloc[idx, :].T.iloc[idx, :]

    return corr_array[idx, :][:, idx]

def cut_bins(x):
    if x <= 5:
        return(5.0)
    elif 5 < x <= 10:
        return 10.0
    elif 10 < x <= 15:
        return 15.0
    elif 15 < x <= 20:
        return 20.0
    elif 20 < x <= 25:
        return 25.0
    elif 25 < x <= 30:
        return 30.0
    elif 30 < x <= 40:
        return 40.0
    elif 40 < x <= 50:
        return 50.0
    elif 50 < x <= 60:
        return 60.0
    else:
        return 70.0


# Getting features as functions for sequence features
def reverse_complement(seq):
    """Returns the reverse complement of the sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    try:
        rc = "".join(complement.get(base, base) for base in reversed(seq))
    except:
        print(seq)
    return rc

def complement(seq):
    """Returns the reverse complement of the sequence."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    try:
        c = "".join(complement.get(base, base) for base in seq)
    except:
        print(seq)
    return c

def get_length(x):
    """ Calculate length."""
    return len(x)

def get_smallerI(x, i):
    """Return true if string x is smaller or equal to i. """
    if len(x) <= i:
        return True
    else:
        return False

def get_higherI(x, i):
    """Return true if string x is smaller or equal to i. """
    if len(x) > i:
        return True
    else:
        return False

def get_countN(x,n):
    """Count the number of nucleotide n in the string."""
    return x.upper().count(n.upper())

def get_Ncomp(x, n):
    if n == 'GC':
        return (get_countN(x,'G') + get_countN(x,'C')) / len(x)
    else:
        return get_countN(x,n) / len(x)

def get_Nrun(x,n):
    """Look for 4 or more consecutive occurences of n in x."""
    my_regex = r"(?i)" + n + "+" + n + n + n
    if bool(re.search(my_regex, x)) == True:
        return True
    else:
        return False

def get_Nrun_max(x,n):
    """Find all consecutive occurences of n and maximum length"""
    my_regex = r"(?i)" + n + n + "+"
    try:
        count = max([len(i) for i in re.findall(my_regex, x)])
        if count > 6:
            return 6
        else:
            return count
    except:
        return 0

def pair_bases_pos(s1, s2, pos):
    """Pairing the n nucleotide after the nicking site with the first nucleotide of the insert"""    
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    s1 = s1 + "T" # spacer + start of scaffold
    base1 = s1[(-3+pos)].upper()
    base2 = s2[pos].upper()
    
    if base2 == complement[base1]:
        return True
    else:
        return False

def score_alignment(x, y):
    """ Use Biopython to score pairwise alignment"""
    # Set up alignment properties
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.2
    # Retrieve score
    score = aligner.align(x.upper(), y.upper()).score
    return score

def length_maxmicrohomology(seq1, seq2):
    t = [[0]*(1+len(seq2)) for i in range(1+len(seq1))]
    l, xl = 0, 0
    
    for x in range(1,1+len(seq1)):
        for y in range(1,1+len(seq2)):
            if seq1[x-1] == seq2[y-1]:
                t[x][y] = t[x-1][y-1] + 1
                if t[x][y]>l:
                    l = t[x][y]
                    xl  = x
            else:
                t[x][y] = 0
    hom = seq1[xl-l: xl]       
    return len(hom)

# Getting features as functions for pairing features
def pairing_bases(seq, spacer, pbs):
    if len(seq) <= 3:
        x = spacer + "GTN&" + reverse_complement(pbs + seq + 'NNN')
    else:
        x = spacer + "GTN&" + reverse_complement(pbs + seq[:3] + 'NNN')
    brackets = str(RNA.cofold(x)[0])
    count = brackets.count("(", 17, 20)
    return count

def scaffold_intact(seq):
    l1 = False
    l2 = False
    scaffold = 'gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc'
    x = scaffold + 'NNN&' + reverse_complement(seq)
    brackets = RNA.cofold(x)[0]
    if brackets[0:30] == "(((((((.((((....))))...)))))))":
        l1 = True
    if brackets[48:76] == "((((....)))).((((((...))))))":
        l2 = True
    return [l1, l2]

def loop1_intact(seq):
    scaffold = 'gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc'
    x = scaffold + 'NNN&' + reverse_complement(seq)
    brackets = RNA.cofold(x)[0]
    if brackets[0:30] == "(((((((.((((....))))...)))))))":
        return True
    else:
        return False

# Getting features as functions for pairing features

def get_tm(x):
    """ Calculates the melting temperature of x based on Biopython's TM_NN."""
    return mt.Tm_NN(x)

def get_vf(x):
    """ Calculates secondary structure approximated by Vienna Fold Energy. """
    vf = RNA.fold(x)
    if vf[1] == np.nan:
        return 0
    else:
        return vf[1]
        
# Generate features
def enhance_feature_sequence(df, seq='insert', rtt='RTT', pbs='PBS', guide = 'spacer'):
    """Calculates relevant sequence features based on insert sequence, RTT, PBS and MMR status.""" 
    # Generate sequences
    df[seq] = df[seq].astype('str')
    df[rtt] = df[rtt].astype('str')
    df[pbs] = df[pbs].astype('str')
    if "insertRTT" not in df.columns:
        df['insertRTT'] = df[seq] + df[rtt]
    # length features
    df['length'] = df[seq].apply(get_length)
    df['length_insertRTT'] = df['insertRTT'].apply(get_length)
    # Nucleotide composition
    df['percA'] = df[seq].apply(lambda x: get_Ncomp(x, 'A'))
    df['percC'] = df[seq].apply(lambda x: get_Ncomp(x, 'C'))
    df['percT'] = df[seq].apply(lambda x: get_Ncomp(x, 'T'))
    df['percG'] = df[seq].apply(lambda x: get_Ncomp(x, 'G'))
    df['percGC'] = df[seq].apply(lambda x: get_Ncomp(x, 'GC'))
    # Find runs
    df['Arun_maxlen'] = df[seq].apply(lambda x: get_Nrun_max(x, 'A'))
    return df

def enhance_feature_pairing(df, seq='insert', rtt='RTT', pbs='PBS', guide = 'spacer'):
    df['pairedbases'] = df.apply(lambda x: pairing_bases(x[seq], x[guide], x[pbs]), axis = 1)
    df['pos1compl'] = df.apply(lambda x: pair_bases_pos(x[guide], x['insertRTT'], 0), axis=1)
    df['loop1_intact'] = df[seq].parallel_apply(loop1_intact)
    return df

def enhance_feature_structure(df):
    df['VF_RTT'] = df['insertRTT'].parallel_apply(get_vf)
    df['VF_insert'] = df['sequence_original'].parallel_apply(get_vf)
    return df

def enhance_feature_structure_z(df1, df2, normcol, colname, on=['length'], how='left', mean_name = 'mean', std_name = 'std'):
    df1 = df1.merge(df2[[mean_name, std_name] + on], on=on, how=how)
    df1[colname] = (df1[normcol] - df1[mean_name]) / df1[std_name]
    df1[colname] = df1[colname].apply(lambda x: 0 if math.isnan(x) else x)
    df1[colname] = df1[colname].apply(lambda x: 0 if math.isinf(x) else x)
    df1 = df1.drop([mean_name, std_name], axis=1)
    return df1

##### models
def call_model(name, x = None, y = None, do_fit = False):

    if name == 'L1':
        model = l1_model(x,y, do_fit = do_fit)
    
    elif name == 'L2':
        model = l2_model(x,y, do_fit = do_fit)
        
    elif name == 'RandomForest':
        model = random_forest(x,y, do_fit = do_fit)
    
    elif name == 'XGBoost':
        model = xgb(x, y, do_fit = do_fit)
    
    elif name == 'MLP':
        model = mlp(x,y, do_fit = do_fit)

    return model

def l1_model(x,y, do_fit=False):
    model = linear_model.Lasso(alpha=0.01)
    if do_fit:
        model.fit(x, y)
    return model

def l2_model(x,y, do_fit=False):
    model = linear_model.Ridge(alpha=0.01)
    if do_fit:
        model.fit(x, y)
    return model

def mlp(x,y, do_fit = False):
    model = MLPRegressor(
        alpha = 1,
        hidden_layer_sizes = (1000, 100))
    
    if do_fit:
        model.fit(x, y)
        
    return model

def xgb(x,y, do_fit = False):
    model = xgboost.XGBRegressor(
        n_estimators=100,
        max_depth=4,
        learning_rate=0.1,
        reg_alpha= 0.00001,
        reg_lambda= 0.1,
        colsample_bytree=1,
        gamma=0.1,
        objective='reg:squarederror',
    )
    
    if do_fit:
        model.fit(x, y)
    
    return model

def random_forest(x,y, do_fit=False):
    model = RandomForestRegressor(
        max_depth=None, 
        random_state=0,
        n_estimators=1000,
        min_samples_leaf= 5,
    )
    if do_fit:
        model.fit(x, y) 
    return model

def DNA(length):
    return ''.join(random.choice('CGTA') for _ in range(length))

# Calculate baseline VF values for inserts that are not at the nicking site
def get_VFmeanstd_splitRTT(df, id_length='length', id_rtt='RTT', rttmap=None):
    VF_baseline = {}
    loi = list(set(df[id_length]))

    #  Generate 1000 random sequences for each length of interest
    for l in loi:
        seqtemp = []
        VF_baseline[l] = {}
        for k in range(1000):
            seqtemp.append(DNA(l))

        # Calculate VF values for all those sequences per target site and append to list
        for t in df[id_rtt].unique():
            if "-" in t:
                # Split
                tsplit = t.split("-")
                VFtemp = Parallel(n_jobs=8)(delayed(get_vf)(tsplit[0] + j+ tsplit[1])
                                            for j in seqtemp)
            else:                
                # For each RTT, create variableseq + RTT and calculate its VF value
                VFtemp = Parallel(n_jobs=8)(delayed(get_vf)(j + t)
                                            for j in seqtemp)
            # Calculate mean and std of the VF values
            mean = np.mean(VFtemp)
            std = np.std(VFtemp)
            # append to dictionary
            VF_baseline[l][t] = [mean, std]
    # Create dataframe
    baseline_df = pd.melt(pd.DataFrame.from_dict(VF_baseline, orient='index').reset_index(),id_vars=["index"]).rename(columns={'index': 'length'})
    baseline_df[['mean', 'std']] = baseline_df.value.values.tolist()
    baseline_df['target'] = baseline_df.variable.map(rttmap)
    return baseline_df

# Calculate baseline VF values
def get_VFmeanstd(df, id_length='length', id_rtt='RTT', rttmap=None):
    VF_baseline = {}
    loi = list(set(df[id_length]))

    #  Generate 1000 random sequences for each length of interest
    for l in loi:
        seqtemp = []
        VF_baseline[l] = {}
        for k in range(1000):
            seqtemp.append(DNA(l))

        # Calculate VF values for all those sequences per target site and append to list
        for t in df[id_rtt].unique():
            # For each RTT, create variableseq + RTT and calculate its VF value
            VFtemp = Parallel(n_jobs=8)(delayed(get_vf)(j + t)
                                        for j in seqtemp)
            # Calculate mean and std of the VF values
            mean = np.mean(VFtemp)
            std = np.std(VFtemp)
            # append to dictionary
            VF_baseline[l][t] = [mean, std]
    # Create dataframe
    baseline_df = pd.melt(pd.DataFrame.from_dict(VF_baseline, orient='index').reset_index(),id_vars=["index"]).rename(columns={'index': 'length'})
    baseline_df[['mean', 'std']] = baseline_df.value.values.tolist()
    baseline_df['target'] = baseline_df.variable.map(rttmap)
    return baseline_df

def zscore(rate, mean, std):
    """Calculates the Z-score from the mean and std."""
    zscore = (rate - mean) / std
    return zscore

def scale_zscore(zscore, mean, std):
    """Calculates the predicited insertion efficiency from the Z-score."""
    zscaled = zscore * std + mean
    return zscaled

#############################################################
###########             DATA IMPORT               ########### 
#############################################################

# Import screening data
measurements = pd.read_csv(screendata)
logging.info(f'Screen data size: {measurements.shape}')

# Import reference files for secondary structure
VF_baseline_length = pd.read_csv(VF_baseline_data_insert, sep = '\t').rename(columns = {'mean':'VF_mean_random','std':'VF_std_random'})  # random sequences across many lengths
VF_baseline_RTT = pd.read_csv(VF_baseline_data_RTT, sep = '\t').rename(columns = {'mean':'VF_mean_insertHA','std':'VF_std_insertHA'}) # random sequences + specific for each homology arm
Tm_baseline_length = pd.read_csv(Tm_baseline_data_insert, sep = '\t').rename(columns = {'mean':'Tm_mean_random','std':'Tm_std_random'})  # random sequences across many lengths
Tm_baseline_RTT = pd.read_csv(Tm_baseline_data_RTT, sep = '\t').rename(columns = {'mean':'Tm_mean_insertHA','std':'Tm_std_insertHA'})  # random sequences across many lengths
logging.info('Baseline structural data imported.')

#############################################################
###########            ADD FEATURES               ########### 
#############################################################

# # # Map experimental fixed factors based on target site
pegRNAmaps = pd.read_csv(pegRNA_mapping).set_index('target').to_dict()
logging.debug(f'pegRNA mapping imported: {pegRNAmaps}')
measurements['PBS'] = measurements.target.map(pegRNAmaps['PBS'])
measurements['RTT'] = measurements.target.map(pegRNAmaps['RTT'])
measurements['spacer'] = measurements.target.map(pegRNAmaps['spacer'])
measurements['mmr'] = np.where(measurements.experiment.str.contains("_HAP1_"), True, False)

# for each target site and cell line, we'll calculate a z-score based on the mean and the standard deviation. Those are screen specific features.
scaling_factor_z = {}
for i in list(set(measurements.axis_name)):
    mean_site = measurements[measurements.axis_name == i]['percIns'].mean()
    std_site = measurements[measurements.axis_name == i]['percIns'].std()
    scaling_factor_z[i] = [mean_site, std_site]
measurements['percIns_z'] = measurements.apply(lambda x: zscore(x['percIns'], scaling_factor_z[x['axis_name']][0], scaling_factor_z[x['axis_name']][1]), axis = 1)

# Add insert sequence specific features # this takes a tiny bit
logging.info('Calculating features for the screen data')
measurements = enhance_feature_sequence(measurements, seq = 'sequence_original')
measurements = enhance_feature_exploration(measurements, seq = 'sequence_original')
measurements['bin'] = measurements['length'].apply(cut_bins)
measurements = enhance_feature_pairing(measurements, seq = 'sequence_original')
measurements = enhance_feature_structure(measurements, seq='sequence_original')
logging.info('Calculating features for the screen data - done')

# Z-score normalization of secondary structure
measurements = enhance_feature_structure_z(measurements, VF_baseline_length, 'VF_insert', 'VF_insert_z', on=['length'], mean_name = 'VF_mean_random', std_name = 'VF_std_random')
measurements = enhance_feature_structure_z(measurements, Tm_baseline_length, 'Tm_insert', 'Tm_insert_z', on=['length'], mean_name = 'Tm_mean_random', std_name = 'Tm_std_random')
measurements = enhance_feature_structure_z(measurements, VF_baseline_RTT, 'VF_RTT', 'VF_RTT_z', on=['length', 'target'], mean_name = 'VF_mean_insertHA', std_name = 'VF_std_insertHA')
measurements = enhance_feature_structure_z(measurements, Tm_baseline_RTT, 'Tm_RTT', 'Tm_RTT_z', on=['length', 'target'], mean_name = 'Tm_mean_insertHA', std_name = 'Tm_std_insertHA')

# One hot encoding of categorical features
measurements = pd.get_dummies(measurements, columns = ['NN', 'N1', 'N2', 'N3', 'N4'], drop_first = False)

# # Save the featurized data
measurements.to_csv(data_featurized, index=False)

#############################################################
###########    Set up train and test data sets    ########### 
#############################################################

# Read in measurements data
# measurements = pd.read_csv(data_featurized)

# Split Test and Train data
traindata, testdata = splitdata(measurements, mainexperiment = 'HEK3_293T_PE2_1', fx = 'sequence_original', fy1 = 'percIns_z', fy2 = 'insbins', seed = setseed,  test_size = 0.3)
traindata.to_csv(datapath_train, index=False)
testdata.to_csv(datapath_test, index=False)

#############################################################
###########               Modelling               ########### 
#############################################################

# Read in train and test data, if continuing from here
# traindata = pd.read_csv(datapath_train)
# testdata = pd.read_csv(datapath_test)

# Splits for kfold splitting for modeleing
Xseq = list(set(traindata['sequence_original']))
split_num = 10
kf = sklearn.model_selection.KFold(n_splits=split_num, random_state=setseed, shuffle=True)

### Feature sets
feat_all = ['length', 'VF_RTT_z',  'percC', 'pairedbases','mmr', 'Arun_maxlen','percT', 'percA', 'percG', 'loop1_intact', 'loop2_intact', 'pos1compl', 
    'percGC', 'smaller4', 'smaller7', 'smaller13', 'longer40', 'N1_A', 'N2_A', 'N3_A', 'N4_A', 'NN_A', 'N1_T', 'N2_T', 'N3_T', 'N4_T', 'NN_T',
    'N1_G', 'N2_G', 'N3_G', 'N4_G', 'NN_G', 'N1_C', 'N2_C', 'N3_C', 'N4_C', 'NN_C', 'Arun', 'Crun', 'Trun', 'Grun', 'Trun_maxlen', 'Crun_maxlen', 'Grun_maxlen',
    'nick2match', 'nick3match', 'align_InsSpc', 'maxmicrohom_InsHA', 'VF_insert', 'VF_RTT', 'VF_insert_z', 'Tm_insert', 'Tm_insert_z']
feat_core = ['length', 'VF_RTT_z', 'mmr', 'percC', 'pairedbases', 'Arun_maxlen',  'percA', 'percT','pos1compl', 'loop1_intact']

### Comparing different models with 10-fold cross-validation --> XGBoost
logging.info(f'Model comaprison is running...')
cross_scores = {}
model_list = ['L1', 'L2','MLP','RandomForest', 'XGBoost'] 
featset = feat_core
# Calculate correlation scores for all models
for m in model_list:
    rscores = []
    # We split the unique sequences into kfolds and then assign their corresponding data
    for train_index, test_index in kf.split(Xseq): 
        # list of test and train sequences
        trainseq = [Xseq[i] for i in train_index]
        testseq = [Xseq[i] for i in test_index]
        # assign data points to test and train
        X_train, X_test = traindata[traindata['sequence_original'].isin(trainseq)][featset], traindata[traindata['sequence_original'].isin(testseq)][featset] 
        y_train, y_test = traindata[traindata['sequence_original'].isin(trainseq)]['percIns_z'], traindata[traindata['sequence_original'].isin(testseq)]['percIns_z']
        # train and test the model
        model = call_model(m, x = X_train, y = y_train, do_fit = True)
        pred = model.predict(X_test)
        corr, _ = stats.pearsonr(pred, y_test.values.ravel())
        rscores.append(corr)
    cross_scores[m]= rscores
pd.DataFrame.from_dict(cross_scores).stack().reset_index().rename(columns = {'level_0': 'model_id', 'level_1': 'model', 0: 'score'}).to_csv(model_crossvalidation_path, index=False)
logging.info(f'Model architectures were computed and output is stored in {model_crossvalidation_path}')

### Feature selection: add sets of features based on discoveries about sequence features, system features
featset_d = [
    ['length'], # length
    ['VF_RTT_z','pairedbases','pos1compl',  'loop1_intact'], # structural features
    ['percC', 'percA', 'percT'], # nucleotide composition
    ['mmr','Arun_maxlen'], # system features (= dropped sequence features)
    ['length', 'percC', 'percA', 'percT', 'VF_RTT_z', 'pairedbases','pos1compl',  'loop1_intact'], # sequence features (= dropped system features)
    ['length', 'VF_RTT_z', 'mmr', 'percC', 'pairedbases', 'Arun_maxlen',  'percA', 'percT','pos1compl', 'loop1_intact'], # all core features
    ['length', 'VF_RTT_z',  'percC', 'pairedbases','mmr', 'Arun_maxlen','percT', 'percA', 'percG', 'loop1_intact', 'loop2_intact', 'pos1compl', 
    'percGC', 'smaller4', 'smaller7', 'smaller13', 'longer40', 'N1_A', 'N2_A', 'N3_A', 'N4_A', 'NN_A', 'N1_T', 'N2_T', 'N3_T', 'N4_T', 'NN_T',
    'N1_G', 'N2_G', 'N3_G', 'N4_G', 'NN_G', 'N1_C', 'N2_C', 'N3_C', 'N4_C', 'NN_C', 'Arun', 'Crun', 'Trun', 'Grun', 'Trun_maxlen', 'Crun_maxlen', 'Grun_maxlen',
    'nick2match', 'nick3match', 'align_InsSpc', 'maxmicrohom_InsHA',  'Tm_insert', 'Tm_insert_z'], # all features, but 'VF_RTT', to avoid that the model can guess the target site directly
]
featnames = ['Length', 'Structure', '%N', 'System', 'Sequence','Model','all']
scores = {k: [] for k in featnames}
for train_index, test_index in kf.split(Xseq):
    # list of test and train sequences
    trainseq = [Xseq[i] for i in train_index]
    testseq = [Xseq[i] for i in test_index]
    for f in range(len(featnames)):
        feat_temp = featset_d[f]
        # assign data points to test and train
        X_train, X_test = traindata[traindata['sequence_original'].isin(trainseq)][feat_temp], traindata[traindata['sequence_original'].isin(testseq)][feat_temp]
        y_train, y_test = traindata[traindata['sequence_original'].isin(trainseq)]['percIns_z'], traindata[traindata['sequence_original'].isin(testseq)]['percIns_z']
        # train and test the model
        model = call_model('XGBoost', x=X_train, y=y_train, do_fit=True)
        pred = model.predict(X_test)
        corr, _ = stats.pearsonr(pred, y_test.values.ravel())
        scores[featnames[f]].append(corr)
cross_scores_R = pd.DataFrame.from_dict(scores, orient = 'index').stack().reset_index().rename(columns = {'level_0': 'feature', 'level_1': 'feature_id', 0: 'score'})
cross_scores_R.to_csv(model_rmvfeatures_path, index=False)
logging.info(f'Evaluating different feature sets is finished')

#Final model for future applications
featset = feat_core
X_train = traindata[featset]
y_train = traindata['percIns_z']
model = call_model('XGBoost', x = X_train, y = y_train, do_fit = True)

## Save model
pickle.dump(model, open(model_path, 'wb'))
logging.info(f'Model is stored in {model_path}')

### Test model on held out dataset
testdata['percIns_z_pred'] = model.predict(testdata[featset])
corr, _ = stats.pearsonr(testdata['percIns_z_pred'], testdata['percIns_z'].values.reshape(-1))
logging.info(f'Pearson R = {corr:.3f}')
pd.DataFrame(testdata.to_csv(model_performance_path, index = False))
logging.info(f'Performance evaluation on testdata is finished and output is stored in {model_performance_path}')


#############################################################
###########        Model characterization         ########### 
#############################################################

#Load model
# model = pickle.load(open(model_path, 'rb'))
featset = feat_core

### Get and plot shap values
shap_values = shap.TreeExplainer(model).shap_values(traindata[featset])
shap.summary_plot(shap_values, traindata[featset], featset, show = False,  plot_size=(5,5))
plt.savefig(model_shap_impact_path, dpi = 300, transparent=True, bbox_inches='tight')
plt.close()
shap.summary_plot(shap_values, features=traindata[featset], feature_names=featset, plot_type='bar',show = False, plot_size=(5,5))   #plot_size=(5,5)
plt.savefig(model_shap_importance_path, dpi = 300, transparent=True, bbox_inches='tight')
plt.close()

#############################################################
###########          Model application            ########### 
#############################################################

## Use for prediction on new target sites: Canadian goose screen
logging.info('Analyzing canadian goose screen data for 6 target sites in HEK3')
# prepare data
data_canadian = pd.read_csv(canadian_path)
scaling_factor_z = {}
for i in list(set(data_canadian.axis_name)):
    mean_site = data_canadian[data_canadian.axis_name == i]['percIns'].mean()
    std_site = data_canadian[data_canadian.axis_name == i]['percIns'].std()
    scaling_factor_z[i] = [mean_site, std_site]
data_canadian['percIns_z'] = data_canadian.apply(lambda x: zscore(x['percIns'], scaling_factor_z[x['axis_name']][0], scaling_factor_z[x['axis_name']][1]), axis = 1)
data_canadian = enhance_feature_sequence(data_canadian, seq = 'sequence_original')
data_canadian = enhance_feature_pairing(data_canadian, seq = 'sequence_original')
data_canadian = enhance_feature_structure(data_canadian)
VF_baseline_canadian = pd.read_csv(VF_baseline_canadian_data_RTT)
data_canadian = enhance_feature_structure_z(data_canadian, VF_baseline_canadian, 'VF_RTT', 'VF_RTT_z', on=['length', 'target'], mean_name = 'mean', std_name = 'std')
data_canadian.to_csv(data_canadian_featurized, index=False)
# Load data without going through feature generation
# data_canadian = pd.read_csv(data_canadian_featurized)
data_canadian['percIns_z_pred'] = model.predict(data_canadian[featset])
corr, _ = stats.pearsonr(data_canadian['percIns_z_pred'], data_canadian['percIns_z'])
data_canadian.to_csv(model_performance_ca_path)
logging.info(f'Performance evaluation (R = {corr:.3f}) on canadian goose testdata is finished and output is stored in {model_performance_ca_path}')

### Use for prediction on new target sites and fusion proteins: Barnacle screen
logging.info('Analyzing barnacle goose screen data for tags in many target sites')
data_barnacle = pd.read_csv(barnacle_path)
data_barnacle = enhance_feature_sequence(data_barnacle, seq = 'sequence_original')
data_barnacle = enhance_feature_pairing(data_barnacle, seq = 'sequence_original')
data_barnacle = enhance_feature_structure(data_barnacle)
VF_baseline_barnacle = get_VFmeanstd_splitRTT(data_barnacle, id_length='length', id_rtt='RTT', rttmap={v: k for k, v in pegRNAmaps['RTT'].items()})
VF_baseline_barnacle.to_csv(VF_baseline_barnacle_data_RTT)
data_barnacle = enhance_feature_structure_z(data_barnacle, VF_baseline_barnacle, 'VF_RTT', 'VF_RTT_z', on=['length', 'target'], mean_name = 'mean', std_name = 'std')
data_barnacle.to_csv(data_barnacle_featurized, index=False)
# Load data without feature generation
data_barnacle = pd.read_csv(data_barnacle_featurized)
data_barnacle['percIns_z_pred'] = model.predict(data_barnacle[featset])
corr, _ = stats.pearsonr(data_barnacle['percIns_z_pred'], data_barnacle['percIns_z'])
# Categorize them into top and bottom performing sequences
cat_ba = data_barnacle.groupby(['target','tag'])['percIns_z_pred'].mean().to_dict()
meanpred_dict = {}
for t in set(data_barnacle.target):
    meanpred_dict[t] = {}
for k,v in cat_ba.items():
    meanpred_dict[k[0]][k[1]] = v
print(cat_ba)
data_barnacle['mean_pred'] = data_barnacle.apply(lambda x: meanpred_dict[x['target']][x['tag']], axis = 1)
data_barnacle['predgroup'] = data_barnacle.apply(lambda x: 'top' if (x['mean_pred'] < x['percIns_z_pred']) else 'bottom', axis=1)
data_barnacle.to_csv(model_performance_ba_path)
logging.info(f'Performance evaluation (R = {corr:.3f}) on barnacle goose testdata is finished and output is stored in {model_performance_ba_path}')


#############################################################
###########          Supplementals                ########### 
#############################################################

###########         Feature exploration

### Feature correlation
corrMatrix = cluster_corr(traindata[feat_all].corr())
corrMatrix.to_csv(path_corrmatrix_all, index=False)
logging.info(f'Feature correlation is finished and output is stored in {path_corrmatrix_all}')


# ###########        Model characterization

### Feature correlation
corrMatrix = cluster_corr(traindata[feat_core].corr())
corrMatrix.to_csv(path_corrmatrix_core, index=False)
logging.info(f'Feature correlation is finished and output is stored in {path_corrmatrix_core}')


### Normalized vs non-normalized VF value on canadian data
featset_d = [['length', 'VF_RTT_z', 'mmr', 'percC', 'pairedbases', 'Arun_maxlen',  'percA', 'percT','pos1compl', 'loop1_intact'], # feat_core
['length', 'VF_RTT', 'mmr', 'percC', 'pairedbases', 'Arun_maxlen',  'percA', 'percT','pos1compl', 'loop1_intact'], # feat_core with VF_RTT
]
# Train the model with VF and without VF on the main training data
X_train = traindata[featset_d[0]]
y_train = traindata['percIns_z']
model_vfz = call_model('XGBoost', x = X_train, y = y_train, do_fit = True)
X_train = traindata[featset_d[1]]
y_train = traindata['percIns_z']
model_novfz = call_model('XGBoost', x = X_train, y = y_train, do_fit = True)
# predict on canadian data with both models
data_barnacle['percIns_z_pred'] = model_vfz.predict(data_barnacle[featset_d[0]])
data_barnacle['percIns_z_pred_novfz'] = model_novfz.predict(data_barnacle[featset_d[1]])
data_barnacle.to_csv(model_vffeatures_path)

## MMR in HAP1
featset_d = [['length', 'VF_RTT_z', 'mmr', 'percC', 'pairedbases', 'Arun_maxlen',  'percA', 'percT','pos1compl', 'loop1_intact'], # feat_core
['length', 'VF_RTT_z', 'percC', 'pairedbases', 'Arun_maxlen',  'percA', 'percT','pos1compl', 'loop1_intact'], # feat_core without mmr
]
# Train the model with VF and without VF on the test dataset for hap1 sequences
X_train = traindata[featset_d[0]]
y_train = traindata['percIns_z']
model_mmr = call_model('XGBoost', x = X_train, y = y_train, do_fit = True)
X_train = traindata[featset_d[1]]
y_train = traindata['percIns_z']
model_nommr = call_model('XGBoost', x = X_train, y = y_train, do_fit = True)
# predict on canadian data with both models
mmr_data = testdata[testdata.cell_line == 'HAP1'].copy()
mmr_data['percIns_z_pred'] = model_mmr.predict(mmr_data[featset_d[0]])
mmr_data['percIns_z_pred_nommr'] = model_nommr.predict(mmr_data[featset_d[1]])
mmr_data.to_csv(model_mmrfeatures_path)

##########        Model application

### Use for sequence padding
logging.info('Analyzing data for short sequences that are padded to increase editing rates')
data_padding = pd.read_csv(padding_path)
data_padding = enhance_feature_sequence(data_padding, seq = 'sequence_original')
data_padding = enhance_feature_pairing(data_padding, seq = 'sequence_original')
data_padding = enhance_feature_structure(data_padding)
VF_baseline_padding = get_VFmeanstd(data_padding, id_length='length', id_rtt = 'RTT', rttmap={v: k for k, v in pegRNAmaps['RTT'].items()})
VF_baseline_padding.to_csv(VF_baseline_padding_data_RTT)
data_padding = enhance_feature_structure_z(data_padding, VF_baseline_padding, 'VF_RTT', 'VF_RTT_z', on=['length', 'target'], mean_name = 'mean', std_name = 'std')
data_padding.to_csv(data_padding_featurized, index=False)
# Load data without feature generation
# data_padding = pd.read_csv(data_padding_featurized)
data_padding['percIns_z_pred'] = model.predict(data_padding[featset])
corr, _ = stats.pearsonr(data_padding['percIns_z_pred'], data_padding['percIns_z'])
data_padding.to_csv(model_performance_pa_path)
logging.info(f'Performance evaluation (R = {corr:.3f}) on padding testdata is finished and output is stored in {model_performance_pa_path}')
