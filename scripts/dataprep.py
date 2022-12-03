import os
import sys
import logging
import math
import re
import RNA
import random
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from joblib import delayed, Parallel
from pandarallel import pandarallel
from collections import defaultdict
from datetime import datetime
from scipy import stats
from scipy.cluster import hierarchy as sch
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio import Align


pandarallel.initialize()

#############################################################
###########             INITIALIZE                ########### 
#############################################################

verbosity = 'info'
workdir = '../'
workdir_input = os.path.join(workdir, 'files/input/')
workdir_output = os.path.join(workdir, 'files/output/') # Change to your prefered output directory if you don't want to overwrite the files here
setseed = 0

# Input
screendata = os.path.join(workdir_input, 'screendata.csv')
pegRNA_mapping = os.path.join(workdir_input, 'Table_S2_pegRNAs.csv')

canadian_path = os.path.join(workdir_input,'data_canadian.csv')
barnacle_path = os.path.join(workdir_input,'data_barnacle.csv')
padding_path = os.path.join(workdir_input,'data_padding.csv')

path_N3 = os.path.join(workdir_input,'HEK293T-T1-eN3-EditScore-Table.csv')
path_N6 = os.path.join(workdir_input,'HEK293T-T1-eN6-EditScore-Table.csv')
path_deeppe = os.path.join(workdir_input,'lib2.csv')

VF_baseline_data_insert = os.path.join(workdir_input, 'VF_baseline_meanstd.tsv')
VF_baseline_data_RTT = os.path.join(workdir_input, 'VF_baseline_RTT_meanstd.tsv')
VF_baseline_data_PBSinsert = os.path.join(workdir_input, 'VF_baseline_PBSinsert_meanstd.csv')
VF_baseline_data_ext = os.path.join(workdir_input, 'VF_baseline_ext_meanstd_2.csv')
Tm_baseline_data_insert = os.path.join(workdir_input,'Tm_baseline_meanstd.tsv')
Tm_baseline_data_RTT = os.path.join(workdir_input, 'Tm_baseline_RTT_meanstd.tsv')
VF_baseline_data_deeppe = os.path.join(workdir_input, 'VF_baseline_deeppe_meanstd.tsv')

# Output
data_featurized = os.path.join(workdir_output, 'Data_features_onehot.csv')
data_canadian_featurized = os.path.join(workdir_output,'data_canadian_featurized.csv')
data_barnacle_featurized = os.path.join(workdir_output,'data_barnacle_featurized.csv')
data_padding_featurized = os.path.join(workdir_output,'data_padding_featurized.csv')
data_choi_featurized_eN3 = os.path.join(workdir_output,'data_choi_featurized_eN3.csv')
data_choi_featurized_eN6 = os.path.join(workdir_output,'data_choi_featurized_eN6.csv')
data_deeppe_featurized = os.path.join(workdir_output,'data_deeppe_featurized.csv')


# Set up logging based on the verbosity level set by the command line arguments:
# logging.basicConfig(format='%(levelname)s: %(message)s', level=verbosity.upper())
logging.getLogger('matplotlib').setLevel(logging.WARNING)

logging.basicConfig(filename= os.path.join(workdir_output,f'log_{datetime.now().strftime("%d%m%Y-%H%M%S")}.txt'),
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)

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

def enhance_feature_exploration(df, seq='insert', ha='HA', pbs='PBS', guide = 'spacer'):
    logging.debug("Adding explorative features")
    # sequences
    df['ext'] = df[pbs] + df[seq] + df[ha]
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
    logging.debug("Adding explorative alignment and structure features.")
    df['nick2match'] = df.apply(lambda x: pair_bases_pos(x[guide], x['RTT'], 1), axis=1)
    df['nick3match'] = df.apply(lambda x: pair_bases_pos(x[guide], x['RTT'], 2), axis=1)
    df['align_InsSpc'] = df.parallel_apply(lambda x: score_alignment(x[seq],x[guide]), axis=1)
    df['maxmicrohom_RTT'] = df.apply(lambda x: length_maxmicrohomology(x[seq], x[ha]),axis=1)
    df['loops_intact'] = df[seq].parallel_apply(scaffold_intact)
    df['loop1_intact'] = df['loops_intact'].apply(lambda x: x[0])
    df['loop2_intact'] = df['loops_intact'].apply(lambda x: x[1])
    # Structure
    df['VF_insert'] = df[seq].parallel_apply(get_vf)
    df['Tm_insert'] = df[seq].apply(get_tm)
    df['VF_ext'] = df['ext'].parallel_apply(get_vf)
    df['VF_PBSinsert'] = df['PBSinsert'].parallel_apply(get_vf)
    df['Tm_RTT'] = df['RTT'].apply(get_tm)
    return df


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

def split_RTT(row, ha, ins):
    if "-" in row[ha]:
        tsplit = row['ha'].split("-")
        rtt = tsplit[0] + row[ins] + tsplit[1]
    else:
        rtt = row[ins] + row[ha]
    return(rtt)

# Generate features
def enhance_feature_sequence(df, seq='insertion', ha='HA', pbs='PBS', guide = 'spacer', split = False):
    """Calculates relevant sequence features based on insert sequence, ha, PBS and MMR status.""" 
    # Generate sequences
    df[seq] = df[seq].astype('str')
    df[ha] = df[ha].astype('str')
    df[pbs] = df[pbs].astype('str')
    if "RTT" not in df.columns:
        if split == False:
            df['RTT'] = df[seq] + df[ha]
        elif split == True:
            df['RTT'] = df.apply(split_RTT(x, ha, seq), axis = 1)
        else:
            error("this doesn't work.")
    # length features
    df['length'] = df[seq].apply(get_length)
    df['length_RTT'] = df['RTT'].apply(get_length)
    # Nucleotide composition
    df['percA'] = df[seq].apply(lambda x: get_Ncomp(x, 'A'))
    df['percC'] = df[seq].apply(lambda x: get_Ncomp(x, 'C'))
    df['percT'] = df[seq].apply(lambda x: get_Ncomp(x, 'T'))
    df['percG'] = df[seq].apply(lambda x: get_Ncomp(x, 'G'))
    df['percGC'] = df[seq].apply(lambda x: get_Ncomp(x, 'GC'))
    # Find runs
    df['Arun_maxlen'] = df[seq].apply(lambda x: get_Nrun_max(x, 'A'))
    return df

def enhance_feature_pairing(df, seq='insert', ha='HA', pbs='PBS', guide = 'spacer'):
    df['pairedbases'] = df.apply(lambda x: pairing_bases(x[seq], x[guide], x[pbs]), axis = 1)
    df['pos1compl'] = df.apply(lambda x: pair_bases_pos(x[guide], x['RTT'], 0), axis=1)
    df['loop1_intact'] = df[seq].parallel_apply(loop1_intact)
    return df

def enhance_feature_structure(df):
    df['VF_RTT'] = df['RTT'].parallel_apply(get_vf)
    df['VF_insert'] = df['insertion'].parallel_apply(get_vf)
    return df

def enhance_feature_structure_z(df1, df2, normcol, colname, on=['length'], how='left', mean_name = 'mean', std_name = 'std'):
    df1 = df1.merge(df2[[mean_name, std_name] + on], on=on, how=how)
    df1[colname] = (df1[normcol] - df1[mean_name]) / df1[std_name]
    df1[colname] = df1[colname].apply(lambda x: 0 if math.isnan(x) else x)
    df1[colname] = df1[colname].apply(lambda x: 0 if math.isinf(x) else x)
    df1 = df1.drop([mean_name, std_name], axis=1)
    return df1

def DNA(length):
    return ''.join(random.choice('CGTA') for _ in range(length))

# Calculate baseline VF values for inserts that are not at the nicking site
def get_VFmeanstd_splitRTT(df, id_length='length', id_part='HA', rttmap=None, part = 'RTT'):
    VF_baseline = {}
    loi = list(set(df[id_length]))

    #  Generate 1000 random sequences for each length of interest
    for l in loi:
        seqtemp = []
        VF_baseline[l] = {}
        for k in range(1000):
            seqtemp.append(DNA(l))

        # Calculate VF values for all those sequences per target site and append to list
        if part == 'RTT':
            rttmap = {v: k for k, v in pegRNAmaps[id_part].items()}
            for t in df[id_part].unique():
                if "-" in t:
                    # Split
                    tsplit = t.split("-")
                    VFtemp = Parallel(n_jobs=8)(delayed(get_vf)(tsplit[0] + j+ tsplit[1]) for j in seqtemp)
                else:                
                    # For each RTT, create variableseq + HA and calculate its VF value
                    VFtemp = Parallel(n_jobs=8)(delayed(get_vf)(j + t) for j in seqtemp)
                # Calculate mean and std of the VF values
                mean = np.mean(VFtemp)
                std = np.std(VFtemp)
                # append to dictionary
                VF_baseline[l][t] = [mean, std]
        elif part == 'PBSinsert':
            rttmap = {v: k for k, v in pegRNAmaps[id_part].items()}
            for t in df[id_part].unique():        
                # For each PBS, create PBS + variableseq and calculate its VF value
                VFtemp = Parallel(n_jobs=8)(delayed(get_vf)(t + j) for j in seqtemp)
                # Calculate mean and std of the VF values
                mean = np.mean(VFtemp)
                std = np.std(VFtemp)
                # append to dictionary
                VF_baseline[l][t] = [mean, std]
        elif part == 'ext':
            rttmap = {v: k for k, v in pegRNAmaps[id_part].items()}
            for t in df[id_part].unique():
                # get matching PBS (Part must be RTT)
                target = rttmap[t]
                pbs = pegRNAmaps['PBS'][target]
                if "-" in t:
                    # For each RTT, create PBS + nicking site + variableseq + rest of HA and calculate its VF value
                    tsplit = t.split("-")
                    VFtemp = Parallel(n_jobs=8)(delayed(get_vf)(pbs + tsplit[0] + j+ tsplit[1]) for j in seqtemp)
                else:                
                    # For each RTT, create PBS + variableseq + HA and calculate its VF value
                    VFtemp = Parallel(n_jobs=8)(delayed(get_vf)(pbs + j + t) for j in seqtemp)
                # Calculate mean and std of the VF values
                mean = np.mean(VFtemp)
                std = np.std(VFtemp)
                # append to dictionary
                VF_baseline[l][t] = [mean, std]
        else:
            print("Part not specified correctly.")

    # Create dataframe
    baseline_df = pd.melt(pd.DataFrame.from_dict(VF_baseline, orient='index').reset_index(),id_vars=["index"]).rename(columns={'index': 'length'})
    baseline_df[['mean', 'std']] = baseline_df.value.values.tolist()
    baseline_df['target'] = baseline_df.variable.map(rttmap)
    return baseline_df

# Calculate baseline VF values
def get_VFmeanstd(df, id_length='length', id_rtt='HA', rttmap=None):
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

def get_VFmeanstd_deeppe(df, id_length='length', id_rtt='HA'):
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
            if std == 0:
                std = 0.001 # 0.001
            # append to dictionary
            VF_baseline[l][t] = [mean, std]
    # Create dataframe
    baseline_df = pd.melt(pd.DataFrame.from_dict(VF_baseline, orient='index').reset_index(),id_vars=["index"]).rename(columns={'index': 'length'})
    baseline_df[['VF_RTT_mean', 'VF_RTT_std']] = baseline_df.value.values.tolist()
    baseline_df['HA'] = baseline_df.variable
    return baseline_df

def get_VF_baseline_PBSinsert(df):
    VF_baseline = {}
    loi = list(set(df['length_PBSinsert']))
    #  Generate 1000 random sequences for each length of interest
    for l in loi:
        seqtemp = []
        VF_baseline[l] = {}
        for k in range(1000):
            seqtemp.append(DNA(l))
        # Calculate VF values for all those sequences per target site and append to list
        for t in df['PBS'].unique():
            # For each PBS, create PBS + variableseq and calculate its VF value
            VFtemp = Parallel(n_jobs=8)(delayed(get_vf)(t + j) for j in seqtemp)
            # Calculate mean and std of the VF values
            mean = np.mean(VFtemp)
            std = np.std(VFtemp)
            # append to dictionary
            VF_baseline[l][t] = [mean, std]
    # Create dataframe
    baseline_df = pd.melt(pd.DataFrame.from_dict(VF_baseline, orient='index').reset_index(),id_vars=["index"]).rename(columns={'index': 'length_PBSinsert','variable':'PBS'})
    baseline_df[['VF_mean_PBSinsert', 'VF_std_PBSinsert']] = baseline_df.value.values.tolist()
    baseline_df.to_csv(VF_baseline_data_PBSinsert, index = False)
    return baseline_df

def get_VF_baseline_ext(df, mapping):
    VF_baseline = {}
    loi = list(set(df['length_ext']))
    #  Generate 1000 random sequences for each length of interest
    for l in loi:
        seqtemp = []
        VF_baseline[l] = {}
        for k in range(1000):
            seqtemp.append(DNA(l))
        # Calculate VF values for all those sequences per target site and append to list
        for t in df['target'].unique():
            ha = mapping['HA'][t]
            pbs = mapping['PBS'][t]
            # For each PBS, create PBS + variableseq + RTT and calculate its VF value
            VFtemp = Parallel(n_jobs=8)(delayed(get_vf)(pbs + j + ha) for j in seqtemp)
            # Calculate mean and std of the VF values
            mean = np.mean(VFtemp)
            std = np.std(VFtemp)
            # append to dictionary
            VF_baseline[l][t] = [mean, std]
    # Create dataframe
    baseline_df = pd.melt(pd.DataFrame.from_dict(VF_baseline, orient='index').reset_index(),id_vars=["index"]).rename(columns={'index': 'length_ext','variable':'target'})
    baseline_df[['VF_mean_ext', 'VF_std_ext']] = baseline_df.value.values.tolist()
    baseline_df.to_csv(VF_baseline_data_ext, index = False)
    return baseline_df

def zscore(rate, mean, std):
    """Calculates the Z-score from the mean and std."""
    zscore = (rate - mean) / std
    return zscore

def scale_zscore(zscore, mean, std):
    """Calculates the predicited insertion efficiency from the Z-score."""
    zscaled = zscore * std + mean
    return zscaled

def get_splitext(x):
    if "-" in x['HA']:
        tsplit = x['HA'].split("-")
        return(tsplit[0] + x['insertion'] + tsplit[1])
    else:            
        return(x['insertion'] + x['HA']) 

def str_contains(subseq, seq):
    if subseq in seq:
        return True
    else:
        return False


#############################################################
###########             Main data               ########### 
#############################################################

# Import screening data
pegRNAmaps = pd.read_csv(pegRNA_mapping).set_index('target').to_dict()
measurements = pd.read_csv(screendata)

measurements['HA'] = measurements.target.map(pegRNAmaps['HA'])
measurements['PBS'] = measurements.target.map(pegRNAmaps['PBS'])
measurements['spacer'] = measurements.target.map(pegRNAmaps['spacer'])
logging.info(f'Screen data size: {measurements.shape}')

# for each target site and cell line, we'll calculate a z-score based on the mean and the standard deviation. Those are screen specific features.
scaling_factor_z = {}
for i in list(set(measurements.axis_name)):
    mean_site = measurements[measurements.axis_name == i]['percIns'].mean()
    std_site = measurements[measurements.axis_name == i]['percIns'].std()
    scaling_factor_z[i] = [mean_site, std_site]

measurements['percIns_z'] = measurements.apply(lambda x: zscore(x['percIns'], scaling_factor_z[x['axis_name']][0], scaling_factor_z[x['axis_name']][1]), axis = 1)
measurements['mmr'] = np.where(measurements.experiment.str.contains("_HAP1_"), True, False)

# Import reference files for secondary structure
VF_baseline_length = pd.read_csv(VF_baseline_data_insert, sep = '\t').rename(columns = {'mean':'VF_mean_random','std':'VF_std_random'})  # random sequences across many lengths
VF_baseline_RTT = pd.read_csv(VF_baseline_data_RTT, sep = '\t').rename(columns = {'mean':'VF_mean_RTT','std':'VF_std_RTT'}) # random sequences + specific for each homology arm
VF_baseline_PBSinsert = pd.read_csv(VF_baseline_data_PBSinsert) # random sequences + specific for each PBS
VF_baseline_ext = pd.read_csv(VF_baseline_data_ext) # random sequences + specific for each homology arm and PBS
Tm_baseline_length = pd.read_csv(Tm_baseline_data_insert, sep = '\t').rename(columns = {'mean':'Tm_mean_random','std':'Tm_std_random'})  # random sequences across many lengths
Tm_baseline_RTT = pd.read_csv(Tm_baseline_data_RTT, sep = '\t').rename(columns = {'mean':'Tm_mean_RTT','std':'Tm_std_RTT'})  # random sequences + specific for each homology arm
logging.info('Baseline structural data imported.')

# add more features
measurements = enhance_feature_sequence(measurements, seq = 'insertion')
measurements = enhance_feature_exploration(measurements, seq = 'insertion')
measurements = enhance_feature_pairing(measurements, seq = 'insertion')
measurements = enhance_feature_structure(measurements)
measurements['bin'] = measurements['length'].apply(cut_bins)

# Z-score normalization of secondary structure
measurements['length_PBSinsert'] = measurements['PBSinsert'].apply(len)
measurements['length_ext'] = measurements['ext'].apply(len)
measurements = enhance_feature_structure_z(measurements, VF_baseline_length, 'VF_insert', 'VF_insert_z', on=['length'], mean_name = 'VF_mean_random', std_name = 'VF_std_random')
measurements = enhance_feature_structure_z(measurements, Tm_baseline_length, 'Tm_insert', 'Tm_insert_z', on=['length'], mean_name = 'Tm_mean_random', std_name = 'Tm_std_random')
measurements = enhance_feature_structure_z(measurements, VF_baseline_RTT, 'VF_RTT', 'VF_RTT_z', on=['length', 'target'], mean_name = 'VF_mean_RTT', std_name = 'VF_std_RTT')
measurements = enhance_feature_structure_z(measurements, Tm_baseline_RTT, 'Tm_RTT', 'Tm_RTT_z', on=['length', 'target'], mean_name = 'Tm_mean_RTT', std_name = 'Tm_std_RTT')
measurements = enhance_feature_structure_z(measurements, VF_baseline_PBSinsert, 'VF_PBSinsert', 'VF_PBSinsert_z', on=['length_PBSinsert', 'PBS'], mean_name = 'VF_mean_PBSinsert', std_name = 'VF_std_PBSinsert')
measurements = enhance_feature_structure_z(measurements, VF_baseline_ext, 'VF_ext', 'VF_ext_z', on=['length_ext', 'target'], mean_name = 'VF_mean_ext', std_name = 'VF_std_ext')

# One hot encoding of categorical features
measurements = pd.get_dummies(measurements, columns = ['NN', 'N1', 'N2', 'N3', 'N4'], drop_first = False)

# Save the featurized data
measurements.to_csv(data_featurized, index=False)


#############################################################
###########          Model applications           ###########
#############################################################

# Use for prediction on new target sites: Canadian goose screen
logging.info('Analyzing canadian goose screen data for 6 target sites in HEK3')
data_canadian = pd.read_csv(canadian_path)
scaling_factor_z = {}
for i in list(set(data_canadian.axis_name)):
    mean_site = data_canadian[data_canadian.axis_name == i]['percIns'].mean()
    std_site = data_canadian[data_canadian.axis_name == i]['percIns'].std()
    scaling_factor_z[i] = [mean_site, std_site]
data_canadian['percIns_z'] = data_canadian.apply(lambda x: zscore(x['percIns'], scaling_factor_z[x['axis_name']][0], scaling_factor_z[x['axis_name']][1]), axis = 1)
data_canadian = enhance_feature_sequence(data_canadian, seq = 'insertion')
data_canadian = enhance_feature_pairing(data_canadian, seq = 'insertion')
data_canadian = enhance_feature_structure(data_canadian)
VF_baseline_canadian = get_VFmeanstd_splitRTT(data_canadian, id_length='length', id_part='HA', rttmap={v: k for k, v in pegRNAmaps['HA'].items()}, part = 'RTT')
data_canadian = enhance_feature_structure_z(data_canadian, VF_baseline_canadian, 'VF_RTT', 'VF_RTT_z', on=['length', 'target'], mean_name = 'mean', std_name = 'std')
data_canadian.to_csv(data_canadian_featurized)

# Use for prediction on new target sites and fusion proteins: Barnacle screen
logging.info('Analyzing barnacle goose screen data for tags in many target sites')
data_barnacle = pd.read_csv(barnacle_path)
data_barnacle = enhance_feature_sequence(data_barnacle, seq = 'insertion')
data_barnacle = enhance_feature_pairing(data_barnacle, seq = 'insertion')
data_barnacle = enhance_feature_structure(data_barnacle)
VF_baseline_barnacle = get_VFmeanstd_splitRTT(data_barnacle, id_length='length', id_part='HA', rttmap={v: k for k, v in pegRNAmaps['HA'].items()}, part = 'RTT')
data_barnacle = enhance_feature_structure_z(data_barnacle, VF_baseline_barnacle, 'VF_RTT', 'VF_RTT_z', on=['length', 'target'], mean_name = 'mean', std_name = 'std')
data_barnacle.to_csv(data_barnacle_featurized)

#############################################################
###########         Validation datasets           ###########
#############################################################

# Data from Choi et al.
logging.info('Analyzing ticker tape screen data')

# for NNNGGA
choidata = pd.read_csv(path_N3)
choidata = choidata[(((choidata[' Rep1_read'] > 20) & (choidata[' Rep2_read'] > 20)) & ((choidata[' Rep3_read'] > 20) & (choidata[' Plasmid_read '] > 30)))] # filter out low covered reads as in publication
choidata['spacer'] = 'GGATGATGGTGAGCACGTGA'
choidata['HA'] = 'TGATGGTGA'
choidata['target'] = 'TAPE1'
choidata['PBS'] = 'GATGGTGAGCACG'
choidata['insertion'] = choidata['Insert']
choidata['mmr'] = False
choidata = enhance_feature_sequence(choidata, seq = 'insertion')
choidata = enhance_feature_pairing(choidata, seq = 'insertion')
choidata = enhance_feature_structure(choidata)
seqtemp = []
for k in range(1000):
    seqtemp.append(DNA(3))
VFtemp = Parallel(n_jobs=8)(delayed(get_vf)(j + 'GGATGATGGTGA' ) for j in seqtemp)
choidata['VF_RTT_mean'] =  np.mean(VFtemp)
choidata['VF_RTT_std'] =  np.std(VFtemp)
choidata['VF_RTT_z'] = choidata.apply(lambda x: (x['VF_RTT'] - x['VF_RTT_mean'])/x['VF_RTT_std'], axis = 1)
choidata['VF_RTT_z'] = choidata['VF_RTT_z'].apply(lambda x: 0 if math.isnan(x) else x)
choidata['ins_freq'] = choidata.apply(lambda x: (x[' Rep1_read']/choidata[' Rep1_read'].sum() + x[' Rep2_read']/choidata[' Rep2_read'].sum() + x[' Rep3_read']/choidata[' Rep3_read'].sum()) / 3, axis = 1)
choidata['plasmid_freq'] = choidata.apply(lambda x: x[' Plasmid_read ']/ choidata[' Plasmid_read '].sum(), axis = 1)
choidata['percIns'] = choidata['ins_freq']/choidata['plasmid_freq'] 
choidata['percIns_z'] = choidata['percIns'].apply(lambda x: (x - choidata['percIns'].mean()) / choidata['percIns'].std() )
choidata.to_csv(data_choi_featurized_eN3, index = False)

# for NNNNNNGGA
choidata = pd.read_csv(path_N6)
choidata = choidata[(((choidata[' Rep1_read'] > 20) & (choidata[' Rep2_read'] > 20)) & ((choidata[' Rep3_read'] > 20) & (choidata[' Plasmid_read '] > 30)))] # filter out low covered reads as in publication
choidata['spacer'] = 'GGATGATGGTGAGCACGTGA'
choidata['HA'] = 'TGATGGTGA'
choidata['target'] = 'TAPE1'
choidata['PBS'] = 'GATGGTGAGCACG'
choidata['insertion'] = choidata['Insert']
choidata['mmr'] = False
choidata = enhance_feature_sequence(choidata, seq = 'insertion')
choidata = enhance_feature_pairing(choidata, seq = 'insertion')
choidata = enhance_feature_structure(choidata)
seqtemp = []
for k in range(1000):
    seqtemp.append(DNA(6))
VFtemp = Parallel(n_jobs=8)(delayed(get_vf)(j + 'GGATGATGGTGA' ) for j in seqtemp)
choidata['VF_RTT_mean'] =  np.mean(VFtemp)
choidata['VF_RTT_std'] =  np.std(VFtemp)
choidata['VF_RTT_z'] = choidata.apply(lambda x: (x['VF_RTT'] - x['VF_RTT_mean'])/x['VF_RTT_std'], axis = 1)
choidata['VF_RTT_z'] = choidata['VF_RTT_z'].apply(lambda x: 0 if math.isnan(x) else x)
choidata['ins_freq'] = choidata.apply(lambda x: (x[' Rep1_read']/choidata[' Rep1_read'].sum() + x[' Rep2_read']/choidata[' Rep2_read'].sum() + x[' Rep3_read']/choidata[' Rep3_read'].sum()) / 3, axis = 1)
choidata['plasmid_freq'] = choidata.apply(lambda x: x[' Plasmid_read ']/ choidata[' Plasmid_read '].sum(), axis = 1)
choidata['percIns'] = choidata['ins_freq']/choidata['plasmid_freq'] 
choidata['percIns_z'] = choidata['percIns'].apply(lambda x: (x - choidata['percIns'].mean()) / choidata['percIns'].std() )
choidata.to_csv(data_choi_featurized_eN6, index = False)

# Data from DeepPE
deeppedata = pd.read_csv(path_deeppe)
deeppedata['percIns'] = deeppedata['Measured PE efficiency']

scaling_factor_z = {}
for i in list(set(deeppedata.guide)):
    mean_site = deeppedata[deeppedata.guide == i]['percIns'].mean()
    std_site = deeppedata[deeppedata.guide == i]['percIns'].std()
    scaling_factor_z[i] = [mean_site, std_site]
deeppedata['scaling'] = deeppedata.guide.map(scaling_factor_z) 
deeppedata= deeppedata[deeppedata.model == "PE_type"] # others change the position for G substitution
deeppedata= deeppedata[deeppedata['extension'].str.contains('g|c|a|t', case=True, na=False)] # find the ones that have an insertion/substitution, but no deletion
deeppedata['split'] = deeppedata.extension.apply(lambda x: re.findall(r'[a-z]+|[A-Z]+', x)) # find the edit
deeppedata = deeppedata[deeppedata.split.apply(len) == 3] # only insertions and substitutions in one location
deeppedata['before_edit_if_insertion'] = deeppedata.split.apply(lambda x: x[0]) + deeppedata.split.apply(lambda x: x[2])
deeppedata['before_edit_if_insertion_rc'] = deeppedata['before_edit_if_insertion'].apply(reverse_complement)
deeppedata['is_ins'] = deeppedata.apply(lambda x: str_contains(x['before_edit_if_insertion_rc'], x['wide_seq']), axis = 1) 
deeppedata = deeppedata[deeppedata.is_ins == True] # filter for insertions
deeppedata = deeppedata.groupby('guide').filter(lambda x: x['is_ins'].count() ==7) # only have target sites with data for all 7 insertions
deeppedata['percIns_z'] = deeppedata.apply(lambda x:  (x['percIns'] - x['scaling'][0]) / x['scaling'][1], axis = 1)
deeppedata['insertion'] = deeppedata.split.apply(lambda x: reverse_complement(x[1]))
deeppedata['PBS'] = deeppedata.split.apply(lambda x: reverse_complement(x[2]))
deeppedata['HA'] = deeppedata.split.apply(lambda x: reverse_complement(x[0]))
deeppedata['RTT'] = deeppedata['insertion'] + deeppedata['HA']
deeppedata['spacer'] = deeppedata['guide']
deeppedata['target'] = deeppedata['guide'] # no naming in dataset
deeppedata['mmr'] = False
deeppedata = enhance_feature_sequence(deeppedata, seq = 'insertion')
deeppedata = enhance_feature_pairing(deeppedata, seq = 'insertion')
deeppedata = enhance_feature_structure(deeppedata)
VF_baseline_deeppe = pd.read_csv(VF_baseline_data_deeppe)
deeppedata = deeppedata.merge(VF_baseline_deeppe, on=['HA', 'length'], how = 'left')
deeppedata['VF_RTT_z'] = deeppedata.apply(lambda x: (x['VF_RTT'] - x['VF_RTT_mean'])/x['VF_RTT_std'], axis = 1)
deeppedata.to_csv(data_deeppe_featurized, index = False)