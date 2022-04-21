import regex as re
from Bio.SeqUtils import MeltingTemp as mt
import RNA
import pandas as pd
import numpy as np
from Bio.Data.IUPACData import ambiguous_dna_values
from itertools import product
import more_itertools
import math
import Bio
from joblib import delayed, Parallel
import random

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
def get_countN(x,n):
    """Count the number of nucleotide n in the string."""
    return x.upper().count(n.upper())
def get_Ncomp(x, n):
    if n == 'GC':
        return (get_countN(x,'G') + get_countN(x,'C')) / len(x)
    else:
        return get_countN(x,n) / len(x)

def get_Nrun_max(x,n):
    """Find all consecutive occurences of n and maximum length"""
    my_regex = r"(?i)" + n + n + "+"
    try:
        return max([len(i) for i in re.findall(my_regex, x)])
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

def pairing_bases(seq, spacer, pbs):
    if len(seq) <= 3:
        x = spacer + "GTN&" + reverse_complement(pbs + seq + 'NNN')
    else:
        x = spacer + "GTN&" + reverse_complement(pbs + seq[:3] + 'NNN')
    brackets = str(RNA.cofold(x)[0])
    count = brackets.count("(", 17, 20)
    return count

def loop1_intact(seq):
    scaffold = 'gttttagagctagaaatagcaagttaaaataaggctagtccgttatcaacttgaaaaagtggcaccgagtcggtgc'
    x = scaffold + 'NNN&' + reverse_complement(seq)
    brackets = RNA.cofold(x)[0]
    if brackets[0:30] == "(((((((.((((....))))...)))))))":
        return True
    else:
        return False

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
        
def DNA(length):
    return ''.join(random.choice('CGTA') for _ in range(length))

# Set up data
def init_df(inserts, spacer, pbs, ha, mmr, mean, std):
    """Generates a pandas dataframe based on the user's input."""    
    allinserts = []
    for insert in inserts: 
        if re.match(".*[RYWSMKHBVDN].*",insert, re.IGNORECASE):
            # If ambigious, find all possibilities and explode dataframe
            insertlist = extend_ambiguous_dna(insert)
            allinserts.extend(insertlist)
        else:
            allinserts.extend([insert])
    # create dataframe with all (exploded) inserts
    allinserts = list(more_itertools.collapse(allinserts))
    df = pd.DataFrame(allinserts, columns =['insert'])
    # add all the batch information
    df = add_peginfo(df, spacer, ha, pbs)
    df = add_batchinfo(df, mmr, mean, std)
    return(df)

def add_batchinfo(df, mmr, mean = None, std = None):
    df['mmr'] = mmr
    df['mean'] = mean
    df['std'] = std
    return (df)

def add_peginfo(df, spacer, ha, pbs):
    df['spacer'] = str(spacer).upper()
    df['HA'] = str(ha).upper()
    df['PBS'] = str(pbs).upper()
    return (df)

def extend_ambiguous_dna(seq):
   """return list of all possible sequences given an ambiguous DNA input"""
   return list(map("".join, product(*map(ambiguous_dna_values.get, seq))))

def extend_nt(df):
    df['seqlist'] = df['insert'].apply(lambda x: extend_ambiguous_dna(x))
    df = df.explode('seqlist')
    df = df.rename(columns = {'insert': 'inputsequence', 'seqlist': 'insert'})
    return df

def aa2dna(aa):
    aminoacids = {
        'A': ['GCA', 'GCC', 'GCG', 'GCT'],'C': ['TGC', 'TGT'],'D': ['GAC', 'GAT'],
        'E': ['GAA', 'GAG'],
        'F': ['TTC', 'TTT'],
        'G': ['GGA', 'GGC', 'GGG', 'GGT'],
        'H': ['CAC', 'CAT'],
        'I': ['ATA', 'ATC', 'ATT'],
        'K': ['AAA', 'AAG'],
        'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
        'M': ['ATG'],
        'N': ['AAC', 'AAT'],
        'P': ['CCA', 'CCC', 'CCG', 'CCT'],
        'Q': ['CAA', 'CAG'],
        'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
        'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
        'T': ['ACA', 'ACC', 'ACG', 'ACT'],
        'V': ['GTA', 'GTC', 'GTG', 'GTT'],
        'W': ['TGG'],
        'Y': ['TAC', 'TAT'],
        '_': ['TAA', 'TAG', 'TGA'],
        '*': ['TAA', 'TAG', 'TGA']
        }
    codonlist = [aminoacids[x] for x in aa]
    dna = list(map("".join, product(*codonlist)))
    return dna

def extend_aa(df):
    df['dnalist'] = df.insert.apply(lambda x: aa2dna(x))
    df = df.explode('dnalist')
    df = df.rename(columns = {'insert': 'protein', 'dnalist': 'insert'})
    return df

def get_pegrna(seq, halen, pbslen, spclen):
    seq = str(seq)
    # find the first brackets {}
    try:
        bracketpos = re.search(r'\{\}', seq).start() # this is the position of the nucleotide after the nick
    except ValueError:
        print("There were no brackets in the sequence to indicate the position of the edit.")
    # Remove the brackets from the sequence
    seq = seq[:bracketpos] + seq[bracketpos+2:]
    # Check if this edit is at the nicking site
    try:
        pampos = re.search(r'.GG', seq[bracketpos+3:bracketpos+6]).start() + bracketpos + 3 # this is the first position of the pam
        # Generate features
        spacer = seq[pampos-spclen:pampos]
        ha = seq[pampos-3:pampos-3+halen]
        pbs = seq[pampos-3-pbslen:pampos-3]
    except:
        raise NotImplementedError('We currently do not accept inserts that are not at the nicking site.')
        # try:
        #     # check if there is a Pam nearby
        #     pampos = bracketpos - 3 - re.search(r'.GG', seq[bracketpos-25:bracketpos+3][::-1]).start()  # take the potential ranges of pam, reverse the sequence to find the one that is closest to the nick
        #     # Generate the features
        #     print("The insert is not at nicking site, but RTT is extended to accomodate this.")
        #     spacer = seq[pampos-spclen:pampos]
        #     rtt = seq[pampos-3:pampos-3+rttlen] # this RTT must include the insert! Which is why this code is not finished yet
        #     pbs = seq[pampos-3-pbslen:pampos-3]
        # except ValueError:
        #     print("Could not find a Pam after indicated insertion position.")

    return spacer, ha, pbs

def get_VF_baseline(df):
    VF_baseline = {}
    loi = list(set(df['length']))
    #  Generate 1000 random sequences for each length of interest
    for l in loi:
        seqtemp = []
        VF_baseline[l] = {}
        for k in range(1000):
            seqtemp.append(DNA(l))
        # Calculate VF values for all those sequences per target site and append to list
        for t in df['HA'].unique():
            # For each RTT, create variableseq + RTT and calculate its VF value
            VFtemp = Parallel(n_jobs=8)(delayed(get_vf)(j + t) for j in seqtemp)
            # Calculate mean and std of the VF values
            mean = np.mean(VFtemp)
            std = np.std(VFtemp)
            # append to dictionary
            VF_baseline[l][t] = [mean, std]
    # Create dataframe
    baseline_df = pd.melt(pd.DataFrame.from_dict(VF_baseline, orient='index').reset_index(),id_vars=["index"]).rename(columns={'index': 'length','variable':'HA'})
    baseline_df[['mean_RTT', 'std_RTT']] = baseline_df.value.values.tolist()
    return baseline_df

def enhance_feature_df(df):
    df['RTT'] = df['insert'] + df['HA']
    # length features
    df['length'] = df['insert'].apply(get_length)
    df['length_RTT'] = df['RTT'].apply(get_length)
    # Nucleotide composition
    df['percA'] = df['insert'].apply(lambda x: get_Ncomp(x, 'A'))
    df['percC'] = df['insert'].apply(lambda x: get_Ncomp(x, 'C'))
    df['percT'] = df['insert'].apply(lambda x: get_Ncomp(x, 'T'))
    # Find runs
    df['Arun_maxlen'] = df['insert'].apply(lambda x: get_Nrun_max(x, 'A'))
    # Pairing
    df['pairedbases'] = df.apply(lambda x: pairing_bases(x['insert'], x['spacer'], x['PBS']), axis = 1)
    df['pos1compl'] = df.apply(lambda x: pair_bases_pos(x['spacer'], x['RTT'], 0), axis=1)
    df['loop1_intact'] = df['insert'].parallel_apply(loop1_intact)
    # Structure
    df['VF_RTT'] = df['RTT'].parallel_apply(get_vf)
    # Normalized structure
    VF_baseline = get_VF_baseline(df)
    df = df.merge(VF_baseline[['mean_RTT', 'std_RTT','length','HA']], on = ['length','HA'], how = 'left')
    df['VF_RTT_z'] = (df['VF_RTT'] - df['mean_RTT']) / df['std_RTT']
    df['VF_RTT_z'] = df['VF_RTT_z'].apply(lambda x: 0 if math.isnan(x) else x)
    df['VF_RTT_z'] = df['VF_RTT_z'].apply(lambda x: 0 if math.isinf(x) else x)
    return df
