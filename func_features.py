import regex as re
from Bio.SeqUtils import MeltingTemp as mt
import RNA
import pandas as pd
import numpy as np
from Bio.Data.IUPACData import ambiguous_dna_values
from itertools import product
import itertools
import more_itertools
from pandarallel import pandarallel


# Getting features as functions
def reverse_complement(seq):
    """ Get the reverse complement for DNA."""
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    rc = "".join(complement.get(base, base) for base in reversed(seq))
    return rc

def get_length(x):
    """ Calculate length."""
    return len(x)

def get_smaller3(x):
    """ Returns True if length shorter than 3 nucleotides. """
    if len(x) <=3:
        return True
    else:
        return False

def get_countN(x,n):
    """ Counts the occurence of variable n in x."""
    return x.upper().count(n.upper())

def get_Nrun(x,n):
    """ Returns True if sequence has run of n in x. """
    my_regex = r"(?i)" + n + "+" + n + n + n
    if bool(re.search(my_regex, x)) == True:
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
        

# Set up data
def init_df(inserts, spacer, pbs, rtt, mmr, mean, std):
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
    # allinserts = list(itertools.chain.from_iterable(allinserts)) # this doesn't work
    allinserts = list(more_itertools.collapse(allinserts))
    df = pd.DataFrame(allinserts, columns =['insert'])
    # add all the batch information
    df = add_peginfo(df, spacer, rtt, pbs)
    df = add_batchinfo(df, mmr, mean, std)
    return(df)

def add_batchinfo(df, mmr, mean = None, std = None):
    df['mmr'] = mmr
    df['mean'] = mean
    df['std'] = std
    return (df)

def add_peginfo(df, spacer, rtt, pbs):
    df['spacer'] = spacer
    df['RTT'] = rtt
    df['PBS'] = pbs
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


def get_pegrna(seq, rttlen, pbslen, spclen):
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
        rtt = seq[pampos-3:pampos-3+rttlen]
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

    return spacer, rtt, pbs



# Generate features
def enhance_feature_df(df, seq = 'insert', ext = 'extension', full = 'full') -> pd.DataFrame:
    """Calculates relevant features based on insert sequence, RTT, PBS and MMR status."""
    pandarallel.initialize(progress_bar=False)
    # Generate sequences
    df[seq] = df[seq].astype('str')
    df['RTT_rc'] = df['RTT'].apply(lambda x: reverse_complement(x))
    df['PBS_rc'] = df['PBS'].apply(lambda x: reverse_complement(x))
    df['insert_rc'] = df[seq].apply(lambda x: reverse_complement(x))
    df[ext] = df['RTT_rc'] + df['insert_rc'] + df['PBS_rc']
    df[full] = df['PBS'] + df[seq] + df['RTT']
    # Length features
    df['length'] = df[seq].apply(get_length)
    df['length_ext'] = df[ext].apply(get_length)
    df['smaller3'] = df[seq].apply(get_smaller3)
    # Bases count
    df['countC'] = df[seq].apply(lambda x: get_countN(x,'C'))
    df['countG'] = df[seq].apply(lambda x: get_countN(x,'G'))
    df['countA'] = df[seq].apply(lambda x: get_countN(x,'A'))
    df['countT'] = df[seq].apply(lambda x: get_countN(x,'T'))
    df['countC_ext'] = df[ext].apply(lambda x: get_countN(x,'C'))
    df['countG_ext'] = df[ext].apply(lambda x: get_countN(x,'G'))
    df['countA_ext'] = df[ext].apply(lambda x: get_countN(x,'A'))
    df['countT_ext'] = df[ext].apply(lambda x: get_countN(x,'T'))
    # Relative content
    df['percC'] = df['countC'] / df['length'] *100
    df['percG'] = df['countG'] / df['length'] *100
    df['percA'] = df['countA'] / df['length'] *100
    df['percT'] = df['countT'] / df['length'] *100
    df['percC_ext'] = df['countC_ext'] / df['length_ext'] *100
    df['percG_ext'] = df['countG_ext'] / df['length_ext'] *100
    df['percA_ext'] = df['countA_ext'] / df['length_ext'] *100
    df['percT_ext'] = df['countT_ext'] / df['length_ext'] *100
    df['percGC'] = (df['countG'] + df['countC'])/df['length'] *100
    df['percGC_ext'] = (df['countG_ext'] + df['countC_ext'])/df['length_ext'] *100
    # Find runs
    df['Arun_ext'] = df[ext].apply(lambda x: get_Nrun(x,'A'))
    df['Crun_ext'] = df[ext].apply(lambda x: get_Nrun(x,'C'))
    df['Trun_ext'] = df[ext].apply(lambda x: get_Nrun(x,'T'))
    df['Grun_ext'] = df[ext].apply(lambda x: get_Nrun(x,'G'))
    # Secondary structure
    df['Tm_NN_ext'] = df[ext].parallel_apply(get_tm)
    df['Tm_NN'] = df[seq].parallel_apply(get_tm)
    df['VF_full'] = df[full].parallel_apply(get_vf)
    df['VF_ext'] = df[ext].parallel_apply(get_vf)
    df['VF_insert'] = df[seq].parallel_apply(get_vf)
    return df

