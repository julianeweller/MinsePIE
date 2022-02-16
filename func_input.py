from ctypes import ArgumentError
import pickle
import regex as re
from Bio.SeqUtils import MeltingTemp as mt
import os, glob
import pandas as pd
import numpy as np
import argparse

from Bio.Data.IUPACData import ambiguous_dna_values
from itertools import product
import itertools
from Bio import SeqIO

# Parsing
def parse_args(defaults, cellchoices):
    # Create the command line interface:
    parser = argparse.ArgumentParser(description='Modeling insertion efficiencies for Prime Editing experiments')

    # Have subparsers for different modes
    subparser = parser.add_subparsers(dest='command')

    ### Single target site mode
    single_parser = subparser.add_parser('single', help='Predict insertion rates for a single target site')

    single_parser.add_argument('-i', '--insert', dest = 'insert', type = validate_seq, nargs='+', help ='Insert seuquence', required=True) # this can be list
    single_parser.add_argument('-p', '--pbs', dest = 'pbs', type = validate_nt, help = 'Primer binding site of pegRNA')
    single_parser.add_argument('-r', '--rtt', dest = 'rtt', type = validate_nt, help = 'Reverse transcriptase template of pegRNA')
    single_parser.add_argument('-g', '--spacer', dest = 'spacer', type = validate_nt, help = 'Spacer of pegRNA')
    single_parser.add_argument('-f', '--fasta', dest = 'fasta', type = validate_fasta, help = 'Target sequence with brackets for place of insertion')
    
    single_parser.add_argument('-rl', '--rttlen', dest = 'rttlen', type = str, default = defaults['rttlen'], help = 'Length of RTT')
    single_parser.add_argument('-pl', '--pbslen', dest = 'pbslen', type = str, default = defaults['pbslen'], help = 'Length of RTT')
    single_parser.add_argument('-gl', '--spclen', dest = 'spclen', type = str, default = defaults['spclen'], help = 'Length of RTT')
    single_parser.add_argument('-im', '--inputmode', dest = 'inputmode', choices = ['dna','protein'], default = defaults['inputmode'], help ='Is your insert sequence given as dna or amino acid sequence?')

    single_parser.add_argument('-m', '--mmr', dest = 'mmr', type = int, help ='MMR status of cell line', required=False)
    single_parser.add_argument('-c', '--cellline', dest = 'cell', choices = cellchoices, default = defaults['cellline'], help ='Choose your cell line for MMR status', required=False)

    single_parser.add_argument('-o', '--outdir', dest = 'outdir', type = dir_path, help ='Path to output directory', required=False)
    single_parser.add_argument('-a', '--mean', dest = 'mean', type = float, default = defaults['mean'], help ='Expected mean editing efficiency for experimental setup', required=False)
    single_parser.add_argument('-s', '--std', dest = 'std', type = float, default = defaults['std'], help ='Expected standard deviation for editing efficiency of experimental setup', required=False)
        
    ### Batch mode for several target sites
    batch_parser = subparser.add_parser('batch', help='Predict insertion rates for a different target sites')

    batch_parser.add_argument('-i', '--insert', dest = 'insert', type = validate_seq, nargs='+', help ='Insert seuquence', required=True) # this can be list
    batch_parser.add_argument('-f', '--fasta', dest = 'fasta', type = validate_fasta, help = 'Target sequence with brackets for place of insertion')

    batch_parser.add_argument('-rl', '--rttlen', dest = 'rttlen', type = str, default = defaults['rttlen'], help = 'Length of RTT')
    batch_parser.add_argument('-pl', '--pbslen', dest = 'pbslen', type = str, default = defaults['pbslen'], help = 'Length of RTT')
    batch_parser.add_argument('-gl', '--spclen', dest = 'spclen', type = str, default = defaults['spclen'], help = 'Length of RTT')
    batch_parser.add_argument('-im', '--inputmode', dest = 'inputmode', choices = ['dna','protein'], default = defaults['inputmode'], help ='Is your insert sequence given as dna or amino acid sequence?')

    batch_parser.add_argument('-m', '--mmr', dest = 'mmr', type = int, help ='MMR status of cell line', required=False)
    batch_parser.add_argument('-c', '--cellline', dest = 'cell', choices = cellchoices, default = defaults['cellline'], help ='Choose your cell line for MMR status', required=False)

    batch_parser.add_argument('-o', '--outdir', dest = 'outdir', type = dir_path, help ='Path to output directory', required=False)
    batch_parser.add_argument('-a', '--mean', dest = 'mean', type = float, default = defaults['mean'], help ='Expected mean editing efficiency for experimental setup', required=False)
    batch_parser.add_argument('-s', '--std', dest = 'std', type = float, default = defaults['std'], help ='Expected standard deviation for editing efficiency of experimental setup', required=False)

    args = parser.parse_args()

    return args

# functions to validate input
def validate_seq(input):
    """ Make sure input is a list and contains valid nucleotides or amino acids."""

    # Everything should be either IUPAC for nucleotides or amino acids
    if all(bool(re.search('[^-\*AGCTGATCRYWSMKHBVDNDEFHIKLMNPQV]', item, re.IGNORECASE)) == False for item in input) == False:
        raise argparse.ArgumentTypeError(f"Your input should only contain IUPAC nucleotides (upper or lower case) or amino acids.")
    else:
        return input

def validate_nt(input):
    """ Make sure input is a list and contains valid nucleotides or amino acids."""

    # Everything should be unique nucleotides
    if all(bool(re.search('[^AGCT\{\}]', item, re.IGNORECASE)) == False for item in input) == False:
        raise argparse.ArgumentTypeError(f"Your input should only contain IUPAC nucleotides (upper or lower case) or amino acids.")
    else:
        return input

def dir_path(string):
    """Check that provided string is a directory."""
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError

def validate_table(file):
    """Validates that privded file is csv or tsv."""
    if not os.path.exists(file):
        raise argparse.ArgumentTypeError(f"{file} does not exist")
    if not (file.endswith('.tsv') or file.endswith('.csv')):
        raise argparse.ArgumentTypeError(f"{file} is not the right format. Supply as {suffixes[0]} or {suffixes[1]}")
    return file

def validate_fasta(input):
    if type(input) == str:
        if (input.endswith(".fasta") or input.endswith(".txt") or input.endswith(".rtf")):
            with open(input) as handle:
                n = 0
                for record in SeqIO.FastaIO.FastaTwoLineIterator(handle):
                    seq = record.seq
                    n += 1
                    if n == 1:
                        return seq
                    else:
                        pass
        else:
            return(validate_nt(input)) # returns input
    else:
        raise argparse.ArgumentTypeError("Please provide sequence string or fasta file.")

# Load
def load_model(modeldir):
    """Loads the models from a directory into a dictionary. Returns dictionary."""
    modellist = [os.path.basename(d) for d in glob.glob(modeldir+  '/*.sav')]
    model_dict = {}
    for model in modellist:
        modelpath = os.path.join(modeldir, model)
        model_temp = pickle.load(open(modelpath, 'rb'))
        model_dict[model] = model_temp
    return model_dict

def read_table(filepath):
    if filepath.endswith('.csv'):
        request = pd.read_csv(filepath)
    elif filepath.endswith('.tsv'):
        request = pd.read_csv(filepath, sep='\t')
    else:
        raise ArgumentError("There was an error reading in the input file. Please check the format.")
    return request