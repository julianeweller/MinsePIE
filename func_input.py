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

from func_features import *
from minsepie import *

# Parsing
def parse_args(defaults, cellchoices):
    # Create the command line interface:
    parser = argparse.ArgumentParser(description='Modeling insertion efficiencies for Prime Editing experiments')

    # Have subparsers for different modes
    subparser = parser.add_subparsers(dest='command')

    ### Single target site mode
    single_parser = subparser.add_parser('single', help='Predict insertion rates for a single target site')

    single_parser.add_argument('-i', '--insert', dest = 'insert', type = val_seq, nargs='+', help ='Insert seuquence', required=True) # this can be list
    single_parser.add_argument('-p', '--pbs', dest = 'pbs', type = val_nt, help = 'Primer binding site of pegRNA')
    single_parser.add_argument('-h', '--ha', dest = 'ha', type = val_nt, help = 'Homologous sequence to DNA in reverse transcriptase template of pegRNA')
    single_parser.add_argument('-g', '--spacer', dest = 'spacer', type = val_nt, help = 'Spacer of pegRNA')
    single_parser.add_argument('-f', '--fasta', dest = 'fasta', type = val_fasta, help = 'Target sequence with brackets for place of insertion')
    
    single_parser.add_argument('-hl', '--halen', dest = 'halen', type = str, default = defaults['halen'], help = 'Length of HA')
    single_parser.add_argument('-pl', '--pbslen', dest = 'pbslen', type = str, default = defaults['pbslen'], help = 'Length of PBS')
    single_parser.add_argument('-gl', '--spclen', dest = 'spclen', type = str, default = defaults['spclen'], help = 'Length of spacer')

    single_parser.add_argument('-im', '--inputmode', dest = 'inputmode', choices = ['dna','protein'], default = defaults['inputmode'], help ='Is your insert sequence given as dna or amino acid sequence?')
    single_parser.add_argument('-m', '--mmr', dest = 'mmr', type = int, help ='MMR status of cell line', required=False)
    single_parser.add_argument('-c', '--cellline', dest = 'cell', choices = cellchoices, default = defaults['cellline'], help ='Choose your cell line for MMR status', required=False)

    single_parser.add_argument('-o', '--outdir', dest = 'outdir', type = val_path, help ='Path to output directory', required=False)
    single_parser.add_argument('-a', '--mean', dest = 'mean', type = float, default = defaults['mean'], help ='Expected mean editing efficiency for experimental setup', required=False)
    single_parser.add_argument('-s', '--std', dest = 'std', type = float, default = defaults['std'], help ='Expected standard deviation for editing efficiency of experimental setup', required=False)
        
    ### Batch mode for several target sites
    batch_parser = subparser.add_parser('batch', help='Predict insertion rates for a different target sites')
    batch_parser.add_argument('-i', '--input', dest = 'dfpath', type = val_table, help ='Path to csv/tsv/txt table with insert sequences and pegRNA features', required=True) # this is the file with all the features

    batch_parser.add_argument('-im', '--inputmode', dest = 'inputmode', choices = ['dna','protein'], default = defaults['inputmode'], help ='Is your insert sequence given as dna or amino acid sequence?')
    batch_parser.add_argument('-m', '--mmr', dest = 'mmr', type = int, help ='MMR status of cell line', required=False)

    batch_parser.add_argument('-o', '--outdir', dest = 'outdir', type = val_path, help ='Path to output directory', required=False)
    batch_parser.add_argument('-long', '--longoutput', dest = 'longoutput', type = val_path, help ='Set to True if full output dataframe is prefered', required=False)
    batch_parser.add_argument('-a', '--mean', dest = 'mean', type = float, default = defaults['mean'], help ='Expected mean editing efficiency for experimental setup', required=False)
    batch_parser.add_argument('-s', '--std', dest = 'std', type = float, default = defaults['std'], help ='Expected standard deviation for editing efficiency of experimental setup', required=False)

    args = parser.parse_args()

    return args


# Cellline MMR status
def load_celllines(file: str, head = 'mmr') -> dict:
    if not os.path.exists(file):
        raise ArgumentError(f"{file} does not exist")
    if file.endswith(".csv"):
        celllinestatus = pd.read_csv(file, index_col = 0).to_dict()[head]
    else:
        raise ArgumentError(f"{file} should be a csv file.")
    return celllinestatus


def read_table(filepath: str):
    if filepath.endswith('.csv'):
        request = pd.read_csv(filepath)
    elif filepath.endswith('.tsv'):
        request = pd.read_csv(filepath, sep='\t')
    elif filepath.endswith('.txt'):
        request = pd.read_csv(filepath, sep='\t')
    else:
        raise ArgumentError("There was an error reading in the input file. Please check the format (either .csv, .tsv or .txt).")
    return request

def load_model(modeldir):
    """Loads the models from a directory into a dictionary. Returns dictionary."""
    modellist = [os.path.basename(d) for d in glob.glob(str(modeldir / '*.sav'))]
    model_dict = {}
    for model in modellist:
        modelpath = modeldir / model
        model_temp = pickle.load(open(modelpath, 'rb'))
        model_dict[model] = model_temp
    return model_dict