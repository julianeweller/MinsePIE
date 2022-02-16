import pickle
import regex as re
from Bio.SeqUtils import MeltingTemp as mt
import RNA
import os
import glob
import pandas as pd
import numpy as np
import argparse
import xgboost as xgb
from pandarallel import pandarallel

from minsepie import *

def validate_table(file):
    if not os.path.exists(file):
        raise argparse.ArgumentTypeError(f"{file} does not exist")
    if not (file.endswith('.tsv') or file.endswith('.csv') or file.endswith('.txt')):
        raise argparse.ArgumentTypeError(f"{file} is not the right format. Supply as {suffixes[0]} or {suffixes[1]}")
    return file

def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError

def add_batchinfo(df, pbs, rtt, mmr, mean = None, std = None):
    df['PBS'] = pbs
    df['RTT'] = rtt
    df['mmr'] = mmr
    df['mean'] = mean
    df['std'] = std
    return (df)


def main_batchinsert():
    packagedir = os.path.dirname(os.path.realpath(__file__))
    modeldir = packagedir + '/models/'

    try:
        assert os.path.exists(modeldir)
    except:
        raise FileNotFoundError(f'Could not find saved models. Looked in: {modeldir}')

    # Create the command line interface:
    parser = argparse.ArgumentParser(description='Modeling insertion efficiencies for Prime Editing experiments')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')    
    # all the inserts will have the same RTT, PBS and MMR status
    required.add_argument('-i', '--input', dest = 'input', type = validate_table, help ='Path to csv or tsv table with insert sequences', required=True)
    required.add_argument('-p', '--pbs', dest = 'pbs', type = validate_nt, help = 'Primer binding site of pegRNA', required=True)
    required.add_argument('-r', '--rtt', dest = 'rtt', type = validate_nt, help = 'Reverse transcriptase template of pegRNA', required=True)
    required.add_argument('-o', '--outdir', dest = 'outdir', type = dir_path, help ='Path to output directory', required=True)
    optional.add_argument('-m', '--mmr', dest = 'mmr', type = int, default = 0,help ='MMR status of cell line. 0: deficient, 1: proficient')
    # All are expected to have the same experimental mean and standard deviation
    optional.add_argument('-a', '--mean', dest = 'mean', type = float, default = np.NAN,help ='Expected mean editing efficiency for experimental setup')
    optional.add_argument('-s', '--std', dest = 'std', type = float, default = np.NAN,help ='Expected standard deviation for editing efficiency of experimental setup')
    
    # Parse the CLI:
    args = parser.parse_args()
    if not args.mmr:
        print("MMR status of cell line is considered as 0 (MMR deficient)")

    # Load the model
    model_dict = load_model(modeldir)
    
    # Bring data into shape
    if args.input.endswith('.csv'):
        request = pd.read_csv(args.input, header = None, names = ['insert'])
    elif args.input.endswith('.tsv'):
        request = pd.read_csv(args.input, header = None, sep='\t', names = ['insert'])
    elif args.input.endswith('.txt'):
        request = pd.read_csv(args.input, header = None, sep='\t', names = ['insert'])
    else:
        print("There was an error reading in the input file. Please check the format.")
    
    request = add_batchinfo(request, args.pbs, args.rtt, args.mmr, args.mean, args.std)
    request = enhance_feature_df(request)

    # Predict
    request = predict(request, model_dict)
    
    # Display maximal inserted sequences
    n3_largest = request.nlargest(3,'zscore')

    # Save results
    outfile = os.path.basename(args.input).split('.')[0] + '_minsepie.csv'
    outpath = os.path.join(args.outdir,outfile)

    if (args.mean is not np.NAN) and (args.std is not np.NAN):
        request[['insert','zscore', 'percIns_predicted']].to_csv(outpath)
        print(n3_largest[['insert','zscore','percIns_predicted']])
    else:
        request[['insert','zscore']].to_csv(outpath)
        print(n3_largest[['insert','zscore']])
    
    print(f'Prediction results are saved as {outpath}.')


if __name__ == '__main__':
    pandarallel.initialize(verbose = 1)
    main_batchinsert()