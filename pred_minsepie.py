from distutils.log import error
import os
import numpy as np
from func_input import *
from func_features import *
from func_score import *
from pandarallel import pandarallel
from datetime import datetime
import pandas as pd
import itertools

import regex as re
from Bio.SeqUtils import MeltingTemp as mt
import RNA
import os, glob
import pandas as pd
import numpy as np
import argparse
import xgboost as xgb
from pandarallel import pandarallel

from Bio.Data.IUPACData import ambiguous_dna_values
from itertools import product
from datetime import datetime


# python pred_minsepie.py single -i ATAACTTCGATAATGTGATGCTATACGAAGTTAT -p CAGACTGAGCACG -r TGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCA -a 4.86 -s 4.28
# python pred_minsepie.py single -i ATAACTT CGATAATGTGATGCT ATACGAAGTTAT -p CAGACTGAGCACG -r TGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCA -a 4.86 -s 4.28
# python pred_minsepie.py single -i ATAACTTC GATAATGTGATG CTATACGAAGTTAT -f AAAAAAAACGTGCAGTCGTCGATGC{}ACTCGGAAACCCGGGTTTAAACCCGGGTTTAAACCGGGTTTAAAC -a 4.86 -s 4.28
# python pred_minsepie.py single -i ATAACTTC GATAATGTGATG CTATACGAAGTTAT -f ./examples/shortsequence.fasta
# python pred_minsepie.py batch -i ATAACTTC GATAATGTGATG CTATACGAAGTTAT -f ./examples/shortsequence.fasta



######## ToDos
# allow for inserts that are not at the nicking site
# batch mode for several target sites (file with all kind of inputs + all bracket locations from fasta)

######## Done
# look up cell lines
# allow for list
# allow for ambigious sequences
# have mutliple inserts at once
# allow for different fasta formats
# raise error when fasta is too short to find the required sequences
# batch import of features with file
# add spacer into input and dataframe

def main():
    packagedir = os.path.dirname(os.path.realpath(__file__))
    modeldir = packagedir + '/models/'

    try:
        assert os.path.exists(modeldir)
    except:
        raise FileNotFoundError(f'Could not find saved models. Looked in: {modeldir}')

    # Set defaults for analysis
    celllinestatus = {  
        # Leukemia cell lines https://academic.oup.com/carcin/article/24/1/31/2608347 https://www.ncbi.nlm.nih.gov/pmc/articles/PMC484197/ 
        'other': 1,
        'HEK293T': 0,
        'HAP1': 1,
        'HAP1dMLH1': 0,
        'Hela': 1,
        'TK6': 1,
        'K562': 1,
        'HL60': 1,
        'Raji': 1,
        'Daudi': 1,
        'SW48':0, 
        'LoVo': 0,
        'Hec-1-A': 0,
        'HCT15': 0,
        'JURKAT': 0,
        'Molt4': 0,
        'CCRF-CEM':0,
        'Nalm6': 0,
        'Reh': 0,
        'SW480':0,
        'HCT116': 0, 
        'PSC': 1, # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4148860/
        'IPSC': 1, # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4148860/
        }
    cellchoices = list(celllinestatus.keys())
    defaults = {'mean' : np.NAN, 'std' : np.NAN, 'cellline' : cellchoices[0], 'mmr': 1, 'rttlen': 30, 'pbslen': 13, 'spclen': 20, 'inputmode': 'dna'}

    # Get input arguments
    args = parse_args(defaults, cellchoices)

    # Define mmr status
    if not args.mmr:
        """If mmr status is given, this will be used. otherwise, if cell line is given, mmr status will be concluded. If neither is given, assume that cell line is mismatch repair proficient"""
        if args.cell:
            try: 
                args.mmr = celllinestatus[args.cell]
            except:
                args.mmr = defaults['mmr']
        else:
            args.mmr = defaults['mmr']

    # Retrieve pegRNA features from input
    if (args.rtt is not None) and (args.pbs is not None): # and (args.spacer is not None)
        pass
    elif args.fasta is not None:
        try:
            args.spacer, args.rtt, args.pbs = get_pegrna(args.fasta, args.rttlen, args.pbslen, args.spclen)
        except:
            raise argparse.ArgumentError("Please check your target site input and define the input position with brackets.")
    else:
        # see if the pegRNA features are in csv file
        if args.command == 'batch':
            # care about those features later
            pass 
        else:
            raise argparse.ArgumentError("Please provide either a target sequence as fasta or the PBS, RTT and spacer sequence.")
    
    # Load the model
    model_dict = load_model(modeldir)
    
    # Create the dataframe
    if args.command == 'single':
        request = init_df(args.insert, args.spacer, args.pbs, args.rtt, args.mmr, args.mean, args.std)
    elif args.command == 'batch':
        request = read_table(args.input)
        # Add the batch information
        if (args.rtt is not None) and (args.pbs is not None) and (args.spacer is not None):
            request = add_batchinfo(request, args.spacer, args.pbs, args.rtt, args.mmr, args.mean, args.std)
        else:
            # check if the information has been in the csv file
            if len(request.columns) > 4:
                # make sure the naming fo columns is correct
                request = request.rename(columns={request.columns[0]: "insert", request.columns[1]: "spacer" , request.columns[2]: "RTT", request.columns[3]: "PBS" }, inplace = True)
                request['mmr'] = args.mmr
                request['mean'] = args.mean
                request['std'] = args.std
            else:
                raise error("Something with the input file seems wrong.")
    else:
        raise argparse.ArgumentError("something went wrong with the input modes.")

    # We have the basic table with sequences now, but if it's an amino acid sequence, we need to convert it to DNA, or if it's DNA we need to accept ambigiuity
    if args.inputmode == 'protein':
        request = extend_aa(request)
    elif args.inputmode == 'dna':
        request = extend_nt(request)
    
    # Calculate features
    request = enhance_feature_df(request)

    # Predict
    request = predict(request, model_dict)

    # Print result
    if len(request.index) > 1:
        # Sort by z-score and return value for all 
        request = request.sort_values(by='zscore', ascending=False)
        # iterate through rows and return value
        if (args.mean is int) and (args.std is int):
            print(request[['insert','zscore', 'percIns_predicted']].head(10))
        else:
            print(request[['insert','zscore']].head(10))
        print("Up to top 10 inserts are printed. If you expect a longer list, please provide an output directory.")
    else:
        zscore = request['zscore'][0]
        scaledz = request['percIns_predicted'][0]
        if (args.mean is int) and (args.std is int):
            print(f'Insertion of {args.insert[0]} \n Z-score: {zscore} \n Scaled score based on provided mean and standard deviation {scaledz}')
        else:
            print(f'Insertion of {args.insert[0]} \n Z-score: {zscore}')
    
    # Save result if outdir given
    if (args.outdir is not None):
        now = datetime.now()
        dt_string = now.strftime("%Y%m%d-%H%M%S")
        outfile = dt_string + '_minsepie_prediction.csv'
        outpath = os.path.join(args.outdir,outfile)
        request[['insert','zscore', 'percIns_predicted']].to_csv(outpath)

if __name__ == '__main__':
    pandarallel.initialize(verbose = 1)
    main()