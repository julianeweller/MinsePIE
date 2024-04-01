from pandarallel import pandarallel
from datetime import datetime
import pandas as pd
from argparse import ArgumentError
import xgboost
from func_features import *
from func_input import *
from func_score import *
import pathlib

# minsepie.predict(['TGTCA'], pbs = 'CAGACTGAGCACG', ha = 'TGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCA', spacer = 'GGCCCAGACTGAGCACGTGA', mmr = 0, outdir = "./")
current_dir = pathlib.Path(__file__).parent.resolve()

def predict(insert, fasta = None, pbs = None, ha = None, spacer = None,  pbslen = 13, halen = 15, spclen = 20, mmr = 0, inputmode = None, cellline = None, outdir = None, mean = None, std = None, model = None):
    pandarallel.initialize()
    """ Predicts editing outcomes for insert sequences based on pegRNA features given individually or determined from fasta sequence. """
    # Clean up variables
    if model == None:
        model = 'MinsePIE_v3.sav'
    if inputmode == None:
        inputmode = 'dna'

    # Retrieve pegRNA features
    if (ha is not None) and ((pbs is not None) and (spacer is not None)):
        pass
    elif fasta is not None:
        try:
            spacer, ha, pbs = get_pegrna(fasta, halen, pbslen, spclen)
        except:
            raise ArgumentError("Please check your target site input and define the input position with brackets.")
    else:
        raise ArgumentError("Please provide either a target sequence as fasta or the PBS, RTT and spacer sequence.")
    
    # Retrieve mmr if mmr is not given
    if mmr == None:
        mmr = cellline2mmr(cellline, current_dir / 'celllines.csv', head = 'mmr')

    # Load the model
    model_dict = load_model(current_dir / 'models')
    
    # Create the dataframe
    request = init_df(insert, spacer, pbs, ha, mmr, mean, std)

    # We have the basic table with sequences now, but if it's an amino acid sequence, we need to convert it to DNA, or if it's DNA we need to accept ambigiuity
    if inputmode == 'protein':
        request = extend_aa(request)
    elif inputmode == 'dna':
        request = extend_nt(request)
    
    # Calculate features
    request = enhance_feature_df(request)

    # Prediction
    request = prediction(request, model_dict, model)

    if not outdir == None:
        outpath = os.path.join(outdir, datetime.now().strftime("%Y%m%d-%H%M%S") + '_minsepie_prediction.csv')
        request[['insert','zscore', 'percIns_predicted']].to_csv(outpath, index = False)

    return request

# minsepie.predict_batch("./examples/his6_codonvariants_batchinput_1.txt", outdir = "./")
# minsepie.predict_batch("./examples/his6_codonvariants_batchinput_2.txt", outdir = "./", longoutput = True)

def predict_batch(dfpath, inputmode = None, mmr = None, outdir = None, mean = None, std = None, model = None, longoutput = None):
    """ Predicts editing outcomes from batch input. Required columns are 'insert', 'spacer', 'pbs' and 'ha'. Default for mmr is mmr-deficient. If you are using cell lines with mmr, please provide mmr = 1. """
    
    batch_df = pd.read_csv(dfpath, sep = '\t')
    
    # Initialize parallel processing
    pandarallel.initialize()
    
    # Clean up variables
    if model == None:
        model = 'MinsePIE_v3.sav'
    if inputmode == None:
        inputmode = 'dna'
    if mmr == None:
        mmr = 0
    if mean == None:
        mean = 0.5
    if std == None:
        std = 0.05
    if longoutput == None:
        longoutput = False

    # Check what variables have been given and are correct
    try:
        batch_df['insert'] = batch_df['insert'].apply(val_nt)
        batch_df['spacer'] = batch_df['spacer'].apply(val_nt)
        batch_df['PBS'] = batch_df['PBS'].apply(val_nt)
        batch_df['HA'] = batch_df['HA'].apply(val_nt)
    except:
        raise ArgumentError("Please provide the correct input sequences in the correct columns.")

    # Retrieve mmr if mmr is not given
    if 'mmr' not in batch_df.columns:
        batch_df['mmr'] = mmr
    else:   
        batch_df['mmr'] = batch_df['mmr'].apply(lambda x: x if x is not None else mmr)

    # Add mean and std if not given
    if 'mean' not in batch_df.columns:
        batch_df['mean'] = mean
    else:
        batch_df['mean'] = batch_df['mean'].apply(lambda x: x if x is not None else mean)
    
    if 'std' not in batch_df.columns:
        batch_df['std'] = std
    else:   
        batch_df['std'] = batch_df['std'].apply(lambda x: x if x is not None else std)

    # Load the model
    model_dict = load_model(current_dir / 'models')
    
    # We have the basic table with sequences now, but if it's an amino acid sequence, we need to convert it to DNA, or if it's DNA we need to accept ambigiuity
    if inputmode == 'protein':
        batch_df = extend_aa(batch_df)
    elif inputmode == 'dna':
        batch_df = extend_nt(batch_df)

    # Calculate features
    request = enhance_feature_df(batch_df)

    # Prediction
    request = prediction(request, model_dict, model)

    if not outdir == None:
        outpath = os.path.join(outdir, datetime.now().strftime("%Y%m%d-%H%M%S") + '_minsepie_prediction.csv')
        request[['insert','PBS','HA','mmr', 'mean', 'std', 'percIns_predicted', 'zscore']].to_csv(outpath, index = False)
        if longoutput == True:
            outpath_long = os.path.join(outdir, datetime.now().strftime("%Y%m%d-%H%M%S") + '_minsepie_prediction_long.csv')
            request.to_csv(outpath_long, index = False)
    return request


def cellline2mmr(cellline: str, file: str, head = 'mmr') -> int:
    cellline_dict = load_celllines(file, head)
    mmr_status = int(cellline_dict[cellline])
    return mmr_status

def val_seq(input: list) -> list:
    """ Raises an error if input is not a list that contains valid nucleotides or amino acid letters."""

    # Everything should be either IUPAC for nucleotides or amino acids
    if all(bool(re.search('[^-\*AGCTGRYWSMKHBVDNDEFHIKLMNPQV]', item, re.IGNORECASE)) == False for item in input) == False:
        raise argparse.ArgumentTypeError(f"Your input should only contain IUPAC nucleotides (upper or lower case) or amino acids.")
    else:
        return input

def val_nt(input: str) -> str:
    """ Raises an error if input is not a valid nucleotide sequence. """

    # Everything should be unique nucleotides
    if all(bool(re.search('[^AGCT\{\}]', item, re.IGNORECASE)) == False for item in input) == False:
        raise argparse.ArgumentTypeError(f"Your input should only contain the IUPAC nucleotides A, C, T and G.")
    else:
        return input

def val_dir(string: str) -> str:
    """Check that provided string is a directory."""
    if os.path.isdir(string):
        return string
    else:
        raise NotADirectoryError

def val_table(file: str) -> str:
    """Validates that provided file endswith with correct file extension."""
    suffixes = ['.csv', '.tsv', 'txt']
    if not os.path.exists(file):
        raise argparse.ArgumentTypeError(f"{file} does not exist")
    if not (file.endswith(suffixes[0]) or (file.endswith(suffixes[1]) or file.endswith(suffixes[2]))):
        raise argparse.ArgumentTypeError(f"{file} is not the right format. Supply as {suffixes[0]} or {suffixes[1]} or {suffixes[2]}")
    return file

def val_fasta(input: str) -> str:
    """Reads in fasta or text file and returns a nucleotide sequence that corresponds to the target site. Only the first record of the file is used."""
    if type(input) == str:
        if (input.endswith(".fasta") or input.endswith(".txt") or input.endswith(".rtf") or input.endswith(".fa")):
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
            return(val_nt(input)) # returns input
    else:
        raise argparse.ArgumentTypeError("Please provide sequence string or fasta file.")

def cellline2mmr(cellline: str, file: str, head = 'mmr') -> int:
    cellline_dict = load_celllines(file, head)
    mmr_status = int(cellline_dict[cellline])
    return mmr_status

def val_path(input: str) -> str:
    if os.path.exists(input):
        return input
    else:
        raise argparse.ArgumentTypeError(f"{input} does not exist. Please provide an existing path.")