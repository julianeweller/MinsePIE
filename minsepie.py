import pickle
import regex as re
from Bio.SeqUtils import MeltingTemp as mt
import RNA
import os, glob
import pandas as pd
import numpy as np
import argparse
import xgboost as xgb
from pandarallel import pandarallel

# Getting features as functions
def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A','a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    rc = "".join(complement.get(base, base) for base in reversed(seq))
    return rc
def get_length(x):
    return len(x)
def get_smaller3(x):
    if len(x) <=3:
        return True
    else:
        return False
def get_countN(x,n):
    return x.upper().count(n.upper())
def get_Nrun(x,n):
    my_regex = r"(?i)" + n + "+" + n + n + n
    if bool(re.search(my_regex, x)) == True:
        return True
    else:
        return False
def get_tm(x):
    return mt.Tm_NN(x)
def get_vf(x):
    vf = RNA.fold(x)
    if vf[1] == np.nan:
        return 0
    else:
        return vf[1]


# Load model
def load_model(modeldir):
    """Loads the models from a directory into a dictionary. Returns dictionary."""
    modellist = [os.path.basename(d) for d in glob.glob(modeldir+  '/*.sav')]
    model_dict = {}
    for model in modellist:
        modelpath = os.path.join(modeldir, model)
        model_temp = pickle.load(open(modelpath, 'rb'))
        model_dict[model] = model_temp
    return model_dict

# Generate features
def create_feature_df(insert, rtt, pbs, mmr):
    """Generates a pandas dataframe based on the user's input."""
    df = pd.DataFrame(data= {'insert' : insert, 'RTT': rtt, 'PBS' : pbs, 'mmr': mmr}, index=[0])
    return df

def enhance_feature_df(df, seq = 'insert', ext = 'extension', full = 'full'):
    """Calculates relevant features based on insert sequence, RTT, PBS and MMR status."""
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

def scale_zscore(zscore, mean, std):
    """Calculates the predicited insertion efficiency from the Z-score."""
    zscaled = zscore * std + mean
    return zscaled

def predict(insert, pbs, rtt, mmr, model_dict, mean = None, std = None, model = 'MinsePIE_v2.sav'):
    """Uses the loaded model to predict insertion efficiency from an input sequence."""
    features = ['length', 'smaller3',  'percGC', 'percG', 'percC', 'percT', 'percA', 
    'countG',  'countT', 'countC', 'countA', 'Crun_ext', 'Arun_ext', 'Trun_ext', 'Grun_ext',
    'VF_insert', 'VF_full', 'VF_ext', 'Tm_NN_ext', 'mmr']

    # make dataframe
    request = create_feature_df(insert, rtt, pbs, mmr)
    # get features
    request = enhance_feature_df(request)
    # choose model
    print(f'Prediction model {model}')
    pred_model = model_dict[model]
    # predict
    zscore = pred_model.predict(xgb.DMatrix(request[features]))
    if (mean is not None) and (std is not None):
        scaled_score = scale_zscore(zscore, mean, std)
    else:
        scaled_score = [np.NAN,np.NAN]
    
    return zscore, scaled_score[0]

def main():
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
    required.add_argument('-i', '--insert', dest = 'insert', type = str, help ='Insert seuquence', required=True)
    required.add_argument('-p', '--pbs', dest = 'pbs', type = str, help = 'Primer binding site of pegRNA', required=True)
    required.add_argument('-r', '--rtt', dest = 'rtt', type = str, help = 'Reverse transcriptase template of pegRNA', required=True)
    optional.add_argument('-m', '--mmr', dest = 'mmr', type = int, default = 0,help ='MMR status of cell line')
    optional.add_argument('-a', '--mean', dest = 'mean', type = float, default = None,help ='Expected mean editing efficiency for experimental setup')
    optional.add_argument('-s', '--std', dest = 'std', type = float, default = None,help ='Expected standard deviation for editing efficiency of experimental setup')
    
    # Parse the CLI:
    args = parser.parse_args()
    if not args.mmr:
        print("MMR status of cell line is considered as 0 (MMR deficient)")

    # Load the model
    model_dict = load_model(modeldir)
    
    # Predict
    zscore, scaled_score = predict(args.insert, args.pbs, args.rtt, args.mmr, model_dict, mean = args.mean, std = args.std)
    
    if (args.mean is not None) and (args.std is not None):
        print(f'Insertion of {args.insert} \n Z-score: {zscore[0]} \n Scaled score based on provided mean and standard deviation {scaled_score}')
    else:
        print(f'Insertion of {args.insert} \n Z-score: {zscore[0]}')

if __name__ == '__main__':
    pandarallel.initialize()
    main()
