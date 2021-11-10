import pickle
import regex as re
from Bio.SeqUtils import MeltingTemp as mt
import RNA
import os, glob
import pandas as pd
import numpy as np
import argparse
import xgboost as xgb



# Load model
def load_model(modeldir):
    """Loads the models from a directory into a dictionary. Returns dictionary."""
    modellist = [os.path.basename(d) for d in glob.glob(modeldir+  '/*')]
    print(modellist)
    model_dict = {}
    for model in modellist:
        modelpath = os.path.join('/Users/jw38/Documents/MinsePIE/models/', model)
        model_temp = pickle.load(open(modelpath, 'rb'))
        model_dict[model] = model_temp
    return model_dict

# Generate features
def create_feature_df(insert, rtt, pbs, mmr):
    """Generates a pandas dataframe based on the user's input."""
    df = pd.DataFrame(data= {'sequence_original' : insert, 'RTT': rtt, 'PBS' : pbs, 'mmr': mmr}, index=[0])
    return df

def enhance_feature_df(df):
    """Calculates relevant features based on insert sequence, RTT, PBS and MMR status."""
    # Create dataframe with relevant features
    df['length'] = df['sequence_original'].apply(lambda x: len(x))
    df['smaller3'] = df['sequence_original'].apply(lambda x: 1 if len(x) <=3 else 0)
    # Bases count
    df.sequence_original = df.sequence_original.astype('str')
    df['countC'] = df['sequence_original'].apply(lambda x: x.count('c' or 'C'))
    df['countG'] = df['sequence_original'].apply(lambda x: x.count('g' or 'G'))
    df['countA'] = df['sequence_original'].apply(lambda x: x.count('a' or 'A'))
    df['countT'] = df['sequence_original'].apply(lambda x: x.count('t' or 'T'))
    # Relative content
    df['percC'] = df['countC'] / df['length']
    df['percG'] = df['countG'] / df['length']
    df['percA'] = df['countA'] / df['length']
    df['percT'] = df['countT'] / df['length']
    df['percGC'] = (df['countG'] + df['countC'])/df['length']
    # Find runs
    df['Arun'] = df['sequence_original'].apply(lambda x: 1 if re.search(r'(?i)a+aaa', x) else 0)
    df['Crun'] = df['sequence_original'].apply(lambda x: 1 if re.search(r'(?i)c+ccc', x) else 0)
    df['Trun'] = df['sequence_original'].apply(lambda x: 1 if re.search(r'(?i)t+ttt', x) else 0)
    df['Grun'] = df['sequence_original'].apply(lambda x: 1 if re.search(r'(?i)g+ggg', x) else 0)
    # Sequence specific features
    df['Tm_NN'] = df['sequence_original'].apply(lambda x: mt.Tm_NN(x))
    df['VF_insert'] = df['sequence_original'].apply(lambda x: (RNA.fold(x)))
    df['VF_insert'] = df['VF_insert'].apply(lambda x: 0 if x[1] == np.nan else x[1])
    df['extension'] = df['PBS'] + df['sequence_original'] + df['RTT']
    df['VF_full'] = df['extension'].apply(lambda x: (RNA.fold(x)))
    df['VF_full'] = df['VF_full'].apply(lambda x: 0 if x[1] == np.nan else x[1])
    return df

def scale_zscore(zscore, mean, std):
    """Calculates the predicited insertion efficiency from the Z-score."""
    zscaled = zscore * std + mean
    return zscaled

def predict(insert, pbs, rtt, mmr, model_dict, mean = None, std = None, model = 'MinsePIE.sav'):
    """Uses the loaded model to predict insertion efficiency from an input sequence."""
    features = [
            'length', 'smaller3',  
            'percGC', 'percG', 'percC', 'percT', 'percA', 
            'countG',  'countT', 'countC', 'countA',
            'Crun', 'Arun', 'Trun', 'Grun',
            'VF_insert','VF_full', 'Tm_NN','mmr']

    # make dataframe
    request = create_feature_df(insert, rtt, pbs, mmr)
    # get features
    request = enhance_feature_df(request)
    # choose model
    pred_model = model_dict[model]
    # predict
    zscore = pred_model.predict(xgb.DMatrix(request[features]))
    if (mean is not None) and (std is not None):
        scaled_score = scale_zscore(zscore, mean, std)
    else:
        scaled_score = np.NAN
    return zscore, scaled_score

def main():
    packagedir = os.path.dirname(os.path.realpath(__file__))
    modeldir = packagedir + '/models/'

    try:
        assert os.path.exists(modeldir)
    except:
        raise FileNotFoundError(f'Could not find saved models. Looked in: {modeldir}')

    # Create the command line interface:
    parser = argparse.ArgumentParser(description='Modeling insertion efficiencies for Prime Editing experiments')
    parser.add_argument('-i', '--insert', dest = 'insert', type = str, help ='Insert seuquence')
    parser.add_argument('-p', '--pbs', dest = 'pbs', type = str, help = 'Primer binding site of pegRNA')
    parser.add_argument('-r', '--rtt', dest = 'rtt', type = str, help = 'Reverse transcriptase template of pegRNA')
    parser.add_argument('-m', '--mmr', dest = 'mmr', type = int, default = 0,help ='MMR status of cell line')
    parser.add_argument('-a', '--mean', dest = 'mean', type = float, default = None,help ='Expected mean editing efficiency for experimental setup')
    parser.add_argument('-s', '--std', dest = 'std', type = float, default = None,help ='Expected standard deviation for editing efficiency of experimental setup')
    
    # Parse the CLI:
    args = parser.parse_args()
    # Load the model
    model_dict = load_model(modeldir)
    # Predict
    zscore, scaled_score = predict(args.insert, args.pbs, args.rtt, args.mmr, model_dict, mean = args.mean, std = args.std)
    
    print(f'Insertion of {args.insert} \n Z-score: {zscore[0]} \n Scaled score based on provided mean and standard deviation {scaled_score}')

if __name__ == '__main__':
    main()