
# import pandas as pd
from argparse import ArgumentError
import xgboost as xgb
from func_input import *
from func_features import *
from pandarallel import pandarallel
from datetime import datetime

def scale_zscore(zscore, mean, std):
    """Calculates the predicited insertion efficiency from the Z-score."""
    zscaled = zscore * std + mean
    return zscaled

def prediction(df, model_dict, model = 'MinsePIE_v2.sav'):
    """Uses the loaded model to predict insertion efficiency from an input sequence."""
    features = ['length', 'smaller3',  'percGC', 'percG', 'percC', 'percT', 'percA', 
    'countG',  'countT', 'countC', 'countA', 'Crun_ext', 'Arun_ext', 'Trun_ext', 'Grun_ext',
    'VF_insert', 'VF_full', 'VF_ext', 'Tm_NN_ext', 'mmr']

    # choose model
    #print(f'Prediction model {model}')
    pred_model = model_dict[model]
    # predict
    df['zscore']= pred_model.predict(xgb.DMatrix(df[features]))
    df['percIns_predicted'] = df['zscore'] * df['std'] + df['mean']

    return df



def predict(insert, fasta = None, pbs = None, rtt = None, spacer = None,  pbslen = 13, rttlen = 15, spclen = 20, mmr = 0, inputmode = None, cellline = None, outdir = None, mean = None, std = None, model = None):
    """ Predicts editing outcomes for insert sequences based on pegRNA features given individually or determined from fasta sequence. """
    # Clean up variables
    if model == None:
        model = 'MinsePIE_v2.sav'
    if inputmode == None:
        inputmode = 'dna'

    # Retrieve pegRNA features
    if (rtt is not None) and ((pbs is not None) and (spacer is not None)):
        pass
    elif fasta is not None:
        try:
            spacer, rtt, pbs = get_pegrna(fasta, rttlen, pbslen, spclen)
        except:
            raise ArgumentError("Please check your target site input and define the input position with brackets.")
    else:
        raise ArgumentError("Please provide either a target sequence as fasta or the PBS, RTT and spacer sequence.")
    
    # Retrieve mmr if mmr is not given
    if mmr == None:
        mmr = cellline2mmr(cellline, './celllines.csv', head = 'mmr')

    # Load the model
    model_dict = load_model('./models/')
    
    # Create the dataframe
    request = init_df(insert, spacer, pbs, rtt, mmr, mean, std)


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

predict(['TGTCA'], pbs = 'CAGACTGAGCACG', rtt = 'TGATGGCAGAGGAAAGGAAGCCCTGCTTCCTCCA', spacer = 'GGCCCAGACTGAGCACGTGA', mmr = 0, outdir = "./")
