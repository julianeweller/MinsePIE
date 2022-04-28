from func_input import *
from func_features import *

def scale_zscore(zscore, mean, std):
    """Calculates the predicited insertion efficiency from the Z-score."""
    zscaled = zscore * std + mean
    return zscaled

def prediction(df, model_dict, model = 'MinsePIE_v3.sav'):
    """Uses the loaded model to predict insertion efficiency from an input sequence."""
    features = ['length', 'VF_RTT_z', 'mmr', 'percC', 'pairedbases', 'Arun_maxlen',  'percA', 'percT','pos1compl', 'loop1_intact']
    # choose model
    pred_model = model_dict[model]
    # predict
    df['zscore']= pred_model.predict(df[features])
    df['percIns_predicted'] = df['zscore'] * df['std'] + df['mean']

    return df

def zscore(rate, mean, std):
    """Calculates the Z-score from the mean and std."""
    zscore = (rate - mean) / std
    return zscore