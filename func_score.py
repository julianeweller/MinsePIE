
# import pandas as pd
import xgboost as xgb

def scale_zscore(zscore, mean, std):
    """Calculates the predicited insertion efficiency from the Z-score."""
    zscaled = zscore * std + mean
    return zscaled

def predict(df, model_dict, model = 'MinsePIE_v2.sav'):
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
