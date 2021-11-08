import pandas as pd
import os
import random
import numpy as np
import sklearn as sk
from sklearn import linear_model
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPRegressor
from scipy import stats
import xgboost as xgb
import matplotlib.pyplot as plt
from collections import defaultdict
import seaborn as sns
from statannot import add_stat_annotation
import shap
from mycolorpy import colorlist as mcp
import pickle

plt.rcParams['pdf.use14corefonts'] = True

def l1_model(x,y, do_fit=False):
    model = linear_model.Lasso(alpha=0.1)
    if do_fit:
        model.fit(x, y)
    return model

def l2_model(x,y, do_fit=False):
    params = {'alpha': 0.03064709917930863, 
              'fit_intercept': True, 
              'normalize': True, 
              'solver': 'sag'
             } # optimized
    model = sk.linear_model.Ridge(**params)
    if do_fit:
        model.fit(x, y)
    return model
    
def random_forest(x,y, do_fit=False):
    model = RandomForestRegressor(
        max_depth=4, 
        random_state=0,
        n_estimators=100,
    )
    if do_fit:
        model.fit(x, y) 
    return model

def grad_boost_tree(x, y, n_trees = 100, min_samples_leaf = 2, max_depth = 4, learning_rate = 0.1, do_fit=False):
    model = sk.ensemble.GradientBoostingRegressor(
        n_estimators = n_trees,
        min_samples_leaf = min_samples_leaf,
        max_depth = max_depth,
        learning_rate = learning_rate,       
        )
    if do_fit:
        model.fit(x, y)
    return model

def xgboost(x, y, do_fit=False, gamma = 0.1, n_trees=100, learning_rate=0.1, max_depth=5, l1_penalty=0.0, l2_penalty=0.375, colsample = 0.95, test_x=None, test_y=None):
    
    xgb_params = {
        'max_depth': max_depth, 
        'eta': learning_rate, 
        'objective':'reg:squarederror', 
        'lambda': l2_penalty, 
        'alpha': l1_penalty,
        'colsample_bytree': colsample,
        'gamma': gamma,
    }
    
    model = xgb.XGBRegressor(**xgb_params)
    xgb_traindata = xgb.DMatrix(x, label=y)
    
    if do_fit:
        model = xgb.train(
        xgb_params, 
        xgb_traindata, 
        n_trees,  
        verbose_eval=False)               
    
    return model

def mlp(x,y, do_fit = False):
    model = MLPRegressor(
        random_state=1, 
        max_iter=500,
        solver = 'lbfgs')
    
    if do_fit:
        model.fit(x, y)
        
    return model


def call_model(name, x = None, y = None, do_fit = False):

    if name == 'L1':
        model = l1_model(x,y)
    
    elif name == 'L2':
        model = l2_model(x,y)
        
    elif name == 'RandomForest':
        model = random_forest(x,y)
    
    elif name == 'GradientBoostedTree':
        model = grad_boost_tree(x,y)
    
    elif name == 'XGBoost':
        model = xgboost(x, y)
    
    elif name == 'MLP':
        model = mlp(x,y)
        
    if do_fit == True:
        model.fit(x, y)
    
    return model


##### Import one hot encoded data
data = pd.read_csv(os.path.join('data_onehotencoded.csv'))
workdir = './'

feat_nommr = [
    'length', 'smaller3',  
    'percGC', 'percG', 'percC', 'percT', 'percA', 
    'countG',  'countT', 'countC', 'countA',
    'Crun', 'Arun', 'Trun', 'Grun',
    'VF_insert','VF_full', 
    'hairpins', 'Tm_NN',
    'N_1_A', 'N_1_C',
    'N_1_G', 'N_1_T', 'N_2_A', 'N_2_C', 'N_2_G', 'N_2_T', 'N_3_A',
    'N_3_C', 'N_3_G', 'N_3_T', 'N_4_A', 'N_4_C', 'N_4_G', 'N_4_T',
    'N_5_A', 'N_5_C', 'N_5_G', 'N_5_T', 'N_6_A', 'N_6_C', 'N_6_G',
    'N_6_T', 'N_7_A', 'N_7_C', 'N_7_G', 'N_7_T', 'N_8_A', 'N_8_C',
    'N_8_G', 'N_8_T', 'N_9_A', 'N_9_C', 'N_9_G', 'N_9_T', 'N_10_A',
    'N_10_C', 'N_10_G', 'N_10_T', 'N_11_A', 'N_11_C', 'N_11_G',
    'N_11_T', 'N_12_A', 'N_12_C', 'N_12_G', 'N_12_T', 'N_13_A',
    'N_13_C', 'N_13_G', 'N_13_T', 'N_14_A', 'N_14_C', 'N_14_G',
    'N_14_T', 'N_15_A', 'N_15_C', 'N_15_G', 'N_15_T', 'N_16_A',
    'N_16_C', 'N_16_G', 'N_16_T', 'N_17_A', 'N_17_C', 'N_17_G',
    'N_17_T', 'N_18_A', 'N_18_C', 'N_18_G', 'N_18_T', 'N_19_A',
    'N_19_C', 'N_19_G', 'N_19_T', 'N_20_A', 'N_20_C', 'N_20_G',
    'N_20_T', 'N_21_A', 'N_21_C', 'N_21_G', 'N_21_T', 'N_22_A',
    'N_22_C', 'N_22_G', 'N_22_T', 'N_23_A', 'N_23_C', 'N_23_G',
    'N_23_T', 'N_24_A', 'N_24_C', 'N_24_G', 'N_24_T', 'N_25_A',
    'N_25_C', 'N_25_G', 'N_25_T', 'N_26_A', 'N_26_C', 'N_26_G',
    'N_26_T', 'N_27_A', 'N_27_C', 'N_27_G', 'N_27_T', 'N_28_A',
    'N_28_C', 'N_28_G', 'N_28_T', 'N_29_A', 'N_29_C', 'N_29_G',
    'N_29_T', 'N_30_A', 'N_30_C', 'N_30_G', 'N_30_T', 'N_31_A',
    'N_31_C', 'N_31_G', 'N_31_T', 'N_32_A', 'N_32_C', 'N_32_G',
    'N_32_T', 'N_33_A', 'N_33_C', 'N_33_G', 'N_33_T', 'N_34_A',
    'N_34_C', 'N_34_G', 'N_34_T', 'N_35_A', 'N_35_C', 'N_35_G',
    'N_35_T', 'N_36_A', 'N_36_C', 'N_36_G', 'N_36_T', 'N_37_A',
    'N_37_C', 'N_37_G', 'N_37_T', 'N_38_A', 'N_38_C', 'N_38_G',
    'N_38_T', 'N_39_A', 'N_39_C', 'N_39_G', 'N_39_T', 'N_40_A',
    'N_40_C', 'N_40_G', 'N_40_T', 'N_41_A', 'N_41_C', 'N_41_G',
    'N_41_T', 'N_42_A', 'N_42_C', 'N_42_G', 'N_42_T', 'N_43_A',
    'N_43_C', 'N_43_G', 'N_43_T', 'N_44_A', 'N_44_C', 'N_44_G',
    'N_44_T', 'N_45_A', 'N_45_C', 'N_45_G', 'N_45_T', 'N_46_A',
    'N_46_C', 'N_46_G', 'N_46_T', 'N_47_A', 'N_47_C', 'N_47_G',
    'N_47_T', 'N_48_A', 'N_48_C', 'N_48_G', 'N_48_T', 'N_49_A',
    'N_49_C', 'N_49_G', 'N_49_T', 'N_50_A', 'N_50_C', 'N_50_G',
    'N_50_T', 'N_51_A', 'N_51_C', 'N_51_G', 'N_51_T', 'N_52_A',
    'N_52_C', 'N_52_G', 'N_52_T', 'N_53_A', 'N_53_C', 'N_53_G',
    'N_53_T', 'N_54_A', 'N_54_C', 'N_54_G', 'N_54_T', 'N_55_A',
    'N_55_C', 'N_55_G', 'N_55_T', 'N_56_A', 'N_56_C', 'N_56_G',
    'N_56_T', 'N_57_A', 'N_57_C', 'N_57_G', 'N_57_T', 'N_58_A',
    'N_58_C', 'N_58_G', 'N_58_T', 'N_59_A', 'N_59_C', 'N_59_G',
    'N_59_T', 'N_60_A', 'N_60_C', 'N_60_G', 'N_60_T', 'N_61_A',
    'N_61_C', 'N_61_G', 'N_61_T', 'N_62_A', 'N_62_G', 'N_63_C',
    'N_63_G', 'N_64_A', 'N_64_G', 'N_65_A', 'N_65_C', 'N_65_G',
    'N_66_A', 'N_67_C', 'N_68_A', 'N_69_T']

feat_mmr = [
    'length', 'smaller3',  
    'percGC', 'percG', 'percC', 'percT', 'percA', 
    'countG',  'countT', 'countC', 'countA',
    'Crun', 'Arun', 'Trun', 'Grun',
    'VF_insert','VF_full', 
    'hairpins', 'Tm_NN',
    'mmr',     
    'N_1_A', 'N_1_C',
    'N_1_G', 'N_1_T', 'N_2_A', 'N_2_C', 'N_2_G', 'N_2_T', 'N_3_A',
    'N_3_C', 'N_3_G', 'N_3_T', 'N_4_A', 'N_4_C', 'N_4_G', 'N_4_T',
    'N_5_A', 'N_5_C', 'N_5_G', 'N_5_T', 'N_6_A', 'N_6_C', 'N_6_G',
    'N_6_T', 'N_7_A', 'N_7_C', 'N_7_G', 'N_7_T', 'N_8_A', 'N_8_C',
    'N_8_G', 'N_8_T', 'N_9_A', 'N_9_C', 'N_9_G', 'N_9_T', 'N_10_A',
    'N_10_C', 'N_10_G', 'N_10_T', 'N_11_A', 'N_11_C', 'N_11_G',
    'N_11_T', 'N_12_A', 'N_12_C', 'N_12_G', 'N_12_T', 'N_13_A',
    'N_13_C', 'N_13_G', 'N_13_T', 'N_14_A', 'N_14_C', 'N_14_G',
    'N_14_T', 'N_15_A', 'N_15_C', 'N_15_G', 'N_15_T', 'N_16_A',
    'N_16_C', 'N_16_G', 'N_16_T', 'N_17_A', 'N_17_C', 'N_17_G',
    'N_17_T', 'N_18_A', 'N_18_C', 'N_18_G', 'N_18_T', 'N_19_A',
    'N_19_C', 'N_19_G', 'N_19_T', 'N_20_A', 'N_20_C', 'N_20_G',
    'N_20_T', 'N_21_A', 'N_21_C', 'N_21_G', 'N_21_T', 'N_22_A',
    'N_22_C', 'N_22_G', 'N_22_T', 'N_23_A', 'N_23_C', 'N_23_G',
    'N_23_T', 'N_24_A', 'N_24_C', 'N_24_G', 'N_24_T', 'N_25_A',
    'N_25_C', 'N_25_G', 'N_25_T', 'N_26_A', 'N_26_C', 'N_26_G',
    'N_26_T', 'N_27_A', 'N_27_C', 'N_27_G', 'N_27_T', 'N_28_A',
    'N_28_C', 'N_28_G', 'N_28_T', 'N_29_A', 'N_29_C', 'N_29_G',
    'N_29_T', 'N_30_A', 'N_30_C', 'N_30_G', 'N_30_T', 'N_31_A',
    'N_31_C', 'N_31_G', 'N_31_T', 'N_32_A', 'N_32_C', 'N_32_G',
    'N_32_T', 'N_33_A', 'N_33_C', 'N_33_G', 'N_33_T', 'N_34_A',
    'N_34_C', 'N_34_G', 'N_34_T', 'N_35_A', 'N_35_C', 'N_35_G',
    'N_35_T', 'N_36_A', 'N_36_C', 'N_36_G', 'N_36_T', 'N_37_A',
    'N_37_C', 'N_37_G', 'N_37_T', 'N_38_A', 'N_38_C', 'N_38_G',
    'N_38_T', 'N_39_A', 'N_39_C', 'N_39_G', 'N_39_T', 'N_40_A',
    'N_40_C', 'N_40_G', 'N_40_T', 'N_41_A', 'N_41_C', 'N_41_G',
    'N_41_T', 'N_42_A', 'N_42_C', 'N_42_G', 'N_42_T', 'N_43_A',
    'N_43_C', 'N_43_G', 'N_43_T', 'N_44_A', 'N_44_C', 'N_44_G',
    'N_44_T', 'N_45_A', 'N_45_C', 'N_45_G', 'N_45_T', 'N_46_A',
    'N_46_C', 'N_46_G', 'N_46_T', 'N_47_A', 'N_47_C', 'N_47_G',
    'N_47_T', 'N_48_A', 'N_48_C', 'N_48_G', 'N_48_T', 'N_49_A',
    'N_49_C', 'N_49_G', 'N_49_T', 'N_50_A', 'N_50_C', 'N_50_G',
    'N_50_T', 'N_51_A', 'N_51_C', 'N_51_G', 'N_51_T', 'N_52_A',
    'N_52_C', 'N_52_G', 'N_52_T', 'N_53_A', 'N_53_C', 'N_53_G',
    'N_53_T', 'N_54_A', 'N_54_C', 'N_54_G', 'N_54_T', 'N_55_A',
    'N_55_C', 'N_55_G', 'N_55_T', 'N_56_A', 'N_56_C', 'N_56_G',
    'N_56_T', 'N_57_A', 'N_57_C', 'N_57_G', 'N_57_T', 'N_58_A',
    'N_58_C', 'N_58_G', 'N_58_T', 'N_59_A', 'N_59_C', 'N_59_G',
    'N_59_T', 'N_60_A', 'N_60_C', 'N_60_G', 'N_60_T', 'N_61_A',
    'N_61_C', 'N_61_G', 'N_61_T', 'N_62_A', 'N_62_G', 'N_63_C',
    'N_63_G', 'N_64_A', 'N_64_G', 'N_65_A', 'N_65_C', 'N_65_G',
    'N_66_A', 'N_67_C', 'N_68_A', 'N_69_T']

model_list = ['L1', 'RandomForest','L2','MLP','XGBoost','XGBoost MMR']

##### Cross-validation
cross_scores = []
n = 10

for model in model_list:
    result = []
    
    if model == 'XGBoost MMR':
        feat_cross = feat_mmr
        model = 'XGBoost'
    else:
        feat_cross = feat_nommr
    
    for i in range(n):        
        # get list of shuffled unique sequence original
        get_unique = list(set(data['sequence_original']))
        get_unique = random.sample(get_unique, len(get_unique))
        get_train = get_unique[:int(len(get_unique)*0.8)]
        get_test = get_unique[int(len(get_unique)*0.8):]

        # Make subdataframe
        unique_train = data.loc[data['sequence_original'].isin(get_train)]
        unique_test = data.loc[data['sequence_original'].isin(get_test)]

        # define X and y
        X_train = unique_train[feat_cross]
        y_train = unique_train[['percIns_z']]
        X_test = unique_test[feat_cross]
        y_test = unique_test[['percIns_z']]

        # train and test the model
        regmodel = call_model(model, x = X_train, y = y_train, do_fit = True)
        predictions = regmodel.predict(X_test)
        
        if model == 'L2':
            predictions = predictions.ravel()
            
        # corr, _ = stats.pearsonr(predictions, y_test.values.reshape(-1))
        corr, _ = stats.pearsonr(predictions, y_test.values.ravel())

        # save the score
        result.append(corr)
    
    # for each model, add result into array
    cross_scores.append(result)
    
# make it a dataframe with columns
cross_scores = pd.DataFrame(data = cross_scores, index = model_list).T

# Train on n-1 target sites, test on 1 dataset
dset_dict_targetsites = {
 'CLYBL_293T_PE2_3': ['EMX1_293T_PE2_3','FANCF_293T_PE2_3','HEK3_293T_PE2_1','FANCF_HAP1_PE2_3','HEK3_HAP1_PE2_2','FANCF_HAP1dMLH1_PE2_4','HEK3_HAP1dMLH1_PE2_4'],
 'EMX1_293T_PE2_3': ['CLYBL_293T_PE2_3','FANCF_293T_PE2_3','HEK3_293T_PE2_1','FANCF_HAP1_PE2_3','HEK3_HAP1_PE2_2','FANCF_HAP1dMLH1_PE2_4','HEK3_HAP1dMLH1_PE2_4'],
 'FANCF_293T_PE2_3': ['CLYBL_293T_PE2_3','EMX1_293T_PE2_3','HEK3_293T_PE2_1','HEK3_HAP1_PE2_2','HEK3_HAP1dMLH1_PE2_4'],
 'HEK3_293T_PE2_1': ['CLYBL_293T_PE2_3','EMX1_293T_PE2_3','FANCF_293T_PE2_3','FANCF_HAP1_PE2_3','FANCF_HAP1dMLH1_PE2_4'],
 'FANCF_HAP1_PE2_3': ['CLYBL_293T_PE2_3','EMX1_293T_PE2_3','HEK3_293T_PE2_1','HEK3_HAP1_PE2_2','HEK3_HAP1dMLH1_PE2_4'],
 'HEK3_HAP1_PE2_2': ['CLYBL_293T_PE2_3','EMX1_293T_PE2_3','FANCF_293T_PE2_3','FANCF_HAP1_PE2_3','FANCF_HAP1dMLH1_PE2_4'],
 'FANCF_HAP1dMLH1_PE2_4': ['CLYBL_293T_PE2_3','EMX1_293T_PE2_3','HEK3_293T_PE2_1','HEK3_HAP1_PE2_2','HEK3_HAP1dMLH1_PE2_4'],
 'HEK3_HAP1dMLH1_PE2_4': ['CLYBL_293T_PE2_3','EMX1_293T_PE2_3','FANCF_293T_PE2_3','FANCF_HAP1_PE2_3','FANCF_HAP1dMLH1_PE2_4']   
}

cross_scores = []
nfold = 10

for dset in dsetOI:
    scores = []
    
    # do 5 iterations
    for i in range(0,nfold):
    
        # get list of unique sequences and split into test and train
        get_unique = list(set(data['sequence_original']))
        get_unique = random.sample(get_unique, len(get_unique))
        
        get_train = get_unique[:int(len(get_unique)*0.8)]
        get_test = get_unique[int(len(get_unique)*0.8):]

        subdata_train = data.loc[data['sequence_original'].isin(get_train)]
        subdata_test = data.loc[data['sequence_original'].isin(get_test)]
        
        # get dataset
        subdata_train = subdata_train.loc[subdata_train['experiment'].isin(dset_dict_targetsites[dset])]
        subdata_test = subdata_test.loc[subdata_test['experiment'] == dset]
        
        X_train = subdata_train[feat_mmr]
        y_train = subdata_train[['percIns_z']]
        
        X_test = subdata_test[feat_mmr]
        y_test = subdata_test[['percIns_z']]

        # Train model
        model = xgboost(X_train, y_train, do_fit = True)

        predictions = model.predict(xgb.DMatrix(X_test))
        corr, _ = stats.pearsonr(predictions, y_test.values.reshape(-1))

        scores.append(corr)
    
    cross_scores.append(scores)
    
cross_scores_n1sites = pd.DataFrame(data = cross_scores, index = dsetOI).T

##### Improving model with mmr
# Use combined model on every target site, but train model either on with or without MMR
n = 10
mmr_scores = []


for i in range(n):
    # get list of shuffled unique sequence original
    get_unique = list(set(data['sequence_original']))
    get_unique = random.sample(get_unique, len(get_unique))
    get_train = get_unique[:int(len(get_unique)*0.8)]
    get_test = get_unique[int(len(get_unique)*0.8):]
    # Make subdataframe
    unique_train = data.loc[data['sequence_original'].isin(get_train)]
    unique_test = data.loc[data['sequence_original'].isin(get_test)]

    # Train model on features with mmr
    X_train = unique_train[feat_mmr]
    y_train = unique_train[['percIns_z']]
    model_mmr = xgboost(X_train,y_train, do_fit = True)

    # Train model on features without mmr
    X_train = unique_train[feat_nommr]
    y_train = unique_train[['percIns_z']]
    model_nommr = xgboost(X_train,y_train, do_fit = True)
    
    for dset in dsetOI:
        # with mmr
        X_test = unique_test.loc[unique_test.experiment == dset][feat_mmr]
        y_test = unique_test.loc[unique_test.experiment == dset][['percIns_z']]

        xgb_test = xgb.DMatrix(X_test, label=y_test)
        predictions = model_mmr.predict(xgb_test)
        corr, _ = stats.pearsonr(predictions, y_test.values.reshape(-1))
        mmr_scores.append([i, dset, corr, 'mmr'])

        # without mmr
        X_test = unique_test.loc[unique_test.experiment == dset][feat_nommr]
        xgb_test = xgb.DMatrix(X_test, label=y_test)
        predictions = model_nommr.predict(xgb_test)
        corr, _ = stats.pearsonr(predictions, y_test.values.reshape(-1))
        mmr_scores.append([i, dset, corr, 'nommr'])
        
mmr_scores_df = pd.DataFrame(mmr_scores, columns = ['rep', 'experiment','R', 'MMR'])

##### Correlation between model and prediction
# Train the model on non-bleeding datasets
# get list of shuffled unique sequence original
get_unique = list(set(data['sequence_original']))
get_unique = random.sample(get_unique, len(get_unique))
get_train = get_unique[:int(len(get_unique)*0.8)]
get_test = get_unique[int(len(get_unique)*0.8):]
# Make subdataframe
unique_train = data.loc[data['sequence_original'].isin(get_train)]
unique_test = data.loc[data['sequence_original'].isin(get_test)]

# Train model on features with mmr
X_train = unique_train[feat_mmr]
y_train = unique_train[['percIns_z']]
X_test = unique_test[feat_mmr]
y_test = unique_test[['percIns_z']]
model = xgboost(X_train,y_train, do_fit = True)

# Evaluate
xgb_test = xgb.DMatrix(X_test, label=y_test)
predictions = model.predict(xgb_test)
corr, _ = stats.pearsonr(predictions, y_test.values.reshape(-1))
r2 = sk.metrics.r2_score(y_test, predictions)

##### SHAP

# Initialize your Jupyter notebook with initjs(), otherwise you will get an error message.
shap.initjs()

# Write in a function
def shap_plot(j):
    explainerModel = shap.TreeExplainer(model)
    shap_values_Model = explainerModel.shap_values(S)
    p = shap.force_plot(explainerModel.expected_value, shap_values_Model[j], S.iloc[[j]])
    return(p)

pd.DataFrame(data=shap_values)

# model.fit(X_train, Y_train)
shap_values = shap.TreeExplainer(model).shap_values(X_train)
shap.summary_plot(shap_values, X_train, feat_mmr_label, show = False, max_display=10)


# model.fit(X_train, Y_train)
shap_values = shap.TreeExplainer(model).shap_values(X_train)
shap.summary_plot(shap_values, X_train, feat_mmr_label, show = False, max_display=10)


# His6 tags
tags =['CACCACCACCACCACCAC', 'CATCACCATCACCATCAC', 'CATCATCATCACCACCAC','CATCATCATCATCATCAT',
       'GTGGTGGTGGTGGTGGTG','GTGGTGGTGATGATGATG','ATGATGATGATGATGATG','GTGATGGTGATGGTGATG']

# train model without his6 tags
data_nohis6 = data.loc[~data.name.str.contains("His6-tag")]

# get list of unique sequences and split into test and train
get_unique = list(set(data_nohis6['sequence_original']))
get_unique = random.sample(get_unique, len(get_unique))

get_train = get_unique[:int(len(get_unique)*0.8)]
get_test = get_unique[int(len(get_unique)*0.8):]

subdata_train = data_nohis6.loc[data['sequence_original'].isin(get_train)]
subdata_test = data_nohis6.loc[data['sequence_original'].isin(get_test)]

X_train = subdata_train[feat_mmr]
y_train = subdata_train[['percIns_z']]

X_test = subdata_test[feat_mmr]
y_test = subdata_test[['percIns_z']]

# Train model
model_nohis6 = xgboost(X_train, y_train, do_fit = True)

predictions = model_nohis6.predict(xgb.DMatrix(X_test))
corr, _ = stats.pearsonr(predictions, y_test.values.reshape(-1))


his6_array_predIns  = {}

for dset in dsetOI:
    
    his6_df = data.loc[(data.name.str.contains("His6")) & (data.experiment == dset)]
    
    mean_exp = data.loc[data.experiment == dset]['percIns'].mean()
    std_exp = data.loc[data.experiment == dset]['percIns'].std()

    # Use the trained model to predict the efficiency for each tag
    his6_df['percIns_z_predicted'] = model_nohis6.predict(xgb.DMatrix(his6_df[feat_mmr]))
    his6_array_predZ.append(his6_df['percIns_z_predicted'].tolist())

    # determine the predicted insertion efficiency from the Z score
    his6_df['percIns_predicted'] = his6_df['percIns_z_predicted'] * std_exp + mean_exp
    insertrates = his6_df['percIns_predicted'].tolist()
    sequences = his6_df['sequence_original'].tolist()

    his6_array_predIns[dset] = pd.DataFrame(insertrates, index = sequences).T

df_his6_predIns = pd.concat([v for k,v in his6_array_predIns.items()]).T
df_his6_predIns.columns = dsetOI

df_his6_predIns

his6_array_measuredIns  = {}

for dset in dsetOI:
    
    his6_df = data.loc[(data.name.str.contains("His6")) & (data.experiment == dset)]
    
    insertrates = his6_df['percIns'].tolist()
    sequences = his6_df['sequence_original'].tolist()

    his6_array_measuredIns[dset] = pd.DataFrame(insertrates, index = sequences).T

df_his6_measuredIns = pd.concat([v for k,v in his6_array_measuredIns.items()]).T
df_his6_measuredIns.columns = dsetOI
