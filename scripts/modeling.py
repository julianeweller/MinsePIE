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
             } 
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

# Import one hot encoded data
data = pd.read_csv('data_onehot.csv')
workdir = './'

# Feature sets
feat_nommr_min = [
    'length', 'smaller3',  
    'percGC', 'percG', 'percC', 'percT', 'percA', 
    'countG',  'countT', 'countC', 'countA',
    'Crun', 'Arun', 'Trun', 'Grun',
    'VF_insert','VF_full', 
    'Tm_NN']

feat_mmr_min = [
    'length', 'smaller3',  
    'percGC', 'percG', 'percC', 'percT', 'percA', 
    'countG',  'countT', 'countC', 'countA',
    'Crun', 'Arun', 'Trun', 'Grun',
    'VF_insert','VF_full', 
     'Tm_NN',
    'mmr']

feat_mmr_min_label = [
    'Length', 'Sequence length <= 3',  
    'GC content [%]', 'G content [%]', 'C content [%]', 'T content [%]', 'A content [%]', 
    'Count G',  'Count T', 'Count C', 'Count A',
    'Has run of C', 'Has run of A', 'Has run of T', 'Has run of G',
    'Secondary structure insert',"Secondary structure 3'extension", 
     'Melting temperature',
    'MMR']

model_list = ['L1', 'RandomForest','L2','MLP','XGBoost']

### Compare models: Cross-validation (Figure 5a)
cross_scores = []
n = 10

for model in model_list:
    result = []
    
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
        X_train = unique_train[feat_nommr_min]
        y_train = unique_train[['percIns_z']]
        X_test = unique_test[feat_nommr_min]
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

### 1 vs 1 comparison (Supplementary Figure 6b)
""" Training the model on one dataset and testing it on another dataset"""
dsetOI = [ 'CLYBL_293T_PE2_3', 'EMX1_293T_PE2_3','FANCF_293T_PE2_3','HEK3_293T_PE2_1', 
          'FANCF_HAP1_PE2_3', 'HEK3_HAP1_PE2_2', 'FANCF_HAP1dMLH1_PE2_4', 'HEK3_HAP1dMLH1_PE2_4']

datas = {}
exps = set(data.experiment)

for i in exps:
    datas[i] = data.loc[data['experiment'] == i] 

crosssite_val_dict = defaultdict(list)    

for dset in dsetOI:
    # make non-overlapping train and test data sets
    # get list of shuffled unique sequence original
    get_unique = list(set(datas[dset]['sequence_original']))
    get_unique = random.sample(get_unique, len(get_unique))
    get_train = get_unique[:int(len(get_unique)*0.8)]
    get_test = get_unique[int(len(get_unique)*0.8):]
    
    # Make subdataframe
    unique_train = datas[dset].loc[datas[dset]['sequence_original'].isin(get_train)]

    # define X and y
    X_train = unique_train[feat_nommr_min]
    y_train = unique_train[['percIns_z']]

    model = xgboost(X_train, y_train, do_fit = True)
    
    for oset in dsetOI:
        unique_test = datas[oset].loc[datas[oset]['sequence_original'].isin(get_test)]
        X_test = unique_test[feat_nommr_min]
        y_test = unique_test[['percIns_z']]

        predictions = model.predict(xgb.DMatrix(X_test))
        r2 = sk.metrics.r2_score(y_test, predictions)
        corr, _ = stats.pearsonr(predictions, y_test.values.reshape(-1))
        
        crosssite_val_dict[dset].append(corr)

crosssite_val_df = pd.DataFrame.from_dict(crosssite_val_dict, orient='index', columns=dsetOI)

### Training on all target sites besides one, testing on a dataset of the target site that was not included in the training (Supplementary Figure 6a)

# Define which dataset is used to train which model
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
        
        X_train = subdata_train[feat_mmr_min]
        y_train = subdata_train[['percIns_z']]
        
        X_test = subdata_test[feat_mmr_min]
        y_test = subdata_test[['percIns_z']]

        # Train model
        model = xgboost(X_train, y_train, do_fit = True)

        predictions = model.predict(xgb.DMatrix(X_test))
        corr, _ = stats.pearsonr(predictions, y_test.values.reshape(-1))

        scores.append(corr)
    
    cross_scores.append(scores)
    
cross_scores_n1sites = pd.DataFrame(data = cross_scores, index = dsetOI).T


### Include MMR as a feature (Figure 5b)
# Use combined model on every target site, but train model either on with or without MMR
n = 5
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
    X_train = unique_train[feat_mmr_min]
    y_train = unique_train[['percIns_z']]
    model_mmr = xgboost(X_train,y_train, do_fit = True)

    # Train model on features without mmr
    X_train = unique_train[feat_nommr_min]
    y_train = unique_train[['percIns_z']]
    model_nommr = xgboost(X_train,y_train, do_fit = True)
    
    for dset in dsetOI:
        # with mmr
        X_test = unique_test.loc[unique_test.experiment == dset][feat_mmr_min]
        y_test = unique_test.loc[unique_test.experiment == dset][['percIns_z']]

        xgb_test = xgb.DMatrix(X_test, label=y_test)
        predictions = model_mmr.predict(xgb_test)
        corr, _ = stats.pearsonr(predictions, y_test.values.reshape(-1))
        mmr_scores.append([i, dset, corr, 'mmr'])

        # without mmr
        X_test = unique_test.loc[unique_test.experiment == dset][feat_nommr_min]
        xgb_test = xgb.DMatrix(X_test, label=y_test)
        predictions = model_nommr.predict(xgb_test)
        corr, _ = stats.pearsonr(predictions, y_test.values.reshape(-1))
        mmr_scores.append([i, dset, corr, 'nommr'])
        
mmr_scores_df = pd.DataFrame(mmr_scores, columns = ['rep', 'experiment','R', 'MMR'])

### Correlation between model prediction and measurement (Figure 5c)
# Train the model on non-bleeding datasets
# get list of shuffled unique sequence original
get_unique = list(set(data_temp['sequence_original']))
get_unique = random.sample(get_unique, len(get_unique))
get_train = get_unique[:int(len(get_unique)*0.8)]
get_test = get_unique[int(len(get_unique)*0.8):]
# Make subdataframe
unique_train = data_temp.loc[data_temp['sequence_original'].isin(get_train)]
unique_test = data_temp.loc[data_temp['sequence_original'].isin(get_test)]

# Train model on features with mmr
X_train = unique_train[feat_mmr_min]
y_train = unique_train[['percIns_z']]
X_test = unique_test[feat_mmr_min]
y_test = unique_test[['percIns_z']]
model = xgboost(X_train,y_train, do_fit = True)

# Evaluate
xgb_test = xgb.DMatrix(X_test, label=y_test)
predictions = model.predict(xgb_test)
corr, _ = stats.pearsonr(predictions, y_test.values.reshape(-1))
r2 = sk.metrics.r2_score(y_test, predictions)
print(f'Pearson R = {corr:.3f} \nPearson R2 = {r2:.3f}')

### SHAP values (Figure 5d)
# Initialize your Jupyter notebook with initjs(), otherwise you will get an error message.
shap.initjs()

# Write in a function
def shap_plot(j):
    explainerModel = shap.TreeExplainer(model)
    shap_values_Model = explainerModel.shap_values(S)
    p = shap.force_plot(explainerModel.expected_value, shap_values_Model[j], S.iloc[[j]])
    return(p)

# model.fit(X_train, Y_train)
shap_values = shap.TreeExplainer(model).shap_values(X_train)
shap.summary_plot(shap_values, X_train, feat_mmr_min_label, show = False, max_display=10)

### His6 tags (Figure 5e and Supplementary Figure 6c)

# train model without his6 tags
data_nohis6 = data.loc[~data.name.str.contains("His6-tag")]

# get list of unique sequences and split into test and train
get_unique = list(set(data_nohis6['sequence_original']))
get_unique = random.sample(get_unique, len(get_unique))

get_train = get_unique[:int(len(get_unique)*0.8)]
get_test = get_unique[int(len(get_unique)*0.8):]

subdata_train = data_nohis6.loc[data['sequence_original'].isin(get_train)]
subdata_test = data_nohis6.loc[data['sequence_original'].isin(get_test)]

X_train = subdata_train[feat_mmr_min]
y_train = subdata_train[['percIns_z']]

X_test = subdata_test[feat_mmr_min]
y_test = subdata_test[['percIns_z']]

# Train model
model_nohis6 = xgboost(X_train, y_train, do_fit = True)

predictions = model_nohis6.predict(xgb.DMatrix(X_test))
corr, _ = stats.pearsonr(predictions, y_test.values.reshape(-1))

# Predict insertion for his tags
his6_array_predIns  = {}

for dset in dsetOI:
    
    his6_df = data.loc[(data.name.str.contains("His6")) & (data.experiment == dset)]
    
    mean_exp = data.loc[data.experiment == dset]['percIns'].mean()
    std_exp = data.loc[data.experiment == dset]['percIns'].std()

    # Use the trained model to predict the efficiency for each tag
    his6_df['percIns_z_predicted'] = model_nohis6.predict(xgb.DMatrix(his6_df[feat_mmr_min]))

    # determine the predicted insertion efficiency from the Z score
    his6_df['percIns_predicted'] = his6_df['percIns_z_predicted'] * std_exp + mean_exp
    insertrates = his6_df['percIns_predicted'].tolist()
    sequences = his6_df['sequence_original'].tolist()

    his6_array_predIns[dset] = pd.DataFrame(insertrates, index = sequences).T

df_his6_predIns = pd.concat([v for k,v in his6_array_predIns.items()]).T
df_his6_predIns.columns = dsetOI

# Get measured insertions for His6 tag
his6_array_measuredIns  = {}

for dset in dsetOI:
    
    his6_df = data.loc[(data.name.str.contains("His6")) & (data.experiment == dset)]
    
    insertrates = his6_df['percIns'].tolist()
    sequences = his6_df['sequence_original'].tolist()

    his6_array_measuredIns[dset] = pd.DataFrame(insertrates, index = sequences).T

df_his6_measuredIns = pd.concat([v for k,v in his6_array_measuredIns.items()]).T
df_his6_measuredIns.columns = dsetOI

# combine data
df_list= []
for i in tags:
    df_list.append(df_his6_predIns[df_his6_predIns.index==i])

df_his6_predIns = pd.concat(df_list)

df_list= []
for i in tags:
    df_list.append(df_his6_measuredIns[df_his6_measuredIns.index==i])

df_his6_measuredIns = pd.concat(df_list)