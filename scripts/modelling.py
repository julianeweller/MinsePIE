
import os
import sys
import logging
import math
import re
import xgboost
import pickle
import sklearn
import RNA
import random
import shap
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from collections import defaultdict
from datetime import datetime

from scipy import stats
from scipy.cluster import hierarchy as sch

from sklearn import linear_model
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import RandomForestRegressor

from Bio.SeqUtils import MeltingTemp as mt
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio import Align

#############################################################
###########             INITIALIZE                ########### 
#############################################################

verbosity = 'info'
workdir = '../'
workdir_input = os.path.join(workdir, 'files/output/')
workdir_output = os.path.join(workdir, 'files/output/') # Change to your prefered output directory if you don't want to overwrite the files here
setseed = 0

# Input
data_featurized = os.path.join(workdir_output, 'Data_features_onehot.csv')
data_canadian_featurized = os.path.join(workdir_input,'data_canadian_featurized.csv')
data_barnacle_featurized = os.path.join(workdir_input,'data_barnacle_featurized.csv')
data_padding_featurized = os.path.join(workdir_input,'data_padding_featurized.csv')
pegRNA_mapping = os.path.join(workdir_input, 'Table_S2_pegRNAs.csv')
data_choi_featurized_eN3 = os.path.join(workdir_input,'data_choi_featurized_eN3.csv')
data_choi_featurized_eN6 = os.path.join(workdir_input,'data_choi_featurized_eN6.csv')
data_deeppe_featurized = os.path.join(workdir_input,'data_deeppe_featurized.csv')

# Output
datapath_train = os.path.join(workdir_output, f'Data_train.csv')
datapath_test = os.path.join(workdir_output, f'Data_test.csv')
model_path = os.path.join(workdir, f'models/minsepie_v3.sav')
path_corrmatrix_all = os.path.join(workdir_output, f'model_features_corrmatrix_all.csv')
path_corrmatrix_core = os.path.join(workdir_output, f'model_features_corrmatrix_core.csv')
model_crossvalidation_path = os.path.join(workdir_output, f'model_architecture_crossvalidation_allfeatures.csv')
model_rmvfeatures_path = os.path.join(workdir_output, f'model_features_sets.csv')
model_sites_crossvalidation_path = os.path.join(workdir_output, f'model_sites_crossvalidation.csv')
model_sites_crosscorr_path = os.path.join(workdir_output, f'model_sites_crosscorr.csv')
model_performance_path = os.path.join(workdir_output,f'model_performance.csv')
model_performance_ca_path = os.path.join(workdir_output,f'model_performance_ca.csv')
model_performance_ba_path = os.path.join(workdir_output,f'model_performance_ba.csv')
model_performance_pa_path = os.path.join(workdir_output,f'model_performance_pa.csv')
model_shap_impact_path = os.path.join(workdir,f'model_shap_impact.pdf')
model_shap_importance_path = os.path.join(workdir,f'model_shap_importance.pdf')
model_performance_choi_eN3_path = os.path.join(workdir_output,'model_performance_choi_eN3.csv')
model_performance_choi_eN6_path = os.path.join(workdir_output,'model_performance_choi_eN6.csv')
model_performance_deeppe_path = os.path.join(workdir_output,'model_performance_deeppe.csv')

# logging
logging.getLogger('matplotlib').setLevel(logging.WARNING)

logging.basicConfig(filename= os.path.join(workdir_output,f'log_{datetime.now().strftime("%d%m%Y-%H%M%S")}.txt'),
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.DEBUG)

#############################################################
###########             Functions                 ########### 
#############################################################
def error(msg, exit_code=1):
    """A function to quit with an error."""
    logging.error(msg)
    sys.exit(exit_code)

def writeDataFile(x, filename, datatype):
    """A function to write data to a given file, with correct logging"""
    try:
        with open(filename, 'w') as output_handle: output_handle.write(x)
    except Exception as err: error('failed to write {} data to {} ({})'.format(datatype, filename, err))

def splitdata_nostrat(data, fx = 'insertion', seed = 0, test_size = 0.3):
    # Get all unique sequences
    seqs = sorted(list(set(data[fx])))
    train_seq, test_seq = train_test_split(seqs, test_size= test_size, random_state=seed)
    logging.info(f"Final Training set sequences: {len(train_seq)}, Test set sequences: {len(test_seq)}")
    
    # Merge data from all target sites depending on sequences
    X_train = data[data[fx].isin(train_seq)]
    X_test = data[data[fx].isin(test_seq)]
    logging.info(f"Original dataset length: {data.shape[0]}, Final Training set:  {X_train.shape[0]}, Test set:  {X_test.shape[0]}")
    
    return(X_train, X_test)

def cluster_corr(corr_array, inplace=False):
    """
    Rearranges the correlation matrix, corr_array, so that groups of highly 
    correlated variables are next to eachother 
    
    Parameters
    ----------
    corr_array : pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix 
        
    Returns
    -------
    pandas.DataFrame or numpy.ndarray
        a NxN correlation matrix with the columns and rows rearranged

    Adapted from https://wil.yegelwel.com/cluster-correlation-matrix/. 
    """
    pairwise_distances = sch.distance.pdist(corr_array)
    linkage = sch.linkage(pairwise_distances, method='complete')
    cluster_distance_threshold = pairwise_distances.max()/2
    idx_to_cluster_array = sch.fcluster(linkage, cluster_distance_threshold, 
                                        criterion='distance')
    idx = np.argsort(idx_to_cluster_array)
    
    if not inplace:
        corr_array = corr_array.copy()
    
    if isinstance(corr_array, pd.DataFrame):
        return corr_array.iloc[idx, :].T.iloc[idx, :]

    return corr_array[idx, :][:, idx]

##### models
def call_model(name, x = None, y = None, do_fit = False):
    if name == 'L1':
        model = l1_model(x,y, do_fit = do_fit)
    elif name == 'L2':
        model = l2_model(x,y, do_fit = do_fit)
    elif name == 'RandomForest':
        model = random_forest(x,y, do_fit = do_fit)
    elif name == 'XGBoost':
        model = xgb(x, y, do_fit = do_fit)
    elif name == 'MLP':
        model = mlp(x,y, do_fit = do_fit)
    return model

def l1_model(x,y, do_fit=False):
    model = linear_model.Lasso(alpha=0.01)
    if do_fit:
        model.fit(x, y)
    return model

def l2_model(x,y, do_fit=False):
    model = linear_model.Ridge(alpha=0.01)
    if do_fit:
        model.fit(x, y)
    return model

def mlp(x,y, do_fit = False):
    model = MLPRegressor(
        alpha = 1,
        hidden_layer_sizes = (1000, 100))
    
    if do_fit:
        model.fit(x, y)
        
    return model

def xgb(x,y, do_fit = False):
    model = xgboost.XGBRegressor(
        n_estimators=100,
        max_depth=4,
        learning_rate=0.1,
        reg_alpha= 0.00001,
        reg_lambda= 0.1,
        colsample_bytree=1,
        gamma=0.1,
        objective='reg:squarederror',
    )
    
    if do_fit:
        model.fit(x, y)
    
    return model

def random_forest(x,y, do_fit=False):
    model = RandomForestRegressor(
        max_depth=None, 
        random_state=0,
        n_estimators=1000,
        min_samples_leaf= 5,
    )
    if do_fit:
        model.fit(x, y) 
    return model

def zscore(rate, mean, std):
    """Calculates the Z-score from the mean and std."""
    zscore = (rate - mean) / std
    return zscore

def scale_zscore(zscore, mean, std):
    """Calculates the predicited insertion efficiency from the Z-score."""
    zscaled = zscore * std + mean
    return zscaled


#############################################################
###########    Set up train and test data sets    ########### 
#############################################################

# Read in measurements data
measurements = pd.read_csv(data_featurized)
traindata, testdata = splitdata_nostrat(measurements, fx = 'insertion', test_size = 0.3, seed = 0)
traindata.to_csv(datapath_train, index = False)
testdata.to_csv(datapath_test, index = False)

# #############################################################
# ###########              Modell setup             ########### 
# #############################################################

### Feature sets
feat_all = ['length', 'VF_RTT_z',  'percC', 'pairedbases','mmr', 'Arun_maxlen','percT', 'percA', 'percG', 'loop1_intact', 'loop2_intact', 'pos1compl', 
    'percGC', 'smaller4', 'smaller7', 'smaller13', 'longer40', 'N1_A', 'N2_A', 'N3_A', 'N4_A', 'NN_A', 'N1_T', 'N2_T', 'N3_T', 'N4_T', 'NN_T',
    'N1_G', 'N2_G', 'N3_G', 'N4_G', 'NN_G', 'N1_C', 'N2_C', 'N3_C', 'N4_C', 'NN_C', 'Arun', 'Crun', 'Trun', 'Grun', 'Trun_maxlen', 'Crun_maxlen', 'Grun_maxlen',
    'nick2match', 'nick3match', 'align_InsSpc', 'maxmicrohom_RTT', 'VF_insert', 'VF_RTT', 'VF_insert_z', 'Tm_insert', 'Tm_insert_z']

featset_d = [
    ['length'], # length
    ['VF_RTT_z','pairedbases','pos1compl',  'loop1_intact'], # structural features
    ['percC', 'percA', 'percT'], # nucleotide composition
    ['mmr','Arun_maxlen'], # system features (= dropped sequence features)
    ['length', 'percC', 'percA', 'percT', 'VF_RTT_z', 'pairedbases','pos1compl',  'loop1_intact'], # sequence features (= dropped system features)
    ['length', 'VF_RTT_z', 'mmr', 'percC', 'pairedbases', 'Arun_maxlen',  'percA', 'percT','pos1compl', 'loop1_intact'], # all core features
    ['length', 'VF_RTT_z', 'mmr', 'percC', 'pairedbases', 'Arun_maxlen',  'percA', 'percT','pos1compl', 'loop1_intact', 'N1_T', 'N1_A','N1_G','N1_C',], # all core features + first NT
    ['length', 'VF_RTT_z',  'percC', 'pairedbases','mmr', 'Arun_maxlen','percT', 'percA', 'percG', 'loop1_intact', 'loop2_intact', 'pos1compl', 
    'percGC', 'smaller4', 'smaller7', 'smaller13', 'longer40', 'N1_A', 'N2_A', 'N3_A', 'N4_A', 'NN_A', 'N1_T', 'N2_T', 'N3_T', 'N4_T', 'NN_T',
    'N1_G', 'N2_G', 'N3_G', 'N4_G', 'NN_G', 'N1_C', 'N2_C', 'N3_C', 'N4_C', 'NN_C', 'Arun', 'Crun', 'Trun', 'Grun', 'Trun_maxlen', 'Crun_maxlen', 'Grun_maxlen',
    'nick2match', 'nick3match', 'align_InsSpc', 'maxmicrohom_RTT',  'Tm_insert', 'Tm_insert_z'], # all features, but 'VF_RTT', to avoid that the model can guess the target site directly (VF_RTT_z is used instead)
]
featnames = ['Length', 'Structure', '%N', 'System', 'Sequence','Model v1','Model v1 + 1st nt','all']

feat_core = ['length', 'VF_RTT_z', 'mmr', 'percC', 'pairedbases', 'Arun_maxlen',  'percA', 'percT','pos1compl', 'loop1_intact']


#############################################################
###########         Find good model               ########### 
#############################################################

#Splits for kfold splitting for modeleing
Xseq = list(set(traindata['insertion']))
split_num = 10
kf = sklearn.model_selection.KFold(n_splits=split_num, random_state=setseed, shuffle=True)

# # Compare differend model architectures with all features
# logging.info(f'Model comparison is running...')
# cross_scores = {}
# model_list = ['L1', 'L2', 'MLP', 'RandomForest', 'XGBoost'] 
# featset =  feat_all
# # Calculate correlation scores for all models
# for m in model_list:
#     rscores = []
#     # We split the unique sequences into kfolds and then assign their corresponding data
#     for train_index, test_index in kf.split(Xseq): 
#         # list of test and train sequences
#         trainseq = [Xseq[i] for i in train_index]
#         testseq = [Xseq[i] for i in test_index]
#         # assign data points to test and train
#         X_train, X_test = traindata[traindata['insertion'].isin(trainseq)][featset], traindata[traindata['insertion'].isin(testseq)][featset] 
#         y_train, y_test = traindata[traindata['insertion'].isin(trainseq)]['percIns_z'], traindata[traindata['insertion'].isin(testseq)]['percIns_z']
#         # train and test the model
#         model = call_model(m, x = X_train, y = y_train, do_fit = True)
#         pred = model.predict(X_test)
#         corr, _ = stats.pearsonr(pred, y_test.values.ravel())
#         rscores.append(corr)
#     cross_scores[m]= rscores
# pd.DataFrame.from_dict(cross_scores).stack().reset_index().rename(columns = {'level_0': 'model_id', 'level_1': 'model', 0: 'score'}).to_csv(model_crossvalidation_path, index=False)
# logging.info(f'Model architectures were computed and output is stored in {model_crossvalidation_path}')


### Feature selection: add sets of features based on discoveries about sequence features, system features
scores = {k: [] for k in featnames}
for train_index, test_index in kf.split(Xseq):
    # list of test and train sequences
    trainseq = [Xseq[i] for i in train_index]
    testseq = [Xseq[i] for i in test_index]
    for f in range(len(featnames)):
        feat_temp = featset_d[f]
        # assign data points to test and train
        X_train, X_test = traindata[traindata['insertion'].isin(trainseq)][feat_temp], traindata[traindata['insertion'].isin(testseq)][feat_temp]
        y_train, y_test = traindata[traindata['insertion'].isin(trainseq)]['percIns_z'], traindata[traindata['insertion'].isin(testseq)]['percIns_z']
        # train and test the model
        model = call_model('XGBoost', x=X_train, y=y_train, do_fit=True)
        pred = model.predict(X_test)
        corr, _ = stats.pearsonr(pred, y_test.values.ravel())
        scores[featnames[f]].append(corr)
cross_scores_R = pd.DataFrame.from_dict(scores, orient = 'index').stack().reset_index().rename(columns = {'level_0': 'feature', 'level_1': 'feature_id', 0: 'score'})
cross_scores_R.to_csv(model_rmvfeatures_path, index=False)
logging.info(f'Evaluating different feature sets is finished')

### Train on one target site, test on another target site
featset = ['length', 'VF_RTT_z', 'mmr', 'percC', 'pairedbases', 'Arun_maxlen',  'percA', 'percT','pos1compl', 'loop1_intact']
feat_crossval_dict = defaultdict(list)
dsetOI =    ['HEK3 HEK293T', 'HEK3 HAP1dMLH1', 'HEK3 HAP1', 'FANCF HEK293T', 'FANCF HAP1dMLH1', 'FANCF HAP1', 'EMX1 HEK293T', 'CLYBL HEK293T'   ]

for dset in dsetOI: # dset is trainset
    Xseq_temp = list(set(traindata[traindata.axis_name.isin([dset])]['insertion']))
    for oset in dsetOI: # oset is testset
        corrscores = []         
        for train_index, test_index in kf.split(Xseq_temp): 
            # list of test and train sequences
            trainseq = [Xseq_temp[i] for i in train_index]
            testseq = [Xseq_temp[i] for i in test_index]
            # assign data points to test and train
            X_train, X_test = traindata[(traindata['insertion'].isin(trainseq)) & (traindata.axis_name.isin([dset]))][featset], traindata[(traindata['insertion'].isin(testseq)) & (traindata.axis_name.isin([oset]))][featset] 
            y_train, y_test = traindata[(traindata['insertion'].isin(trainseq)) & (traindata.axis_name.isin([dset]))]['percIns_z'], traindata[(traindata['insertion'].isin(testseq)) & (traindata.axis_name.isin([oset]))]['percIns_z']
            # train and test the model
            model = call_model('XGBoost', x = X_train, y = y_train, do_fit = True)
            pred = model.predict(X_test)
            try:
                corr, _ = stats.pearsonr(pred, y_test.values.ravel())
            except:
                corr = np.NaN
            corrscores.append(corr)
        feat_crossval_dict[dset].append(np.mean(corrscores)) # key is trainset, value is list of mean crossscores for the test sets in the order of oset
crosssite_val_df = pd.DataFrame.from_dict(feat_crossval_dict, orient='index', columns=dsetOI) # row name is training set, column name is test data
site2site = []
for s in crosssite_val_df.columns: 
    for t in crosssite_val_df.columns: 
        site2site.append([s, t, crosssite_val_df[s][t]])
pd.DataFrame(site2site, columns = ['trainset', 'testset', 'value']).to_csv(model_sites_crossvalidation_path, index=False)
logging.info(f'1 to 1 dataset crossvalidation is finished and output is stored in {model_sites_crossvalidation_path}')

# Also get the correlation between dsets
feat_crossval_dict = defaultdict(list)
for dset in dsetOI: # dset is trainset
    for oset in dsetOI: # oset is testset
        dset_temp = pd.merge(traindata[traindata.axis_name.isin([dset])][['insertion', 'percIns_z']], traindata[traindata.axis_name.isin([oset])][['insertion', 'percIns_z']], how='inner', on=['insertion'])
        corr, _ = stats.pearsonr(dset_temp['percIns_z_y'], dset_temp['percIns_z_x'])
        feat_crossval_dict[dset].append(corr) # key is trainset, value is list of mean crossscores for the test sets in the order of oset
crosssite_val_df = pd.DataFrame.from_dict(feat_crossval_dict, orient='index', columns=dsetOI) # row name is training set, column name is test data
site2site = []
for s in crosssite_val_df.columns: 
    for t in crosssite_val_df.columns: 
        site2site.append([s, t, crosssite_val_df[s][t]])
pd.DataFrame(site2site, columns = ['trainset', 'testset', 'value']).to_csv(model_sites_crosscorr_path, index=False)
logging.info(f'1 to 1 dataset crossvalidation is finished and output is stored in {model_sites_crosscorr_path}')



# #############################################################
# ###########              Final Model              ########### 
# #############################################################

# Final model for future applications
X_train = traindata[feat_core]
y_train = traindata['percIns_z']
logging.info(f'Xtrain shape: {X_train.shape}')
logging.info(f'ytrain shape: {y_train.shape}')

model = call_model('XGBoost', x = X_train, y = y_train, do_fit = True)

# Save model
pickle.dump(model, open(model_path, 'wb'))
logging.info(f'Model is stored in {model_path}')

# ### Test model on held out dataset
testdata_model = testdata.copy().reset_index(drop = True)
testdata_model['percIns_z_pred'] = model.predict(testdata_model[feat_core])
corr, _ = stats.pearsonr(testdata_model['percIns_z_pred'], testdata_model['percIns_z'].values.reshape(-1))
logging.info(f'Pearson R = {corr:.3f}')
pd.DataFrame(testdata_model.to_csv(model_performance_path, index = False))
logging.info(f'Performance evaluation on testdata is finished and output is stored in {model_performance_path}')


#############################################################
###########        Model characterization         ########### 
#############################################################

#Load model
model = pickle.load(open(model_path, 'rb'))
featset = feat_core

## Get and plot shap values
shap_values = shap.TreeExplainer(model).shap_values(traindata[featset])
plt.figure()
shap.summary_plot(shap_values, traindata[featset], featset, show = False,  plot_size=(5,5))
plt.savefig(model_shap_impact_path, dpi = 300, transparent=True, bbox_inches='tight')
plt.close()
shap.summary_plot(shap_values, features=traindata[featset], feature_names=featset, plot_type='bar',show = False, plot_size=(5,5))
plt.savefig(model_shap_importance_path, dpi = 300, transparent=True, bbox_inches='tight')
plt.close()

#############################################################
###########          Model application            ###########  
#############################################################

# Use for prediction on new target sites: Canadian goose screen
logging.info('Analyzing canadian goose screen data for 6 target sites in HEK3')
data_canadian = pd.read_csv(data_canadian_featurized)
data_canadian['percIns_z_pred'] = model.predict(data_canadian[featset])
corr, _ = stats.pearsonr(data_canadian['percIns_z_pred'], data_canadian['percIns_z'])
data_canadian.to_csv(model_performance_ca_path)
logging.info(f'Performance evaluation (R = {corr:.3f}) on canadian goose testdata is finished and output is stored in {model_performance_ca_path}')

# Use for prediction on new target sites and fusion proteins: Barnacle screen
logging.info('Analyzing barnacle goose screen data for tags in many target sites')
data_barnacle = pd.read_csv(data_barnacle_featurized)
data_barnacle['percIns_z_pred'] = model.predict(data_barnacle[featset])
corr, _ = stats.pearsonr(data_barnacle['percIns_z_pred'], data_barnacle['percIns_z'])
logging.info(f"Correlation of prediction and data for barnachle: {corr}")
data_barnacle.to_csv(model_performance_ba_path)
logging.info(f'Performance evaluation (R = {corr:.3f}) on barnacle goose testdata is finished and output is stored in {model_performance_ba_path}')


#############################################################
###########          Others                ########### 
#############################################################

##########         Feature correlation
corrMatrix = cluster_corr(traindata[feat_all].corr())
corrMatrix.to_csv(path_corrmatrix_all, index=False)
logging.info(f'Feature correlation is finished and output is stored in {path_corrmatrix_all}')

corrMatrix = cluster_corr(traindata[feat_core].corr())
corrMatrix.to_csv(path_corrmatrix_core, index=False)
logging.info(f'Feature correlation is finished and output is stored in {path_corrmatrix_core}')

#############################################################
###########         Validation datasets           ###########
#############################################################

# Choi data
choidata_eN3 = pd.read_csv(data_choi_featurized_eN3)
choidata_eN3['percIns_z_pred'] = model.predict(choidata_eN3[feat_core])
choidata_eN3.to_csv(model_performance_choi_eN3_path, index = False)

choidata_eN6 = pd.read_csv(data_choi_featurized_eN6)
choidata_eN6['percIns_z_pred'] = model.predict(choidata_eN6[feat_core])
choidata_eN6.to_csv(model_performance_choi_eN6_path, index = False)

# DeepPE data
deeppedata = pd.read_csv(data_deeppe_featurized)
deeppedata['percIns_z_pred'] = model.predict(deeppedata[feat_core])


deeppedata.to_csv(model_performance_deeppe_path, index = False)
