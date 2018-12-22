import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold, cross_val_score
from sklearn.metrics import accuracy_score
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler, LabelBinarizer
from sklearn.ensemble import RandomForestClassifier
from sklearn import linear_model
import sklearn.metrics
from scipy import stats
from sklearn.metrics import roc_curve, precision_recall_curve, auc, make_scorer, recall_score, accuracy_score, precision_score, confusion_matrix
from sklearn.decomposition import PCA
import gzip
import random

plt.style.use("ggplot")

# random state
rand_state = 233

# cutoffs
lower_bound = 2000
upper_bound = 5000000

# name of the csv file
file = gzip.open('./data/output_features_07_26_no_genelists.csv.gz')
df = pd.read_csv(file, dtype=str)

# drop columns that don't shouldn't be features
df.drop(['genes_in_proximity','chr', 'start', 'end', 'Unnamed: 0', 'drop'], axis=1, inplace=True)
df.drop(['repeat_Other', 'repeat_Unknown'], axis=1, inplace=True)
df = df.astype(float)

# cutoffs
df = df[df["size"] > lower_bound]
df = df[df["size"] < upper_bound]

# fix pathogenicity number
df['pathogenicity'] = df['pathogenicity'].replace(-1, 0)

# split up columns
x_labels = df['pathogenicity'].values
df.drop(['pathogenicity'], axis=1, inplace=True)

# create density columns
df['gene_density'] = df['number_of_genes_in_proximity'] / df['size'] * 100000
cols = [c for c in df.columns if c.lower()[:6] == 'repeat']
for col in cols:
    df[col + '_density'] = df[col] / df['size'] * 100000

cols = [c for c in df.columns if c.lower()[:4] == 'bioc' or c.lower()[:4] == 'kegg' or c.lower()[:5] == 'react']
df = df.drop(cols,axis=1)
to_be_scaled = ['repeat_LTR_density','repeat_SINE_density','repeat_LINE_density','repeat_Low complexity_density','repeat_Transposable element_density',
'repeat_Simple repeat_density','repeat_LINE','repeat_Segmental duplication_density','size','gene_density','mpo_multi_entrez_to_num_phenotypes','mpo_multi_entrez_to_num_phenotypes_using_thresh','mpo_behavior/neurological phenotype',
'mpo_growth/size/body region phenotype','pli_pli_0.0_to_0.1','pli_pli_0.9_to_1.0','pli_pli_0.8_to_0.9',
'pli_pli_0.3_to_0.4','gain_loss','omim_num_diseases', 'number_of_genes_in_proximity']
#to_be_scaled = ['size','gene_density', 'number_of_genes_in_proximity',    'gain_loss',   'mpo_multi_entrez_to_num_phenotypes',   'mpo_behavior/neurological phenotype', 'repeat_LTR_density',    'omim_num_diseases',  'mpo_mortality/aging','pli_pli_0.3_to_0.4','pli_pli_0.0_to_0.1', 'pli_pli_0.9_to_1.0']

# to_be_scaled = df.columns


df= df[to_be_scaled]

param_test1 = {
    'max_depth':range(3,10,2),
    'min_child_weight':range(1,6,2)
}
param_test1 = {
    'max_depth':range(3,10,2)
}
scorers = {
    'precision_score': make_scorer(precision_score),
    'recall_score': make_scorer(recall_score),
    'accuracy_score': make_scorer(accuracy_score),
    'average_precision_score': make_scorer(sklearn.metrics.average_precision_score),
    'roc_auc_score': make_scorer(sklearn.metrics.roc_auc_score)
}
def adjusted_classes(y_scores, t):
    """
    This function adjusts class predictions based on the prediction threshold (t).
    Will only work for binary classification problems.
    """
    return [1 if y >= t else 0 for y in y_scores]

def precision_recall_threshold(p, r, thresholds, y_true, y_score, t=0.5):
    """
    plots the precision recall curve and shows the current value for each
    by identifying the classifier's threshold (t).
    """
    
    # generate new class predictions based on the adjusted_classes
    # function above and view the resulting confusion matrix.
    y_pred_adj = adjusted_classes(y_score, t)
    print(pd.DataFrame(confusion_matrix(y_true, y_pred_adj),
                       columns=['pred_neg', 'pred_pos'], 
                       index=['neg', 'pos']))


def newt_average_precision_score(y_true, y_score, average="macro",
                        sample_weight=None):
    def _binary_uninterpolated_average_precision(
            y_true, y_score, sample_weight=None):
        precision, recall, thresholds = precision_recall_curve(
            y_true, y_score, sample_weight=sample_weight)
        precision_recall_threshold(precision,recall,thresholds, y_true, y_score)
        # Return the step function integral
        # The following works because the last entry of precision is
        # guaranteed to be 1, as returned by precision_recall_curve
        return -np.sum(np.diff(recall) * np.array(precision)[:-1])

    return sklearn.metrics.base._average_binary_score(_binary_uninterpolated_average_precision,
                                 y_true, y_score, average,
                                 sample_weight=sample_weight)


best_scorer = make_scorer(newt_average_precision_score)

gsearch1 = GridSearchCV(estimator = XGBClassifier( learning_rate =0.1, max_depth=5,
 min_child_weight=1, gamma=0, subsample=0.8, colsample_bytree=0.8,
 objective= 'binary:logistic', nthread=1, scale_pos_weight=5, seed=rand_state), 
 param_grid = param_test1, scoring=best_scorer, n_jobs=4,iid=False, refit=best_scorer, cv=StratifiedKFold(n_splits=5,shuffle=True, random_state=rand_state))
gsearch1.fit(df, x_labels)
gsearch1.fit(df, x_labels)
print(gsearch1.grid_scores_)
print(gsearch1.best_params_, gsearch1.best_score_, gsearch1.best_estimator_)

exit()
# >>> print(gsearch1.best_params_, gsearch1.best_score_, gsearch1.best_estimator_)
# {'max_depth': 9, 'min_child_weight': 3}
param_test2 = {
    'max_depth':[8,9,10],
    'min_child_weight':[2,3,4]
}

gsearch2 = GridSearchCV(estimator = XGBClassifier( learning_rate =0.1, max_depth=5,
 min_child_weight=1, gamma=0, subsample=0.8, colsample_bytree=0.8,
 objective= 'binary:logistic', nthread=1, scale_pos_weight=5, seed=rand_state), 
 param_grid = param_test2, scoring='recall',n_jobs=4,iid=False, refit='recall', cv=StratifiedKFold(n_splits=5,shuffle=True, random_state=rand_state))
gsearch2.fit(df, x_labels)
print(gsearch2.grid_scores_)
print(gsearch2.best_params_, gsearch2.best_score_, gsearch2.best_estimator_)

# >>> print(gsearch2.best_params_, gsearch2.best_score_, gsearch2.best_estimator_)
# {'max_depth': 9, 'min_child_weight': 3}


param_test3 = {
    'gamma':[i/10.0 for i in range(0,5)]
}
gsearch3 = GridSearchCV(estimator = XGBClassifier( learning_rate =0.1, max_depth=9,
 min_child_weight=3, gamma=0, subsample=0.8, colsample_bytree=0.8,
 objective= 'binary:logistic', nthread=1, scale_pos_weight=5, seed=rand_state), 
 param_grid = param_test3, scoring='recall',n_jobs=4,iid=False, refit='recall', cv=StratifiedKFold(n_splits=5,shuffle=True, random_state=rand_state))
gsearch3.fit(df, x_labels)
print(gsearch3.grid_scores_)
print(gsearch3.best_params_, gsearch3.best_score_, gsearch3.best_estimator_)

# >>> print(gsearch3.best_params_, gsearch3.best_score_, gsearch3.best_estimator_)
# {'gamma': 0.0}


param_test4 = {
    'subsample':[i/10.0 for i in range(6,10)],
    'colsample_bytree':[i/10.0 for i in range(6,10)]
}
gsearch4 = GridSearchCV(estimator = XGBClassifier( learning_rate =0.1, max_depth=9,
 min_child_weight=3, gamma=0, subsample=0.8, colsample_bytree=0.8,
 objective= 'binary:logistic', nthread=1, scale_pos_weight=5, seed=rand_state), 
 param_grid = param_test4, scoring='recall',n_jobs=4,iid=False, refit='recall', cv=StratifiedKFold(n_splits=5,shuffle=True, random_state=rand_state))
gsearch4.fit(df, x_labels)
print(gsearch4.grid_scores_)
print(gsearch4.best_params_, gsearch4.best_score_, gsearch4.best_estimator_)


# >>> print(gsearch4.best_params_, gsearch4.best_score_, gsearch4.best_estimator_)
# {'colsample_bytree': 0.8, 'subsample': 0.8} 

# XGBClassifier(base_score=0.5, booster='gbtree', colsample_bylevel=1,
#        colsample_bytree=0.8, gamma=0, learning_rate=0.1, max_delta_step=0,
#        max_depth=9, min_child_weight=3, missing=None, n_estimators=100,
#        n_jobs=1, nthread=1, objective='binary:logistic', random_state=0,
#        reg_alpha=0, reg_lambda=1, scale_pos_weight=5, seed=233,
#        silent=True, subsample=0.8)




