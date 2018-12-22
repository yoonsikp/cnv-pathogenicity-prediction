import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from xgboost import XGBClassifier
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold
from sklearn.metrics import accuracy_score
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler, LabelBinarizer
from sklearn.ensemble import RandomForestClassifier
from sklearn import linear_model
from sklearn.metrics import roc_curve, precision_recall_curve, auc, make_scorer, recall_score, accuracy_score, precision_score, confusion_matrix
from sklearn.decomposition import PCA
import gzip
import sklearn.manifold
from MulticoreTSNE import MulticoreTSNE as TSNE
from scipy import stats
import os
import tensorflow as tf
from tensorflow.contrib.tensorboard.plugins import projector
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import random
## Get working directory
PATH = os.getcwd()

## Path to save the embedding and checkpoints generated
LOG_DIR = PATH + '/project-tensorboard/log-1/'

plt.style.use("ggplot")

# random state
rand_state = 293

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


df = df.sample(frac=0.1, random_state = rand_state)
df = df.reset_index(drop=True)

# split up columns
x_labels = df['pathogenicity'].values
# create density columns
df['gene_density'] = df['number_of_genes_in_proximity'] / df['size'] * 100000
cols = [c for c in df.columns if c.lower()[:6] == 'repeat']
for col in cols:
    df[col + '_density'] = df[col] / df['size'] * 100000

to_be_scaled = ['repeat_LTR_density','repeat_SINE_density','repeat_LINE_density','repeat_Low complexity_density','repeat_Transposable element_density',
'repeat_Simple repeat_density','repeat_LINE','repeat_Segmental duplication_density','size','gene_density','mpo_multi_entrez_to_num_phenotypes',
'mpo_multi_entrez_to_num_phenotypes_using_thresh','mpo_behavior/neurological phenotype',
'mpo_growth/size/body region phenotype','pli_pli_0.0_to_0.1','pli_pli_0.9_to_1.0','pli_pli_0.8_to_0.9',
'pli_pli_0.3_to_0.4','gain_loss','omim_num_diseases', 'number_of_genes_in_proximity', 'pathogenicity']
df= df[to_be_scaled]
df['size_scaled'], lmbda = stats.boxcox(df['size'].copy())

out_df = df.copy()
df.drop(['size'], axis=1, inplace=True)
df.drop(['pathogenicity'], axis=1, inplace=True)

to_be_jittered = ['mpo_behavior/neurological phenotype',
'mpo_growth/size/body region phenotype','pli_pli_0.0_to_0.1','pli_pli_0.9_to_1.0','pli_pli_0.8_to_0.9',
'pli_pli_0.3_to_0.4', 'omim_num_diseases']
for col in to_be_jittered:
    for i in range(len(out_df[col])):
        out_df[col].iloc[i] = out_df[col].iloc[i]*random.uniform(0.99, 1.01)

out_df.to_csv("all_metadata.csv")
import csv
file = open("all_metadata.csv", 'r')
wfile = open("all_metadata.tsv", 'w')
csvreader = csv.reader(file)
csvwriter = csv.writer(wfile, delimiter='\t')
for row in csvreader:
    csvwriter.writerow(row)
wfile.close()
metadata = os.path.join(PATH, 'all_metadata.tsv')

# 26500
# df= df[to_be_scaled][:10000]
# x_labels = x_labels[:][:10000]
# df.loc[:,'size'], lmbda = stats.boxcox(df['size'].copy())

if True:
    scaler = StandardScaler()
    scaler.fit(df)
    newdf = pd.DataFrame(scaler.transform(df), columns = df.columns)
else:
    newdf = df

tf_data = tf.Variable(newdf)

with tf.Session() as sess:
    saver = tf.train.Saver([tf_data])
    sess.run(tf_data.initializer)
    saver.save(sess, os.path.join(LOG_DIR, 'tf_data.ckpt'))
    config = projector.ProjectorConfig()
    # One can add multiple embeddings.
    embedding = config.embeddings.add()
    embedding.tensor_name = tf_data.name
    # Link this tensor to its metadata(Labels) file
    embedding.metadata_path = metadata
    # Saves a config file that TensorBoard will read during startup.
    projector.visualize_embeddings(tf.summary.FileWriter(LOG_DIR), config)
