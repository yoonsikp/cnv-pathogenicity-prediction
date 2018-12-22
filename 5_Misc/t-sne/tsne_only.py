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


df = df.sample(frac=0.1, random_state = rand_state)
df = df.reset_index(drop=True)

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
'repeat_Simple repeat_density','repeat_LINE','repeat_Segmental duplication_density','size','gene_density','mpo_multi_entrez_to_num_phenotypes',
'mpo_multi_entrez_to_num_phenotypes_using_thresh','mpo_behavior/neurological phenotype',
'mpo_growth/size/body region phenotype','pli_pli_0.0_to_0.1','pli_pli_0.9_to_1.0','pli_pli_0.8_to_0.9',
'pli_pli_0.3_to_0.4','gain_loss','omim_num_diseases', 'number_of_genes_in_proximity']
df= df[to_be_scaled]


# 26500
# df= df[to_be_scaled][:10000]
# x_labels = x_labels[:][:10000]
# df.loc[:,'size'], lmbda = stats.boxcox(df['size'].copy())

if True:
    df.loc[:,'size'], lmbda = stats.boxcox(df['size'].copy())
    scaler = StandardScaler()
    scaler.fit(df)
    newdf = pd.DataFrame(scaler.transform(df), columns = df.columns)
else:
    newdf = df


exit()
# tsne
tsne = TSNE(n_jobs=4, verbose = 1, perplexity=750, learning_rate=100, n_iter= 800)
tsne_train = tsne.fit_transform(newdf)

finalDf = pd.concat([pd.DataFrame(tsne_train, columns=['tsne1','tsne2']), pd.DataFrame(x_labels, columns=['pathogenicity'])], axis = 1)
finalDf = pd.concat([finalDf, newdf['size']], axis = 1)
finalDf = pd.concat([finalDf, newdf['repeat_LINE_density']], axis = 1)
# plot
fig = plt.figure(figsize = ( 8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('tSNE X', fontsize = 15)
ax.set_ylabel('tSNE Y', fontsize = 15)
ax.set_title('tSNE of select features', fontsize = 20)

# pathogenic coloring
for target, color in zip([0,1], ['b', '#FF4500']):
    indicesToKeep = finalDf['pathogenicity'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'tsne1'] , finalDf.loc[indicesToKeep, 'tsne2'], c = color, s = 4, alpha=0.4)

ax.legend(['non-pathogenic', 'pathogenic'])
ax.grid()


plt.show()
exit()
# general coloring
cmap = matplotlib.cm.get_cmap('plasma')
# normalize = matplotlib.colors.LogNorm(vmin=min(finalDf['size']), vmax=max(finalDf['size']))
normalize = matplotlib.colors.Normalize(vmin=min(finalDf['repeat_LINE_density']), vmax=max(finalDf['repeat_LINE_density']))

colors = [cmap(normalize(value)) for value in finalDf['repeat_LINE_density']]
ax.scatter(finalDf['tsne1'] , finalDf['tsne2'], c = colors, s = 4, alpha=0.4)
cax, _ = matplotlib.colorbar.make_axes(ax)
cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize)


ax.set_facecolor('#FFFFFF')
ax.legend(['non-pathogenic', 'pathogenic'])
ax.grid()
plt.show()