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
from MulticoreTSNE import MulticoreTSNE as TSNE
# pca
finalDf = pd.concat([principalDf_train, pd.DataFrame(y_train, columns=['pathogenicity'])], axis = 1)
fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)
targets = [0,1]
colors = ['r', 'g']
for target, color in zip(targets,colors):
    indicesToKeep = finalDf['pathogenicity'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'pathway_PCA 1'] , finalDf.loc[indicesToKeep, 'pathway_PCA 2'], c = color, s = 2, alpha=0.4)

ax.legend(['non-pathogenic', 'pathogenic'])
ax.grid()

# tsne
tsne = TSNE(n_jobs=4, verbose = 1, perplexity=40)
tsne_train = tsne.fit_transform(principalDf_train)

finalDf = pd.concat([pd.DataFrame(tsne_train, columns=['tsne1','tsne2']), pd.DataFrame(y_train, columns=['pathogenicity'])], axis = 1)
finalDf = pd.concat([finalDf, X_pre_train['size']], axis = 1)
finalDf = pd.concat([finalDf, X_pre_train['repeat_LINE_density']], axis = 1)
# plot
fig = plt.figure(figsize = ( 8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('tSNE X', fontsize = 15)
ax.set_ylabel('tSNE Y', fontsize = 15)
ax.set_title('tSNE of select features', fontsize = 20)
targets = [0,1]

# pathogenic coloring
for target, color in zip(targets,['b', '#FF4500']):
    indicesToKeep = finalDf['pathogenicity'] == target
    ax.scatter(finalDf.loc[indicesToKeep, 'tsne1'] , finalDf.loc[indicesToKeep, 'tsne2'], c = color, s = 4, alpha=0.4)

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