import csv
import pandas as pd
from numpy import *
def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

# random state
rand_state = 23

# cutoffs
lower_bound = 2000
upper_bound = 5000000


# name of the csv file
df = pd.read_csv('merged_plus_features.csv')
df.drop(['genes_in_proximity', 'start', 'end'], axis=1, inplace=True)
df.drop(['chr'], axis=1, inplace=True)
df = df.astype(float)
df = df[df["size"] > lower_bound]
df = df[df["size"] < upper_bound]

#split up columns
x_labels = df['pathogenicity'].values
df.drop(['pathogenicity'], axis=1, inplace=True)
#print(df["gene_density"].values)
#print(x_labels)
for column in df.columns:
	print(column)
	print(corrcoef(df[column],x_labels)[1,0])









