# -*- coding: utf-8 -*-
"""new_feature_creation.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/11UcvGyRWukMPwtrXWxkJbaQeYuf3UekE

# Set up Environment
Update Python and install prerequisite software. Most of this code can be skipped if you are running the Notebook locally, simply make sure you have all the right libraries. I use the following code cell to mount a folder from Google Drive.
"""

!apt-get install -y -qq software-properties-common python-software-properties module-init-tools > /dev/null
!add-apt-repository -y ppa:alessandro-strada/ppa 2>&1 > /dev/null
!apt-get update -qq 2>&1 > /dev/null
!apt-get -y install -qq google-drive-ocamlfuse fuse > /dev/null
from google.colab import auth
auth.authenticate_user()
from oauth2client.client import GoogleCredentials
creds = GoogleCredentials.get_application_default()
import getpass
!google-drive-ocamlfuse -headless -id={creds.client_id} -secret={creds.client_secret} < /dev/null 2>&1 | grep URL
vcode = getpass.getpass()
!echo {vcode} | google-drive-ocamlfuse -headless -id={creds.client_id} -secret={creds.client_secret}

!mkdir -p drive
!google-drive-ocamlfuse drive
!pip install -q keras
!pip install -q intervaltree

"""Change to the Google Drive learning directory, or your directory with the libraries and data files"""

import os
os.chdir("/content/drive/learning")

"""# Import Pandas and Numpy"""

import pandas
import numpy

"""# Set up filenames

Define input and ouput filenames for the files
"""

input_filename = "merged.csv"
output_filename = "output_features_07_26_everything.csv"
output_file = open(output_filename, "w", newline = '')

"""# Load CSV File"""

df = pandas.read_csv(input_filename, dtype=str, keep_default_na=False)

"""# Define Helper Function for splitting Entrez IDs"""

def split_entrez_ids(ids):
    entrez_ids = []
    if ids.strip() != '':
        for entrez_id in ids.strip().split(';'):
            entrez_ids.append(int(entrez_id))
    return entrez_ids

"""# Usage
Simply click the run button for every feature you would like added to the CSV. Once you are done adding the features, go to the last line, called Save CSV, and click the run button for that.

# Keep track of features added
"""

list_of_features_added = []

"""# MPO Feature"""

# import the mpo library from the ./libraries folder
from libraries import mpo

#the mpo database is in the ./libraries folder
mpo_class = mpo.Gene_to_MPO_Features(
    "./libraries/MPO_topPh_GP20180216.tsv")

# keep track of feature name for column naming
feature_name = 'mpo'

# warn if feature is already added!
if feature_name in list_of_features_added:
    print(__name__ + ": Warning: feature '" +
          feature_name + "' already added to dataframe")
list_of_features_added.append(feature_name)

# create dictionary of columns -> list of values for each row
final_mpo_dict = {}

# These lists represent one column each

# multi_entrez_to_num_phenotypes is the total sum of the phenotypes matched
list_multi_entrez_to_num_phenotypes = []

# multi_entrez_to_num_phenotypes_using_thresh is a more conservative sum
# where a gene that matches the same phenotype multitple times is counted only
# once.
list_multi_entrez_to_num_phenotypes_using_thresh = []

# create a key for each column in the dictionary
for column in mpo_class.get_all_phenotypes():
    final_mpo_dict[feature_name + '_' + column] = []

# iterate through every row
for ids in df['genes_in_proximity']:
    list_multi_entrez_to_num_phenotypes.append(
        mpo_class.multi_entrez_to_num_phenotypes(split_entrez_ids(ids)))
    list_multi_entrez_to_num_phenotypes_using_thresh.append(
        mpo_class.multi_entrez_to_num_phenotypes_using_thresh(
            split_entrez_ids(ids)))
    for key, value in mpo_class.multi_entrez_to_phenotypes_using_thresh(
        split_entrez_ids(ids)).items():
        final_mpo_dict[feature_name + '_' + key].append(value)

# add all columns from dictionary
df = pandas.concat([df, pandas.DataFrame.from_dict(final_mpo_dict)], axis=1)

# add new column
df[feature_name + "_multi_entrez_to_num_phenotypes"] = pandas.DataFrame(
    list_multi_entrez_to_num_phenotypes)

# add new column
df[feature_name +
   "_multi_entrez_to_num_phenotypes_using_thresh"] = pandas.DataFrame(
    list_multi_entrez_to_num_phenotypes_using_thresh)

"""# OMIM Feature"""

from libraries import omim
omim_class = omim.Gene_to_OMIM_Features(
    "./libraries/OMIMdiseasePhenotypeOnly_gene_info_GP20180213.tsv")
feature_name = 'omim'
if feature_name in list_of_features_added:
    print(__name__ + ": Warning: feature '" +
          feature_name + "' already added to dataframe")

list_of_features_added.append(feature_name)

# only one column, how many genes are listed in OMIM?
list_omim_num_diseases = []
for ids in df['genes_in_proximity']:
    list_omim_num_diseases.append(
        omim_class.multi_entrez_to_num_diseases(split_entrez_ids(ids)))
df[feature_name + "_num_diseases"] = pandas.DataFrame(list_omim_num_diseases)

"""# pLI Feature"""

# we need two libraries, one to convert from the entrez ID to the Gene Symbol
# and the other to actually give us the pLI of the gene(s)
from libraries import entrez_symbol
entrez_symbol_class = entrez_symbol.Entrez_Symbol_Lookup(
    "./libraries/map.gsymbol.to.enzid.tsv")
from libraries import gene_pli
pli_class = gene_pli.Gene_to_pLI_Features(
    "./libraries/pLI_EXac_broad_institure_2016_03.txt")
feature_name = 'pli'
if feature_name in list_of_features_added:
    print(__name__ + ": Warning: feature '"+
          feature_name + "' already added to dataframe")
list_of_features_added.append(feature_name)

# the pLI values are binned to be quite granular
columns = ["pli_0.0_to_0.1", "pli_0.1_to_0.2", "pli_0.2_to_0.3",
           "pli_0.3_to_0.4", "pli_0.4_to_0.5", "pli_0.5_to_0.6",
           "pli_0.6_to_0.7", "pli_0.7_to_0.8", "pli_0.8_to_0.9",
           "pli_0.9_to_1.0"]

# once agains we keep a dictionary where each key/column -> list of values
final_pli_dict = {}
for column in columns:
    final_pli_dict[feature_name + '_' + column] = []
    
# iterate through every row
for ids in df['genes_in_proximity']:
    # create another temporary dictionary just for the current row,
    # so that every key/column -> a single value
    row_pli_dict = {}
    # initialize all values to 0
    for column in columns:
        row_pli_dict[column] = 0
    for entrez_id in split_entrez_ids(ids):
        # there is more than one symbol returned for a given Entrez ID ...
        list_symbols = entrez_symbol_class.entrez_id_to_symbols(entrez_id)
        max_pli = -2
        # find the max pLI of all the "same" gene symbols
        for symbol in list_symbols:
            max_pli = max(max_pli, pli_class.gene_symbol_to_pLI(symbol))
        # the library returns '-1' if the gene is not found, so if none of the
        # genes are found or no entrez ID translation existed, then max_pli < 0
        if max_pli < 0:
            continue
        elif max_pli <= 0.1:
            row_pli_dict["pli_0.0_to_0.1"] += 1
        elif max_pli <= 0.2:
            row_pli_dict["pli_0.1_to_0.2"] += 1
        elif max_pli <= 0.3:
            row_pli_dict["pli_0.3_to_0.4"] += 1
        elif max_pli <= 0.4:
            row_pli_dict["pli_0.3_to_0.4"] += 1
        elif max_pli <= 0.5:
            row_pli_dict["pli_0.4_to_0.5"] += 1
        elif max_pli <= 0.6:
            row_pli_dict["pli_0.5_to_0.6"] += 1
        elif max_pli <= 0.7:
            row_pli_dict["pli_0.6_to_0.7"] += 1
        elif max_pli <= 0.8:
            row_pli_dict["pli_0.7_to_0.8"] += 1
        elif max_pli <= 0.9:
            row_pli_dict["pli_0.8_to_0.9"] += 1
        elif max_pli <= 1.0:
            row_pli_dict["pli_0.9_to_1.0"] += 1
    # finally add the current row dictionary to the final dictionary
    for key, value in row_pli_dict.items():
        final_pli_dict[feature_name + '_' + key].append(value)
df = pandas.concat([df, pandas.DataFrame.from_dict(final_pli_dict)], axis=1)

"""# Repetitive Elements
Warning: creating this feature may take around an hour.
Even loading the library takes a few minutes
"""

from libraries import repeat_elem
repeat_class = repeat_elem.Gene_Interval_to_Repetitive_Elements(
    "./libraries/RLCRs_DNN-CNV.txt")
feature_name = 'repeat'
if feature_name in list_of_features_added:
    print(__name__ + ": Warning: feature '" +
          feature_name + "' already added to dataframe")
list_of_features_added.append(feature_name)

final_repeat_dict = {}
# make a column for every repetitive element, e.g. "LINE", "SINE"
for key in repeat_class.set_of_element_types:
    final_repeat_dict[feature_name + '_' + key] = []
for index, row in df.iterrows():
    # this process takes so long, we need a progress update
    if index % 100 == 0:
        print(" cnv:" + str(index), end = '')
        
    # add to the dictionary for every repetitive element
    for key, value in repeat_class.get_all_intersecting_elements(
                    row['chr'], int(row['start']), int(row['end'])).items():
        final_repeat_dict[feature_name + '_' + key].append(value)
df = pandas.concat([df, pandas.DataFrame.from_dict(final_repeat_dict)], axis=1)

"""# Pathways"""

from libraries import pathways
pathways_class = pathways.Gene_to_Pathways(
    "./libraries/allSizes_GOincludingIEA_pathways_20180213.GMT")
feature_name = 'pathways'
if feature_name in list_of_features_added:
    print(__name__ + ": Warning: feature '" +
          feature_name + "' already added to dataframe")
list_of_features_added.append(feature_name)

# Get all the Pathways Databases once
databases = pathways_class.get_all_databases()
# NCI is no longer supported and GO has too many genesets
databases.remove("NCI")
databases.remove("GO")

hugedf = pandas.DataFrame(numpy.zeros((len(df), 1)), columns = ['drop'])
for database in databases:
    for line_num in pathways_class._dict_databases_to_lines[database]:
        hugedf[database + '_' + str(line_num)] = 0
    for index, row in df.iterrows():
        for key, descrip, value in pathways_class.multi_entrez_to_gene_sets(
                              split_entrez_ids(row['genes_in_proximity']),
                              database, filter_min = 30, filter_max = 1000):
            hugedf.at[index, database + '_' + str(key)] = str(value)
hugedf.drop(['drop'], axis=1)
df = pandas.concat([df, hugedf], axis=1)

"""# Save CSV"""

df.to_csv(output_filename)