This folder contains the code used to create features based on both the Entrez IDs that were determined in a previous step, as well as the chromosome/start/end locations.

### List of Files
* `feature_creation.ipynb`

    This is an iPython notebook containing the full code required to extract features from the CNV .csv file. The iPython notebook explains what is going on. It is highly recommended to run the .ipynb notebook using something like Jupyter, since each code cell/feature can then be selectively run.

* `feature_creation.py`

    The raw Python code from the notebook is still included in the folder as .py file.

* `output\output_features_07_26_everything.csv.gz`

    The output from feature extraction, including the very large gene_lists feature. The reason it is gzipped is because the file size is so big (242 MB)

* `output\output_features_07_26_no_genelists.csv.gz`

    The output from feature extraction, without the gene lists feature, so this is the file that you probably want. This is only a few MB in size.

### Libraries
The following libraries have documentation, and are the "meat" of the project:

* `libraries\entrez_symbol.py`

    Given either a certain human entrez ID or Gene Symbol, this module will translate to the other gene identifier. This module depends on `libraries\map.gsymbol.to.enzid.tsv`

* `libraries\gene_pli.py`

    Given a gene symbol, this module will return the pLI. This module is extendable to return other columns from the pLI file. This module depends on `libraries\pLI_EXac_broad_institure_2016_03.txt`

* `libraries\mpo.py`

    Given a certain human entrez ID, this module will return MPO features. This module also supports lookups with multiple human entrez IDs. This module depends on `libraries\MPO_topPh_GP20180216.tsv`

* `libraries\omim.py`

    Given a certain human entrez ID, this module will return OMIM features. This module also supports lookups with multiple human entrez IDs. This module depends on `libraries\OMIMdiseasePhenotypeOnly_gene_info_GP20180213.tsv`

* `libraries\pathways.py`

    Given a certain human entrez ID, the module will return all pathways. This module also supports returning pathways for multiple entrez ids, and also supports filtering the pathways based on gene numbers. This module depends on `libraries\allSizes_GOincludingIEA_pathways_20180213.GMT`

* `libraries\repeat_elem.py`

    Given a certain chromosome, start, and end location, this module will return repetitive element features. The repetitive element file is by default assumed to be compressed. This module depends on `libraries\RLCRs_DNN-CNV.txt.gz`
