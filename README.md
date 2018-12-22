# Predicting the Pathogenicity of Copy Number Variations

Copy number variations (CNVs) describe a subset of the wide variety of genetic modifications that occur in humans. However, it remains difficult for researchers to predict the effects a CNV will have on an individual. CNVs exhibit a spectrum of phenotypic effects ranging from benign to pathogenic to even beneficial. This project aims to detect pathogenic CNVs, while safely discarding CNVs that are confidently predicted to be benign.

This repository contains the code and datasets required to replicate the results of the project. Furthermore, the libraries used for Feature Extraction can be repurposed for any project involving regions of genetic data aligned to the hg19 reference genome. Every top level folder contains a descriptive `README`. The following links are example notebooks from the project.

[Feature Extraction](2_Feature_Extraction/feature_creation.ipynb)

[Model Training (XGBoost)](3_Model_Training/classifier_k_fold_xgboost_all_features_and_feature_importance.ipynb)

[Model Testing (XGBoost)](4_Model_Testing/4_Model_Testing/classifier_k_fold_xgboost_all_features_and_feature_importance.ipynb)

### Presentation

[Final Presentation Excerpt](https://docs.google.com/viewer?url=https://github.com/yoonsikp/cnv-pathogenicity-prediction/raw/master/BCB330_Presentation_Excerpt.pdf)

### Requirements
This project depends on Python 3. The Python 3 libraries needed are listed in `requirements.txt`. 
