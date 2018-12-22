This folder contains the code used to train the classifiers. Each classifier is documented and commented inside an iPython notebook. Please open the notebooks either through GitHub, Google Colab, or Jupyter local installation

### List of Notebooks

* `classifier_k_fold_logit_all_features.ipynb`

  5-fold cross validation for logistic regression, using all features

* `classifier_k_fold_neural_all_features_final.ipynb`

  5-fold cross validation for neural networks, using all features

* `classifier_k_fold_neural_select_features_final.ipynb`

  5-fold cross validation for neural networks, using selected features

* `classifier_k_fold_xgboost_all_features_and_feature_importance.ipynb`

  5-fold cross validation for gradient boosted trees, using all features

* `classifier_k_fold_xgboost_select_features_and_feature_importance.ipynb`

  5-fold cross validation for gradient boosted trees, using selected features

### List of Files

* `classifier_k_fold_xgboost_select_features_grid_search.py`

  This file also uses a 5-fold cross validation. However, the goal is to use grid search to determine the optimal hyper paramaters for the `xgboost_select_features` model. The model did not improve by using grid search. The documentation is a bit spotty on this one, and it's here for reference.
