
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.metrics import make_scorer, average_precision_score, roc_auc_score
from sklearn.feature_selection import SelectKBest, f_classif

import numpy as np
import pandas as pd
from scipy.stats import uniform


#ADD COVARIATES STUFF

class BinarySearch():

    def __init__(self, min_value, max_value, increment = 1):
        self.min_value = min_value
        self.max_value = max_value
        self.increment = increment
        self.update_midpoint()

    def update_midpoint(self):
        self.midpoint = (self.min_value - self.max_value) / 2

    def __next__(self, too_high):

        if too_high is None:
            return self.midpoint
        
        if too_high == True:
            self.min_value = self.midpoint + self.increment
        elif too_high == False:
            self.max_value = self.midpoint - self.increment
        else:
            raise AssertionError('too_high argument must be either True, False, or None')

        self.update_midpoint()
        return self.midpoint


LR_model_kwargs = {

}

def midpoint(low, high):
        return (high - low) / 2

def select_best_datasets_LR_model(reg_potential_data, labels, epsilon = 1e-7, num_datasets = 10, max_iters = 100, penalty_min = -1, penalty_max = 10):
    
    penalty = midpoint(penalty_min, penalty_max)

    for _ in range(max_iters):
        #define a model using the fixed parameters chosen in LR_model_kwargs
        coefs = LogisticRegression(C=2**penalty, **LR_model_kwargs).fit(reg_potential_data, labels).coef_
        #get upweighted datasets
        best_datasets = coefs > epsilon
        num_datasets_selected = best_datasets.sum()
        if num_datasets_selected == num_datasets:
            break
        else:
            if num_datasets_selected > num_datasets:
                penalty_max = midpoint
            else:
                penalty_min = midpoint
            penalty = midpoint(penalty_min, penalty_max)

    return best_datasets

    

            










