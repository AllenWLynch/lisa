"""
This file defines the ML models used in Lisa. The class hierarchy is shown below (* = applied directly by Lisa):

BaseEstimator --> EstimatorInterface --> SampleSelectionModel --> LR_BinarySearch_SampleSelectionModel*
                            |
                            |                   
                            --> ChromatinModel --> LR_ChromatinModel*

The EstimatorInterface enforces "fit" and "get_info" functions, and the SampleSelectionModel and ChromatinModel enforce functions
for the selection of discriminative datasets, and the training of an accurate chromatin model, respectively.

Users may extend the SampleSelectionModel or ChromatinModel class to define new behavior for Lisa.
"""

from sklearn.base import BaseEstimator
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler 
from sklearn.feature_selection import SelectKBest, f_classif
import numpy as np
from sklearn.exceptions import ConvergenceWarning
import warnings
warnings.filterwarnings("ignore", category=ConvergenceWarning)
from lisa.core.assays import get_deltaRP_activation, transform_RP

class EstimatorInterface(BaseEstimator):
    #Wraps the sklearn Estimator class, and enforces that the user implement a "fit" and "get_info" method
    def fit(self):
        raise NotImplementedError("Must implement \"fit\" method of ChromatinModel")
    def get_info(self):
        raise NotImplementedError("Must implement \"get_info\" method of ChromatinModel")


class SampleSelectionModel(EstimatorInterface):
    """
    Extends the EstimatorInterface class, and enforces "get_selected_datasets" and "get_num_selected_datasets" methods
    get_selected_datasets should return a boolean index of the "selected status of each dataset in the data matrix
    """
    def get_selected_datasets(self):
        raise NotImplementedError()
    def get_num_selected_datasets(self):
        raise NotImplementedError()

class LR_BinarySearch_SampleSelectionModel(SampleSelectionModel):
    """
    Uses logistic regression to select n best discriminatory datasets. Uses binary search to tune the regularization parameter.
    max_iters: number of binary_search iterations to run, max
    num_datasets_selected: find best n datasets
    epsilon: tolerance. If dataset's regression coef is > epsilon, it is considered "selected"
    penalty: l1, use l1 penalty to choose datasets
    tol: tolerance, stop condition for LR model
    penalty_range: high and low boundaries for the C regularization parameter of the LR model.
        C = 2**-penalty, and C is inverse of regularization. Therefor for range [-1, 10], the minimum C = 2**1 = 2, very slightly regularized,
        and C = 2**-10, or very highly regularized. Binary search finds the optimal parameter in this space
    """
    def __init__(self, num_anova_features, num_datasets_selected, max_iters = 50, 
        epsilon = 1e-7, penalty = 'l1', tol = 0.01, penalty_range = (-1, 10)):
        self.max_iters = max_iters
        self.num_datasets_selected = num_datasets_selected
        self.epsilon = epsilon
        self.penalty = penalty
        self.tol = tol
        self.penalty_range = penalty_range
        self.num_anova_features = num_anova_features

    #binary search to close in on optimal penalty value
    def _binary_search(self, X, y, low, high, iter_num = 0):
        penalty = (high - low) / 2 + low
        self.model.C = 2**-penalty
        self.model.fit(X,y)
        #get upweighted datasets
        num_datasets_selected = self.get_num_selected_datasets()

        #break if the desired amount of datasets were used or max iters reached
        if num_datasets_selected == self.num_datasets_selected or iter_num == self.max_iters: 
            return self

        return self._binary_search(X, y, penalty, high, iter_num + 1) if num_datasets_selected > self.num_datasets_selected \
            else self._binary_search(X, y, low, penalty, iter_num + 1)


    #returns binary index of selected datasets
    def get_selected_datasets(self):
        return np.squeeze(np.abs(self.model.coef_) > self.epsilon)

    #returns the number of selected datasets
    def get_num_selected_datasets(self):
        return self.get_selected_datasets().sum()

    #instantiates a new LR model, then tunes the C parameter to get n datasets
    def fit(self, rp_matrix, labels):

        # Normalize with mean centering and log transformation
        index_array = np.arange(rp_matrix.shape[1])

        X = StandardScaler(with_std = False).fit_transform( transform_RP(rp_matrix) )

        #narrow features with anova selection
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            anova_featues = SelectKBest(f_classif, k=self.num_anova_features).fit(X, labels).get_support()
        #get indices mapping anova features to original features
        index_array = index_array[anova_featues]
        #subset selected features
        X = X[:, anova_featues]

        #instantiate LR search model
        self.model = LogisticRegression(penalty = self.penalty, tol = self.tol, solver = 'liblinear', random_state = 777)
        #use binary search technique to find datasets
        self._binary_search(X, labels, *self.penalty_range)
        #get boolean mask for selected datasets
        lr_selections = self.get_selected_datasets()
        #map LR selections back to original features as boolean mask
        dataset_selections = index_array[lr_selections]

        return dataset_selections

    #gets information on parameters chosen by user, and parameters of LR model to save to json logs
    def get_info(self):
        return dict(
            search_params = self.get_params(),
            search_model_params = self.model.get_params(),
            dataset_coefs = list(np.squeeze(self.model.coef_)),
        )


class ChromatinModel(EstimatorInterface):
    """
    Enforces that a chromtin model class implements get_deltaRP_activation
    """
    def get_deltaRP_activation(self, X):
        raise NotImplementedError()


class LR_ChromatinModel(ChromatinModel):

    def __init__(self, param_grid, kfold_cv = 4, scoring = 'roc_auc', penalty = 'l1'):
        self.kfold_cv = kfold_cv
        self.param_grid = param_grid
        self.scoring = scoring
        self.penalty = penalty

    def fit(self, rp_matrix, labels, n_jobs = 1):

        self.rp_0 = rp_matrix

        self.normalizer = StandardScaler(with_std = False)

        X0 = self.normalizer.fit_transform( transform_RP(self.rp_0) )

        #define a classifier for query vs background genes
        classifier = LogisticRegression(penalty = self.penalty, solver = 'lbfgs' if self.penalty == 'l2' else 'liblinear', random_state = 777)
        #optimize model with grid search
        self.grid_search = GridSearchCV(classifier, param_grid = self.param_grid, scoring = self.scoring, cv = self.kfold_cv, n_jobs = 1)\
            .fit(X0, labels)

        self.model = self.grid_search.best_estimator_

        return self

    def get_deltaRP_activation(self, rp_knockout):
        """
        rp_knockout: is a datacube of shape (genes, samples, TFs),
        this method must implement a transformation into a genes x TFs matrix, sumarrizing the dataset-axis effects
        """
        #subtract define deltaX to be the log2 of the fraction of knocked-out RP
        #deltaX = np.log2(self.rp_0[:,:,np.newaxis] - rp_knockout + PSEUDOCOUNT) - np.log2(self.rp_0[:,:,np.newaxis] + PSEUDOCOUNT)
        
        deltaX = get_deltaRP_activation(self.rp_0[:,:,np.newaxis], rp_knockout)
        
        #flip sign so that more knockout = more deltaR
        return deltaX.transpose(0,2,1).dot(self.model.coef_.reshape(-1))


    def get_info(self):
        return dict(
            model_params = self.model.get_params(),
            coefs = list(np.squeeze(self.model.coef_)),
            auc_score = self.grid_search.best_score_
        )