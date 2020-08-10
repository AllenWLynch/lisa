from sklearn.base import BaseEstimator
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LogisticRegression


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
    def __init__(self, max_iters = 50, num_datasets_selected = 10, 
        epsilon = 1e-7, penalty = 'l1', tol = 0.01, penalty_range = (-1, 10)):
        self.max_iters = max_iters
        self.num_datasets_selected = num_datasets_selected
        self.epsilon = epsilon
        self.penalty = penalty
        self.tol = tol
        self.penalty_range = penalty_range

    #binary search to close in on optimal penalty value
    def _binary_search(self, X, y, low, high, iter_num = 0):
        penalty = (high - low) / 2
        self.model.C = 2**-penalty
        self.model.fit(X,y)
        coefs = self.model.coef_
        #get upweighted datasets
        num_datasets_selected = self.get_num_selected_datasets()
        #break if the desired amount of datasets were used or max iters reached
        if num_datasets_selected == num_datasets_desired or iters == max_iters: 
            return self
        #if too many datasets used, shift window to increase penalty
        if num_datasets_selected > num_datasets:
            #start new search with moved window
            return self._binary_search(X, y, penalty, high, iter_num + 1)
        else:
            #decrease penalty if too few datasets used
            return self._binary_search(X, y, low, penalty, iter_num + 1)

    #returns binary index of selected datasets
    def get_selected_datasets(self):
        return self.model.ceof_ > self.epsilon

    #returns the number of selected datasets
    def get_num_selected_datasets(self):
        return self.get_selected_datasets.sum()

    #instantiates a new LR model, then tunes the C parameter to get n datasets
    def fit(self, X, y):
        self.model = LogisticRegression(penalty = penalty, tol = tol)
        return self._binary_search(X, y, *self.penalty_range)

    #gets information on parameters chosen by user, and parameters of LR model to save to json logs
    def get_info(self):
        return dict(
            search_params = self.get_params(),
            search_model_params = self.model.get_params(),
            dataset_coefs = self.model.coef_
        )


class ChromatinModel(EstimatorInterface):
    """
    Enforces that a chromtin model class implements get_deltaRP_activation
    """
    def get_deltaRP_activation(self, X):
        raise NotImplementedError()


class DeltaRegLogisticRegression(LogisticRegression, ChromatinModel):
    """
    Extends sklearn LR model with the get_deltaRP_activation function
    """
    def get_deltaRP_activation(self, X):
        """
        The regulatory score of a gene is R(X), where X is the regulatory potential of that gene in each of n selected datasets
        The delta regulatory score is the change in R(X) after insilico deletion of a TF, and the effects of ISD are Z.
        So deltaR = R(X - Z) - R(X).
        For logistic regression R(X) = c * X + b,
        so deltaR = c (X - Z) + b - (c * X + b) = cX - cZ + b - cX - b = - cZ
        """
        np.dot(X, self.model.best_estimator_.coef_.reshape((-1, 1))) # X : (genes x datasets) * C (datasets  x 1) = (genes x 1)


class LR_ChromatinModel(ChromatinModel):

    def __init__(self, param_grid, kfold_cv = 4, scoring = 'auc', **regression_kwargs):
        self.regression_kwargs = regression_kwargs
        self.kfold_cv = kfold_cv
        self.param_grid = param_grid
        self.scoring = scoring

    def fit(self, X, n_jobs = 1):
        #define a classifier for query vs background genes
        classifier = DeltaRegLogisticRegression(self.regression_kwargs)
        #optimize re
        self.grid_search = GridSearchCV(classifier, param_grid = self.param_grid, scoring = self.scoring, self.cv = kfold_cv, n_jobs = 1)
        self.model = self.grid_search.best_estimator_

    def get_deltaRP_activation(self, X):
        return self.model.get_deltaRP_activation(X)

    def get_info(self):
        return dict(
            regression_params = self.regression_kwargs,
            model_params = self.model.get_params(),
            coefs = self.model.coef_
        )
    
