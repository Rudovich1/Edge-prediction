import pandas as pd
from sklearn.linear_model import LogisticRegression
import pickle
import core.logistic_regression.metric as metric
import matplotlib.pyplot as plt 

def init_data(data: pd.DataFrame):
    y = data['predict'].astype(int)
    X = data.drop(columns=['predict', 'time_init'])
    return X, y
    
def get_fit_logistic_regression(X: pd.DataFrame, y: pd.Series):
    logreg = LogisticRegression()
    logreg.fit(X, y)
    return logreg

def save_file(data, save_path: str):
    if (type(data) == plt.Figure):
        plt.savefig(save_path)
    else:
        with open(save_path, 'wb') as file:
            pickle.dump(data, file)
    
def load_logistic_regression(load_path: str):
    with open(load_path, 'rb') as file:
        logreg = pickle.load(file)
    return logreg

def score_logistic_regression(logreg: LogisticRegression, X: pd.DataFrame, y: pd.DataFrame):
    return metric.report(logreg, X, y), metric.roc_auc_curve(logreg, X, y)

def predict_logistic_regression(logreg: LogisticRegression, edge_features: pd.DataFrame):
    return bool(logreg.predict(edge_features)[0])

