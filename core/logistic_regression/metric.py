from sklearn.metrics import classification_report
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
import matplotlib.pyplot as plt 
from sklearn.linear_model import LogisticRegression
import pandas as pd

def report(logreg: LogisticRegression, X: pd.DataFrame, y: pd.DataFrame):
    y_pred = logreg.predict(X)
    return classification_report(y, y_pred, zero_division=1)
    
def roc_auc_curve(logreg: LogisticRegression, X: pd.DataFrame, y: pd.DataFrame):
    logit_roc_auc = roc_auc_score(y, logreg.predict(X))
    fpr, tpr, thresholds = roc_curve(y, logreg.predict_proba(X)[:,1])
    fig = plt.figure()
    plt.rc("font", size=14)
    plt.plot(fpr, tpr, label='Logistic Regression (area = %0.2f)' % logit_roc_auc)
    plt.plot([0, 1], [0, 1],'r--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic')
    plt.legend(loc="lower right")
    return fig