import core.logistic_regression.logistic_regression as lg
import pandas as pd
from sklearn.model_selection import train_test_split

def init_logistic_regression(data_path: str, save_path: str, train_size=0.2):
    data = pd.read_csv(data_path)
    X, y = lg.init_data(data)
    X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=train_size)
    logreg = lg.get_fit_logistic_regression(X_train, y_train)
    lg.save_file(logreg, save_path)  
    return lg.score_logistic_regression(logreg, X_test, y_test)
    
def edge_predict(model_path: str, edge_features_path: str):
    logreg = lg.load_logistic_regression(model_path)
    edge_data = pd.read_csv(edge_features_path)
    edge_X, _ = lg.init_data(edge_data)
    return lg.predict_logistic_regression(logreg, edge_X)

def edge_set_score(model_path: str, edge_set_features_path: str):
    logreg = lg.load_logistic_regression(model_path)
    edges_data = pd.read_csv(edge_set_features_path)
    edges_X, edges_y = lg.init_data(edges_data)
    return lg.score_logistic_regression(logreg, edges_X, edges_y)
    
def save_file(data, save_path: str):
    lg.save_file(data, save_path)
    
def raw_data_processing(row_data_path: str, data_path: str):
    file = open(row_data_path)
    res = []
    n = 0
    edges = set()
    h = file.readline()
    while(len(h) > 0):
        h = h.split()
        if len(h) == 3:
            a, b, time = h[0], h[1], h[2]
        else:
            a, b, time = h[0], h[1], h[3]
        a = int(a)
        b = int(b)
        time = float(time)
        if (a > b):
            a, b = b, a
        if (a == b or (a*1000000000) + b in edges or time == 0):
            h = file.readline()
            continue
        n = max(n, b)
        edges.add((a*1000000000) + b)
        res.append([a, b, time])
        h = file.readline()
        
    file.close()
    
    res = sorted(res, key = lambda x: x[1])
    res = sorted(res, key = lambda x: x[0])
    
    file = open(data_path, "w")
    file.write(f"{n} {len(res)}\n")
    
    for i in res:
        file.write(f"{i[0]} {i[1]} {i[2]}\n")
        
                
