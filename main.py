import core.logistic_regression.interface as inter
import subprocess


def test():
    extraction_exe = "core/feature_extraction/extraction.exe"
    edge_features_exe = "core/feature_extraction/edge_features.exe"
    random_edge_features_exe = "core/feature_extraction/random_edge_features.exe"
    file_name = "test"
    
    dir_path = f"data/{file_name}"    
    subprocess.call(["mkdir", dir_path])
    
    edges_path = dir_path + "/row_data.txt"
    features_path = dir_path + "/features.txt"
    data_size = "100"
    feature_set = "f"

    res = subprocess.call([extraction_exe, edges_path, features_path, data_size, feature_set])
 
    if (res == 0):
        print("\nSuccessful extraction!\n")
    else:
        print(f"\nExtraction error: {res}\n")
        return
    
    model_path = dir_path + "/model.pkl"
    model_res, model_fig = inter.init_logistic_regression(features_path, model_path)
    
    model_fig_path = dir_path + "/fit_fig.png"
    inter.save_file(model_fig, model_fig_path)
    
    print(f"\nModel accuracy:\n {model_res}\n")
    
    edge = ["2", "4"]
    t1 = "1"
    t2 = "20"
    
    edge_features_path = dir_path + f"/edge_{edge[0]}_{edge[1]}_{t1}_{t2}.txt"
    
    res = subprocess.call([edge_features_exe, edges_path, edge_features_path, feature_set, edge[0], edge[1], t1, t2])
    
    if (res == 0):
        print("\nSuccessful extraction!\n")
    else:
        print(f"\nExtraction error: {res}\n")
        return
        
    predict = inter.edge_predict(model_path, edge_features_path)
    
    if (predict):
        print(f"\nEdge [{edge[0]}, {edge[1]}] will appear in graph in {t2}\n")
    else:
        print(f"\nEdge [{edge[0]}, {edge[1]}] will not appear in graph in {t2}\n")
        
    dataset_size = "45"
    dataset_path = dir_path + f"/feature_set_{dataset_size}_{t1}_{t2}.txt"
    
    res = subprocess.run([random_edge_features_exe, edges_path, dataset_path, feature_set, dataset_size, t1, t2]).returncode
    
    if (res == 0):
        print("Successful extraction!")
    else:
        print(f"Extraction error: {res}")
        return
        
    model_score, model_fig = inter.edge_set_score(model_path, dataset_path)
    model_fig_path = dir_path + f"/feature_set_{data_size}_{t1}_{t2}_fig.png"
    inter.save_file(model_fig, model_fig_path)
    
    print (f"\nScore'\n: {model_score}\n")

def __main__():
    
    extraction_exe = "core/feature_extraction/extraction.exe"
    extraction_DEBUG_exe = "core/feature_extraction/extraction_DEBUG.exe"
    edge_features_exe = "core/feature_extraction/edge_features.exe"
    edge_feature_set_exe = "core/feature_extraction/edge_feature_set.exe"
    
    while (True):
        print(
"""
_________________________________

0 - Close
1 - Raw data processing
2 - Fit model
3 - Predict edge
4 - Accuracy of model on random edges
_________________________________

""")
        operation = int(input())
        try:
            if (operation == 0):
                break
                
            if (operation == 1):
                file_name = input("Model name: ")
                row_data_path = input("Row data path: ")
                
                dir_path = f"data/{file_name}"
                subprocess.call(["mkdir", dir_path])
                data_path = dir_path + "/edges_data.txt"
                
                inter.raw_data_processing(row_data_path, data_path)
                
            if (operation == 2):
                file_name = input("Model name: ")
                data_size = input("Number of edges to process: ")
                percent_init = input("Percent: ")
                feature_set = input("Full \\ Lite \\ Predict [f \\ l \\ p]: ")
                print()
                
                dir_path = f"data/{file_name}"
                subprocess.call(["mkdir", dir_path])
                
                edges_path = dir_path + "/edges_data.txt"
                features_path = dir_path + "/features.txt"
                
                res = subprocess.call([extraction_exe, edges_path, features_path, data_size, feature_set, percent_init])
                    
                if (res == 0):
                    print("\nSuccessful extraction!\n")
                else:
                    print(f"\nExtraction error: {res}\n")
                    continue
                
                model_path = dir_path + "/model.pkl"
                model_res, model_fig = inter.init_logistic_regression(features_path, model_path)
                
                fig_path = dir_path + "/fit_fig.png"
                inter.save_file(model_fig, fig_path)
                
                print(f"\nModel accuracy:'\n {model_res}\n")
                
            if (operation == 3):
                file_name = input("Model name: ")
                edge = input("Edge: ").split()
                t1 = input("Time 1: ")
                t2 = input("Time 2: ")
                feature_set = input("Full \\ Lite [f \\ l]: ")
                
                dir_path = f"data/{file_name}"
                edges_path = dir_path + "/edges_data.txt"
                edge_features_path = dir_path + f"/edge_{edge[0]}_{edge[1]}_{t1}_{t2}.txt"
                
                res = subprocess.call([edge_features_exe, edges_path, edge_features_path, feature_set, edge[0], edge[1], t1, t2])
                
                if (res == 0):
                    print("\nSuccessful extraction!\n")
                else:
                    print(f"\nExtraction error: {res}\n")
                    continue
                
                model_path = dir_path + "/model.pkl"
                predict = inter.edge_predict(model_path, edge_features_path)
                                
                if (predict):
                    print(f"\nEdge [{edge[0]}, {edge[1]}] will appear in graph in {t2}\n")
                else:
                    print(f"\nEdge [{edge[0]}, {edge[1]}] will not appear in graph in {t2}\n")
                    
            if (operation == 4):
                file_name = input("Model name: ")
                size = input("Number of edges to process: ")
                percent_init = input("Percent 1: ")
                percent_predict = input("Percent 2: ")
                feature_set = input("Full \\ Lite \\ Predict [f \\ l \\ p]: ")
                
                dir_path = f"data/{file_name}"
                edges_path = dir_path + "/edges_data.txt"
                dataset_path = dir_path + f"/feature_set_{percent_init}_{percent_predict}.txt"
                
                res = subprocess.call([edge_feature_set_exe, edges_path, dataset_path, feature_set, size, percent_init, percent_predict])
                
                if (res == 0):
                    print("Successful extraction!")
                else:
                    print(f"Extraction error: {res}")
                    continue
                
                model_path = dir_path + "/model.pkl"
                model_score, model_fig = inter.edge_set_score(model_path, dataset_path)
                
                fig_path = dir_path + f"/feature_set_{percent_init}_{percent_predict}_fig.png"
                inter.save_file(model_fig, fig_path)
                
                print (f"\nScore: {model_score}\n")
        except:
            continue

# test()

__main__()