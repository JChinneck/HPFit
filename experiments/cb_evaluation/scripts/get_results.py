import pandas as pd
import numpy as np
import scipy.stats
import glob
from sklearn.linear_model import LinearRegression
from sklearn import metrics
loc="/home/jpbrooks/HPFit/experiments/cb_evaluation/results"

folnames = ["bm-like", "bm-nox", "clustered_outliers", "olive", "rvd-like", "unclustered_outliers", "bm", "rvd", "clustered_outliers_small"]
reg_names = ["bm-like", "bm-nox", "olive", "rvd-like", "bm", "rvd"]
gen_names = ["clustered_outliers", "unclustered_outliers", "clustered_outliers_small"]
experiments=["evaluation"]
for experiment in experiments:
    writer = pd.ExcelWriter(loc + "/" + experiment + ".xlsx", mode="w")
    
    for folname in folnames:
        results={}
        if folname in gen_names:
            non_outlier = pd.read_csv("../data/" + folname + "/non_outlier_euclid_error.csv",
                                      index_col=0) # may need to add comma at end of first line
        for fname in glob.glob(loc+"/"+experiment+"/"+folname+"/*.csv"):
            f=open(fname, "r")
            line=f.readline()
            line1=line.split(",") 
            dataset=line1[0]
              
            if len(line1) > 5:
                formulation=line1[6]
                if not dataset in results:
                    results[dataset] = {}
                results[dataset][formulation + " Dnon"]=float(line1[7])
                if formulation in ["alg3PCA", "cb"]:
                    results[dataset][formulation + " runtime"]=float(line1[9])
                    results[dataset][formulation + " tsestar"]=float(line1[10])
                elif formulation in ["arob", "bb", "hbreg", "lm", "lts", "lqs"]:
                    results[dataset][formulation + " runtime"]=float(line1[11])
                    results[dataset][formulation + " gamma"]=float(line1[8])
                    results[dataset][formulation + " tsestar"]=float(line1[12])
                if formulation=="cb":
                    results[dataset]["i"]=int(line1[1])
                    results[dataset]["m"]=int(line1[2])
                    results[dataset]["n"]=int(line1[3])
                    results[dataset]["m_normal"]=int(line1[4])
                    results[dataset]["q"]=int(line1[5])
                    # calculate bnd2
                    if folname in reg_names:
                        my_data = pd.read_csv(dataset,
                                              header=None) 
                        m_normal = results[dataset]["m_normal"]
                        my_data = my_data.iloc[range(m_normal),:]
                        my_lm = LinearRegression().fit(my_data.iloc[:,1:], my_data.iloc[:,0])
                        lm_pred = my_lm.predict(my_data.iloc[:,1:])
                        results[dataset]["bnd2"] = float(m_normal)*metrics.mean_squared_error(my_data.iloc[:,0], lm_pred)
                    else:
                        i = results[dataset]["i"] 
                        results[dataset]["bnd2"] = non_outlier.loc[i,"non_outlier_sq_error"]            
            
        results_df = pd.DataFrame(results).transpose()
        #if "alg3PCA runtime" not in results_df:
        #    results_df["alg3PCA runtime"] = np.nan
        #    results_df["alg3PCA gamma"] = np.nan
        #    results_df["alg3PCA Dnon"] = np.nan
        #    results_df["alg3PCA LTS"] = np.nan
        results_df = results_df[["i","m","n","m_normal","q",
                             "alg3PCA runtime", "arob runtime","bb runtime","hbreg runtime", "lm runtime","lts runtime","lqs runtime", "cb runtime",
                             "alg3PCA tsestar", "arob tsestar","bb tsestar","hbreg tsestar", "lm tsestar","lts tsestar","lqs tsestar", "cb tsestar", "bnd2"]]
        hbreg_fail = []
        for dataset in results.keys():
            if ("lm runtime" in results[dataset]) and ("cb runtime" not in results[dataset]):
                hbreg_fail.append(dataset)
    
        fail_file = open(folname+"hbreg_fail.in", "w")
        for dataset in hbreg_fail:
            fail_file.write(dataset + "\n")
        fail_file.close()
           
    
        methods= ["alg3PCA", "arob","bb","hbreg", "lm","lts","lqs", "cb"]
        ltscols = [m + " tsestar" for m in methods]
        results_df[ltscols] = results_df[ltscols].astype(float) 
        results_df["best"] = results_df[ltscols].min(axis=1)

        ratiocols = [m + " Ratio" for m in methods]
        results_df[ratiocols] = results_df[ltscols].div(results_df["best"], axis=0)
        results_df.sort_values(by=["i"], inplace=True)
        #new_row = [np.nan for i in range(len(results_df.columns))]
        #gmeans = scipy.stats.gmean(results_df[ratiocols],axis=0)
        #print(gmeans)
        #new_row[-len(methods):] = gmeans
        #results_df.loc[len(results_df.index)] = new_row
        results_df.to_excel(writer, sheet_name=folname, float_format="%f")
        #print(results_df)
    
    writer.close()



    
