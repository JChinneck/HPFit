import pandas as pd
import numpy as np
import scipy.stats
import glob
loc="/home/jpbrooks/HPFit/experiments/cb_evaluation/results"

experiments=["evaluation"]
for experiment in experiments:
    writer = pd.ExcelWriter(loc + "/" + experiment + ".xlsx", mode="w")
    
    folnames = ["bm-like", "bm-nox", "clustered_outliers", "olive", "rvd-like", "unclustered_outliers"]
    for folname in folnames:
        results={}
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
        results_df = pd.DataFrame(results).transpose()
        #if "alg3PCA runtime" not in results_df:
        #    results_df["alg3PCA runtime"] = np.nan
        #    results_df["alg3PCA gamma"] = np.nan
        #    results_df["alg3PCA Dnon"] = np.nan
        #    results_df["alg3PCA LTS"] = np.nan
        results_df = results_df[["i","m","n","m_normal","q",
                             "alg3PCA runtime", "arob runtime","bb runtime","hbreg runtime", "lm runtime","lts runtime","lqs runtime", "cb runtime",
                             "alg3PCA tsestar", "arob tsestar","bb tsestar","hbreg tsestar", "lm tsestar","lts tsestar","lqs tsestar", "cb tsestar" ]]
        hbreg_fail = []
        for dataset in results.keys():
            if "hbreg runtime" not in results[dataset]:
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
        new_row = [np.nan for i in range(len(results_df.columns))]
        gmeans = scipy.stats.gmean(results_df[ratiocols],axis=0)
        print(gmeans)
        new_row[-len(methods):] = gmeans
        results_df.loc[len(results_df.index)] = new_row
        results_df.to_excel(writer, sheet_name=folname, float_format="%f")
        #print(results_df)
    
    writer.close()



    
