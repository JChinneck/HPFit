import pandas as pd
import numpy as np
import scipy.stats
import glob
loc="/home/jpbrooks/hyperplane_fitting/rbm_evaluation/results"



#experiments=["evaluation","alg3test"]
experiments=["evaluation"]
for experiment in experiments:
    writer = pd.ExcelWriter(loc + "/" + experiment + ".xlsx", mode="w")
    
    if experiment == "evaluation":
        folnames = ["evaluation", "bm-like", "bm-nox", "olive", "rvd-like", "evaluation_2", "one_pt_per_cluster"]
    if experiment == "alg3test":
        folnames = ["evaluation"]
    for folname in folnames:
        results={}
        for fname in glob.glob(loc+"/"+experiment+"/"+folname+"/*.csv"):
            f=open(fname, "r")
            line=f.readline()
            line1=line.split(",") 
            dataset=line1[0]
            if len(line1) > 5:
                formulation=line1[6]
                if formulation in ["arob","bb","hbreg", "lm","lts","lqs"]:
                    dataset = dataset.replace("results/evaluation", "data")
                if not dataset in results:
                    results[dataset] = {}
                results[dataset][formulation + " Dnon"]=float(line1[7])
                if formulation in ["alg3", "alg3LP", "alg3PCA", "cb"]:
                    results[dataset][formulation + " runtime"]=float(line1[9])
                    results[dataset][formulation + " gamma"]=float(line1[8])
                    results[dataset][formulation + " LTS"]=float(line1[10])
                    #print(dataset,formulation, results[dataset]["alg3PCA LTS"])
                elif formulation in ["arob", "bb", "hbreg", "lm", "lts", "lqs"]:
                    results[dataset][formulation + " runtime"]=float(line1[11])
                    results[dataset][formulation + " gamma"]=float(line1[8])
                    results[dataset][formulation + " LTS"]=float(line1[12])
                elif formulation == "rbm-mio3":
                    results[dataset][formulation + " runtime"]=float(line1[8])
                    results[dataset][formulation + " status"]=line1[9]
                    results[dataset][formulation + " gamma"]=float(line1[10])
                    results[dataset][formulation + " bestbound"]=float(line1[11])  
    
                    try:
                        results[dataset][formulation + " RBM Dnon"]=float(line1[24])
                        results[dataset][formulation + " Step 1 Dnon"]=float(line1[12])
                        results[dataset][formulation + " Step 2 Dnon"]=float(line1[14])
                        results[dataset][formulation + " Step 3 Dnon"]=float(line1[16])
    
                        results[dataset][formulation + " RBM numClose"]=float(line1[20])
                        results[dataset][formulation + " Step 1 numClose"]=float(line1[21])
                        results[dataset][formulation + " Step 2 numClose"]=float(line1[22])
                        results[dataset][formulation + " Step 3 numClose"]=float(line1[23])

                        results[dataset][formulation + " RBM LTS"]=float(line1[25])
                        results[dataset][formulation + " Step 1 LTS"]=float(line1[26])
                        results[dataset][formulation + " Step 2 LTS"]=float(line1[27])
                        results[dataset][formulation + " Step 3 LTS"]=float(line1[28])
                        results[dataset][formulation + " LTS"]=float(line1[28])
    
                        results[dataset][formulation + " num outliers in q"]=float(line1[18])
                    except:
                        x = 1
                elif formulation == "rbm":
                    results[dataset][formulation + " runtime"]=float(line1[8])
                    results[dataset][formulation + " gamma"]=float(line1[9])
                    results[dataset][formulation + " LTS"]=float(line1[10])
     
                if (formulation=="rbm") or (formulation=="alg3LP"):
                    results[dataset]["i"]=int(line1[1])
                    results[dataset]["m"]=int(line1[2])
                    results[dataset]["n"]=int(line1[3])
                    results[dataset]["m_normal"]=int(line1[4])
                    results[dataset]["q"]=int(line1[5])
        results_df = pd.DataFrame(results).transpose()
        if "alg3PCA runtime" not in results_df:
            results_df["alg3PCA runtime"] = np.nan
            results_df["alg3PCA gamma"] = np.nan
            results_df["alg3PCA Dnon"] = np.nan
            results_df["alg3PCA LTS"] = np.nan
        if experiment == "evaluation":
            results_df = results_df[["i","m","n","m_normal","q",
                                 "alg3PCA runtime", "arob runtime","bb runtime","hbreg runtime", "lm runtime","lts runtime","lqs runtime", "rbm runtime","rbm-mio3 runtime", "cb runtime",
                                 "rbm-mio3 status",
                                 "alg3PCA gamma", "arob gamma","bb gamma","hbreg gamma", "lm gamma","lts gamma","lqs gamma", "rbm gamma","rbm-mio3 gamma",
                                 "rbm-mio3 bestbound",
                                 "rbm-mio3 num outliers in q",
                                 "rbm-mio3 RBM numClose", "rbm-mio3 Step 1 numClose","rbm-mio3 Step 2 numClose","rbm-mio3 Step 3 numClose",
                                 "rbm-mio3 RBM Dnon", "rbm-mio3 Step 1 Dnon","rbm-mio3 Step 2 Dnon","rbm-mio3 Step 3 Dnon",
                                 "alg3PCA Dnon", "arob Dnon","bb Dnon","hbreg Dnon", "lm Dnon","lts Dnon","lqs Dnon", "rbm Dnon","rbm-mio3 Dnon", "cb Dnon",
                                 "rbm-mio3 RBM LTS", "rbm-mio3 Step 1 LTS","rbm-mio3 Step 2 LTS","rbm-mio3 Step 3 LTS",
                                 "alg3PCA LTS", "arob LTS","bb LTS","hbreg LTS", "lm LTS","lts LTS","lqs LTS", "rbm LTS","rbm-mio3 LTS", "cb LTS"
]]
        if experiment == "alg3test":
            print(results_df)
            results_df = results_df[["i","m","n","m_normal","q",
                                 "alg3LP runtime", "alg3PCA runtime",
                                 "alg3LP gamma", "alg3PCA gamma",
                                 "alg3LP Dnon", "alg3PCA Dnon"]]
        hbreg_fail = []
        for dataset in results.keys():
            if "rbm runtime" not in results[dataset]:
                hbreg_fail.append(dataset)
    
        fail_file = open(folname+"hbreg_fail.in", "w")
        for dataset in hbreg_fail:
            fail_file.write(dataset + "\n")
        fail_file.close()
           
    
        if experiment == "evaluation":
            methods= ["alg3PCA", "arob","bb","hbreg", "lm","lts","lqs", "rbm", "rbm-mio3", "cb"]
            ltscols = [m + " LTS" for m in methods]
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



    
