import pandas as pd
import numpy as np
import scipy.stats
import glob
loc="/home/jpbrooks/hyperplane_fitting/mio_evaluation/results"



experiment="comparison"
writer = pd.ExcelWriter(loc + "/" + experiment + ".xlsx", mode="w")

folnames = ["no_outliers", "one_clust_vary_dist", "olive"]

for folname in folnames:
    results={}
    for fname in glob.glob(loc+"/"+experiment+"/"+folname+"/m*.csv"):
        f=open(fname, "r")
        line=f.readline()
        line1=line.split(",") 
        dataset=line1[0]
        if not dataset in results:
            results[dataset] = {}
        formulation=line1[6]
        results[dataset][formulation + " Dnon"]=float(line1[7])
        results[dataset][formulation + " runtime"]=float(line1[8])
        #if folname == "olive":
        results[dataset][formulation + " status"]=line1[9]
        results[dataset][formulation + " gamma"]=float(line1[10])
        results[dataset][formulation + " bestbound"]=float(line1[11])
        if formulation=="mio3":
            results[dataset][formulation + " Step 1 Dnon"]=float(line1[12])
            results[dataset][formulation + " Step 2 Dnon"]=float(line1[14])
            results[dataset][formulation + " Step 3 Dnon"]=float(line1[16])
            results[dataset][formulation + " num outliers in q"]=float(line1[18])
        else:
            results[dataset][formulation + " num outliers in q"]=float(line1[12])
        #else:
        #    line2=f.readline()
        #    line2=line2.split(",")
        #    results[dataset][formulation + " gamma"]=float(line2[1])
        #    results[dataset][formulation + " bestbound"]=float(line2[2])
        #    if formulation=="mio3":
        #        results[dataset][formulation + " Step 1 Dnon"]=float(line2[3])
        #        results[dataset][formulation + " Step 2 Dnon"]=float(line2[5])
        #        results[dataset][formulation + " Step 3 Dnon"]=float(line2[7])
        if formulation=="mio-bm":
            results[dataset]["i"]=int(line1[1])
            results[dataset]["m"]=int(line1[2])
            results[dataset]["n"]=int(line1[3])
            results[dataset]["m_normal"]=int(line1[4])
            results[dataset]["q"]=int(line1[5])
    results_df = pd.DataFrame(results).transpose()
    #print(results_df)
    results_df = results_df[["i","m","n","m_normal","q",
                             "mio-bm runtime","mio1 runtime","mio3 runtime",
                             "mio-bm status","mio1 status","mio3 status",
                             "mio-bm gamma","mio1 gamma","mio3 gamma",
                             "mio-bm bestbound","mio1 bestbound","mio3 bestbound",
                             "mio-bm num outliers in q","mio1 num outliers in q","mio3 num outliers in q",
                             "mio3 Step 1 Dnon","mio3 Step 2 Dnon","mio3 Step 3 Dnon",
                             "mio-bm Dnon","mio1 Dnon","mio3 Dnon"]]
    
    results_df[["mio-bm Dnon", "mio1 Dnon", "mio3 Dnon"]] = results_df[["mio-bm Dnon", "mio1 Dnon", "mio3 Dnon"]].astype(float) 
    results_df["best"] = results_df[["mio-bm Dnon", "mio1 Dnon", "mio3 Dnon"]].min(axis=1)
    results_df["mio-bm Ratio"] = results_df["mio-bm Dnon"]/results_df["best"]
    results_df["mio1 Ratio"]   = results_df["mio1 Dnon"]/results_df["best"]
    results_df["mio3 Ratio"]   = results_df["mio3 Dnon"]/results_df["best"]
    results_df.sort_values(by=["i"], inplace=True)
    #print(results_df["mio-bm Dnon"])
    new_row = [np.nan for i in range(len(results_df.columns))]
    gmeans = scipy.stats.gmean(results_df[["mio-bm Ratio", "mio1 Ratio", "mio3 Ratio"]],axis=0)
    new_row[-3:] = gmeans
    #print(results_df.shape)
    #print(new_row)
    results_df.loc[len(results_df.index)] = new_row
    results_df.to_excel(writer, sheet_name=folname, float_format="%f")
    #print(results_df)

writer.close()

    
experiment="check_q"
writer = pd.ExcelWriter(loc + "/" + experiment + ".xlsx", mode="w")

folnames = ["no_outliers", "one_clust_vary_dist", "olive"]

for folname in folnames:
    results={}
    for fname in glob.glob(loc+"/"+experiment+"/"+folname+"/m*.csv"):
        f=open(fname, "r")
        line=f.readline()
        line1=line.split(",") 
        dataset=line1[0]
        if not dataset in results:
            results[dataset] = {}
        formulation=line1[6]
        results[dataset][formulation + " Dnon"]=float(line1[7])
        results[dataset][formulation + " runtime"]=float(line1[8])
        results[dataset][formulation + " status"]=line1[9]
        results[dataset][formulation + " gamma"]=float(line1[10])
        results[dataset][formulation + " bestbound"]=float(line1[11])
        results[dataset][formulation + " num outliers in q"]=float(line1[12])
        results[dataset]["i"]=int(line1[1])
        results[dataset]["m"]=int(line1[2])
        results[dataset]["n"]=int(line1[3])
        results[dataset]["m_normal"]=int(line1[4])
        results[dataset]["q"]=int(line1[5])
    results_df = pd.DataFrame(results).transpose()
    results_df = results_df[["i","m","n","m_normal","q",
                             "mio1 runtime",
                             "mio1 status",
                             "mio1 gamma",
                             "mio1 bestbound",
                             "mio1 num outliers in q",
                             "mio1 Dnon"]]
    results_df[["mio1 Dnon"]] = results_df[["mio1 Dnon"]].astype(float) 
    results_df.sort_values(by=["i"], inplace=True)
    results_df.to_excel(writer, sheet_name=folname, float_format="%f")

writer.close()

 



   
        

