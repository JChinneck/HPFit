import pandas as pd
import numpy as np
import scipy.stats
import glob
loc="/home/jpbrooks/HPFit/experiments/mio_evaluation/results"


experiment="heuristics"
writer = pd.ExcelWriter(loc + "/" + experiment + ".xlsx", mode="w")

folnames = ["olive", "bm", "rvd", "clustered_outliers_small", "bm_small"]

for folname in folnames:
    results={}
    for fname in glob.glob(loc+"/"+experiment+"/"+folname+"/*.csv"):
        f=open(fname, "r")
        line=f.readline()
        line1=line.split(",") 
        dataset=line1[0]
        if not dataset in results:
            results[dataset] = {}
        formulation=line1[6]
        if formulation=="mio-bm-first":
            results[dataset]["i"]=int(line1[1])
            results[dataset]["m"]=int(line1[2])
            results[dataset]["n"]=int(line1[3])
            results[dataset]["m_normal"]=int(line1[4])
            results[dataset]["q"]=int(line1[5])
            #results[dataset][formulation + " Dnon"]=float(line1[7])
            results[dataset][formulation + " runtime"]=float(line1[8])
            results[dataset][formulation + " status"]=line1[9]
            results[dataset][formulation + " gamma"]=float(line1[10])
            results[dataset][formulation + " bestbound"]=float(line1[11])
            results[dataset][formulation + " tse"]=float(line1[13])
            results[dataset][formulation + " tsestar"]=float(line1[14])

        #if folname == "olive":
        if formulation=="mio1-first":
            results[dataset][formulation + " runtime"]=float(line1[8])
            results[dataset][formulation + " status"]=line1[9]
            results[dataset][formulation + " gamma"]=float(line1[10])
            results[dataset][formulation + " bestbound"]=float(line1[11])
            results[dataset][formulation + " tse"]=float(line1[13])
            results[dataset][formulation + " tsestar"]=float(line1[14])
        if formulation=="alg3PCA":
            results[dataset]["alg3 gamma"]=float(line1[8])
            results[dataset]["alg3 runtime"]=float(line1[9])
        if formulation=="cbq": 
            results[dataset]["cbq gamma"]=float(line1[7])
            results[dataset]["cbq runtime"]=float(line1[9])
        if formulation=="lqs": 
            results[dataset]["lqs gamma"]=float(line1[8])
            results[dataset]["lqs runtime"]=float(line1[11])
        elif formulation in ["mio-bm", "mio1", "alg3-mio-bm", "alg3-mio1", "lqs-mio-bm", "lqs-mio1", "cbq-mio-bm", "cbq-mio1"]:
            results[dataset][formulation + " runtime"]=float(line1[8])
            results[dataset][formulation + " gamma"]=float(line1[10])
  

    results_df = pd.DataFrame(results).transpose()
    print(results_df.columns)
    results_df = results_df[["i","m","n","m_normal","q",
                             "mio-bm-first runtime","mio1-first runtime","alg3 runtime","lqs runtime","cbq runtime",
                             "mio-bm runtime", "mio1 runtime", "alg3-mio-bm runtime", "alg3-mio1 runtime", "lqs-mio-bm runtime", "lqs-mio1 runtime", "cbq-mio-bm runtime", "cbq-mio1 runtime", 
                             "mio-bm-first gamma","mio1-first gamma","alg3 gamma","lqs gamma","cbq gamma",
                             "mio-bm gamma", "mio1 gamma", "alg3-mio-bm gamma", "alg3-mio1 gamma", "lqs-mio-bm gamma", "lqs-mio1 gamma", "cbq-mio-bm gamma", "cbq-mio1 gamma"
                             ]]
    results_df.sort_values(by=["i"], inplace=True)
    results_df.to_excel(writer, sheet_name=folname, float_format="%f")

writer.close()

experiment="comparison"
writer = pd.ExcelWriter(loc + "/" + experiment + ".xlsx", mode="w")

folnames = ["olive", "bm", "rvd", "clustered_outliers_small", "bm_small"]

for folname in folnames:
    results={}
    for fname in glob.glob(loc+"/"+experiment+"/"+folname+"/*.csv"):
        f=open(fname, "r")
        line=f.readline()
        line1=line.split(",") 
        dataset=line1[0]
        if not dataset in results:
            results[dataset] = {}
        try:
            formulation=line1[6]
            if formulation=="mio-bm": 
                results[dataset]["i"]=int(line1[1])
                results[dataset]["m"]=int(line1[2])
                results[dataset]["n"]=int(line1[3])
                results[dataset]["m_normal"]=int(line1[4])
                results[dataset]["q"]=int(line1[5])
                #results[dataset][formulation + " Dnon"]=float(line1[7])
                results[dataset][formulation + " runtime"]=float(line1[8])
                results[dataset][formulation + " status"]=line1[9]
                results[dataset][formulation + " gamma"]=float(line1[10])
                results[dataset][formulation + " bestbound"]=float(line1[11])
                results[dataset][formulation + " tse"]=float(line1[13])
                results[dataset][formulation + " tsestar"]=float(line1[14])

            if formulation=="alg3-mio-bm" or formulation=="lqs-mio-bm" or formulation=="cbq-mio-bm" or formulation=="alg3-mio1" or formulation=="lqs-mio1" or formulation=="cbq-mio1" or formulation=="mio1":  
                results[dataset][formulation + " runtime"]=float(line1[8])
                results[dataset][formulation + " status"]=line1[9]
                results[dataset][formulation + " gamma"]=float(line1[10])
                results[dataset][formulation + " bestbound"]=float(line1[11])
                results[dataset][formulation + " tse"]=float(line1[13])
                results[dataset][formulation + " tsestar"]=float(line1[14])
                results[dataset][formulation + " timelimit"]=float(line1[15])

        except:
            print(fname)


    results_df = pd.DataFrame(results).transpose()
    print(results_df.columns)
    results_df = results_df[["i","m","n","m_normal","q",
                             "mio-bm bestbound","mio1 bestbound","alg3-mio-bm bestbound","alg3-mio1 bestbound","lqs-mio-bm bestbound", "lqs-mio1 bestbound", "cbq-mio-bm bestbound", "cbq-mio1 bestbound",
                             "mio-bm status","mio1 status","alg3-mio-bm status","alg3-mio1 status","lqs-mio-bm status", "lqs-mio1 status", "cbq-mio-bm status", "cbq-mio1 status",
                             "mio-bm runtime","mio1 runtime","alg3-mio-bm runtime","alg3-mio1 runtime","lqs-mio-bm runtime", "lqs-mio1 runtime", "cbq-mio-bm runtime", "cbq-mio1 runtime",
                             "alg3-mio-bm timelimit","alg3-mio1 timelimit","lqs-mio-bm timelimit", "lqs-mio1 timelimit", "cbq-mio-bm timelimit", "cbq-mio1 timelimit",
                             "mio-bm gamma","mio1 gamma","alg3-mio-bm gamma","alg3-mio1 gamma","lqs-mio-bm gamma", "lqs-mio1 gamma", "cbq-mio-bm gamma", "cbq-mio1 gamma"
                             ]]
    results_df.sort_values(by=["i"], inplace=True)
    results_df.to_excel(writer, sheet_name=folname, float_format="%f")

writer.close()


experiment="comparison"
writer = pd.ExcelWriter(loc + "/" + experiment + "_tse.xlsx", mode="w")

folnames = ["olive", "bm", "rvd", "clustered_outliers_small", "bm_small", "bm-like", "rvd-like"]

for folname in folnames:
    results={}
    for fname in glob.glob(loc+"/"+experiment+"/"+folname+"/*.csv"):
        f=open(fname, "r")
        line=f.readline()
        line1=line.split(",") 
        dataset=line1[0]
        if not dataset in results:
            results[dataset] = {}
        try:
            formulation=line1[6]
            if formulation=="mio1": 
                results[dataset]["i"]=int(line1[1])
                results[dataset]["m"]=int(line1[2])
                results[dataset]["n"]=int(line1[3])
                results[dataset]["m_normal"]=int(line1[4])
                results[dataset]["q"]=int(line1[5])
                #results[dataset][formulation + " Dnon"]=float(line1[7])
                results[dataset][formulation + " runtime"]=float(line1[8])
                results[dataset][formulation + " status"]=line1[9]
                results[dataset][formulation + " gamma"]=float(line1[10])
                results[dataset][formulation + " bestbound"]=float(line1[11])
                results[dataset][formulation + " tse"]=float(line1[13])
                results[dataset][formulation + " tsestar"]=float(line1[14])
                results[dataset][formulation + " timelimit"]=float(line1[15])

        except:
            print(fname)


    results_df = pd.DataFrame(results).transpose()
    print(results_df.columns)
    results_df = results_df[["i","m","n","m_normal","q",
                             "mio1 bestbound",
                             "mio1 status",
                             "mio1 runtime",
                             "mio1 timelimit",
                             "mio1 tsestar",
                             "mio1 tse",
                             "mio1 gamma"
                             ]]
    results_df.sort_values(by=["i"], inplace=True)
    results_df.to_excel(writer, sheet_name=folname, float_format="%f")

writer.close()
    
#experiment="check_q"
#writer = pd.ExcelWriter(loc + "/" + experiment + ".xlsx", mode="w")
#
#folnames = ["olive"]
#
#for folname in folnames:
#    results={}
#    for fname in glob.glob(loc+"/"+experiment+"/"+folname+"/m*.csv"):
#        f=open(fname, "r")
#        line=f.readline()
#        line1=line.split(",") 
#        dataset=line1[0]
#        if not dataset in results:
#            results[dataset] = {}
#        formulation=line1[6]
#        results[dataset][formulation + " Dnon"]=float(line1[7])
#        results[dataset][formulation + " runtime"]=float(line1[8])
#        results[dataset][formulation + " status"]=line1[9]
#        results[dataset][formulation + " gamma"]=float(line1[10])
#        results[dataset][formulation + " bestbound"]=float(line1[11])
#        results[dataset][formulation + " num outliers in q"]=float(line1[12])
#        results[dataset]["i"]=int(line1[1])
#        results[dataset]["m"]=int(line1[2])
#        results[dataset]["n"]=int(line1[3])
#        results[dataset]["m_normal"]=int(line1[4])
#        results[dataset]["q"]=int(line1[5])
#    results_df = pd.DataFrame(results).transpose()
#    results_df = results_df[["i","m","n","m_normal","q",
#                             "mio1 runtime",
#                             "mio1 status",
#                             "mio1 gamma",
#                             "mio1 bestbound",
#                             "mio1 num outliers in q",
#                             "mio1 Dnon"]]
#    results_df[["mio1 Dnon"]] = results_df[["mio1 Dnon"]].astype(float) 
#    results_df.sort_values(by=["i"], inplace=True)
#    results_df.to_excel(writer, sheet_name=folname, float_format="%f")
#
#writer.close()

 



   
        

