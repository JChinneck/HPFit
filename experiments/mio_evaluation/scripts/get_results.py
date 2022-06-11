import pandas as pd
import numpy as np
import scipy.stats
import glob
loc="/home/jpbrooks/HPFit/experiments/mio_evaluation/results"



experiment="comparison"
writer = pd.ExcelWriter(loc + "/" + experiment + ".xlsx", mode="w")

folnames = ["olive", "bm", "rvd"]

for folname in folnames:
    results={}
    for fname in glob.glob(loc+"/"+experiment+"/"+folname+"/*.csv"):
        f=open(fname, "r")
        line=f.readline()
        line1=line.split(",") 
        dataset=line1[0]
        if not dataset in results:
            results[dataset] = {}
        print(line1, fname)
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

        #if folname == "olive":
        if formulation=="mio3":
            results[dataset]["mio1 runtime"]=float(line1[8])
            results[dataset]["mio3 runtime"]=float(line1[28])
            results[dataset]["mio1 status"]=line1[9]
            results[dataset]["mio1 bestbound"]=float(line1[11])
            results[dataset]["mio1 gamma"]=float(line1[10])
            results[dataset]["mio2 gamma"]=float(line1[26])
            results[dataset]["mio3 gamma"]=float(line1[27])
            results[dataset]["mio1 tse"]=float(line1[20])
            results[dataset]["mio2 tse"]=float(line1[21])
            results[dataset]["mio3 tse"]=float(line1[22])
            results[dataset]["mio1 tsestar"]=float(line1[23])
            results[dataset]["mio2 tsestar"]=float(line1[24])
            results[dataset]["mio3 tsestar"]=float(line1[25])
        if formulation=="cbmio3": 
            results[dataset]["cbmio2 runtime"]=float(line1[7])
            results[dataset]["cbmio2 status"]=line1[8]
            results[dataset]["cbmio2 bestbound"]=float(line1[10])

            results[dataset]["cbmio1 tsestar"]=float(line1[11])
            results[dataset]["cbmio1 tse"]=float(line1[12])
            results[dataset]["cbmio1 gamma"]=float(line1[13])
            results[dataset]["cbmio2 tsestar"]=float(line1[14])
            results[dataset]["cbmio2 tse"]=float(line1[15])
            results[dataset]["cbmio2 gamma"]=float(line1[16])
            results[dataset]["cbmio3 tsestar"]=float(line1[17])
            results[dataset]["cbmio3 tse"]=float(line1[18])
            results[dataset]["cbmio3 gamma"]=float(line1[19])
            results[dataset]["cbmio4 tsestar"]=float(line1[20])
            results[dataset]["cbmio4 tse"]=float(line1[21])
            results[dataset]["cbmio4 gamma"]=float(line1[22])

            results[dataset]["cbmio1 runtime"]=float(line1[28])
            results[dataset]["cbmio3 runtime"]=float(line1[29])


    results_df = pd.DataFrame(results).transpose()
    #print(results_df.columns)
    results_df = results_df[["i","m","n","m_normal","q",
                             "mio-bm runtime","mio1 runtime","mio3 runtime", "cbmio2 runtime", "cbmio1 runtime",
                             "mio-bm status","mio1 status", "cbmio2 status",
                             "mio-bm bestbound","mio1 bestbound", "cbmio2 bestbound",
                             "mio-bm gamma","mio1 gamma","mio2 gamma","mio3 gamma","cbmio1 gamma","cbmio2 gamma","cbmio3 gamma","cbmio4 gamma",
                             "mio-bm tse","mio1 tse","mio2 tse","mio3 tse","cbmio1 tse","cbmio2 tse","cbmio3 tse","cbmio4 tse",
                             "mio-bm tsestar","mio1 tsestar","mio2 tsestar","mio3 tsestar","cbmio1 tsestar","cbmio2 tsestar","cbmio3 tsestar","cbmio4 tsestar"
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

 



   
        

