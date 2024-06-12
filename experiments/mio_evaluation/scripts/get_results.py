import pandas as pd
import numpy as np
import scipy.stats
import glob
loc="/lustre/home/jpbrooks/HPFit/experiments/mio_evaluation/results"


experiment="miostarts"
writer = pd.ExcelWriter(loc + "/" + experiment + ".xlsx", mode="w")

folnames = ["olive", "bm", "rvd", "clustered_outliers_small", "bm_small", "bm-like", "rvd-like"]
#folnames = ["bm"]

for folname in folnames:
    results={}
    for fname in glob.glob(loc+"/"+experiment+"/"+folname+"/*.csv"):
        f=open(fname, "r")
        line=f.readline()
        line1=line.split(",") 
        dataset=line1[0]
        #print(fname)
        if not dataset in results:
            results[dataset] = {}
        print(dataset) 
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
            results[dataset][formulation + " gammaHeur"]=float(line1[16])
            results[dataset][formulation + " tsestarHeur"]=float(line1[17])
            results[dataset][formulation + " gamma60"]=float(line1[18])
            results[dataset][formulation + " tsestar60"]=float(line1[19])
            results[dataset][formulation + " gamma3600"]=float(line1[20])
            results[dataset][formulation + " tsestar3600"]=float(line1[21])
            results[dataset][formulation + " bestbound"]=float(line1[11])
            results[dataset][formulation + " tse"]=float(line1[13])
            results[dataset][formulation + " tsestar"]=float(line1[14])
            results[dataset][formulation + " timelimit"]=float(line1[15])

        #if folname == "olive":
        else: 
            results[dataset][formulation + " runtime"]=float(line1[8])
            results[dataset][formulation + " status"]=line1[9]
            results[dataset][formulation + " gamma"]=float(line1[10])
            results[dataset][formulation + " bestbound"]=float(line1[11])
            results[dataset][formulation + " tse"]=float(line1[13])
            results[dataset][formulation + " tsestar"]=float(line1[14])
            results[dataset][formulation + " gammaHeur"]=float(line1[16])
            results[dataset][formulation + " tsestarHeur"]=float(line1[17])
            results[dataset][formulation + " gamma60"]=float(line1[18])
            results[dataset][formulation + " tsestar60"]=float(line1[19])
            results[dataset][formulation + " gamma3600"]=float(line1[20])
            results[dataset][formulation + " tsestar3600"]=float(line1[21])
            results[dataset][formulation + " bestbound"]=float(line1[11])
            results[dataset][formulation + " timelimit"]=float(line1[15])
  

    results_df = pd.DataFrame(results).transpose()
    print(results_df.columns)
    # need to add mio-bm-first and mio1-first
    results_df = results_df[["i","m","n","m_normal","q",
                             "mio-bm timelimit", "mio-bm runtime", "mio-bm bestbound", "mio-bm status", "mio-bm gammaHeur", "mio-bm tsestarHeur", "mio-bm gamma60", "mio-bm tsestar60", "mio-bm gamma3600", "mio-bm tsestar3600", 
                             "lqs-mio-bm timelimit", "lqs-mio-bm runtime", "lqs-mio-bm bestbound", "lqs-mio-bm status", "lqs-mio-bm gammaHeur", "lqs-mio-bm tsestarHeur", "lqs-mio-bm gamma60", "lqs-mio-bm tsestar60", "lqs-mio-bm gamma3600", "lqs-mio-bm tsestar3600", 
                             "alg3-mio-bm timelimit", "alg3-mio-bm runtime", "alg3-mio-bm bestbound", "alg3-mio-bm status", "alg3-mio-bm gammaHeur", "alg3-mio-bm tsestarHeur", "alg3-mio-bm gamma60", "alg3-mio-bm tsestar60", "alg3-mio-bm gamma3600", "alg3-mio-bm tsestar3600", 
                             "cbq-mio-bm timelimit", "cbq-mio-bm runtime", "cbq-mio-bm bestbound", "cbq-mio-bm status", "cbq-mio-bm gammaHeur", "cbq-mio-bm tsestarHeur", "cbq-mio-bm gamma60", "cbq-mio-bm tsestar60", "cbq-mio-bm gamma3600", "cbq-mio-bm tsestar3600"
                             , 
                             "mio1 timelimit", "mio1 runtime", "mio1 bestbound", "mio1 status", "mio1 gammaHeur", "mio1 tsestarHeur", "mio1 gamma60", "mio1 tsestar60", "mio1 gamma3600", "mio1 tsestar3600", 
                             "lqs-mio1 timelimit", "lqs-mio1 runtime", "lqs-mio1 bestbound", "lqs-mio1 status", "lqs-mio1 gammaHeur", "lqs-mio1 tsestarHeur", "lqs-mio1 gamma60", "lqs-mio1 tsestar60", "lqs-mio1 gamma3600", "lqs-mio1 tsestar3600", 
                             "alg3-mio1 timelimit", "alg3-mio1 runtime", "alg3-mio1 bestbound", "alg3-mio1 status", "alg3-mio1 gammaHeur", "alg3-mio1 tsestarHeur", "alg3-mio1 gamma60", "alg3-mio1 tsestar60", "alg3-mio1 gamma3600", "alg3-mio1 tsestar3600", 
                             "cbq-mio1 timelimit", "cbq-mio1 runtime", "cbq-mio1 bestbound", "cbq-mio1 status", "cbq-mio1 gammaHeur", "cbq-mio1 tsestarHeur", "cbq-mio1 gamma60", "cbq-mio1 tsestar60", "cbq-mio1 gamma3600", "cbq-mio1 tsestar3600",
                             "mio-bm-first runtime", "mio-bm-first gamma60", "mio-bm-first tsestar60", 
                             "mio1-first runtime", "mio1-first gamma60", "mio1-first tsestar60"
                             ]]
    results_df.sort_values(by=["i"], inplace=True)
    results_df.to_excel(writer, sheet_name=folname, float_format="%f")

writer.close()

