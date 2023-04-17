from A00omics2pathlist import omics2pathlist
from A01resampling import resample
from A02prediction import prediction
from A03supp_age import ageTransfer
from A04stage1 import stage1
from A05stage2 import stage2, stage2pediction
from xmlrpc.client import boolean
from itertools import repeat
import pandas as pd
import numpy as np
from typing import List, Optional
import time
from functools import reduce
 
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, r
from rpy2.robjects.conversion import localconverter

import yaml

startread = time.time()
readRDS = r['readRDS']
colnames  = r["colnames"]
intersect = r["intersect"]
Reduce = r["Reduce"]
List = r["list"]

# import data paths 
with open('./config.yml', 'r') as file:
     root = yaml.safe_load(file)

phase3 = root["pathway"]["phase3"]
pathCommon = root["pathway"]["common"]
phase3Result = root["pathway"]["phase3Result"]


autoMat = phase3.format("gwMat.rds")
covariateData =  phase3.format("covMat.rds")
CpGList = phase3.format("CpGList.xlsx")
golist = pathCommon.format("goBP202209.rds")
cpgAnno = pathCommon.format("geneIDmatV1Dataframe.rds")
chrXY = pathCommon.format("cpgAnno.rds")

# step0
# prepare data 
methylData = readRDS(autoMat)
covariateData = readRDS(covariateData)
CpGList = pd.read_excel(CpGList, index_col=0).CpG.values.tolist() 
golist = readRDS(golist)
cpgAnno = readRDS(cpgAnno)
chrXY = readRDS(chrXY)

with localconverter(ro.default_converter + pandas2ri.converter):
    methylData = ro.conversion.rpy2py(methylData)
    cpgAnno = ro.conversion.rpy2py(cpgAnno)
    covariateData = ro.conversion.rpy2py(covariateData)
    chrXY = ro.conversion.rpy2py(chrXY)

chrXY = chrXY[~chrXY['CHR'].isin(["X", "Y"])][["entrezID", "CHR"]].drop_duplicates()
cpgAnno = cpgAnno.merge(chrXY, how="left", on =["entrezID"]).dropna().drop(columns = ["CHR"])
age = covariateData[["Age", "Label"]]

age['Age'] = pd.to_numeric(age['Age'], errors='coerce')
age = age[age['Age'].notna()].astype({'Age': float})

methylData = methylData[CpGList].join(age, how= "left")
methylData = methylData[methylData["Label"].astype(str).eq("Control")]

methylData = methylData.drop(columns=["Label"])
methylData['Age'] = methylData['Age'].apply(lambda x: ageTransfer(x))
print(methylData)

golist = dict(zip(golist.names, map(list,list(golist))))
omics2pathlist = omics2pathlist(methylData, golist, cpgAnno)
print(len(omics2pathlist))
endread = time.time()
print("Data ready in {} seconds".format(endread - startread))

# call pathway Age model
sampleMode = "CrossValidation"
predictionMode = "GradientBoosting"
nfold = 5
reconData = True
outerPath = phase3Result.format("resultControl/dataForStage2Regression_fold{}.csv")
path = phase3Result.format("resultControl/data4Stage2.csv")
hyperParam = { 
   "max_features": "sqrt",
   "random_state": 20,
    }
# step1
# pathway Age model stage1
# for test only ues 5 pathways!!
stage1(
    omics2pathlist[:5],
    methylData,
    sampleMode,
    nfold,
    predictionMode,
    hyperParam,
    outerPath,
    path,
    reconData,
)

# step2.A
# pathway Age model stage2 

filePath =  phase3Result.format('resultControl/dataForStage2Regression_fold{}.csv')
savePath = phase3Result.format('resultControl/predictionAge.csv')
stage2(
    filePath,
    sampleMode,
    nfold,
    predictionMode,
    hyperParam,
    age,
    savePath,
)

# step2.B
# pathway Age model stage2 

# predictFilePath = phase3Result.format("resultCase/dataForStage2Regression_fold{}.csv")
# modelFilePath = phase3Result.format("resultControl/data4Stage2.csv")
# savePath =  phase3Result.format('resultControl/phase3ControlModelPhase3CasePrediction.csv')

# stage2pediction(
#     predictFilePath,
#     modelFilePath,
#     sampleMode,
#     nfold,
#     predictionMode,
#     hyperParam,
#     savePath,
# )

endrun = time.time()
print("Data run time {}".format(endrun - endread))