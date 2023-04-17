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
import json
import yaml
import os
from pathlib import Path

import warnings
warnings.simplefilter('ignore')

def metaData():
    """
    return:
        cpgAnno: a dataframe maping CpGs to genes
        golist: a dictionary of genes on each pathway      
    """
    # import data paths 
    with open('./config.yml', 'r') as file:
        root = yaml.safe_load(file)
    metaData = root["pathway"]["metaData"]
    golist = metaData.format("golist.json")
    cpgAnno = metaData.format("cpgAnno.csv")
    # import data
    cpgAnno = pd.read_csv(cpgAnno, index_col=0)
    cpgAnno['entrezID'] = cpgAnno['entrezID'].astype(str)
    # Opening JSON file
    with open(golist) as json_file:
        golist = json.load(json_file)
    return cpgAnno, golist


def pathwayAge(
    methylData: pd.DataFrame,
    resultName: str,
    methylTestData: Optional[pd.DataFrame] = pd.DataFrame(),
    restrictUp: Optional[int] = 200,
    restrictDown: Optional[int] = 10,
    minPathSize: Optional[int] = 5,
    nfold: Optional[int] = 5,
    randomState: Optional[int] = 6677,
    predictionMode: Optional[str] = "GradientBoosting",
    reconData: Optional[boolean] = False,
    tuneHyperParam: Optional[boolean] = False,
    hyperParam: Optional[dict]= None,
    cores: Optional[int]= 5,
):
    """



    """
    cpgAnno, golist= metaData()
    methylData["Age"] = methylData["Age"].apply(lambda x: ageTransfer(x))
    age = methylData[["Age"]]
    methyl2PathList = omics2pathlist(
        methylData, 
        golist,
        cpgAnno, 
        restrictUp = restrictUp,
        restrictDown= restrictDown,
        minPathSize = minPathSize,
        )

    if methylTestData.empty:
        reconData = True
        print("no specified test data, parameter reconData FORCE to be TURE!")
        _, cvList = stage1(
            methyl2PathList[:5],
            age,
            nfold=nfold,
            randomState = randomState,
            predictionMode=predictionMode,
            reconData = reconData,
            hyperParam = hyperParam,
            tuneHyperParam = tuneHyperParam,
            cores = cores,
        )

        stage2(
            cvList = cvList,
            age = age,
            resultName = resultName,
            nfold = nfold,
            randomState = randomState,
            predictionMode = predictionMode,
            tuneHyperParam = tuneHyperParam,
            hyperParam = hyperParam,
        )
    else:
        data4Stage2, _ = stage1(
                methyl2PathList[:5],
                age,
                nfold=nfold,
                randomState = randomState,
                predictionMode = predictionMode,
                reconData = False,
                hyperParam = hyperParam,
                tuneHyperParam = tuneHyperParam,
                cores = cores,
            ) 

        methylTestData["Age"] = methylTestData["Age"].apply(lambda x: ageTransfer(x))
        testAge = methylTestData[["Age"]]
        methyl2PathTestList = omics2pathlist(
            methylTestData, 
            golist,
            cpgAnno, 
            restrictUp = restrictUp,
            restrictDown= restrictDown,
            minPathSize = minPathSize,
        )  

        _, cvList = stage1(
            methyl2PathTestList[:5],
            age = testAge,
            nfold=nfold,
            randomState = randomState,
            predictionMode = predictionMode,
            reconData = True,
            hyperParam = hyperParam,
            tuneHyperParam = tuneHyperParam,
            cores = cores,
            methylTestData = methylTestData
        )

        stage2pediction(
            predict = cvList,
            model = data4Stage2, 
            resultName = resultName,
            nfold = nfold,
            randomState = randomState,
            predictionMode = predictionMode,
            tuneHyperParam = tuneHyperParam,
            hyperParam = hyperParam,
        )


# to test pathwayAge

methylData = "./data/methlyData.csv"
covariateData =  "./data/covariateData.csv"
predictionMode = "SVR"
nfold = 3
randomState = 00
reconData = True
tuneHyperParam = True
cores = 10
### GB ###
# hyperParam = { 
#  'learning_rate': np.arange(0.05,0.21,0.05),
# }
### SVM ###
hyperParam = {
    "kernel": ["rbf", "poly"],
    "C": [20, 10, 1, 0.1],
    # "epsilon": [0.001, 0.01, 0.1, 1]
}

methylData = pd.read_csv(methylData, index_col="CpG")
methylData = methylData.T
testData = methylData.iloc[:100,:10000]
covariateData = pd.read_csv(covariateData, index_col ="Sample")
age = covariateData["Age"]
methylDataAge = testData.join(age)
methylTestData = methylData.iloc[101:201,:10000]
methylTestData = methylTestData.join(age)
resultName = "Emmaaaaa"

startread = time.time()

pathwayAge(
    methylData = methylDataAge,
    resultName = resultName,
    minPathSize = 10,
    nfold = nfold,
    tuneHyperParam  = tuneHyperParam,
    hyperParam = hyperParam,
    cores= cores,
    methylTestData = methylTestData,
    predictionMode = predictionMode,
)

endDeal = time.time()
print("Data dealed in {} seconds".format(endDeal - startread))


# hyperParam = { 
# 'learning_rate': np.arange(0.05,0.21,0.05),
# }

# cwd = os.getcwd()
# folderPath = cwd + "/" + resultFolderName
# Path(folderPath).mkdir(parents=True, exist_ok=True)
# path = folderPath + "/{}"
# outerPath = folderPath + "/dataForStage2Regression_fold{}.csv"

# hyperParam = { 
#    "max_features": "sqrt",
#    "random_state": 20,
#     }