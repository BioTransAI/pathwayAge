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
    """ pathwayAge for Age prediction
   this estimator builds a two-stage model for biological age prediction.
   It allows for the optimization of arbitrary differentiable regression model.
   In each stage a regression optimization is allowed.
   If methylTestData not provided, pathwayAge will predict suject age in methylTestData 
   by Cross Validation. 
   Parameters
   ----------
   methylData: methylation matrix with Age info;
   resultName: output file name;
   methylTestData: test data, default empty dataframe;
   restrictUp: the upper limit of number of genes in one pathway, default 200;
   restrictDown: the lower limit of number of genes in one pathway, default 10;
   minPathSize: CpG sites minimal number in each pathway, default 5;
   nfold: N-fold of cross validation, default 5;
   randomState: seed for random number generator, default 6677;
   predictionMode: regression model, default GradientBoosting;
   reconData: the way of generate the stage1 data(the age score of each BP pathway).
        #TODO
        default False
   tuneHyperParam: whether tuning the hyperparamter of each regression stage or not.
        default False;
    
   hyperParam: hyperparameter of regression model, default None;
   cores: number of processors on the machine, default 5.
   return
   ----------

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
            methyl2PathList,
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
                methyl2PathList,
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
            methyl2PathTestList,
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
