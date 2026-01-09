import concurrent.futures
from xmlrpc.client import boolean
from itertools import repeat
import pandas as pd
from typing import List, Optional
import time
from functools import reduce
import os
from pathlib import Path
from A03data4stage2 import dataForStage2, reconTrainTestData, dataForStage2Prediction


def stage1(
omics2pathlist,
age, 
resultName, 
nfold: int,
randomState: int,
predictionMode: str,
reconData: boolean,
tuneHyperParam: boolean,
hyperParam: dict,
cores: int,
omics2pathTestList: Optional[List] = [], 
):

    """stage1
    
        Generate the intermediate data for the next step called "stage2" parallelly. 

    Parameter
    ----------
        omics2pathlist: a list of dataframe, each dataframe contains CpGs of one pathway;
        age: a dataframe contains age of training datset;
        nfold: Number of folds, default=5;
        randomState: utilizing Customized hyperparameters if TURE;
        predictionMode: methods to generate the model, 'Ridge', 'SVR' or 'GradientBoosting';
        reconData: Reconstructing data if TRUE;
        tuneHyperParam: utilizing Customized hyperparameters if TURE;
        hyperParam: Customized hyperparameters;
        i: The ith fold in K-Folds;
        cores: set multiple CPU cores, default 5.
        omics2pathTestList: testing data,  a list of dataframe, each dataframe contains CpGs of one pathway;

    Return
    ----------
        a dataframe. 
    """ 
    cvList = []
    if reconData: 
        for i in range(nfold):
            with concurrent.futures.ProcessPoolExecutor(max_workers = cores) as executor:
                predictionMeanList = list(executor.map(dataForStage2,
                                                omics2pathlist,
                                                repeat(predictionMode),
                                                repeat(tuneHyperParam),
                                                repeat(hyperParam),
                                                repeat(i),
                                                repeat(reconData),
                                                repeat(nfold),
                                                repeat(randomState)
                                                ))

            data4Stage2 = reduce(lambda df1,df2: pd.concat([df1, df2], axis=1), predictionMeanList)
            data4Stage2 = data4Stage2.join(age[["Age"]])
            # print(data4Stage2)
            # data4Stage2.to_csv(outerPath.format(i))  
            cvList.append(data4Stage2)
        print("Outer fold labels remain unseen in training!")
        return data4Stage2, cvList

    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
            predictionListTrain = list(executor.map(dataForStage2,
                                            omics2pathlist,
                                            repeat(predictionMode),
                                            repeat(tuneHyperParam),
                                            repeat(hyperParam),
                                            repeat(None),
                                            repeat(None),
                                            repeat(nfold),
                                            repeat(randomState)
                                            ))

        data4Stage2Train = reduce(lambda df1,df2: pd.concat([df1, df2], axis=1), predictionListTrain)
        data4Stage2Train = data4Stage2Train.join(omics2pathlist[0][["Age"]])
        print(data4Stage2Train.shape)
        print("all labels are used for training without cross-validation!")
        
        dataList = list(zip(omics2pathlist, omics2pathTestList))
        with concurrent.futures.ProcessPoolExecutor(max_workers=50) as executor:
            predictionListTest = list(executor.map(dataForStage2Prediction,
                                            dataList,
                                            repeat(predictionMode),
                                            repeat(tuneHyperParam),
                                            repeat(hyperParam),
                                            ))
        data4Stage2Test = reduce(lambda df1,df2: pd.concat([df1, df2], axis=1), predictionListTest)
        data4Stage2Test = data4Stage2Test.join(age[["Age"]])

        print(data4Stage2Test.shape)
        print("data4Stage2 testing dataset!")
        
        return data4Stage2Train, data4Stage2Test
    