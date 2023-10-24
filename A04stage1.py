import concurrent.futures
from xmlrpc.client import boolean
from itertools import repeat
import pandas as pd
from typing import List, Optional
import time
from functools import reduce
import os
from pathlib import Path
from A03data4stage2 import dataForStage2, reconTrainTestData


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
methylTestData: Optional[pd.DataFrame] = pd.DataFrame(), 
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
        methylTestData: testing data, default empty dataframe.

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
            print(data4Stage2)
            # data4Stage2.to_csv(outerPath.format(i))  
            cvList.append(data4Stage2)
        print("data4Stage2 datasets keeping the outer data unseen complete!")

    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
            predictionMeanList = list(executor.map(dataForStage2,
                                            omics2pathlist,
                                            repeat(predictionMode),
                                            repeat(tuneHyperParam),
                                            repeat(hyperParam),
                                            repeat(None),
                                            repeat(None),
                                            repeat(nfold),
                                            repeat(randomState)
                                            ))

        data4Stage2 = reduce(lambda df1,df2: pd.concat([df1, df2], axis=1), predictionMeanList)
        data4Stage2 = data4Stage2.join(age[["Age"]])
        print(data4Stage2)
        data4Stage2.to_csv("{}data4Stage2.csv".format(resultName))
        print("data4Stage2 datasets all data been seen complete!")
    
    return data4Stage2, cvList