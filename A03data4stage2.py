
from xmlrpc.client import boolean
from A00omics2pathlist import omics2pathlist
from A01resampling import resample
from A02prediction import prediction
from A03supp_age import ageTransfer 
import concurrent.futures
from itertools import repeat
import pandas as pd
from typing import List, Optional
import time
from functools import reduce



def reconTrainTestData(
    data: list,
    nfold: int,
    randomState : int
):
    """
    reconstruct the training and testing data. 
    keep the outter test data unseen. 
    use the inner test data to train the model.
    return: list of training and test datasets.
    """
    dataReconList= []
    [trainData, testData] = data
    data = resample(
        trainData = trainData, 
        nfold = nfold,
        randomState = randomState,
        ).resample()
    for fold in range(nfold):
        dataReconList.append([data[fold][0], 
                              pd.concat([data[fold][1], testData])
                              ]) 
    return dataReconList

def dataForStage2(
    data: pd.DataFrame,
    predictionMode: str,
    tuneHyperParam: boolean,
    hyperParam: dict,
    i,
    reconData,
    nfold,
    randomState,
):
  """
  generate the data for stage2. 
  data: a dataframe contains CpG in one pathway; d
  predictionMode: can be swithed for different tasks, here we ues "GradientBoosting";
  i: ith fold of cross validation;
  return: a dataframe -> prediction mean based on nfold inner training data models.
  """
  data = resample(
            trainData = data, 
            nfold = nfold,
            randomState = randomState).resample()
  if reconData:
    data = reconTrainTestData(data[i], nfold, randomState)
  predictionList = []
  for elem in data:
    [trainData, testData] = elem
    result = prediction(predictionMode, trainData, testData, tuneHyperParam, hyperParam).predictionMode()   
    predictionList.append(result)
  predictionAll = reduce(lambda df1,df2: pd.concat([df1, df2], axis=1), predictionList)
  predictionAll["mean"] = predictionAll.mean(axis = 1, numeric_only=True, skipna=True)
  predictionMean = predictionAll[["mean"]].rename(columns = {"mean": predictionAll.columns[0]})
  return predictionMean

