
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
import numpy as np



def reconTrainTestData(
    data: list,
    nfold: int,
    randomState : int
):


    """reconTrainTestData
    
      Data reassembly processing. 
      Reconstructing training/testing data into sepcfic form.
      Details please see the tutorial documents.

    Parameter
    ----------
      data: a list of training and testing datasets;
      nfold: Number of folds, default=5;
      randomState: random_state affects the ordering of the indices, which controls the randomness of each fold.
                   defult = 6677;

    Return
    ----------
      a list of dataframe. 
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


  """dataForStage2
  
    Generate the intermediate data for the next step called "stage2". 

  Parameter
  ----------
    data: a dataframe contains CpG in one pathways;
    predictionMode: methods to generate the model, 'Ridge', 'SVR' or 'GradientBoosting';
    tuneHyperParam: utilizing Customized hyperparameters if TURE;
    hyperParam: Customized hyperparameters;
    i: The ith fold in K-Folds;
    reconData: Reconstructing data if TRUE;
    nfold: Number of folds, default=5;
    randomState: random_state affects the ordering of the indices, which controls the randomness of each fold,
                 defult = 6677;

  Returnƒ
  ----------
    a dataframe. 
  """ 

  indexOrder = pd.DataFrame(
    index=data.index, 
    columns=["rank"], 
    data = np.arange(data.shape[0]))
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
  predictionMean = indexOrder.join(predictionMean).sort_values(by=["rank"]).drop(columns=["rank"])

  return predictionMean

