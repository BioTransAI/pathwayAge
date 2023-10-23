from typing import List, Optional
import pandas as pd
import math

class resample(object):
  def __init__(
      self,
      trainData: pd.DataFrame,
      nfold: int,
      randomState: int,
      sampleMode: Optional[str] = "CrossValidation",
  ):
    """resample
      
      for "CrossValidation": K-Folds splitters.
      Provides train/test indices to split data in train/test sets. 
      Split dataset into k consecutive folds (with fixed shuffling by default).
      Each fold is then used once as a validation while the k - 1 remaining folds form the training set.
      This Function could be extended, for example, "bootstraping" by chaning the sampleMode. 
      
      Parameters
    ----------
      trainData: methylation matrix;
      nfold: Number of folds, default=5; 
      randomState: random_state affects the ordering of the indices, which controls the randomness of each fold.
                   defult = 6677;
      sampleMode: CrossValidation
    Return
    ----------
      a list of dataframe. 
    """ 
    self.sampleMode = sampleMode
    self.trainData = trainData
    self.nfold = nfold
    self.randomState = randomState
  def resample(self):
      method=getattr(self, self.sampleMode, lambda :'Invalid sample mode')
      return method()
    
  def CrossValidation(self):
    dataCVList = []
    foldSize = math.ceil(len(self.trainData.index)/self.nfold)
    trainData = self.trainData.sample(frac=1, random_state=self.randomState)
    listFold = [list(trainData.index)[i:i + foldSize] for i in range(0, len(trainData.index), foldSize)]
    if  len(listFold) + 1 == self.nfold:
      foldSize = math.floor(len(self.trainData.index)/self.nfold)
      listFold = [list(trainData.index)[i:i + foldSize] for i in range(0, foldSize*(self.nfold -1), foldSize)]
      listFold.append(list(trainData.index)[foldSize*(self.nfold -1): len(trainData.index)])
    for fold in listFold:
      testDataCV = trainData.loc[fold]
      trainDataCV = trainData[~trainData.index.isin(fold)]
      dataCVList.append([trainDataCV, testDataCV])
    return dataCVList



