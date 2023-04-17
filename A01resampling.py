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
    listFold = [list(self.trainData.index)[i:i + foldSize] for i in range(0, len(self.trainData.index), foldSize)]
    for fold in listFold:
      trainData = self.trainData.sample(frac=1, random_state=self.randomState)
      testDataCV = trainData.loc[fold]
      trainDataCV = trainData[~trainData.index.isin(fold)]
      dataCVList.append([trainDataCV, testDataCV])
    return dataCVList



