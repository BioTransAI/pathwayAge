from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import RandomizedSearchCV
from typing import List, Optional
from xmlrpc.client import boolean
import pandas as pd
import numpy as np
from sklearn.svm import SVR


class prediction(object):
    def __init__(
        self,
        predictionModule: str,
        trainData: pd.DataFrame,
        testData: pd.DataFrame,
        tuneHyperParam: boolean,
        hyperParam: dict,
        showTuningDetails: Optional[boolean] = True
    ):
      self.predictionModule = predictionModule
      self.trainData = trainData
      self.testData = testData
      self.hyperParam = hyperParam
      self.tuneHyperParam = tuneHyperParam
      self.showTuningDetails = showTuningDetails
    def predictionMode(self):
      method=getattr(self, self.predictionModule, lambda :'Invalid prediction mode')
      return method()

    def GradientBoosting(self):
        trainData = self.trainData.drop(columns=["GOName","Age"])
        if self.tuneHyperParam: 
          if self.hyperParam:
            model = GradientBoostingRegressor(**self.hyperParam, max_features= "sqrt", random_state = 20)
            Rsearch = RandomizedSearchCV(estimator = model, 
                                  param_distributions = self.hyperParam, 
                                  scoring='neg_mean_absolute_error', 
                                  n_iter=100, 
                                  random_state=0,
                                  n_jobs=30, 
                                  cv=5)
            Rsearch.fit(trainData, self.trainData["Age"])
            if self.showTuningDetails:
              print(" Results from random Search " )
              print(Rsearch.best_params_, Rsearch.best_score_)
            # use the tuning hyperParameter to train the model
            bestParams = Rsearch.best_params_ 
            GB = GradientBoostingRegressor(**bestParams)
          else:
            print("hyperParam type not valid, please check again!")
        else:
          if self.hyperParam:
            hyperParam = self.hyperParam
          else:
            hyperParam = {"random_state": 20}
          GB = GradientBoostingRegressor(**hyperParam)
        GBMod = GB.fit(trainData, self.trainData["Age"])
        result = GBMod.predict(self.testData.drop(columns=["GOName","Age"]))
        return pd.DataFrame(index = self.testData.index, 
          columns = [self.testData["GOName"][0]], 
          data = result)
  
    def SVR(self):
      trainData = self.trainData.drop(columns=["GOName","Age"])
      if self.tuneHyperParam: 
        if self.hyperParam:
          print(**self.hyperParam)
          model = SVR(**self.hyperParam)
          Rsearch = RandomizedSearchCV(estimator = model, 
                              param_distributions = self.hyperParam, 
                              scoring='neg_mean_absolute_error', 
                              n_iter=18, 
                              random_state=0,
                              n_jobs=30, 
                              cv=5)
          Rsearch.fit(trainData, self.trainData["Age"])
          if self.showTuningDetails:
              print(" Results from Grid Search " )
              print(Rsearch.best_params_, Rsearch.best_score_)
          bestParams = Rsearch.best_params_ 
          regr = SVR(**bestParams)
        else:
          print("hyperParam type not valid, please check again!")
      else:
        if self.hyperParam:
          hyperParam = self.hyperParam
        else:
          hyperParam = {"kernel":"poly"}
        regr = SVR(**hyperParam)
      SVRmodel = regr.fit(trainData, self.trainData["Age"])
      result = SVRmodel.predict(self.testData.drop(columns=["GOName","Age"]))
      return pd.DataFrame(index = self.testData.index, 
                    columns = [self.testData["GOName"][0]], 
                    data = result)

class prediction4stage2(object):
    def __init__(
        self,
        data: pd.DataFrame,
        predictionModule: str,
        tuneHyperParam: boolean,
        hyperParam: dict,
        showTuningDetails: Optional[boolean] = True
    ):
      self.predictionModule = predictionModule
      self.data = data
      self.tuneHyperParam = tuneHyperParam
      self.hyperParam = hyperParam
      self.showTuningDetails = showTuningDetails
    
    def predictionMode(self):
      method=getattr(self, self.predictionModule, lambda :'Invalid prediction mode')
      return method()
    
    def GradientBoosting(self):
        [trainData, testData] = self.data
        if self.tuneHyperParam: 
            if self.hyperParam:
              model = GradientBoostingRegressor(**self.hyperParam, max_features= "sqrt", random_state = 20)
              Rsearch = RandomizedSearchCV(estimator = model, 
                                    param_distributions = self.hyperParam, 
                                    scoring='neg_mean_absolute_error', 
                                    n_iter=100, 
                                    random_state=0,
                                    n_jobs=30, 
                                    cv=5)
              Rsearch.fit(trainData, trainData["Age"])
              if self.showTuningDetails:
                print(" Results from random Search " )
                print(Rsearch.best_params_, Rsearch.best_score_)
              # use the tuning hyperParameter to train the model
              bestParams = Rsearch.best_params_ 
              GB = GradientBoostingRegressor(**bestParams)
            else:
              print("hyperParam type not valid, please check again!")
        else:
          if self.hyperParam:
            hyperParam = self.hyperParam
          else:
            hyperParam = {"random_state": 20}
          GB = GradientBoostingRegressor(**hyperParam)
        GBMod = GB.fit(trainData.drop(columns=["Age"]), trainData["Age"])
        result = GBMod.predict(testData.drop(columns=["Age"]))
        prediction = pd.DataFrame(index = testData.index)
        prediction["prediction"] = result
        return prediction

    def SVR(self):
        [trainData, testData] = self.data
        if self.tuneHyperParam: 
            if self.hyperParam:
              model = SVR(**self.hyperParam)
              Rsearch = RandomizedSearchCV(estimator = model, 
                                    param_distributions = self.hyperParam, 
                                    scoring='neg_mean_absolute_error', 
                                    n_iter=100, 
                                    random_state=0,
                                    n_jobs=30, 
                                    cv=5)
              Rsearch.fit(trainData, trainData["Age"])
              if self.showTuningDetails:
                print(" Results from random Search " )
                print(Rsearch.best_params_, Rsearch.best_score_)
              # use the tuning hyperParameter to train the model
              bestParams = Rsearch.best_params_ 
              regr = SVR(**bestParams)
            else:
              print("hyperParam type not valid, please check again!")
        else:
          if self.hyperParam:
            hyperParam = self.hyperParam
          else:
            hyperParam = {"kernel":"poly"}
          regr = SVR(**hyperParam)
        SVRmodel = regr.fit(trainData.drop(columns=["Age"]), trainData["Age"])
        result = SVRmodel.predict(testData.drop(columns=["Age"]))
        prediction = pd.DataFrame(index = testData.index)
        prediction["prediction"] = result
        return prediction


