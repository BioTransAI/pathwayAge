
from A01resampling import resample
from A02prediction import prediction4stage2
from A03supp_age import antiAgeTransfer, ageTransfer
import pandas as pd
from functools import reduce
from xmlrpc.client import boolean
from typing import List, Optional

def stage2(
    cvList: list,
    age: pd.DataFrame,
    resultName: str,
    nfold: int,
    randomState: int,
    predictionMode: str,
    tuneHyperParam: boolean,
    hyperParam: dict,
): 
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
    resultList = []
    for i in range(nfold):
        dataForStage2 = cvList[i]
        data = resample(
            trainData = dataForStage2,
            nfold = nfold,
            randomState =randomState).resample()
        result = prediction4stage2(
            data[i], 
            predictionMode, 
            tuneHyperParam,
            hyperParam).predictionMode()
        resultList.append(result)
    prediction = reduce(lambda df1,df2: pd.concat([df1, df2], axis=0), resultList)
    prediction = prediction.astype(float)
    prediction = prediction.join(age[["Age"]])
    prediction = prediction.applymap(lambda x: antiAgeTransfer(x))

    print(prediction)
  
    mean_squared_error = mean_squared_error(prediction["Age"].values, prediction["prediction"].values)
    mean_absolute_error = mean_absolute_error(prediction["Age"].values, prediction["prediction"].values)
    r2_score = r2_score(prediction["Age"].values, prediction["prediction"].values)
    dataCorr = prediction.drop(columns=["Age"]).corrwith(prediction['Age'], method='pearson')
    print("mean_squared_error: {}".format(mean_squared_error))
    print("mean_absolute_error: {}".format(mean_absolute_error))
    print("r2_score: {}".format(r2_score))
    print("correlation: {}".format(dataCorr))
    prediction.to_csv("{}.csv".format(resultName))


def stage2pediction(
    predict: list,
    model: pd.DataFrame, 
    resultName: str,
    nfold: int,
    randomState: int,
    predictionMode: str,
    tuneHyperParam: boolean,
    hyperParam: dict,
):
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
    resultList = []
    testDataList = []
    for i in range(nfold):
        predictionDataForStage2 = predict[i]
        predictionData = resample(
            trainData = predictionDataForStage2,
            nfold = nfold,
            randomState= randomState).resample()
        [_, testData] = predictionData[i]
        testDataList.append(testData)
    predictTestData = reduce(lambda df1,df2: pd.concat([df1, df2], axis=0), testDataList)
    for i in range(nfold):
        dataForStage2 = model
        data = resample(
            trainData = dataForStage2,
            nfold = nfold,
            randomState = randomState).resample()
        [innertrainData, _] = data[i]
        data = [innertrainData, predictTestData]
        result = prediction4stage2(
            data, 
            predictionMode, 
            tuneHyperParam, 
            hyperParam).predictionMode()
        resultList.append(result)

    prediction = reduce(lambda df1,df2: pd.concat([df1, df2], axis=1), resultList)
    prediction["mean"] = prediction.mean(axis = 1, numeric_only=True, skipna=True)
    prediction = prediction.astype(float)
    prediction = prediction.join(predict[0][["Age"]])
    print("prediction", prediction)
    prediction = prediction[["mean", "Age"]].applymap(lambda x: antiAgeTransfer(x))
    prediction = prediction.rename(columns={"mean": "prediction"})
    print("prediction", prediction)

    mean_squared_error = mean_squared_error(prediction["Age"].values, prediction["prediction"].values)
    mean_absolute_error = mean_absolute_error(prediction["Age"].values, prediction["prediction"].values)
    r2_score = r2_score(prediction["Age"].values, prediction["prediction"].values)
    dataCorr = prediction[["prediction"]].corrwith(prediction['Age'], method='pearson')
    print("mean_squared_error: {}".format(mean_squared_error))
    print("mean_absolute_error: {}".format(mean_absolute_error))
    print("r2_score: {}".format(r2_score))
    print("correlation: {}".format(dataCorr))
    prediction.to_csv("{}.csv".format(resultName))

