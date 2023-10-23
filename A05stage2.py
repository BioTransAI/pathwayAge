
from A01resampling import resample
from A02prediction import prediction4stage2
from A03supp_age import antiAgeTransfer
import pandas as pd
from functools import reduce
from xmlrpc.client import boolean


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

    """stage2
    
        Finial step for the biological age prediction.
        utilizing the intermediate data from stage1 building predictor.

    Parameter
    ----------
        cvList: a list of dataframe, each dataframe contains CpGs of one pathway;
        age: a dataframe contains age of training datset;
        resultName: the file name of biological age prediction;
        nfold: Number of folds, default=5;
        randomState: random_state affects the ordering of the indices, which controls the randomness of each fold.
                     defult = 6677;
        predictionMode: methods to generate the model, 'Ridge', 'SVR' or 'GradientBoosting';
        tuneHyperParam: utilizing Customized hyperparameters if TURE;
        hyperParam: Customized hyperparameters;

    Return
    ----------
        a dataframe. 
    """ 
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
    prediction = prediction[["mean", "Age"]].applymap(lambda x: antiAgeTransfer(x))
    prediction = prediction.rename(columns={"mean": "prediction"})
    print("prediction", prediction)
    prediction.to_csv("{}.csv".format(resultName))

