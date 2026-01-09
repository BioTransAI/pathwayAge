
from A01resampling import resample
from A02prediction import prediction4stage2
from A03supp_age import antiAgeTransfer
import pandas as pd
from functools import reduce
from xmlrpc.client import boolean
import yaml
from typing import List, Optional


with open('./config.yml', 'r') as file:
        root = yaml.safe_load(file)

resultPath = root["pathway"]["result"]


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
    TrainDataList = []
    for i in range(nfold):
        dataForStage2 = cvList[i]
        data = resample(
            trainData = dataForStage2,
            nfold = nfold,
            randomState =randomState).resample()
        [_, testData] = data[i]
        TrainDataList.append(testData)
        result = prediction4stage2(
            data[i], 
            predictionMode, 
            tuneHyperParam,
            hyperParam).predictionMode()
        resultList.append(result)
    TrainData = reduce(lambda df1,df2: pd.concat([df1, df2], axis=0), TrainDataList)
    prediction = reduce(lambda df1,df2: pd.concat([df1, df2], axis=0), resultList)
    prediction = prediction.astype(float)
    prediction = prediction.join(age[["Age"]])
    prediction = prediction.applymap(lambda x: antiAgeTransfer(x))
    print(prediction)
    # prediction.to_csv("{}.csv".format(resultName))
    TrainData.to_csv(resultPath.format(resultName) + "TrainingCV_Data4Stage2.csv")
    prediction.to_csv(resultPath.format(resultName) + "TrainingCV_Prediction.csv")



def stage2pediction(
    predict: pd.DataFrame,
    model: pd.DataFrame, 
    resultName: str,
    # nfold: int,
    # randomState: int,
    predictionMode: str,
    tuneHyperParam: boolean,
    hyperParam: dict,
):

    data = [model[list(predict.columns)], predict]
    prediction = prediction4stage2(
        predictionModule = predictionMode, 
        data = data, 
        tuneHyperParam = tuneHyperParam, 
        hyperParam =hyperParam).predictionMode()
    prediction = prediction.join(predict[["Age"]])
    prediction = prediction.applymap(lambda x: antiAgeTransfer(x))

    print("prediction", prediction)
    # prediction.to_csv("{}.csv".format(resultName))
    prediction.to_csv(resultPath.format(resultName) + "Testing_Prediction.csv")




