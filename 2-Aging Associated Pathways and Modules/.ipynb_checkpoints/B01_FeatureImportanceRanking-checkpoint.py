
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import yaml
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, r
import rpy2.robjects as robjects
import json

def linearRegressionScore(X, y):
    # using linear regression score function to compute the r2
    from sklearn.linear_model import LinearRegression
    model = LinearRegression()
    model.fit(X, y)
    r_squared = model.score(X, y)
    return r_squared

def LR(data):
    ### normalize the confounders, where u is the mean of the training samples
    dataCopy = data.copy()
    try:
        dataCopyX = dataCopy.drop(columns = ["prediction", "label"])
    except:
        dataCopyX = dataCopy.drop(columns = ["prediction"])
    model = make_pipeline(StandardScaler(), LinearRegression())
    # print(dataCopyX.columns)
    residual = model.fit(dataCopyX, dataCopy["prediction"])
    dataCopy["LRPrediction"] = residual.predict(dataCopyX)
    return dataCopy 

def LRWithCol(data, col):
    ### normalize the confounders, where u is the mean of the training samples
    dataCopy = data.copy()
    dataCopyX = dataCopy.drop(columns = [col])
    model = make_pipeline(StandardScaler(), LinearRegression())
    # print(dataCopyX.columns)
    residual = model.fit(dataCopyX, dataCopy[col])
    dataCopy["LRPrediction"] = residual.predict(dataCopyX)
    return dataCopy 

def evaluation(predictionAge):
    from sklearn.metrics import mean_absolute_error, mean_squared_error
    mean_squared_error = mean_squared_error(predictionAge["Age"].values, predictionAge["prediction"].values)
    mean_absolute_error = mean_absolute_error(predictionAge["Age"].values, predictionAge["prediction"].values)
    r2_score = linearRegressionScore(predictionAge[["Age"]], predictionAge[["prediction"]])
    dataCorr = predictionAge.drop(columns=["Age"]).corrwith(predictionAge['Age'], method='pearson')
    return mean_squared_error, mean_absolute_error, r2_score, dataCorr.values[0]

def AgeAccelerationResidual(
    control,
    case,
    method, #AA or IEAA or Basic or Full 
    ):
    model = make_pipeline(StandardScaler(), LinearRegression())
    control = control.astype(float)
    case = case.astype(float)
    dataX = control.drop(columns = ["prediction", "LRPrediction"])
    caseX = case.drop(columns = ["prediction", "LRPrediction"])  
    # print(dataX.columns)
    controlModel = model.fit(dataX, control["prediction"])
    case[method] = case["prediction"] - controlModel.predict(caseX) 
    control[method] = control["prediction"] - control["LRPrediction"]
    return control, case


def AgeAcc(predictionResult, cofunder): 

    label = cofunder[["Label"]].astype(int)
    prediction = label.join(predictionResult).dropna()
    controlOriginal = prediction[prediction.Label.eq(0)].drop(columns="Label")
    caseOriginal = prediction[prediction.Label.eq(1)].drop(columns="Label")
    caseOriginalShape = caseOriginal.shape[0]
    if caseOriginalShape == 0 :
        caseOriginal = controlOriginal

    cofunder = cofunder.drop(columns = ["Age", "Label", "Cohort"])
    Full = cofunder
    
    controlFull = controlOriginal.join(Full)
    controlFull = LR(controlFull)
    caseFull = controlFull

    caseFull = caseOriginal.join(Full)
    caseFull = LR(caseFull)
    
    controlFull, caseFull = AgeAccelerationResidual(controlFull, caseFull, "Full")

    controlFull["Tag"] = "Control"
    caseFull["Tag"] = "Case"

    ageGapFullAdjusted = pd.concat([controlFull, caseFull])

    if caseOriginalShape == 0 :
        ageGapFullAdjusted = controlFull

    return ageGapFullAdjusted[["Full", "Tag"]]

def AgeAccPerGO(data4Stage2, cofunder):
    AgeGapList = []

    for col in data4Stage2.columns:
        Full = cofunder.drop(columns = ["Label", "Cohort"])
        # print("confounder: ", Full.columns)
        controlFull = data4Stage2[[col]].join(Full).dropna()   
        controlFull = LRWithCol(controlFull, col)
        controlFull = controlFull.astype(float)
        controlFull["Full"] = controlFull[col] - controlFull["LRPrediction"]
        controlFull = controlFull[["Full"]].rename(columns={"Full": col})
        # print(col)
        AgeGapList.append(controlFull)
    AgeGapGO = pd.concat(AgeGapList, axis=1)
    AgeGapGO.to_csv("../Demo Results/AgeAccPerGO.csv")

    return AgeGapGO

def AgeAccCorrWithGO():
    AgeAccControl = AgeAcc()
    AgeGapGO = AgeAccPerGO()

    control = AgeAccControl[AgeAccControl["Tag"].eq("Control")][["Full"]].join(AgeGapGO)
    correlationControl = control.drop(columns= ["Full"]).corrwith(control["Full"])
    correlationControl = pd.DataFrame(data=correlationControl.values, 
                                    index=correlationControl.index,
                                    columns=["Rho"])


    goDescription = pd.read_csv('../Demo Meta Data/GO_Description.csv', index_col = 0)
    goDescription = dict(zip(goDescription.names, map(list,list(goDescription))))
    correlationControl["Description"] = [goDescription[keys][0] for keys in correlationControl.index]
    #check the correlation -+ 
    correlationControl['Direction'] = np.where(correlationControl["Rho"]<0, "Negtive", "Positive")
    correlationControl["RhoAbs"] = correlationControl["Rho"].abs()
    correlationControl = correlationControl.sort_values(by=["RhoAbs"], ascending=False)
    print(correlationControl.nlargest(20,"RhoAbs"))
    print(correlationControl[correlationControl.Direction.eq("Positive")].describe())
    print(correlationControl[correlationControl.Direction.eq("Negtive")].describe())

    return correlationControl





