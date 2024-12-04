
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

def evaluationControl():

    predictionResult = "./Result/PredictionAge.csv"
    predictionResult = pd.read_csv(predictionResult, index_col=0)

    cofunder = "./Demo Meta Data/CovariateData.csv"
    cofunder = pd.read_csv(cofunder, index_col="Sample")

    label = cofunder[["Label"]].astype(int)
    prediction = label.join(predictionResult)
    controlOriginal = prediction[prediction.Label.eq(0)].drop(columns="Label")
    mean_squared_error, mean_absolute_error, r2_score, dataCorr = evaluation(controlOriginal)
    data = [mean_squared_error,mean_absolute_error,r2_score, dataCorr]
    Control = pd.DataFrame(data=data, columns = ["Training"], index=["MSE", "MAE", "R2", "Rho"])
    Control["Tag"] = "Control" 

    return Control

def AgeAcc(): 

    cofunder = "./Demo Meta Data/CovariateData.csv"
    cofunder = pd.read_csv(cofunder, index_col="Sample")
    cofunder = cofunder.drop(columns = ["Age", "Label"])

    predictionResult = "./Result/PredictionAge.csv"
    predictionResult = pd.read_csv(predictionResult, index_col=0)

    label = cofunder[["Label"]].astype(int)
    prediction = label.join(predictionResult)
    controlOriginal = prediction[prediction.Label.eq(0)].drop(columns="Label")

    Full = cofunder
    basic = cofunder[["Female", "CohortTag"]]
    cellComposition = ['CD8.naive', 'CD8pCD28nCD45RAn', 'PlasmaBlast', 'CD4T', 'NK', 'Mono', 'Gran', 'CohortTag']
    cellCount = cofunder[cellComposition]

    controlAA = LR(controlOriginal)
    caseAA = controlAA

    controlCellCount = controlOriginal.join(cellCount)
    controlIEAA = LR(controlCellCount)
    caseIEAA = controlIEAA

    controlBasic = controlOriginal.join(basic)
    controlBasic = LR(controlBasic)
    caseBasic = controlBasic

    controlFull = controlOriginal.join(Full)
    controlFull = LR(controlFull)
    caseFull = controlFull


    controlAA, caseAA = AgeAccelerationResidual(controlAA, caseAA, "AA")
    controlBasic, caseBasic = AgeAccelerationResidual(controlBasic, caseBasic, "Basic")
    controlIEAA, caseIEAA = AgeAccelerationResidual(controlIEAA, caseIEAA, "IEAA")
    controlFull, caseFull = AgeAccelerationResidual(controlFull, caseFull, "Full")

    controlAA["Basic"] = controlBasic["Basic"]
    caseAA["Basic"] = caseBasic["Basic"]
    controlAA["IEAA"] = controlIEAA["IEAA"]
    caseAA["IEAA"] = caseIEAA["IEAA"]
    controlAA["Full"] = controlFull["Full"]
    caseAA["Full"] = caseFull["Full"]

    caseAA["Tag"] = "Case"
    controlAA["Tag"] = "Control"

    ageGapFullAdjusted = controlAA
    print(ageGapFullAdjusted)

    return ageGapFullAdjusted

def AgeAccPerGO():
    AgeGapList = []
    data4Stage2 = './Result/Data4Stage2Sub.csv'
    data4Stage2 = pd.read_csv(data4Stage2, index_col=0)   
    data4Stage2 = data4Stage2.drop(columns= ["Age"]) 
    print(data4Stage2)  

    cofunder = "./Demo Meta Data/CovariateData.csv"
    cofunder = pd.read_csv(cofunder, index_col="Sample") 

    for col in data4Stage2.columns:
        Full = cofunder
        controlFull = data4Stage2[[col]].join(Full).dropna()   
        controlFull = LRWithCol(controlFull, col)
        controlFull = controlFull.astype(float)
        controlFull["Full"] = controlFull[col] - controlFull["LRPrediction"]
        controlFull = controlFull[["Full"]].rename(columns={"Full": col})
        print(col)
        AgeGapList.append(controlFull)
    AgeGapGO = pd.concat(AgeGapList, axis=1)
    AgeGapGO.to_csv("./Result/AgeAccPerGOSub.csv")
    print(AgeGapGO)
    return AgeGapGO

def AgeAccCorrWithGO():
    AgeAccControl = AgeAcc()
    AgeGapGO = AgeAccPerGO()

    control = AgeAccControl[AgeAccControl["Tag"].eq("Control")][["Full"]].join(AgeGapGO)
    correlationControl = control.drop(columns= ["Full"]).corrwith(control["Full"])
    correlationControl = pd.DataFrame(data=correlationControl.values, 
                                    index=correlationControl.index,
                                    columns=["Rho"])

    readRDS = r['readRDS']
    goDescription = readRDS('./Demo Meta Data/GO_Description.rds')
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

def GORankingPlot():
    # data for plot Bubble plot
    golist = "./Demo Meta Data/golist.json"
    with open(golist) as json_file:
        golist = json.load(json_file)

    correlationControl = AgeAccCorrWithGO()
    correlationControlTop20 = correlationControl.nlargest(20,"RhoAbs")
    correlationControlTop20["GeneCount"] = [len(value) for key, value in golist.items() for index in correlationControlTop20.index if (key == index)]
    plotTag = pd.read_excel("./Demo Meta Data/plotTag.xlsx", index_col=0)
    correlationControlTop20 = correlationControlTop20.join(plotTag)
    print(correlationControlTop20)
    # correlationControlTop20.to_excel("Temp.xlsx")
    coloMap = {
        "developmental process":"#7FB3D5", 
        "metabolic process":"#BC8B56", 
        "response to stimulus": "#F6E758", 
        "cellular process" :"#6D936E",
        "biological regulation" :"#949b9b",
        "localization": "#8c7e78",
        "multicellular organismal process": "#3c4444"
        }

    correlationControlTop20["color"] = correlationControlTop20["Tag"].replace(coloMap)
    correlationControlTop20 = correlationControlTop20.sort_values(by = ["RhoAbs"], ascending= True)
    correlationControlTop20.to_csv("./Result/GORankSub.csv", index=False)



