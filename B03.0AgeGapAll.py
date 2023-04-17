
from scipy.stats import wilcoxon, mannwhitneyu
import pandas as pd
import numpy as np
import time
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri, r
from rpy2.robjects.conversion import localconverter
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
import yaml
from functools import reduce

startread = time.time()

def LR(data):
    """
    linear regression cofunders on biological age
    """
    ### normalize the confounders, where u is the mean of the training samples
    dataCopy = data.copy()
    dataCopyX = dataCopy.drop(columns = ["prediction"])
    model = make_pipeline(StandardScaler(), LinearRegression())
    residual = model.fit(dataCopyX, dataCopy["prediction"])
    dataCopy["LRPrediction"] = residual.predict(dataCopyX)
    return dataCopy 

def AgeAccelerationResidual(
 control,
 case,
 method, #AA or IEAA or Basic or Full 
):
    """
    compute the residuals
    """
    model = make_pipeline(StandardScaler(), LinearRegression())
    control = control.astype(float)
    case = case.astype(float)
    dataX = control.drop(columns = ["prediction", "LRPrediction"])
    caseX = case.drop(columns = ["prediction", "LRPrediction"])  
    print(caseX.columns)
    controlModel = model.fit(dataX, control["prediction"])
    case[method] = case["prediction"] - controlModel.predict(caseX) 
    control[method] = control["prediction"] - control["LRPrediction"]
    return control, case


#################################################################################################################################################################
#################################################################################################################################################################
#################################################################################################################################################################

readRDS = r['readRDS']

# import data paths 
with open('./config.yml', 'r') as file:
     root = yaml.safe_load(file)

pathCommon = root["pathway"]["common"]
phase3 = root["pathway"]["phase3"]
phase3Result = root["pathway"]["phase3Result"]


dataDict = {
    # could be extended 
    "Phase3":phase3, 
    } 
phaseDict = {
    "phase3": phase3Result
    }

for phase, result in phaseDict.items():
    DataList = []
    for data, path in dataDict.items():
        # if traning dataset
        control = result.format('resultControl/predictionAge.csv')
        controlOriginal = pd.read_csv(control, index_col = 0)
        # else 
        # control = result.format('resultGBControl/{}ControlModel{}ControlPrediction.csv'.format(phase, data))
        # controlOriginal = pd.read_csv(control, index_col = 0)[["mean", "Age"]].rename(columns= {"mean": "prediction"})
        case = result.format('resultControl/{}ControlModel{}CasePrediction.csv'.format(phase, data))
        caseOriginal = pd.read_csv(case, index_col = 0)[["mean", "Age"]].rename(columns= {"mean": "prediction"})
        controlOriginal = controlOriginal[controlOriginal["Age"].le(60)]
        caseOriginal = caseOriginal[caseOriginal["Age"].le(60)]
        
        covariateData =  path.format("covariateData.rds")    
        confounders = readRDS(covariateData)

        with localconverter(ro.default_converter + pandas2ri.converter):
            confounders = ro.conversion.rpy2py(confounders)
        confounders = confounders.drop(columns = ["Age", "label"])

        Full = confounders
        basic = confounders[["Female"]]
        cellComposition = ['CD8.naive', 'CD8pCD28nCD45RAn', 'PlasmaBlast', 'CD4T', 'NK', 'Mono', 'Gran', "Age"]
        cellCount = confounders[cellComposition]
   
        controlAA = LR(controlOriginal)
        caseAA = LR(caseOriginal)

        controlBasic = controlOriginal.join(basic)
        caseBasic = caseOriginal.join(basic)
        controlBasic = LR(controlBasic)
        caseBasic = LR(caseBasic)

        controlCellCount = controlOriginal.join(cellCount).drop(columns=["Age"])
        caseCellCount = caseOriginal.join(cellCount).drop(columns=["Age"])
        controlIEAA = LR(controlCellCount)
        caseIEAA = LR(caseCellCount)

        controlFull = controlOriginal.join(Full)
        caseFull = caseOriginal.join(Full)
        controlFull = LR(controlFull)
        caseFull = LR(caseFull)

        controlAA, caseAA = AgeAccelerationResidual(controlAA, caseAA, "AA")
        controlBasic, caseBasic = AgeAccelerationResidual(controlBasic, caseBasic, "Basic")
        controlIEAA, caseIEAA = AgeAccelerationResidual(controlIEAA, caseIEAA, "IEAA")
        controlFull, caseFull = AgeAccelerationResidual(controlFull, caseFull, "Full")
      
        controlEEAA = controlAA
        caseEEAA = caseAA
        controlEEAA["AA"] = controlAA["AA"]
        caseEEAA["AA"] = caseAA["AA"]
        controlEEAA["IEAA"] = controlIEAA["IEAA"]
        caseEEAA["IEAA"] = caseIEAA["IEAA"]
        controlEEAA["Basic"] = controlBasic["Basic"]
        caseEEAA["Basic"] = caseBasic["Basic"]
        controlEEAA["Full"] = controlFull["Full"]
        caseEEAA["Full"] = caseFull["Full"]

        caseEEAA["Tag"] = "Case"
        controlEEAA["Tag"] = "Control"
        deltaAge = pd.concat([controlEEAA, caseEEAA,])
        deltaAge = deltaAge.drop(columns= ["LRPrediction"])
        print(deltaAge)
        colunmNames = ["AA", "IEAA",  "Basic", "Full"]
        pList = []
        meanControllist = []
        meanCaseList = []
        for col in colunmNames:
            _, p = mannwhitneyu(controlEEAA[col], caseEEAA[col])
            meanControl = round(controlEEAA[col].mean(),3)
            meanCase = round(caseEEAA[col].mean(),3)
            pList.append(p)
            meanControllist.append(meanControl)
            meanCaseList.append(meanCase)

        Data = pd.DataFrame(data=[pList, meanControllist, meanCaseList], columns=colunmNames, index=["P", "Control Mean", "Case Mean"]).round(3)
        Data["Data"] = data
        DataList.append(Data)
    result = reduce(lambda df1,df2: pd.concat([df1, df2]),  DataList)
    result.columns = pd.MultiIndex.from_product([["BioMM{}".format(phase)], result.columns])
    print(result)
