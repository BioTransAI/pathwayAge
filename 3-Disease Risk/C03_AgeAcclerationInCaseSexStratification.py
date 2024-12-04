import pandas as pd
from functools import reduce
from scipy.stats import mannwhitneyu
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression

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


# import Data
def loadingData():

    GOPredictionFemale = pd.read_csv("./Result/PredictionGO_Female.csv", index_col = 0)
    GOPredictionMale = pd.read_csv("./Result/PredictionGO_Male.csv", index_col = 0)
    KEGGPredictionFemale = pd.read_csv("./Result/PredictionKEGG_Female.csv", index_col = 0)
    KEGGPredictionMale = pd.read_csv("./Result/PredictionKEGG_Male.csv", index_col = 0)

    GOPredictionFemale = GOPredictionFemale.rename(columns= {"prediction": "PathwayAge_Female(GO)"})
    GOPredictionMale = GOPredictionMale.rename(columns= {"prediction": "PathwayAge_Male(GO)"})
    KEGGPredictionFemale = KEGGPredictionFemale.rename(columns= {"prediction": "PathwayAge_Female(KEGG)"})
    KEGGPredictionMale = KEGGPredictionMale.rename(columns= {"prediction": "PathwayAge_Male(KEGG)"})

    resultList = [
        GOPredictionFemale,
        GOPredictionMale,
        KEGGPredictionFemale,
        KEGGPredictionMale  
    ]
    pathwayAge = reduce(lambda df1,df2: pd.concat([df1.join(df2.drop(columns = ["Age"]))], axis = 1),  resultList)
    pathwayAge = pathwayAge[pathwayAge.Age.ge(18) & pathwayAge.Age.le(80)]
    print(pathwayAge)

    confunder = pd.read_csv("./Demo Meta Data/CovariateData.csv")
    pathwayAgePredcition = pathwayAge.join(confunder[["Label", "Cohort"]])
    pathwayAgePredcition["Label"] = pathwayAgePredcition["Label"].apply(pd.to_numeric, errors='coerce')
    return pathwayAgePredcition


# Test whether the case has a significant effect on age acceleration.
def testAgeAccInCase():
    methodList = ["PathwayAge_Female(GO)", "PathwayAge_Male(GO)", "PathwayAge_Female(KEGG)", "PathwayAge_Male(KEGG)"]
    resultList = []
    deltaAgeList = []
    pathwayAgePredcition = loadingData()
    for method in methodList:
        TestSub = pathwayAgePredcition.rename(columns = {method: "prediction"})
        TestSub = TestSub[["prediction", "Age", "Label","Cohort"]]
        DataList = []
        for cohort in pathwayAgePredcition.Cohort.unique():
            prediction = TestSub[TestSub.Cohort.eq(cohort)]
            controlOriginal = prediction[prediction.Label.eq(0)].drop(columns= ["Label", "Cohort"]).dropna()
            caseOriginal = prediction[prediction.Label.eq(1)].drop(columns= ["Label", "Cohort"]).dropna()
            
            if len(prediction.Label.unique()) > 1:
                confunderTest = cofunder[[
                    'Smoking', "Age", "Batch",
                    'CD8.naive','CD8pCD28nCD45RAn', 'PlasmaBlast', 'CD4T', 'NK', 'Mono', 'Gran',
                    'population_0', 'population_1', 'population_2', 'population_3', 'population_4', 'population_5',
                    'population_6', 'population_7', 'population_8', 'population_9'
                    ]]
                cofunder  = confunderTest.drop(columns = ["Age"])
                Full = cofunder
                cellComposition = ['CD8.naive', 'CD8pCD28nCD45RAn', 'PlasmaBlast', 'CD4T', 'NK', 'Mono', 'Gran', "Batch"]

                if cohort in [ 
                    "GSE87640_monocytes",  
                    "GSE87640_CD8", 
                    "GSE87640_CD4",  
                    "GSE71955_CD4 T cells",
                    "GSE71955_CD8 T cells"
                ]:
            
                    controlFull = controlOriginal.join(Full.drop(columns = cellComposition)).dropna()
                    caseFull = caseOriginal.join(Full.drop(columns = cellComposition)).dropna()
                    controlFull = LR(controlFull)
                    caseFull = LR(caseFull) 
    
        
                else:
                    controlFull = controlOriginal.join(Full).dropna()
                    caseFull = caseOriginal.join(Full).dropna()
                    controlFull = LR(controlFull)
                    caseFull = LR(caseFull)
            
                controlFull, caseFull = AgeAccelerationResidual(controlFull, caseFull, "Full")
            
                caseFull["Tag"] = "Case"
                controlFull["Tag"] = "Control"
        
                deltaAge = pd.concat([controlFull[["Full", "Tag"]], caseFull[["Full", "Tag"]]])
                deltaAge["Cohort"] = cohort
                deltaAge["Method"] = method
                deltaAgeList.append(deltaAge)
        
                colunmNames = ["Full"]
                pList = []
                meanControllist = []
                meanCaseList= []
                
                for col in colunmNames:
                    _, p = mannwhitneyu(controlFull[col].dropna(), caseFull[col].dropna())
                    meanControl = round(controlFull[col].median(),3)
                    meanCase = round(caseFull[col].median(),3)
                    differ = meanCase - meanControl
                    pList.append(p)
                    meanControllist.append(meanControl)
                    meanCaseList.append(differ)
        
                Data = pd.DataFrame(data=[pList, meanCaseList], columns=colunmNames, index=["P", "Differ"]).round(3)
                Data["Data"] = cohort
                DataList.append(Data)
        result = reduce(lambda df1,df2: pd.concat([df1, df2]),  DataList)
        result[method] = result.Full
        resultList.append(result[[method, "Data"]])
    AllClock = reduce(lambda df1,df2: pd.concat([df1.drop(columns = ["Data"]), df2], axis = 1),  resultList)
    print(AllClock)
    print(deltaAgeList.shape)

    return AllClock, deltaAgeList