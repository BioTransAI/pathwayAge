
from A05stage2 import stage2pediction
import pandas as pd
import random
from scipy.stats import mannwhitneyu
from functools import reduce
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline

sampleMode = "CrossValidation"
predictionMode = "Ridge" 
hyperParam = {
    'alphas': [0.001, 0.01, 0.1, 0.5, 1, 10, 100, 1000,]
}
Save = "./Result/CaseTop20&Permutation/{}"


def loadingData():
    AgingTop20 = pd.read_excel("./Result/Aging_top_20_GO.xlsx", index_col = 0)
    AgingTop20["Disease"] = AgingTop20["Cohort"].map(disease_Map)

    AgingGO = pd.read_csv("./Result/Data4Stage2Sub.csv", index_col = 0)
    GOs = list(AgingGO.drop(columns = list(AgingTop20.index) + ["Age"]))
    return GOs

def random_sample_with_different_seeds(data_list, num_elements=20, iterations=100):
    results = []
    for i in range(iterations):
        seed = i + 1  
        random.seed(seed)
        sample = random.sample(data_list, num_elements)
        sample.append("Age")
        results.append(sample)
    return results

def RandomGOPathwahAgePrediction():
    GOs = loadingData()
    samples = random_sample_with_different_seeds(GOs)
    PredictionList = []
    for i in range(100):
        print(i)
        predictFilePath = "./Result/Data4Stage2TestinSub.csv"
        modelFilePath = "./Result/Data4Stage2Sub.csv"
        savePath = Save.format("AgingTop20Case_Random_GO_"  + str(i) + "_Prediction.csv")
        Prediction = stage2pediction(
            predictFilePath = predictFilePath,
            modelFilePath = modelFilePath,
            predictionMode = predictionMode,
            hyperParam= hyperParam,
            savePath = savePath,
            GOList = samples[i]
        )
        PredictionList.append(Prediction)
    print("Done")
    return PredictionList

def Cohen(control, case, col):
    # Compute means and standard deviations
    mean_group1 = control[col].mean()
    mean_group2 = case[col].mean()
    std_group1 = control[col].std()
    std_group2 = case[col].std()
    
    # Compute Cohen's d
    pooled_std = np.sqrt(((len(control) - 1) * std_group1**2 + (len(case) - 1) * std_group2**2) / (len(control) + len(case) - 2))
    cohens_d = (mean_group2 - mean_group1) / pooled_std
    
    # print("Cohen's d:", cohens_d)
    return cohens_d

def effectSizeRandomGOPathwahAge():
    PredictionList = RandomGOPathwahAgePrediction()
    permutationTest = []
    iteration = 0 
    confunder = pd.read_csv("./Demo Meta Data/CovariateData.csv")
    for df in PredictionList:
        PathwayAge = df.join(confunder[["Label", "Cohort"]])
        DataList = []
        for cohort in disease_Map.keys:
            prediction = PathwayAge[PathwayAge.Cohort.eq(cohort)]
            controlOriginal = prediction[prediction.Label.eq(0)].drop(columns= ["Label", "Cohort"]).dropna()
            caseOriginal = prediction[prediction.Label.eq(1)].drop(columns= ["Label", "Cohort"]).dropna()
            
            if len(prediction.Label.unique()) > 1:
                confunder = pd.read_csv("./Demo Meta Data/CovariateData.csv")
                confunderTest = confunder[[
                    'Smoking', 'Female', "Age", "Batch",
                    'CD8.naive','CD8pCD28nCD45RAn', 'PlasmaBlast', 'CD4T', 'NK', 'Mono', 'Gran',
                    'population_0', 'population_1', 'population_2', 'population_3', 'population_4', 'population_5',
                    'population_6', 'population_7', 'population_8', 'population_9'
                    ]]
                confunderTest  = confunderTest.drop(columns = ["Age"])
                Full = confunderTest
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
                
                colunmNames = ["Full"]
                effectSizeList = []
                pList = []
                meanControllist = []
                meanCaseList= []
                
                for col in colunmNames:
                    _, p = mannwhitneyu(controlFull[col].dropna(), caseFull[col].dropna())
                    effectSize = Cohen(controlFull, caseFull, col)
                    meanControl = round(controlFull[col].mean(),3)
                    meanCase = round(caseFull[col].mean(),3)
                    pList.append(p)
                    effectSizeList.append(effectSize)
                    meanControllist.append(meanControl)
                    meanCaseList.append(meanCase)
        
                Data = pd.DataFrame(data=[effectSizeList, pList, meanCaseList], columns=colunmNames, index=["effectSize", "P", "Case Mean"])
                Data["Data"] = cohort
                DataList.append(Data)
        iteration +=1 
        result = reduce(lambda df1,df2: pd.concat([df1, df2]),  DataList)
        result = result.drop(columns = ["AA", "IEAA", "Basic"])
        result["Index"] = result.index
        result["Index_new"] = result.Data + "_"+ result.Index
        result.index = result.Index_new
        result = result.drop(columns =["Index_new", "Index"] )
        result = result.rename(columns = {"Full": "PermutationTest_" + str(iteration)})
        permutationTest.append(result)
    permutation = reduce(lambda df1,df2: df1.drop(columns = ["Data"]).join(df2),  permutationTest)
    permutation["Disease"] = permutation.Data.replace(disease_Map)

    return permutation

def effectSizeTopGOPathwahAge():
    confunder = pd.read_csv("./Demo Meta Data/CovariateData.csv")
    PathwayAge = pd.read_csv("./Result/PredictionAge.csv")
    PathwayAge = PathwayAge.join(confunder[["Label", "Cohort"]])
    DataList = []
    for cohort in disease_Map.keys:
        prediction = PathwayAge[PathwayAge.Cohort.eq(cohort)]
        controlOriginal = prediction[prediction.Label.eq(0)].drop(columns= ["Label", "Cohort"]).dropna()
        caseOriginal = prediction[prediction.Label.eq(1)].drop(columns= ["Label", "Cohort"]).dropna()
        
        if len(prediction.Label.unique()) > 1:
            confunder = pd.read_csv("./Demo Meta Data/CovariateData.csv")
            confunderTest = confunder[[
                'Smoking', 'Female', "Age", "Batch",
                'CD8.naive','CD8pCD28nCD45RAn', 'PlasmaBlast', 'CD4T', 'NK', 'Mono', 'Gran',
                'population_0', 'population_1', 'population_2', 'population_3', 'population_4', 'population_5',
                'population_6', 'population_7', 'population_8', 'population_9'
                ]]
            confunderTest  = confunderTest.drop(columns = ["Age"])
            Full = confunderTest
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
            
            colunmNames = ["Full"]
            effectSizeList = []
            pList = []
            meanControllist = []
            meanCaseList= []
            
            for col in colunmNames:
                _, p = mannwhitneyu(controlFull[col].dropna(), caseFull[col].dropna())
                effectSize = Cohen(controlFull, caseFull, col)
                meanControl = round(controlFull[col].mean(),3)
                meanCase = round(caseFull[col].mean(),3)
                pList.append(p)
                effectSizeList.append(effectSize)
                meanControllist.append(meanControl)
                meanCaseList.append(meanCase)
    
            Data = pd.DataFrame(data=[effectSizeList, pList, meanCaseList], columns=colunmNames, index=["effectSize", "P", "Case Mean"])
            Data["Data"] = cohort
            DataList.append(Data)
    result = reduce(lambda df1,df2: pd.concat([df1, df2]),  DataList)
    result["Index"] = result.index
    result["Index_new"] = result.Data + "_"+ result.Index
    result.index = result.Index_new
    result = result.drop(columns =["Index_new", "Index"] )

    return result

def PermutationTest():
    Aging_EffectSize_Permutation = effectSizeRandomGOPathwahAge()
    Aging_EffectSize = effectSizeTopGOPathwahAge()
    Aging_EffectSize_Permutation = Aging_EffectSize_Permutation.join(Aging_EffectSize[["Full", "Index"]])
    Aging_EffectSize_Permutation = Aging_EffectSize_Permutation[Aging_EffectSize_Permutation.Index.eq("effectSize")]
    Aging_EffectSize_Permutation = Aging_EffectSize_Permutation.drop(columns = ["Data", "Disease", "Index"])
    Aging_EffectSize_Permutation['count_less_than_target'] = Aging_EffectSize_Permutation.apply(lambda row: (row < row['Full']).sum() - 1, axis=1)
    Aging_EffectSize_Permutation['Permutation_P'] = Aging_EffectSize_Permutation['count_less_than_target'].apply(Permutation_P)
    Aging_EffectSize_Permutation = Aging_EffectSize_Permutation[["Full", "P", "count_less_than_target","count_less_than_target", "Permutation_P"]].rename(columns = {"Full": "EffectSize"}).round(3)
    return Aging_EffectSize_Permutation

def Permutation_P(f, k = 100):
    return (f+1)/(k+1)

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

disease_Map = {
    "GSE80417": "SCZ1",
    "GSE84727":	"SCZ2",
    "GSE152026": "SCZ3",
    "GSE152027": "SCZ4",
    "GSE144858": "AD",
    "GSE111629": "PD1",
    "GSE72774": "PD2",
    "GSE145361": "PD3",
    "GSE77696": "HIV",
    "GSE71955": "GD",
    "GSE87640": "CD",
    "GSE139404": "CRC",
    "GSE183040": "PC",
    "GSE222595": "OB",
   }