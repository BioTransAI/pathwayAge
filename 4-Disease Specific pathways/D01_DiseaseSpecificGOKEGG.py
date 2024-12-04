
import pandas as pd 
import numpy as np
from scipy.stats import  mannwhitneyu
from functools import reduce


def Z_scroe(x, y, U_Test):
    # Sizes of the two samples
    n_x = len(x)
    n_y = len(y)
    # Mean (μ_U) and standard deviation (σ_U) of U under the null hypothesis
    mean_U = n_x * n_y / 2
    std_U = np.sqrt(n_x * n_y * (n_x + n_y + 1) / 12)
    # Z-score calculation
    Z = (U_Test - mean_U) / std_U
    return Z

def loadingData():
    controlRanking = pd.read_csv("./Result/AgingTrainingRanking.csv", index_col = 0 )
    GO = pd.read_csv("./Result/AgeAccPerGO.csv", index_col = 0)
    GO = GO.drop(columns = ["Tag", "Cohort"]).astype(float).join(GO[["Tag", "Cohort"]])
    KEGG = pd.read_csv("./Result/AgeAccPerKEGG.csv", index_col = 0)
    KEGG = KEGG.drop(columns = ["Tag", "Cohort"]).astype(float).join(GO[["Tag", "Cohort"]])
    return controlRanking, GO, KEGG

def SignificatePathway():
    DiseaseListGO = []
    DiseaseListKEGG = []
    TestList = list(disease_Map.keys())
    controlRanking, GO, KEGG = loadingData()
    for cohort in TestList:
        SubGO = GO[GO.Cohort.eq(cohort)].dropna()
        if SubGO.shape[0] > 1:
            U_TestList = []
            pList = []
            meanControllist = []
            meanCaseList = []
            medianControllist = []
            medianCaseList = []
            colunmNames = list(SubGO.drop(columns = ["Tag",	"Cohort"]).columns)
            for col in colunmNames:
                control = SubGO[SubGO.Tag.eq("Control")][col].fillna(0)
                case = SubGO[SubGO.Tag.eq("Case")][col].fillna(0)
                U_Test, p = mannwhitneyu(control, case)
                ZScore = Z_scroe(control, case, U_Test)
                meanControl = round(control.mean(),3)
                meanCase = round(case.mean(),3)
                medianControl = round(control.median(),3)
                medianCase = round(case.median(),3)
                U_TestList.append(ZScore)
                pList.append(p)
                meanControllist.append(meanControl)
                meanCaseList.append(meanCase)
                medianControllist.append(medianControl)
                medianCaseList.append(medianCase)
            Data = pd.DataFrame(data=[pList, U_TestList, meanControllist, meanCaseList, medianControllist, medianCaseList], 
                                columns=colunmNames, index=["P", "ZScore", "Control Mean", "Case Mean", "Control Median", "Case Median"])
            Disease = controlRanking.join(Data.T)
            Disease["Cohort"] = cohort
            DiseaseListGO.append(Disease)
        
    for cohort in TestList:
        SubKEGG = KEGG[KEGG.Cohort.eq(cohort)].dropna()
        if SubKEGG.shape[0] > 1:
            U_TestList = []
            pList = []
            meanControllist = []
            meanCaseList = []
            medianControllist = []
            medianCaseList = []
            colunmNames = list(SubKEGG.drop(columns = ["Tag",	"Cohort"]).columns)
            for col in colunmNames:
                control = SubKEGG[SubKEGG.Tag.eq("Control")][col].fillna(0)
                case = SubKEGG[SubKEGG.Tag.eq("Case")][col].fillna(0)
                U_Test, p = mannwhitneyu(control, case)
                ZScore = Z_scroe(control, case, U_Test)
                meanControl = round(control.mean(),3)
                meanCase = round(case.mean(),3)
                medianControl = round(control.median(),3)
                medianCase = round(case.median(),3)
                U_TestList.append(ZScore)
                pList.append(p)
                meanControllist.append(meanControl)
                meanCaseList.append(meanCase)
                medianControllist.append(medianControl)
                medianCaseList.append(medianCase)
            Data = pd.DataFrame(data=[pList, U_TestList, meanControllist, meanCaseList, medianControllist, medianCaseList], 
                                columns=colunmNames, index=["P", "ZScore", "Control Mean", "Case Mean", "Control Median", "Case Median"])
            Disease = controlRanking.join(Data.T)
            Disease["Cohort"] = cohort
            DiseaseListKEGG.append(Disease)

    GOData = reduce(lambda df1,df2: pd.concat([df1, df2], axis=0), DiseaseListGO)
    KEGGData = reduce(lambda df1,df2: pd.concat([df1, df2], axis=0), DiseaseListKEGG)
    
    return GOData, KEGGData



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
