# to test pathwayAge
from pathwayAge import pathwayAge
import pandas as pd 
import time

# 1. no test data

# methylData: pd.DataFrame,
# methylData = "./data/methlyData.csv"
methylData = "/home/lipan231/code/pathwayAgeCluster/data/BrainMethylData.csv"
# covariateData =  "./data/covariateData.csv"
covariateData =  "/home/lipan231/code/pathwayAgeCluster/data/BrainCovariateData.csv"
methylData = pd.read_csv(methylData, index_col="CpG")
methylData = methylData.T
covariateData = pd.read_csv(covariateData, index_col ="Sample")
age = covariateData["Age"]
methylDataAge = methylData.iloc[:100,:10000]

print(methylDataAge)

# resultName: str,
resultName = "test"

# methylTestData: Optional[pd.DataFrame] = pd.DataFrame(),
methylTestData = methylData.iloc[101:201,:10000]
methylTestData = methylTestData

# restrictUp: Optional[int] = 200,
restrictUp = 20

# restrictDown: Optional[int] = 10,
restrictDown = 10

# minPathSize: Optional[int] = 5,
minPathSize = 10

# nfold: Optional[int] = 5,
nfold = 6

# randomState: Optional[int] = 6677,
randomState = 00

# predictionMode: Optional[str] = "GradientBoosting",
predictionMode = "Ridge"

# reconData: Optional[boolean] = False,
reconData = True

# tuneHyperParam: Optional[boolean] = False,
tuneHyperParam = True

# hyperParam: Optional[dict]= None,
    # if tuneHyperParam == False
### GB ###
# hyperParam = {
#     "n_estimators": 50,
#     "max_depth": 4,
#     "min_samples_split": 5,
#     "learning_rate": 0.01,
#     "loss": "squared_error",
# }
# ### SVM ###
# hyperParam = {
#     "kernel": ["rbf", "poly"],
#     "C": [20, 10, 1, 0.1],
#     # "epsilon": [0.001, 0.01, 0.1, 1]
#         }

hyperParam = {'alphas': [0.01, 0.1, 1, 10]}

# cores: Optional[int]= 5,
cores = 10

startread = time.time()

pathwayAge(
    methylData = methylDataAge,
    resultName = resultName,
    minPathSize = 10,
    nfold = nfold,
    tuneHyperParam  = tuneHyperParam,
    hyperParam = hyperParam,
    cores= cores,
    methylTestData = methylTestData,
    predictionMode = predictionMode,
)

endDeal = time.time()
print("Data dealed in {} seconds".format(endDeal - startread))