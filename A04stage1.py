import concurrent.futures
from xmlrpc.client import boolean
from itertools import repeat
import pandas as pd
from typing import List, Optional
import time
from functools import reduce
import os
from pathlib import Path
from A03data4stage2 import dataForStage2, reconTrainTestData


def stage1(
omics2pathlist,
age, 
nfold: int,
randomState: int,
predictionMode: str,
reconData: boolean,
tuneHyperParam: boolean,
hyperParam: dict,
cores: int,
methylTestData: Optional[pd.DataFrame] = pd.DataFrame(), 
):
    cvList = []
    if reconData: 
        for i in range(nfold):
            with concurrent.futures.ProcessPoolExecutor(max_workers = cores) as executor:
                predictionMeanList = list(executor.map(dataForStage2,
                                                omics2pathlist,
                                                repeat(predictionMode),
                                                repeat(tuneHyperParam),
                                                repeat(hyperParam),
                                                repeat(i),
                                                repeat(reconData),
                                                repeat(nfold),
                                                repeat(randomState)
                                                ))

            data4Stage2 = reduce(lambda df1,df2: pd.concat([df1, df2], axis=1), predictionMeanList)
            data4Stage2 = data4Stage2.join(age[["Age"]])
            print(data4Stage2)
            # data4Stage2.to_csv(outerPath.format(i))  
            cvList.append(data4Stage2)
        print("data4Stage2 datasets keeping the outer data unseen complete!")

    if  methylTestData.empty:
        with concurrent.futures.ProcessPoolExecutor(max_workers=cores) as executor:
            predictionMeanList = list(executor.map(dataForStage2,
                                            omics2pathlist,
                                            repeat(predictionMode),
                                            repeat(tuneHyperParam),
                                            repeat(hyperParam),
                                            repeat(None),
                                            repeat(None),
                                            repeat(nfold),
                                            repeat(randomState)
                                            ))

        data4Stage2 = reduce(lambda df1,df2: pd.concat([df1, df2], axis=1), predictionMeanList)
        data4Stage2 = data4Stage2.join(age[["Age"]])
        print(data4Stage2)
        # data4Stage2.to_csv(path.format("data4Stage2.csv"))
        print("data4Stage2 datasets all data been seen complete!")
    else: 
        data4Stage2 = pd.DataFrame()
    
    return data4Stage2, cvList