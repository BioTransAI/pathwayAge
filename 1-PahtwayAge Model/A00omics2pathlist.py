import pandas as pd
from typing import List, Optional
import concurrent.futures

def omics2pathlist(
    data: pd.DataFrame,
    pathDict: dict,
    featureAnno: pd.DataFrame,
    restrictUp: int,
    restrictDown: int,
    minPathSize: int
):
  """omics2pathlist
    create a list of dataframe, each dataframe contains CpGs of one pathway.
  Parameters
  ----------
    data: methylation matrix with Age info;
    pathDict: dictionary key is the pathway, value is the genes associated with the pathway;
    featureAnno: contains gene "entrezID" and the CpG sites "ID" aroud its gene;
    restrictUp: the upper limit of number of genes in one pathway;
    restrictDown: the lower limit of number of genes in one pathway;
    minPathSize: CpG sites minimal number in each pathway;
  Return
  ----------
    a list of dataframe. 
  """ 
  ### todo
  ### check the type of the input data
  assert "Age" in data.columns, "No 'Age' column in 'methlyData'!"
  assert "ID" in featureAnno.columns, "No 'ID' column in 'CpGAnno'!"
  assert "entrezID" in featureAnno.columns, "No 'entrezID' column in 'CpGAnno'! "
  methylPerPathList = []
  for (key, value) in pathDict.items():
     if restrictDown<= len(value)<= restrictUp:
       cpgIDPerPath = featureAnno.loc[featureAnno["entrezID"].isin(value)]
       methylPerPath = data[data.columns[data.columns.isin(cpgIDPerPath["ID"].values)]].copy()
       if len(methylPerPath.columns) > minPathSize:  
          methylPerPath["GOName"] = key
          methylPerPath["Age"] = data["Age"].values
          methylPerPathList.append(methylPerPath)
  # print("methyl2PathList data is ready! {} pathways in total".format(len(methylPerPathList)))
  return methylPerPathList

def pathSummary(
    data: list
):
  print("Remove nbumber of too small siza path, the retained pathway {}"
        .format(len(data))
       )
  ### -2 because ["GOName", "Age"]
  summary = pd.DataFrame([len(item.columns)-2 for item in data])
  return summary[0].describe().apply("{0:.1f}".format)



