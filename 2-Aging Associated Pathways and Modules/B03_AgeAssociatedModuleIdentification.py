
import numpy as np
import pandas as pd 
from functools import reduce
import json
from flask import Flask, render_template

# step0 import data
root = "../Demo Results/GOCluster/{}"
clusterRandIndexPath = root.format("moduleRandIndex.csv")
GOAncestorPath = "./Demo Meta Data/GO_Ancestor.csv"
clusterRandIndex = pd.read_csv(clusterRandIndexPath,  index_col=0)
GOAncestor = pd.read_csv(GOAncestorPath)

print(GOAncestor)

# step1 find the cutHeight and minCluster size
selectedParameter = clusterRandIndex.nlargest(1, "randIndex")
cutHeight = selectedParameter["cutHeight"].values[0]
minClusterSize = selectedParameter["minClusterSize"].values[0]
print(cutHeight, minClusterSize)

# step2 Cluster with fixed cutHeight and minCluster size 
ClusterWithFixedParaPath = root.format("_trueModuleCluster_" + str(cutHeight) + "_" + str(minClusterSize) + "_.csv")
cluster = pd.read_csv(ClusterWithFixedParaPath, index_col=0).rename(columns={"Label": "Cluster"})
cluster.GO = cluster.GO.str.replace(".", ":")
# step3 select cluster with the threhold at 60              
dataList = []

for clusterIndex in cluster.Cluster.unique():
    subCluster = cluster[cluster.Cluster.eq(clusterIndex)]
    data = GOAncestor[GOAncestor.GO.isin(subCluster.GO)]
    print(data)
    data_pivot = data.pivot(index='GO', columns='Ancestor', values='GO')
    percent = 100 - (data_pivot.isnull().sum() * 100 / len(data_pivot))
    missing_value_df = pd.DataFrame({'column_name': data_pivot.columns,
                                    'percent': percent})
    data = missing_value_df.nlargest(1, "percent")[["percent"]]
    data["Cluster"] = clusterIndex
    print(data)
    dataList.append(data)
data = reduce(lambda df1,df2: pd.concat([df1, df2]),  dataList).sort_values(by= ["Cluster"])
data["Count"] = cluster.groupby(["Cluster"]).agg("count")["GO"].values
data = data.sort_values(by=["percent"])
selectedCluster = data[data.percent.gt(60) & data.Count.le(100)]
print(selectedCluster)

# step4 find the GOs In the above cluster
listGO = []
for Cluster in selectedCluster.Cluster.unique():
    result =  GOAncestor[GOAncestor.GO.isin(cluster[cluster.Cluster.eq(Cluster)].GO)].drop_duplicates(subset=["Name"])
    result["Cluster"] = Cluster
    # print(selectedCluster[selectedCluster.Cluster.eq(Cluster)].index.values)
    result["Type"] = str(selectedCluster[selectedCluster.Cluster.eq(Cluster)].index[0])
    listGO.append(result)
GOs = reduce(lambda df1,df2: pd.concat([df1, df2]),  listGO).sort_values(by= ["Cluster"])
GOs = GOs.drop(columns=[ "Ancestor", "AncestorGO"])
print(GOs)

GOs.to_csv("../Demo Meta Data/GOWGCNACluster.csv")