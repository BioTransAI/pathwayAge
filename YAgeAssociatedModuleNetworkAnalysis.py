
import numpy as np
import pandas as pd 
from functools import reduce
import json
from flask import Flask, render_template

# step0 import data
root = "./discovery3KGOClusterSub/{}"
clusterRandIndexPath = root.format("moduleRandIndex.csv")
GOAncestorPath = "./analysis/GO_Ancestor.csv"
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
# step3 select cluster with the threhold at 70              
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
selectedCluster = data[data.percent.gt(50)]
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

# step5 calculate the correlation between the GOs
AgeGapPerGOPath = "discovery3KAgeGapPerGOSub.csv"
AgeGapPerGO = pd.read_csv(AgeGapPerGOPath, index_col=0)
AgeGapPerGO = AgeGapPerGO[list(GOs.GO.values)]
heatmap = AgeGapPerGO.corr()
upper_triangle_values =  heatmap.mask(np.triu(np.ones(heatmap.shape)).astype(bool)).fillna(0)
print(upper_triangle_values)

# step6 prepare data for plot
nodes = GOs[["GO", "Cluster"]].rename(columns={"GO": "id", "Cluster": "group"})
print(GOs.Cluster.unique())
nodes = nodes.to_dict(orient='records')
link = []
for index, row in upper_triangle_values.iterrows():
    for column in upper_triangle_values.columns:
        if row[column] < -0.25 or row[column] > 0.25:
            dict = {}
            # print(f"Row {index}, Column {column}: Value is not equal to 1 ({row[column]})")
            dict["source"] = index 
            dict["target"] = column
            dict["value"] =  row[column]
            link.append(dict)
result = {}
result["nodes"] = nodes
result["links"] = link
file_name = "./dataSub.json"

with open(file_name, "w") as json_file:
    json.dump(result, json_file)

# step7 plot in D3
app = Flask(__name__)
@app.route('/')
def index():
    return render_template('./index.html')
if __name__ == '__main__':
    app.run()