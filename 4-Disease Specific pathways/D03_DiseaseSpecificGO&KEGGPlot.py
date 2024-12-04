import pandas as pd 
from D01_DiseaseSpecificGOKEGG import SignificatePathway
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import leaves_list
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import matplotlib.pyplot as plt


def Data4Plot():
    GOData, _ = SignificatePathway()
    GOData['Differ'] = GOData["Case Median"] - GOData["Control Median"]
    DataSub = GOData[GOData.P.lt(0.05) & GOData.Differ.gt(0)]
    top_n = DataSub.groupby('Cohort', group_keys=False).apply(lambda x: x.nlargest(3000, 'Case Median'))
    top_n["Disease"] = top_n.Cohort.replace(disease_Map)
    top_n["Index"] = top_n.index 
    CaseScore = pd.pivot_table(top_n, values='Differ', index='Index', columns='Cohort')
    CaseScore = CaseScore.reset_index()
    CaseScore.index = CaseScore["Index"]
    CaseScore = CaseScore.drop(columns = ["Index"])
    return CaseScore

def SpecificPathwayPlot():

    CaseScore = Data4Plot()
    scaler = StandardScaler()
    result = scaler.fit_transform(CaseScore)
    CaseScore_scale = pd.DataFrame(result, columns=CaseScore.columns, index = CaseScore.index)

    linked = linkage(CaseScore_scale, method='ward')

    ordered_indices = leaves_list(linked)
    Ranking_df = CaseScore_scale.iloc[ordered_indices, :]
    CaseScore_scale.index.name = None

    sns.set_theme()
    g = sns.clustermap(Ranking_df, center=0, cmap="vlag", 
                    method='ward', 
                    dendrogram_ratio=(.1, .2),
                    cbar_pos=(.95, .3, .03, .2),
                    figsize=(10, 12),
                    )
    g.ax_heatmap.set_yticks([])
    plt.title('AgeAcc(GO)')
    plt.savefig('./Result/AgeAcc(GO).png', bbox_inches='tight')  # Save as PNG with high resolution
    plt.show()

disease_Map = {
    "GSE80417": "SCZ1",
    "GSE84727": "SCZ2",
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