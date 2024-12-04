import pandas as pd 
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt

def getOverlapPathways():
    Trans_GO = pd.read_csv("./Result/Transcription_GO_Ranking.csv", index_col = 0) 
    Methly_GO = pd.read_csv("./Result/Methylation_GO_Ranking.csv", index_col = 0 )
    # alternatively KEGG
    # Trans_KEGG = pd.read_csv("./Result/Transcription_All_KEGG.csv", index_col = 0) 
    # Methly_KEGG = pd.read_csv("./Result/Methylation_All_KEGG.csv", index_col = 0) 

    Trans_GO = Trans_GO.rename(columns = {"Rho" : "Trans_Rho"})
    Trans_Methly_GO = Methly_GO.join(Trans_GO.Trans_Rho)
    GOoverlap = list(set(Methly_GO.nsmallest(100, "Ranking").index) & set(Trans_GO.nlargest(100, "Rho").index))

    Trans_Methly_GO["Tag"] = np.where(Trans_Methly_GO.index.isin(list(Methly_GO.nsmallest(100, "Ranking").index)), "Top 100 GO in Methylome only", "")
    Trans_Methly_GO["Tag"] = np.where(Trans_Methly_GO.index.isin(list(Trans_GO.nlargest(100, "Trans_Rho").index)), "Top 100 GO in Transcriptome only", Trans_Methly_GO["Tag"])
    Trans_Methly_GO["Tag"] = np.where(Trans_Methly_GO.index.isin(GOoverlap), "Overlap of Top 100 GO terms", Trans_Methly_GO["Tag"])

    Trans_Methly_GO_clean = Trans_Methly_GO.dropna()
    Trans_Methly_GO_clean = Trans_Methly_GO_clean.rename(columns = {"Rho": "Rho in Methylome", "Trans_Rho": "Rho in Transcriptome"})
    return Trans_Methly_GO_clean


def plotOverlapPathways():
    Trans_Methly_GO_clean = getOverlapPathways()
    custom_palette = { "Overlap of Top 100 GO terms": '#E75454', "": '#B0B0B0', "Top 100 GO in Methylome only": "#007B9F", "Top 100 GO in Transcriptome only": "#5B9A6E"}
    f, ax = plt.subplots(figsize=(6.5, 6.5))
    sns.scatterplot(x="Rho in Methylome", 
                    y="Rho in Transcriptome", 
                    hue="Tag", 
                    palette=custom_palette, 
                    sizes=(1, 8), 
                    linewidth=0,
                    data=Trans_Methly_GO_clean, 
                    ax=ax
                )
    plt.annotate('Brain Development', xy=(0.596394, 0.31759), xytext=(0.596394 + 0.05, 0.31759), arrowprops=dict(facecolor='black', shrink=0.05))
    plt.annotate('Osteoblast Differentiation', xy=(0.592747, 0.287680), xytext=(0.592747 + 0.05, 0.287680), arrowprops=dict(facecolor='black', shrink=0.05))
    plt.annotate('Actin Cytoskeleton Organization', xy=(0.575973, 0.381161), xytext=(0.575973 + 0.05, 0.381161), arrowprops=dict(facecolor='black', shrink=0.05))
    legend = ax.legend(title="",  loc='upper left', bbox_to_anchor=(0, 1.1), fontsize='small')
    ax.grid(False)
    sns.despine()  
    ax.set_ylim([0, 0.5])
    path_p1 = "./GO_Top100_Overlap.png"
    plt.savefig(path_p1, dpi=600, bbox_inches='tight')  # Save as PNG with high resolution
    plt.show()