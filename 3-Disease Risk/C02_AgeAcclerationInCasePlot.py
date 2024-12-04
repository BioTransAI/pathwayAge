import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from C01_AgeAcclerationInCase import testAgeAccInCase

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

A = testAgeAccInCase()
A.rename( columns =  {
    "DNAmAge": "Horvath", 
    "DNAmAgeHannum": "Hannum", 
    "DNAmPhenoAge": "PhenoAge", 
    "DNAmGrimAgeBasedOnRealAge": "GrimAge",
    "Epigenetic Age (Zhang)": "Zhang",
}, inplace=True)

order = list(disease_Map.values())
A["Index"] = A.index

CaseMean_A = A[A.Index.str.contains("Differ")]
CaseMean_A['Disease'] = pd.Categorical(CaseMean_A['Disease'], categories=order, ordered=True)
CaseMean_A = CaseMean_A.sort_values('Disease')
CaseMean_A.index = CaseMean_A.Disease
CaseMean_A = CaseMean_A.drop(columns = ["Data", "Disease", "Index"])

def Color(value):
    if value < 0:
        return 0
    elif  value == 0:
        return 1
    elif value > 0 :
        return 2

CaseMean_A = CaseMean_A.applymap(Color)
CaseMean_A_T = CaseMean_A.T

Pvalue_A = A[A.Index.str.contains("P")]
Pvalue_A['Disease'] = pd.Categorical(Pvalue_A['Disease'], categories=order, ordered=True)
Pvalue_A = Pvalue_A.sort_values('Disease')
Pvalue_A.index = Pvalue_A.Disease
Pvalue_A = Pvalue_A.drop(columns = ["Data", "Disease", "Index"])
Pvalue_A_T = Pvalue_A.T

def significance_stars(p):
    if p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return ''
        
Pvalue_A_T = Pvalue_A_T.applymap(significance_stars)
Pvalue_A_T
f, ax = plt.subplots(figsize=(12, 7))

# Generate a custom diverging colormap
cmap = sns.color_palette(["#6496ae", "gray", "#c4593c"])

# Draw the heatmap with the mask and correct aspect ratio
sns.heatmap(
    CaseMean_A_T, 
    cmap=cmap, 
    # annot=Pvalue_sorted, 
    cbar=False, 
    square=True, 
    linewidths=1, 
    cbar_kws={"shrink": .15}
)

# Overlay the significance stars
for i in range(Pvalue_A_T.shape[0]):
    for j in range(Pvalue_A_T.shape[1]):
        ax.text(j + 0.5, i + 0.5, Pvalue_A_T.iloc[i, j],
                ha='center', va='center', color='black', fontsize=12)
        
plt.title("")
plt.xticks(rotation=0) 
plt.savefig('./Result/DiseaseRisk.png', dpi=300, bbox_inches='tight')  # Save as PNG with high resolution

# Optionally, you can show the plot (if needed)
plt.show()
