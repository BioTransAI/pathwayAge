import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from C03_AgeAcclerationInCase import testAgeAccInCase
from functools import reduce


def AgeAcclerationInCaseSexStratificationPlot():
    Pvalue_sortedT, deltaAgeList = testAgeAccInCase()
    deltaAgeAll = reduce(lambda df1,df2: pd.concat([df1, df2], axis = 0),  deltaAgeList)
    deltaAgeAll["Disease"] = deltaAgeAll.Cohort.replace(disease_Map)
    deltaAgeAll["Method_temp"]= np.where(deltaAgeAll['Method'].str.contains("GO", case=False), "GO", "KEGG")
    deltaAgeAll.sort_values(by = ["Disease", "Method_temp"])

    Female = deltaAgeAll[deltaAgeAll.Method.isin(["PathwayAge_Female(KEGG)", "PathwayAge_Female(GO)"])]
    FemaleSCZ_ND = Female[Female.Disease.isin(['SCZ1', 'SCZ2', 'SCZ3', 'SCZ4', 'AD', 'PD1', 'PD2', 'PD3'])]
    FemaleSCZ_ND["Disease"] = FemaleSCZ_ND.Disease + "(" +FemaleSCZ_ND.Method_temp+ ")"

    FemaleSCZ_ND['Disease'] = pd.Categorical(FemaleSCZ_ND['Disease'], categories=SCZ_ND, ordered=True)
    FemaleSCZ_ND = FemaleSCZ_ND.sort_values(by='Disease')


    plt.figure(figsize=(20, 6))
    plt.ylim(-10, 10)
    sns.boxplot(x='Disease', y='Full', hue='Tag', data=FemaleSCZ_ND, showfliers=False, palette=['#616A6B', '#CB4335'], width=0.5)

    # Grouped boxplot
    Pvalue_sortedT_Female = Pvalue_sortedT[['SCZ1', 'SCZ2', 'SCZ3', 'SCZ4', 'AD', 'PD1', 'PD2', 'PD3']].T[["PathwayAge_Female(GO)", "PathwayAge_Female(KEGG)"]]
    Pvalue_sortedT_Female = Pvalue_sortedT_Female.to_numpy()
    flattened = Pvalue_sortedT_Female.flatten()
    flattened 
    # Annotate the P-values on the plot and add an underline
    for i, p_val in enumerate(flattened):
        x_coord = i
        y_coord = 1.05 * 10  # Y position of the P-value text
        plt.text(x_coord, y_coord, p_val, ha='center', fontsize=20)
        # Add an underline below the P-value text
        plt.plot([x_coord - 0.1, x_coord + 0.1], [y_coord - 0.02, y_coord - 0.02], color='black', lw=1.5)

    plt.ylabel("AgeAcc")
    plt.xlabel("")
    ax = plt.gca()  # Get current axes
    for spine in ax.spines.values():
        spine.set_visible(False)
    plt.legend(loc='center left', bbox_to_anchor=(0.97, 0.5))

    plt.savefig("./Restult/Disease_Female_SCZ_ND.png", dpi=600, bbox_inches='tight')  # Save as PNG with high resolution
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

SCZ_ND =[
    'SCZ1(GO)', 
    'SCZ1(KEGG)', 
    'SCZ2(GO)',
    'SCZ2(KEGG)', 
    'SCZ3(GO)',
    'SCZ3(KEGG)', 
    'SCZ4(GO)', 
    'SCZ4(KEGG)', 
    'AD(GO)',
    'AD(KEGG)', 
    'PD1(GO)',
    'PD1(KEGG)',
    'PD2(GO)', 
    'PD2(KEGG)',
    'PD3(GO)', 
    'PD3(KEGG)',  
]
