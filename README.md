PathwayAge
==============

PathwayAge is a two-stage AI model specially designed to predict biological age and analyze age-related pathways, offering both predictive accuracy and explainability.

 <!-- ![](Image/pathwayAgeProcess.png) -->

Installation
------------------- 


1. clone the repo:
 
```bash

git clone git@github.com:BioTransAI/pathwayAge.git

```
2. environments:
  - Make sure the  [Anaconda](https://www.anaconda.com/products/individual)  or [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/) is already installed.  <br>
  - Move the file "pyBioTrans.yaml" under the path ~/miniconda3/bin/. <br>

  ```bash
  
    mv pyBioTrans.yml ~/miniconda3/bin/
    cd ~/miniconda3/bin/
    conda env create -f pyBioTrans.yml
    conda activate pyBioTrans

  ```


Example Usage:
------------------- 

1. To predict biological age, please call the "PathwayAge" function.

  - Supplying both the training dataset and the testing dataset:

  ```python

    from pathwayAge import pathwayAge
    
    pathwayAge(
      # import Training data
      methylData = methylTrainData,
      # import testing data
      methylTestData = methylTestData,
      # name the file for the prediction results
      resultName = resultName,
      # Select one machine learning method
      predictionMode = predictionMode,
    )
  ```

  - Using only the training dataset, the model will automatically perform nested cross validation to prevent data leakage into the testing data labels (Age) in stage1.

  <!-- ![](Image/crossValidation.png) -->

  ```python

    from pathwayAge import pathwayAge
    
    pathwayAge(
      # import Training data
      methylData = methylTrainData,
      # name the file for the prediction results
      resultName = resultName,
      # select one machine learning method
      predictionMode = predictionMode,
    )
  ```

2. The confounder-adjusted pathway feature importance ranking.

    The Python script will rank the adjusted age-related biological pathways and present the results using an R script.

    ```bash 
    
    # Execute the Python script to obtain feature ranking
    python3 XcofounderAdjustedPathwayFeatureImportanceRanking.py

    # Run the R script to visualize the ranking
    R XcofounderAdjustedPathwayFeatureImportanceRankingPlot.R

    ```

3. Age-Associated Modules Identification and Network Analysis

    First, identify the Age-Associated Modules using WGCNA with prior knowledge of pathways.
    Second, analyze the relationship between the modules, between the modules and aging, then visualize the connections using a user-friendly D3 graph.


  ```bash 
    
    # Execute the R script to discover the modules.
    R YAgeAssociatedModuleIdentification.R

    # Run the Python script to analyze the modules and generate visual results.
    python3 YAgeAssociatedModuleNetworkAnalysis.py

  ```


Tutorial
------------------- 
- [Data input and output](tutorials/DataFormat.md): How to preprocess your input data for PathwayAge model and what output data format you will get.
- [Quick Start](tutorials/QuickStart.ipynb): What data and parameters need to be load into PathwayAge and predict biological Age.
- [confounder Adjusted Pathway Feature Importance Ranking](tutorials/FeatureRanking.ipynb): How to ranking the adjusted Age-related biological pathways.
- [Age Associated Module Identification](tutorials/ModuleIdentification.R): How to iidentify the Age-Associated Modules using WGCNA with prior knowledge of pathways.
- [Age Associated Module Network Analysis](tutorials/NetworkAnalysis.ipynb): How to analyze the relationship between the modules, between the modules and aging, then visualize the connections using a user-friendly D3 graph.

Citation
------------------- 