PathwayAge
==============

PathwayAge is a two-stage AI module designed for predicting biological age and analyzing age-related pathways.

 ![](Image/pathwayAgeProcess.png)

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

1. PathwayAge Biological age presdiction   

  - Provide both the training dataset and the testing dataset:

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

  - Using only the training dataset, the model will automatically perform cross-validation-style predictions in each outer cross-validation loop to prevent  data leakage into the testing data labels (Age).

  ![](Image/crossValidation.png)

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

2. The co-founder adjusted the pathway feature importance ranking

    The python script will ranking the adjusted Age-related biological pathways and display the results in by R srcip.

  ```bash 
   
   # Execute the Python script to obtain feature ranking
   python3 XcofounderAdjustedPathwayFeatureImportanceRanking.py

   # Run the R script to visualize the ranking
   R XcofounderAdjustedPathwayFeatureImportanceRankingPlot.R

  ```

3. Identification of Age-Associated Modules and Network Analysis

    First, identify the Age-Associated Modules using WGCNA with prior knowledge of pathways.
    Second, analyze the relationship between the modules, between the modules and aging, then visualize the connections using a user-friendly D3 graph.


  ```bash 
    
    # Execute the R script to discover the modules.
    R YAgeAssociatedModuleIdentification.R

    # Run the Python script to analyze the modules and generate visual results.
    python3 YAgeAssociatedModuleNetworkAnalysis.py

  ```


Citation
------------------- 

 
