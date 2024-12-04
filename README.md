PathwayAge
==============

PathwayAge is a two-stage AI model specially designed to predict biological age and analyze age-related pathways, offering both predictive accuracy and explainability.

 <!-- ![](Image/pathwayAgeProcess.png) -->

Repository Structure
------------------- 
- `1-PathwayAge Model/`: This directory contains the code for constructing the PathwayAge model and predicting biological age.
- `2-Aging Associated Pathways and Modules/`: This directory explains how to identify age-associated  modules using WGCNA or GO terms.
- `3-Disease Risk/`: This directory evaluates the differences in AgeAcc between patients and healthy individuals.
- `4-Disease Specific pathways/`: This directory analyzes differential pathways across different disease states, including sex-stratified differences within each disease. It also includes a permutation test for the top pathways.
- `5-Reproduction in Transcriptomics/`: This directory contains an analysis using transcriptomics to validate the reproducibility of the top methylation pathways.

Installation
------------------- 


1. clone the repo:
 
```bash

git clone git@github.com:BioTransAI/pathwayAge_NA.git

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

 To predict biological age, please call the "PathwayAge" function.

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
