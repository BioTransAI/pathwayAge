PathwayAge
==============

PathwayAge is a two-stage AI module designed for the prediction of biological age and the analysis of age-related pathways.

Installation
------------------- 


1. clone the repo:
 
```python{install}

git clone git@github.com:BioTransAI/pathwayAge.git

```
2. environments:
  1. Make sure the  [Anaconda](https://www.anaconda.com/products/individual)  or [Miniconda](https://docs.conda.io/projects/miniconda/en/latest/) is aready installed.  <br>
  2. Move the file "pyBioTrans.yaml" under the path ~/miniconda3/bin/.
  ```r{install}
    mv pyBioTrans.yml ~/miniconda3/bin/
    cd ~/miniconda3/bin/
    conda env create -f pyBioTrans.yml

        install.packages("BiocManager")

    BiocManager::install("methylclock")
    ```

 
```python{env}

git clone git@github.com:BioTransAI/pathwayAge.git

```

`ddBD(tr, outgroup, root.time = 1, measure = c("SSE", "KL"))`

  `tr` is an object of class "phylo" specifying the relative times from RelTime. 
	
  `outgroup` is a vector of character specifying all outgroup tips
  
  `root.time` is the age of root in the current time unit. The default value is 1. 
  
   `measure` is the method for selecting the initial values in grid search. The best initial values can be selected by minimzing the sum of squared errors (SSE) or Kullback-Leibler divergence (KL). The  default is method is SSE.
 	
Users need to provide a relative timetree inferred by RelTime without any calibrations in order to run `ddBD`. To get the relative timetree, one can use MEGA X (https://www.megasoftware.net).  

This program will produce the parameters for birth-death model (i.e., birth rate, death rate, and sampling fraction). All three parameters are automatically estimated simultaneously. The program that only estimates birth rate and death rate with a user-specified sampling fraction rate will come soon. Currently, the estimated parameters can only be directly used in MCMCTree for dating analysis.

Note that the program currently works well with R version 3.6.x. R version 4.0.x gives an error. I will fix it soon.


Directory Structure
------------------- 

"code" directory contains `ddBD` R function.

"simulation" directory contains all simulated data, estimated birth-death parameter values and MCMCTree results obtained using the true (BD_10_5_0.99), default flat (BD_2_2_0.1), and inferred (BD_inf) tree prior. 

"empirical-data" directory contains the empirical data, estimated birth-death parameter values and MCMCTree results. 


If you have more questions, please email qiqing.tao@temple.edu.


Citation
------------------- 

If you use ddBD, please cite:
Tao Q, Barba-Montoya J, and Kumar S. 2021. Data-driven speciation tree prior for better species divergence times in calibration-poor molecular phylogenies. bioRxiv. doi:10.1101/2021.03.27.437326. (In press by Bioinformatics)

If you use RelTime from MEGA X, please also cite:
Kumar S, Stecher G, Li M, Knyaz C, and Tamura K. 2018. MEGA X: Molecular Evolutionary Genetics Analysis across Computing Platforms. Mol. Biol. Evol. 35:1547-1549.
 
