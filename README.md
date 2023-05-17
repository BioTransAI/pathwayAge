# pathwayAge

PathwayAge is a 2-stage biological age predictor that first builds separate machine learning models for methylation sites mapping to each of the BP pathways, yielding 1 machine-learning model per pathway (first-stage). This procedure compresses data from individual methylation sites into a pathway-level feature. Then, a second-stage algorithm integrates these pathway-level features into a systems-level regression.
 
## Installation:

1. Clone the repo:
```bash
git clone https://github.com/BioTransAI/pathwayAge.git
```
2. Environment:
make sure [Anaconda](https://docs.anaconda.com/free/anaconda/install/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) is already installed. Then move the file "pyBioTrans.yml" under the ~/miniconda3/bin/. 
```bash
mv pyBioTrans.yml ~/miniconda3/bin/
conda env create -f pyBioTrans.yml
```

## Installation:
Example usage:
```python
import pathwayAge from pathwayAge

pathwayAge(
    methylData = methylDataAge,
    methylTestData = methylTestData,
    resultName = "bioAgePrediction",
)

```
