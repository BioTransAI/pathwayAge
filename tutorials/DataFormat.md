# Data input and output

Recognizing data is the first and essential step for the PathwayAge model, 
as well as for any other models. Here are some criteria for pre-processing 
the methylation data and show you the structure of the input, internal and 
output datasets. Both input and output datasets are saved as CSV files.

* [Data Input](#input_data)
  * [MethylData format](#MethylData_format)
    * [Handling Missing Values in Methylation Sites](#MissValue)
    * [Methylation Site Alignment](#Alignment)
  * [confunders fornmat](#confunder)

* [Data output](#dataOutput)
  * [Prediction format](#Prediction)
  * [data4stage2 format](#data4stage2)
## <a name="input_data"></a>Input data format

To run the 'pathwayAge' function, methylation data should be 
provided as the input data. For later analysis, the confunders
info(covariance data) is also required.

### <a name="MethylData_format"></a>MethylData format

For the input methylation matrix, each row contains information 
about a methylation site, and each column contains the information 
of all methylation sites for that sample. It is important to note 
that there is a crucial row named 'Age' attached to the matrix above.

#### <a name="MissValue"></a>Handling Missing Values in Methylation Sites

If the missing rate of a CpG value in the samples is less than 5%, 
then fill the NaN with the median. Otherwise, remove that CpG row.

#### <a name="Alignment"></a>Methylation Site Alignment

After processing the missing values, the training and testing datasets
should retain the same CpGs. Simply obtain the overlapping CpGs from 
all datasets.

<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>CpG</th>
      <th>GSM2333901</th>
      <th>GSM2333902</th>
      <th>GSM2333903</th>
      <th>GSM2333904</th>
      <th>...</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Age</td>
      <th>72.000</th>
      <th>55.000</th>
      <th>23.000</th>
      <th>86.000</th>
      <th>...</th>
    </tr>
    <tr>
      <td>cg05352250</td>
      <th>0.961</th>
      <th>0.951</th>
      <th>0.967</th>
      <th>0.967</th>
      <th>...</th>
    </tr>
    <tr>
      <td>cg16882684</td>
      <th>0.722</th>
      <th>0.598</th>
      <th>0.680</th>
      <th>0.653</th>
      <th>...</th>
    </tr>
    <tr>
      <td>...</td>
      <th>...</th>
      <th>...</th>
      <th>...</th>
      <th>...</th>
      <th>...</th>
    </tr>
  </tbody>
</table>
</div>


### <a name="confunder"></a>confunders fornmat

Covariate data is a user-designed matrix that includes confounding variables 
related to a specific disease or variables of interest. 

<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>Sample</th>
      <th>Age</th>
      <th>Label</th>
      <th>Female</th>
      <th>...</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>GSM2333901</td>
      <th>43</th>
      <th>0</th>
      <th>0</th>
      <th>...</th>
    </tr>
    <tr>
      <td>GSM2333902</td>
      <th>55</th>
      <th>0</th>
      <th>1</th>
      <th>...</th>
    </tr>
    <tr>
      <td>GSM2333903</td>
      <th>26</th>
      <th>0</th>
      <th>0</th>
      <th>...</th>
    </tr>
    <tr>
      <td>...</td>
      <th>...</th>
      <th>...</th>
      <th>...</th>
      <th>...</th>
    </tr>
  </tbody>
</table>
</div>

## <a name="dataOutput"></a>Output data format

After processing the biological age predictor, two files are generated: 
one containing the final age prediction results and the other containing 
internal datasets called "data4Stage2". 

### <a name="Prediction"></a>Prediction format
Prediction contains the biologically predicted age computed by 'pathwayAge'

<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>Sample</th>
      <th>Prediction</th>
      <th>Age</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>GSM2333901</td>
      <th>56.469699</th>	
      <th>55.0</th>
    </tr>
    <tr>
      <td>GSM2333903</td>
      <th>56.469699</th>
      <th>55.0</th>
    </tr>
    <tr>
      <td>GSM2333903</td>
      <th>22.008102</th>		
      <th>23.0</th>	
    </tr>
    <tr>
      <td>...</td>
      <th>...</th>
      <th>...</th>
    </tr>
  </tbody>
</table>
</div>

### <a name="data4Stage2"></a>data4Stage2 fornmat
'data4Stage2' is the outcome of condensing CpGs into pathways

<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th> </th>
      <th>GO:0000002</th>
      <th>GO:0000012</th>
      <th>GO:0000027</th>
      <th>...</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>7786915049_R03C02</td>
      <th>1.306</th>	
      <th>1.166</th>
      <th>1.069</th>
      <th>...</th>
    </tr>
    <tr>
      <td>7471147041_R05C01</td>
      <th>1.011</th>
      <th>1.767</th>
      <th>1.233</th>
      <th>...</th>
    </tr>
    <tr>
      <td>7507867089_R02C02</td>				
      <th>1.575</th>		
      <th>1.172</th>
      <th>1.511</th>
      <th>...</th>
    </tr>
    <tr>
      <td>...</td>
      <th>...</th>
      <th>...</th>
      <th>...</th>
    </tr>
  </tbody>
</table>
</div>

