# Data input and output

Recognizing data is the first and essential step for the PathwayAge model, 
as well as for any other models. Here are some criteria for pre-processing 
the methylation data and show you the structure of the input, internal and 
output datasets. Both input and output datasets are saved as CSV files.

* [Data Input](#input_data)
  * [MethylData format](#MethylData_format)
    * [Handling Missing Values in Methylation Sites](#MissValue)
    * [Methylation Site Alignment](#Alignment)
  * [cofunders fornmat](#cofunder)

* [Data output](#dataOutput)
  * [Predicton format](#MethylData_format)
  * [cofunders fornmat](#separate_mat)

## <a name="input_data"></a>Input data format

To run the 'pathwayAge' function, methylation data should be 
provided as the input data. For later analysis, the cofunders 
(covariance data) should be used as input data.

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


### <a name="cofunder"></a>cofunders fornmat

Covariate data is a user-designed matrix that includes confounding variables 
related to a specific disease or variables of interest. The following tabel
displays the components I used in this project which 

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

## <a name="dataOutput"></a>Output data format

After processing the biological age predictor, two files are generated: 
one containing the final age prediction results and the other containing 
internal datasets called "data4Stage2".

### <a name="sample_meta"></a>Sample metadata

The sample metadata is a table which contains additional information about each sample, such as timepoint or genotype.
Each row should represent a sample and each column should represent a metadata feature, where the first columns contains the same sample identifier that was used in the gene expression matrix
The rows should be in the same order as the rows of the gene expression matrix, or
the user can specify `order=False`.

<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>Sample_id</th>
      <th>Age</th>
      <th>Tissue</th>
      <th>Sex</th>
      <th>Genotype</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>sample_11615</td>
      <td>4mon</td>
      <td>Cortex</td>
      <td>Female</td>
      <td>5xFADHEMI</td>
    </tr>
    <tr>
      <td>sample_11616</td>
      <td>4mon</td>
      <td>Cortex</td>
      <td>Female</td>
      <td>5xFADWT</td>
    </tr>
  </tbody>
</table>
</div>

### <a name="params"></a>Other parameters
These are other parameters that can be specified.

* **name**: Name of the WGCNA used to visualize data (default: `WGCNA`)

* **save**: Whether to save the results of important steps or not (If you want to set it
`True` you should have a write access on the output directory)

* **outputPath**: Where to save your data, otherwise it will be stored in the same directory as the code.

* **TPMcutoff**: TPM cutoff for removing genes

* **networkType** : Type of network to generate ({`unsigned`, `signed` and `signed hybrid`}, default: `signed hybrid`)

* **adjacencyType**: Type of adjacency matrix to use ({`unsigned`, `signed` and `signed hybrid`}, default: `signed hybrid`)

* **TOMType**: Type of topological overlap matrix(TOM) to use ({`unsigned`, `signed`}, default: `signed`)

For depth-in documentation on these parameters see [here](https://mortazavilab.github.io/PyWGCNA/html/PyWGCNA.html).

## <a name="data_preproc"></a>Data cleaning and preprocessing

PyWGCNA can clean the input data according to the following criteria:
1. Remove genes without any expression more than `TPMcutoff` value (default one) across all samples.
2. Find genes and samples `goodSamplesGenes()` function to find genes and samples with too many missing values.
3. Cluster the samples (uses [hierarchical clustering](https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html#module-scipy.cluster.hierarchy)
from [scipy](https://scipy.org/)) to see if there are any obvious outliers. The user can define value the height by specifying the `cut` value. By default, no samples are removed by hierarchical clustering
