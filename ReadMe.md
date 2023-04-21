# pepDESC

pepDESC is a useful tool for differential expression analysis of single-cell MS-based proteomics data. 

## Features

* pepDESC could analysis label-free quantification MS-based proteomics data generated by various MS machine and search engines. It was designed for the single-cell proteomics data but it also has good performance for low-input and regular-size proteomics data.
* pepDESC tolerates data with large fraction of missing values without data imputation.
* pepDESC filtered out contaminants peptides. It also has the ability to identify falsely-match contaminant peptides and removed them if the average retention time and m/z of peptides were provided.
* pepDESC deciding the differentially expressed proteins by scoring each protein.

## Installation

To install the pepDESC package, start R (version 4.0) and enter:

```R
devtools::install_github('pepDESC')
```

## Prepare your data

pepDESC reads the peptide-level quantification result generated by various search engine, as long as transformed into a readable data frame. Each row of the data is a peptide. The column of the data should incorporate the following  keywords.

* **Accession**: the UNIPROT protein accession of the peptide. If you would like to remove contamination peptides, you need to include a certain contamination indicator in this column.
* **Description**: a brief description of the protein, better include the name of the protein. If  you would like to remove the Keratin peptides or the Trypsin peptides from the data, you need to include the keyword "Keratin" or "Trypsin" respectively in this column.
* **Sequence**: the sequence of each peptide
* **Modifications**: the modifications, if exist, of each peptide
* **Abundance**: the abundance of each peptides in different cells. 

The Initiation function will report an error if any of the mentioned keywords were missing in the column names

* **RT**: the average retention time of the peptide
* **m/z**: the average m/z of the peptide

This two information could be omitted if the contamination match is not desired.

An example of the prepared data frame would be like the *test_raw.rda* in the *data* folder.

> pepDESC could not handle shared peptides. Please remove the shared peptides in advance if they were included in the search engine result.

## Run your analysis

To run a general pepDESC analysis of the test data, start R (version 4.0) and enter:

```R
load(pepDESC)

# This function make your data frame ready for analysis
Peptide<-Initiation(test_raw)

# This function filtered the peptides not apllied for downstream analysis
Peptide<-FilterData(Peptide,type=test_type,'A','B')

# This function removed outlier samples in your data
Peptide<-FilterSample(Peptide)

# This function applied normalization to samples
Peptide<-Normalization(Peptide)

# This function calculate the pairwise ratio of peptides between two groups
Peptide<-PeptideRatio(Peptide)

# This function calculate the final DE-score of each protein
Result<-ScoreProtein(Peptide)

# Finally, you find the differentially expressed proteins by a certain threshold of DE-score
Up<-rownames(Result)[Result$score>0.3]
Down<-rownames(Result)[Result$score<(-0.3)]
```

## Build your own pepDESC workflow

Several parameters could be adjusted in the pepDESC to make it applicable to your data.

### FilterData

* The **FilterDat**a function remove the keratin peptides and the trypsin peptides by default. If that is not desired, change the rmKRT and the **rmTRY** to FALSE. 

  ```R
  Peptide<-FilterData(Peptide,type=test_type,'A','B',rmKRT=FALSE,rmTRY=FALSE)
  ```

* Change the accession keywords of the contaminants protein in the **FilterData** function. For example, if your data includes reverse peptides, which are marked with protein accession starting with "REV", include it in the **Contamination** argument.

  ```R
  Peptide<-FilterData(Peptide,type=test_type,'A','B',Contamination='CON|REV')
  ```

* Change the maximum permitted missing value fraction by the **NAThres** argument in the **FilterData** function considering your data size. By default, the **FilterData** function allows 60% of missing values to the greatest extent. 

### Normalization

* This function normalize the samples by the median values. However, when a sample contains a large fraction of missing values, where the median value would be 0, this function would normalize the samples by the average abundances of the quantified peptides. You could also select the average-value normalization by setting the **method** argument in the **Normalization** function to "Mean"

  ```R
  Peptide<-Normalization(Peptide,method='Mean')
  ```

* If the difference of the size of two batches of samples is nonnegligible, change the **group** argument in the **Normalization** function to be TRUE:

  ```R
  Peptide<-Normalization(Peptide,group=TRUE)
  ```

### PeptideRatio

* The function take the first appeared sample type as the numerator of the calculation. You could change this by changing the **numerator** argument in this function.

  ```R
  Peptide<-PeptideRatio(Peptide,numerator='B')
  ```

### ScoreProtein

* Change the significance for deciding a differentially expressed peptide by the **pvaluecutoff** argument in the **ScoreProtein** function. 

  ```R
  Result<-ScoreProtein(Peptide,pvaluecutoff=0.99)
  ```

**Or freely re-combine the functions provided by pepDESC for your own analysis.**

## Contributing

Contribute this project by [opening an issue](https://github.com/dionezhang/pepDESC/issues).
