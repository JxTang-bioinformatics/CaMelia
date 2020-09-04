# CaMelia: an accurate imputation method for single-cell methylomes

**CaMelia (CAtboost-based method for predicting MEthyLatIon stAtes), a computational imputation method based on the CatBoost gradient boosting model for predicting single-cell methylation states to address the sparsity of the data.** CaMelia imputed the missing methylation states in a cell by borrowing information of the same locus in other cells, based on the locally paired similarity of methylation patterns around the target locus between cells or between cell and bulk data. **The trained CaMelia model can be used for different downstream analyses, including to impute low-coverage methylation profiles for sets of cells or aid the identification of cell types or cell sub-populations.**

The overview of CaMelia model:
![image](https://github.com/JxTang-bioinformatics/PretiMeth/raw/master/images/Diagram_of_PretiMeth.png)

# Before starting

The code of PretiMeth has been tested in Windows with Python 2.7 (and 3.6).

You will need to download the source code and locate them in the same directory.

1) You can download pre-trained models:
GSE65364:
HCCs:
HepG2:
mESC:
GSE56879:
2i:
Serum:
MII:

2) For training models for your own data, you need to download necessary files as fallowing:
   i)  Feature extraction:
   ii) Model training:
   

# Getting started
1) Store the raw data from cells of the same cell type or bulk data of this cell type into a file with the following columns:

* Chromosome (with chr)
* Position of the CpG site on the chromosome starting with one
* Binary methylation state of the CpG sites in cell 1(0=unmethylation, 1=methylated)
* ...
* Binary methylation state of the CpG sites in cell n
* Binary methylation state of the CpG sites in bulk data

Example:

.. code::

chrom  location  cell1  cell2  ...  celln bulk
chr1   136409    1.0    NA     ...  0.0   1.0
chr1   136423    NA     1.0    ...  NA    0.0
chr1   136425    0.0    NA     ...  1.0   1.0
       ...                     ...
chrY   28748270  NA     1.0    ...  1.0   1.0
chrY   28748361  0.0    NA     ...  0.0   0.0
chrY   28773349  NA     NA     ...  0.0   0.0









