# CaMelia: an accurate imputation method for single-cell methylomes

**CaMelia (CAtboost-based method for predicting MEthyLatIon stAtes), a computational imputation method based on the CatBoost gradient boosting model for predicting single-cell methylation states to address the sparsity of the data.** CaMelia imputed the missing methylation states in a cell by borrowing information of the same locus in other cells, based on the locally paired similarity of methylation patterns around the target locus between cells or between cell and bulk data. 
**The trained CaMelia model can be used for different downstream analyses, including to impute low-coverage methylation profiles for sets of cells or aid the identification of cell types or cell sub-populations.**

![image](https://github.com/JxTang-bioinformatics/CaMelia/blob/master/images/forgithub-01.png)

**The overview of CaMelia model.**

# Before starting

The code of CaMelia has been tested in Windows with Python 2.7 (and 3.6).

You will need to download the source code and locate them in the same directory.

For training models for your own data, you need to download necessary files as fallowing:

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

```
chrom  location  cell1  cell2  ...  celln bulk
chr1   136409    1.0    NA     ...  0.0   1.0
chr1   136423    NA     1.0    ...  NA    0.0
chr1   136425    0.0    NA     ...  1.0   1.0
       ...                     ...
chrY   28748270  NA     1.0    ...  1.0   1.0
chrY   28748361  0.0    NA     ...  0.0   0.0
chrY   28773349  NA     NA     ...  0.0   0.0
```

2) Run ``XXX.py`` and ``XXX.py`` to create the input data for CaMelia:

3) Run ``XXX.py`` to train CaMelia and evaluate model performances:

4) Use ``XXX.py`` to impute methylation profiles:


# Examples:



# Trained models:

GSE65364:

HCCs:

HepG2:

mESC:

GSE56879:

2i:

Serum:

MII:

# Contact

* Jianxiong Tang
* jxtang@std.uestc.edu.cn 



