# CaMelia: an accurate imputation method for single-cell methylomes

**CaMelia (CAtboost-based method for predicting MEthyLatIon stAtes), a computational imputation method based on the CatBoost gradient boosting model for predicting single-cell methylation states to address the sparsity of the data. CaMelia imputed the missing methylation states in a cell by borrowing information of the same locus in other cells, based on the locally paired similarity of methylation patterns around the target locus between cells or between cell and bulk data. The trained CaMelia model can be used for different downstream analyses, including to impute low-coverage methylation profiles for sets of cells or aid the identification of cell types or cell sub-populations.

The overview of CaMelia model:
![image](https://github.com/JxTang-bioinformatics/PretiMeth/raw/master/images/Diagram_of_PretiMeth.png)

