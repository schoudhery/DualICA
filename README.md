# DualICA

Dual ICA script:
> python3 dual_ICA.py lfcs_data categories_data NCIC NGIC

    * lfcs data : genes x conditions and each value is the LFC of the condition compared to DMSO.
    * categories data : a metadata file with atelast two colums- Name and Classification, where name is the condition names in lfcs data and classification is the class of the condition name. 
        # NOTE: categories data file is used ONLY to label the output file for the user and not used in the calculation of the blocks. 
    * NCIC : number of compoenents to extract from conditions -> signal matrix of this has genes on the rows
    * NGIC : number of compoenents to extract from genes -> signal matrix of this has conditions on the rows


Linear Regression:
> python3 linear_reg.py G_matrix C_Matrix lfcs_data

    * G_matrix : The signal matrix of the ICA where genes are the rows
    * C_matrix : The signal matrix of the ICA where conditions  are the rows
    * lfcs data : genes x conditions and each value is the LFC of the condition compared to DMSO.
    * categories data : a metadata file with atelast two colums- Name and Classification, where name is the condition names in lfcs data and classification is the class of the condition name. 
        # NOTE: categories data file is used ONLY to label the output file for the user and not used in the calculation of the blocks. 


Outer Product
> python3 linear_reg.py G_matrix C_Matrix lfcs_data

    * G_matrix : The signal matrix of the ICA where genes are the rows
    * C_matrix : The signal matrix of the ICA where conditions  are the rows
    * lfcs data : genes x conditions and each value is the LFC of the condition compared to DMSO.
    * categories data : a metadata file with atelast two colums- Name and Classification, where name is the condition names in lfcs data and classification is the class of the condition name. 
        # NOTE: categories data file is used ONLY to label the output file for the user and not used in the calculation of the blocks. 

