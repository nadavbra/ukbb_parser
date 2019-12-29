What is ukbb_parser?
===============

ukbb_parser is a Python module for loading phenotypic and genetic data from the `UK Biobank <https://www.ukbiobank.ac.uk/>`_ (UKBB) through either programmatic (Python) or command-line interface. 

Key features:

* Loading any set of UKBB fields into a clean dataframe or CSV file.

* Implementing standard filtrations over the loaded dataset of samples (e.g. restricting ethnicity or avoiding kinship).

* Automatically extracting covariates commonly used in GWAS (sex, age, genetic principal components, batch and assessment canter).

* Extracting and parsing phenotypes derived from ICD-10 codes.

* Loading of genetic data, in both PLINK BED format (for raw markers) or the BGEN format (for imputed variants).

Once configured and provided with the paths of the required UKBB files on your filesystem, ukbb_parser will automatically access the raw data files required for each operation.
