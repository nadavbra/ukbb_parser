What is ukbb_parser?
===============

ukbb_parser is a Python module for loading phenotypic and genetic data from the `UK Biobank <https://www.ukbiobank.ac.uk/>`_ (UKBB) through either programmatic (Python) or command-line interface. 

Key features:

* Loading any set of UKBB fields into a clean dataframe or CSV file.

* Implementing standard filtrations over the loaded dataset of samples (e.g. restricting ethnicity or avoiding kinship).

* Automatically extracting covariates commonly used in GWAS (sex, age, genetic principal components, batch and assessment canter).

* Extracting and parsing phenotypes derived from ICD-10 codes.

* Loading genetic data in both PLINK BED format (for raw markers) or BGEN format (for imputed variants).

Once configured and provided with the paths of the UKBB files on your filesystem, ukbb_parser will automatically access the raw data files required for each operation.


Usage
=====

Python API
-----------

A quick example:

.. code-block:: python
    
    from ukbb_parser import create_dataset
    PHENOTYPES = [
        ('height', 50, 'continuous'),
        ('diastolic_blood_preasure', 4079, 'continuous'),
        ('red_blood_cell_count', 30010, 'continuous'),
    ]
    eid, fields, covariates = create_dataset(PHENOTYPES)

After executing the above code, you will get a dataframe (``fields``) of the loaded samples (one per row) with three columns, corresponding to the three requested fields: ``height`` (defined by `UKBB field #50 <http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=50>`_), ``diastolic_blood_preasure`` (defined by `UKBB field #4079 <http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=4079>`_) and ``red_blood_cell_count`` (defined by `UKBB field #30010 <http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=30010>`_).

A more elaborate example is available in the `ukbb_parser demo notebook <https://github.com/nadavbra/ukbb_parser/blob/master/ukbb_parser%20demo.ipynb>`_. For more details on all the different options of the ``create_dataset`` function (and their default - but no necessarily desired - behavior), see the help message of the function.
    
