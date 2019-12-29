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

``create_ICD10_dataset`` is a similar function that would provide you with ICD-10 derived phenotypes, on top of the "regular" UKBB fields requested. ICD-10 codes are derived from the following UKBB fields: `main diagnoses (41202) <http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=41202>`_, `secondary diagnoses (41204) <http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=41204>`_, `cancer type (40006) <http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=40006>`_, `primary cause of death (40001) <http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=40001>`_, `secondary cause of death (40002) <http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=40002>`_ and `external causes (41201) <http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=41201>`_. Usage of this function is also exemplified in the same notebook.

The  `ukbb_parser demo notebook <https://github.com/nadavbra/ukbb_parser/blob/master/ukbb_parser%20demo.ipynb>`_ illustrates most of the core functionality of the Python API. In addition, you may consult with the help message of the module's functions.

Command-line API
-----------

If you are not a fan of Python, and just want to use ukbb_parser to extract some phenotypes from the UKBB and process them elsewhere (e.g. using another language such as R, or using a dedicated software such as PLINK), then you may prefer to work with the command-line API.

The main command provided by ukbb_parser is ``create_ukbb_phenotype_dataset``. A basic usage of this command goes as follows:

.. code-block:: cshell

    create_ukbb_phenotype_dataset --phenotype-specs-file=/path/to/phenotype_specs.py --output-dataset-file=/path/to/created_ukbb_dataset.csv
    
As an example, you can use `examples/phenotype_specs.py <https://github.com/nadavbra/ukbb_parser/blob/master/examples/phenotype_specs.py>`_ which specifies 49 prominent phenotypes from various sources:

.. code-block:: cshell

    cd /tmp
    wget https://raw.githubusercontent.com/nadavbra/ukbb_parser/master/examples/phenotype_specs.py
    create_ukbb_phenotype_dataset --phenotype-specs-file=./phenotype_specs.py --output-dataset-file=./ukbb_dataset.csv --nrows=10000

After a short execution, this should produce a CSV file (ukbb_dataset.csv) with 8,012 rows (of the filtered samples) and 226 columns (that include the specified phenotypes, and additional covariates and indices). We limit the analysis to 10,000 samples (using ``--nrows=10000``) only because we are impatient to see that it actually works. Of course you will want to remove this option when extracting the full dataset to be used in a research.

The main thing to understand about this command is how to specify the requested phenotypes. If you open `examples/phenotype_specs.py <https://github.com/nadavbra/ukbb_parser/blob/master/examples/phenotype_specs.py>`_ you will see that it actually follows a rather straightforward format, namely a list dictionaries, each representing a distinct phenotype to be parsed:

.. code-block:: python

    specs = [
        {
        # Specification of phenotype 1...
        },
        {
        # Specification of phenotype 2...
        },
        # ...
    ]
    
Each phenotype specification requires two settings:

* **name**: The name of the phenotype (to comprise the name of the column in the output CSV)
* **source**: ``'field'`` if taken directly from a UKBB field, or ``'ICD-10'`` if taken from a set of ICD-10 codes, or ``'aggregation'`` if defined by some aggregation function of other sub-specifications (more details on each of these three options below).

*field* specifications should also provide the following two settings:

* **field_id**: The UKBB field ID (as listed in the `UKBB's showcase Crystal system <http://biobank.ctsu.ox.ac.uk/crystal/>`_).
* **field_type**: The type of the field, which could be either ``'continuous'`` (to be parsed as real-valued continuous numbers, which are to be the maximum of all the samples associated with each sample, in case of multiple versions/entries of the field), or ``'binary'`` (to be parsed as 0 or 1), or ``'set'`` (to provide all the values associated with each sample), or ``'function'`` (to provide any Python function to parse the raw values of the field). More details about these options is available in the documentation ukbb_parser.create_dataset. 

*ICD-10* specifications should also provide a **codings** setting, that should contain a list of ICD-10 codings to be considered as part of this phenotype (this will result a binary phenotype with positive (1) values corresponding to each sample that has at least one of the listed ICD-10 codes).

Finally, *aggregation* specifications should provide a **subspecs** setting that lists the sub-specifications to be aggregated, and an **aggregation_function** setting that provided a Python function to perform the aggregation. Each of the sub-specifications is defined exactly like a root specification (i.e. it can be either a *field* or an *ICD-10* or even another *aggregation* specification). The aggregation function should receive the series/dataframe values created for each of the sub-specifications (each as a separate argument) and return the resulting series/dataframe for this field.

Each specification may also provide a **sex_filter** setting, in case that the phenotype is only relevant to Males (``'M'``) or females (``'F'``).

If you want to extract a dataset to work on using the ``create_ukbb_phenotype_dataset`` command, it is highly advisable that you spend a few minutes going throguh its help message (simply run ``create_ukbb_phenotype_dataset --help``) so that you won't be surprised by any unintended consequences of its default behavior and so that you will be able to make the most out of it.

ukbb_parser also provides a ``create_ukbb_genotype_spec_file`` command, that lists all the genotyping files provided by the UKBB. For example:

.. code-block:: cshell

    create_ukbb_genotype_spec_file --genotyping-type=raw --output-file=/tmp/ukbb_raw_marker_genotyping_spec.csv
    
Or:

.. code-block:: cshell

    create_ukbb_genotype_spec_file --genotyping-type=imputation --output-file=/tmp/ukbb_imputation_genotyping_spec.csv
    

Installation
===============

Just run:

.. code-block:: cshell

    git clone https://github.com/nadavbra/ukbb_parser.git /tmp/ukbb_parser_src
    cd /tmp/ukbb_parser_src
    git submodule update --init --recursive
    python setup.py install


Python dependencies
-----------------

(requires Python 3)

* numpy
* scipy
* pandas
* matplotlib
* biopython
* statsmodels


Configuration
===============

Following installtion, you will have to create and configure a ``.ukbb_paths`` file in your home directory (i.e. ``~/.ukbb_paths``). This settings file specifies the paths of all the necessary UKBB data files on your filesystem.







    

    
