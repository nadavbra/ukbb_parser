from setuptools import setup
        
def readme():
    with open('README.rst', 'r') as f:
        return f.read()

setup(
    name = 'ukbb_parser',
    version = '1.0.2',
    description = 'A Python module for loading phenotypic and genetic data from the UK Biobank.',
    long_description = readme(),
    url = 'https://github.com/nadavbra/ukbb_parser',
    author = 'Nadav Brandes',
    author_email  ='nadav.brandes@mail.huji.ac.il',
    license = 'MIT',
    packages = ['ukbb_parser', 'ukbb_parser.shared_utils'],
    scripts = [
        'bin/create_ukbb_phenotype_dataset',
        'bin/create_ukbb_genotype_spec_file',
    ],
    install_requires = [
        'numpy',
        'scipy',
        'pandas',
        'matplotlib',
        'biopython',
        'statsmodels',
    ],
)
