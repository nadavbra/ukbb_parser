'''
The ukbb_parser Python module assumes the existence of a ~/.ukbb_paths.py settings file, that would specify the paths of all the necessary UKBB data
files on your filesystem.
See the detailed instruction above each of the required settings below in order to learn how to obtain the data files and configure them.
'''

import os

# The following two variables are not really part of the settings, and thus aren't necessary, but they can help you organize your settings.
BASE_DIR = '/path/to/ukbiobank/files/on/your/filesystem'
GENETICS_DIR = os.path.join(BASE_DIR, 'genetics')

# This is the primary CSV file with all the provided fields, according to your application. 
# Once your UKBB applicaiton is approved (through the UKBB AMS system), you should receive an email titled "UK Biobank Application [...], data now
# available for Basket [...]" (or something along those lines). In this email, you should receive an access key, and detailed instructions how to download
# your approved data. At the end of the process, you should have generated the file ukb<application_id>.csv, whose path should be set here.
PHENOTYPE_DATASET_CSV_FILE_PATH = os.path.join(BASE_DIR, 'ukb12345.csv')
PHENOTYPE_DATASET_CSV_FILE_ENCODING = 'latin-1'

# Every once in a while, the UKBB should send you an email titled "UK Biobank Application [...], Participant Withdrawal Notification [...]" (or something
# along those lines). In this email they attach a CSV file with the sample IDs (eid) of UKBB participants who have requested to withdraw from the dataset.
# As a UKBB user, you are obligated to abide to these requests and exclude those samples from your analysis.
# In order to manage these withdrawals, simply create a folder and put all the receieved CSV files (attached to these emails) within this directory. Set
# this setting to the path of that directory and ukbb_parser will automatically exclude the samples from any analysis.
PARTICIPANT_WITHDRAWALS_DIR = os.path.join(BASE_DIR, 'participant_withdrawals')

# FAM files are part of the PLINK BED format. In the context of UKBB, it is provided for the genotyping of ~800K raw markers in each sample using the Axiome
# SNP array. ukbb_parser requires the FAM file for two purposes: i) to extract additional covariates (e.g. the genetic batch of each sample), and ii) when
# loading raw-marker genetic data (using the pandas-plink module).
# Since the FAM file is the same for all chromosomes, it's sufficient to download this file only once. It can be downloaded through the following
# command line:
# ukbgene cal -m -c1
# where ukbgene can be downloaded from: http://biobank.ctsu.ox.ac.uk/crystal/download.cgi?id=665&ty=ut
FAM_FILE_PATH = os.path.join(GENETICS_DIR, 'ukb12345_cal_chr1_v2_s123456.fam')

# I you wish to access the ~800,000 raw markers genotyped for each of the UKBB samples (through the ukbb_parser.get_chrom_raw_marker_data function), then
# you will have to download the EGAD00010001226 directory. If you don't wish to access this genetic data, then you can comment-out this setting.
# The EGAD00010001226 directory can be downloaded using, for example, the EGA Download Client (https://ega-archive.org/download/using-ega-download-client)
# or through FTP / AsperaConnect clients (contact the UK Biobank Access Team if you need help downloading this data). 
# Within EGAD00010001226/001/, you should have (among other files, which are not handled by ukbb_parser) files of the form ukb_cal_chr#_v2.bed and 
# ukb_snp_chr#_v2.bim (one file of each of the two forms per chromosome).
CALL_DIR = os.getenv('CALL_DIR', os.path.join(GENETICS_DIR, 'EGAD00010001226/001'))

# If you wish the access the ~100M variants imputed, based on the ~800K called variants, for each of the UKBB samples (through the 
# ukbb_parser.get_chrom_imputation_data function), then you will have to download the EGAD00010001474 directory. If you don't wish to access this genetic
# data, then you can comment-out this setting.
# The EGAD00010001474 directory can be downloaded using, for example, the EGA Download Client (https://ega-archive.org/download/using-ega-download-client)
# or through FTP / AsperaConnect clients (contact the UK Biobank Access Team if you need help downloading this data).
# Within EGAD00010001474/, you should have files of the form ukb_imp_chr#_v3.bgen and ukb_imp_chr#_v3.bgen.bgi (one file of each of the two forms per 
# chromosome). On top of these files, you should also add link files of the form ukb#_imp_chr#_v3.sample (that are specific to your application). These
# .sample files (one per chromosome), can be downloaded using the command:
# ukbgene imp -m -c#
# (where # is the chromosome name, e.g. -c1 for chr1)
# ukbgene can be downloaded from: http://biobank.ctsu.ox.ac.uk/crystal/download.cgi?id=665&ty=ut
IMPUTATION_V3_DIR = os.getenv('IMPUTATION_V3_DIR', os.path.join(GENETICS_DIR, 'EGAD00010001474'))

# Since the name of the .sample files in <IMPUTATION_V3_DIR> (see above) are specific per application, you will have to provide their tamplate name here
# (i.e. just replace 12345 with your specific application ID).
IMPUTATION_V3_SAMPLE_FILE_NAME_TEMPLATE = 'ukb12345_imp_chr%s_v3.sample'

# If you wish to filter samples by kinship (which is the default behavior of ukbb_parser's functions), then you will have to download the kinship table
# by running the following command line:
# ukbgene rel
# where ukbgene can be downloaded from: http://biobank.ctsu.ox.ac.uk/crystal/download.cgi?id=665&ty=ut
# Afterwards, set this setting to the path of the downloaded file. 
KINSHIP_TABLE_FILE_PATH = os.path.join(GENETICS_DIR, 'ukb12345_rel_s654321.dat')

# If you wish to include genotyping metadata as part of the parsed covariates (which the default, unless you specify use_genotyping_metadata = False, e.g.
# through the parse_dataset_covariates_kwargs argument), then you will also have to download this file.
# Files in EGAD00010001225 can be downloaded using, for example, the EGA Download Client (https://ega-archive.org/download/using-ega-download-client)
# or through FTP / AsperaConnect clients (contact the UK Biobank Access Team if you need help downloading this data).
SAMPLE_GENOTYPING_METADATA_FILE_PATH = os.path.join(GENETICS_DIR, 'EGAD00010001225/001/ukb_sqc_v2.txt')

# If you wish to include genotyping metadata as part of the parsed covariates (which is the default; see previous setting), then you will have to download
# this file from: https://biobank.ctsu.ox.ac.uk/crystal/refer.cgi?id=531.
GENETICS_DESCRIPTION_FILE_PATH = os.path.join(BASE_DIR, 'ukb_genetic_file_description.txt')

# Create a folder for field codings, and download into it the following codings:
# - coding10.tsv (assessment centers; required only when include_assessment_centers = True, which is the default) from:
# http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=10
# - coding19.tsv (ICD-10 tree; required only when parsing a dataset derived from ICD-10 codes) from: http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=19
CODINGS_DIR = os.path.join(BASE_DIR, 'codings')

# If you need to work with raw-marker variants (see the CALL_DIR setting above), then you need to also provide a directory in which you have writing
# permissions, so that ukbb_parser will be able to create (rather small) temporary files (mostly links) used by pandas-plink. 
TMP_DIR_FOR_RAW_MARKER_LINKS = '/tmp'