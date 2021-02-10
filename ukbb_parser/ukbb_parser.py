import os
import imp
from getpass import getuser

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from .shared_utils.util import log, create_time_measure_if_verbose, safe_symlink, swap_series_index_and_value, value_df_to_hot_encoding_df, \
        resolve_dummy_variable_trap, get_row_last_values

UKBB_PATHS_SETTINGS_FILE_PATH = os.path.expanduser('~/.ukbb_paths.py')

# See: http://biobank.ctsu.ox.ac.uk/crystal
COVARIATE_FIELDS = [
    # (field_name, field_id, field_type)
    ('sex', 31, 'binary'),
    ('year_of_birth', 34, 'continuous'),
    ('assessment_centers', 54, 'ignore'),
]

# See: http://biobank.ctsu.ox.ac.uk/crystal
ICD10_FIELDS = [
    # (field_name, field_id, field_type)
    ('main diagnoses', 41202, 'raw'),
    ('secondary diagnoses', 41204, 'raw'),
    ('cancer type', 40006, 'raw'),
    ('primary cause of death', 40001, 'raw'),
    ('secondary cause of death', 40002, 'raw'),
    ('external causes', 41201, 'raw'),
]

# See: http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=21000
ETHNIC_BACKGROUND_FIELD_ID = 21000
# See: http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=1001
WHITE_ETHNIC_BACKGROUND_CODES = {1, 1001, 1002, 1003}

# See: http://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=22006
GENETIC_ETHNIC_GROUPING_FIELD_ID = 22006
# See: http://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=1002
CAUCASIAN_GENETIC_ETHNIC_GROUPING_CODE = 1

def get_chrom_raw_marker_data(chrom):
    
    '''
    A helper function to read the UKBB's raw-marker genetic data of a given chromosome using the read_plink function in the
    pandas_plink module (https://pypi.org/project/pandas-plink/). Obviously, using this function requires this module to be installed.
    The function assumes the following paths: <CALL_DIR>/ukb_snp_chr<CHR>_v2.bim, <CALL_DIR>/ukb_cal_chr<CHR>_v2.bed and
    <FAM_FILE_PATH>.
    @param chrom (str): The name of the chromosome to load the data for (could be: '1', '2', ..., '22', 'X', 'Y', 'XY', 'MT').
    @return: The outputs returned by the pandas_plink.read_plink function (bim, fam, G).
    '''
    
    try:
        from pandas_plink import read_plink
    except ImportError:
        raise ImportError('Failed importing pandas_plink.read_plink. Make sure pandas-plink is installed. See: https://pypi.org/project/pandas-plink/.')
    
    _create_chrom_raw_marker_links(chrom)
    return read_plink(_get_chrom_raw_marker_links_path_prefix(chrom))

def get_chrom_imputation_data(chrom):

    '''
    A helper function to read the UKBB's imputation genetic data of a given chromosome using the bgen_parser module
    (https://github.com/nadavbra/bgen_parser). Obviously, using this function requires this module to be installed.
    The function assumes the following paths: <IMPUTATION_V3_DIR>/ukb_imp_chr<CHR>_v3.bgen, <IMPUTATION_V3_DIR>/ukb_imp_chr<CHR>_v3.bgi and
    <IMPUTATION_V3_DIR>/<IMPUTATION_V3_SAMPLE_FILE_NAME_TEMPLATE> % <CHR>.
    @param chrom (str): The name of the chromosome to load the data for (could be: '1', '2', ..., '22', 'X', 'Y', 'XY', 'MT').
    @return: A BgenParser object.
    '''

    try:
        from bgen_parser import BgenParser
    except ImportError:
        raise ImportError('Failed importing bgen_parser.BgenParser. Make sure bgen_parser is installed. See: https://github.com/nadavbra/bgen_parser.')
    
    bgen_file_path = os.path.join(ukbb_paths.IMPUTATION_V3_DIR, 'ukb_imp_chr%s_v3.bgen' % chrom)
    bgi_file_path = os.path.join(ukbb_paths.IMPUTATION_V3_DIR, 'ukb_imp_chr%s_v3.bgen.bgi' % chrom)
    sample_file_path = os.path.join(ukbb_paths.IMPUTATION_V3_DIR, ukbb_paths.IMPUTATION_V3_SAMPLE_FILE_NAME_TEMPLATE % chrom)
    
    return BgenParser(bgen_file_path, bgi_file_path, sample_file_path)
    
def create_dataset(fields_specs, nrows = None, only_caucasians = True, no_kinship = True, parse_dataset_covariates_kwargs = dict(), \
        verbose = True):
    
    '''
    Creates a UKBB dataset for requested fields.
    @param fields_specs (list of tuples): The specification of the requested fields. Each field is a tuple of the form (field_name, field_id,
    field_type) where:
    - field_name (str): The name of the column to assign to the parsed values of the field (could be any user-specified name).
    - field_id (int): The ID of the field according to the UKBB (see http://biobank.ctsu.ox.ac.uk/crystal/).
    - field_tpye (str or function): The type of the field, which could take any of the values provided to the field_type argument of the function
    parse_field_data (i.e. 'binary', 'continuous', 'set' or a function; see the documentation of parse_field_data), or 'raw' to include in the
    results all the columns in the raw_dataset associated with the field (potentially more than one column per field), or 'ignore' to ignore the
    field altogether and not include it in the results.
    @param nrows (int or None): An upper limit to the number of read samples (or None to read all rows).
    @param only_caucasians (bool): Whether to filter out non-Caucasian samples.
    @param no_kinship (bool): Whether to filter out related individuals, leaving only one sample from each kinship group (see the function
    remove_kinships).
    @param parse_dataset_covariates_kwargs (dictionary): Additional kwargs to affect the construction of the covariates. These parameters are
    passed to the function parse_dataset_covariates.
    @param verbose (bool): Whether to log details of the operation of this function.
    @return:
    1. eid (pd.Series): The eid is the primary index used by the UKBB to uniquely identify samples (also used across the genetic datasets).
    2. fields (pd.DataFrame): The parsed fields requested by the fields_specs argument (see the function parse_fields).
    3. covariates (pd.DataFrame): The parsed covariates (see parse_dataset_covariates).
    Importantly, all three returned dataframes/series have a matching index, which indicates the row number of each sample within the primary
    phenotypic CSV file (PHENOTYPE_DATASET_CSV_FILE_PATH).
    '''
    
    # Read the raw dataset.
    
    field_ids = _determine_required_field_ids(fields_specs, only_caucasians)
    raw_dataset = read_raw_dataset(field_ids, nrows = nrows, verbose = verbose)
    
    if verbose:
        log('Read a dataset of %d samples.' % len(raw_dataset))
    
    # Remove withdrawn samples.
    
    withdrawn_eids = get_withdrawn_eids()
    withdrawal_mask = raw_dataset['eid'].isin(withdrawn_eids)
    raw_dataset = raw_dataset[~withdrawal_mask]
    
    if verbose:
        log(('Knowing of %d samples who have wished to withdraw, %d of them are in the loaded dataset. Filtering out these records, the dataset ' + \
                'has reduced from %d to %d samples.') % (len(withdrawn_eids), withdrawal_mask.sum(), len(withdrawal_mask), len(raw_dataset)))
            
    # Remove non-Caucasians (if needed).
    if only_caucasians:
        raw_dataset = _filter_caucasians(raw_dataset, verbose)
    
    # Parse the covariate.
    covariates = parse_dataset_covariates(raw_dataset, verbose = verbose, **parse_dataset_covariates_kwargs)
    raw_dataset = raw_dataset.loc[covariates.index]
    
    # Remove kinships (if needed). It's important that the removal of related individuals will be the very last filtration step (otherwise we
    # may end up filtering out more samples than actually necessary).
    if no_kinship:
        no_kinship_mask = remove_kinships(raw_dataset['eid'], verbose = True)
        raw_dataset = raw_dataset[no_kinship_mask]
        covariates = covariates[no_kinship_mask]
    
    # Parsing the fields and returning the results
    eid = raw_dataset['eid']
    fields = parse_fields(raw_dataset, fields_specs, verbose = verbose)
    return eid, fields, covariates

def create_ICD10_dataset(additional_fields_specs = [], nrows = None, only_caucasians = True, no_kinship = True, \
        parse_dataset_covariates_kwargs = dict(), filter_samples_without_ICD10 = True, desired_ICD10_codes = None, \
        filter_empty_nodes = True, verbose = True):
    
    '''
    Creates a UKBB dataset in the form of an ICD-10 tree, for additioank requested fields (if provided).
    @param additional_fields_specs (list of tuples): If provided, will return, in addition to the ICD-10 tree, a dataframe with the additional
    requested fields (see create_dataset).
    @param nrows (int or None): An upper limit to the number of read samples (or None to read all rows).
    @param only_caucasians (bool): Whether to filter out non-Caucasian samples.
    @param no_kinship (bool): Whether to filter out related individuals, leaving only one sample from each kinship group (see the function
    remove_kinships).
    @param parse_dataset_covariates_kwargs (dictionary): Additional kwargs to affect the construction of the covariates. These parameters are
    passed to the function parse_dataset_covariates.
    @param filter_samples_without_ICD10 (bool): Whether to filter out samples without any ICD-10 code.
    @param desired_ICD10_codes (set of strings or None): If not None, will construct the ICD-10 tree using only nodes descending from these
    ICD-10 codes (making the construction go faster).
    @param filter_empty_nodes (bool): Whether to trim the ICD-10 tree of nodes without any samples.
    @param verbose (bool): Whether to log details of the operation of this function.
    @return:
    1. eid (pd.Series): The eid is the primary index used by the UKBB to uniquely identify samples (also used across the genetic datasets).
    2. ICD10_tree (pd.DataFrame): The created ICD-10 tree, where each row represents an ICD-10 node. The samples associated with each node
    are present in the 'samples' column, which contains, for each ICD-10 node, a set of strings indicating the sample IDs. These sample IDs
    correspond to the index of the three other dataframes/series (eid, additional_fields and covariates). 
    3. additional_fields (pd.DataFrame): The parsed fields requested by the additional_fields_specs argument (see the function parse_fields). If
    no additional fields were requested, it will be an empty dataframe without any columns.
    4. covariates (pd.DataFrame): The parsed covariates (see parse_dataset_covariates).
    5. any_ICD10_mask (pd.Series of booleans): Whether each sample as at least one ICD-10 coding.
    Note that eid, additional_fields, covariates and any_ICD10_mask have a matching index, which indicates the row number of each sample within
    the primary phenotypic CSV file (PHENOTYPE_DATASET_CSV_FILE_PATH). ICD10_tree, on the other hand, is indexed by nodes.
    '''
        
    # Create the basic dataset.
    all_fields_specs = additional_fields_specs + ICD10_FIELDS
    eid, all_fields, covariates = create_dataset(all_fields_specs, nrows = nrows, only_caucasians = only_caucasians, no_kinship = no_kinship, \
            parse_dataset_covariates_kwargs = parse_dataset_covariates_kwargs, verbose = verbose)
    
    # Split the created dataset into raw_ICD10_dataset and additional_fields.
    ICD10_columns = {column for _, field_id, _ in ICD10_FIELDS for column in _get_field_columns(all_fields.columns, field_id)}
    raw_ICD10_dataset = all_fields[list(ICD10_columns)]
    additional_fields = all_fields[[column for column in all_fields.columns if column not in ICD10_columns]]
    any_ICD10_mask = pd.notnull(raw_ICD10_dataset).any(axis = 1)
      
    # Filter samples without ICD-10 codes (if needed).
    if filter_samples_without_ICD10:
        log('Filtering out %d samples wihtout any ICD-10 codes, keeping %d of %d samples.' % ((~any_ICD10_mask).sum(), \
                any_ICD10_mask.sum(), len(any_ICD10_mask)))
        eid = eid[any_ICD10_mask]
        raw_ICD10_dataset = raw_ICD10_dataset[any_ICD10_mask]
        additional_fields = additional_fields[any_ICD10_mask]
        covariates = covariates[any_ICD10_mask]
        any_ICD10_mask = any_ICD10_mask[any_ICD10_mask]
                    
    # Parse the raw_ICD10_dataset into an ICD10_tree.
    with create_time_measure_if_verbose('Parsing the read dataset into an ICD-10 tree...', verbose):
        ICD10_tree = parse_raw_ICD10_dataset_into_tree(raw_ICD10_dataset, desired_ICD10_codes = desired_ICD10_codes, \
                filter_empty_nodes = filter_empty_nodes, verbose = verbose)
    
    # Return the results.
    return eid, ICD10_tree, additional_fields, covariates, any_ICD10_mask
    
def pivot_ICD10_tree_to_dataset(ICD10_tree, samples_index):
    
    '''
    Pivots a dataframe representing an ICD-10 tree (where each row represents an ICD-10 nodes) into a sample-centric dataframe, where each
    row represents a sample.
    @param ICD10_tree (pd.DataFrame): The ICD-10 tree to pivot (see create_ICD10_dataset). The dataframe representing the ICD-10 tree (with 
    each node representing an ICD-10 node) should have two columns:
    1) 'meaning' (each entry a string): the disease name of each ICD-10 code.
    2) 'samples' (each entry a set of strings): the sample IDs associated with each ICD-10 code, where a sample ID should correspond to the
    provided samples_index (will be used as the index of the resulted dataframe).
    @param samples_index (pd.Index): The index of the samples to use, which will become the index of the returned dataframe. It should
    correspond to the sample IDs in ICD10_tree['samples'].
    @return: A pd.DataFrame with the pivoted dataset. The index of the dataframe will correspond to samples (it will be samples_index), and
    the columns will correspond to ICD-10 conditions (according to ICD10_tree['meaning']). The values of the dataset will be binary (0 or 1),
    describing which sample is associated with which ICD-10 codes.
    '''
    
    ICD10_dataset = pd.DataFrame(0, index = samples_index, columns = ICD10_tree['meaning'])
    samples_series = pd.Series(samples_index, index = samples_index)
    
    for _, ICD10_node in ICD10_tree.iterrows():
        ICD10_dataset.loc[samples_series.isin(ICD10_node['samples']), ICD10_node['meaning']] = 1
        
    return ICD10_dataset
    
def filter_ICD10_tree(ICD10_tree, desired_ICD10_codes, verbose = True):
    
    '''
    Filters an ICD-10 tree, trimming all the nodes that are not descendent from a given set of ICD-10 codes.
    @param ICD10_tree (pd.DataFrame): The tree to trim (see construct_ICD10_tree).
    @param desired_ICD10_codes (set of strings): The ICD-10 codes to keep.
    @param verbose (bool): Whether to log details of the operation of this function.
    @return: A pd.DataFrame with the kept sub-tree (i.e. with only the relevant subset of the rows of the original dataframe, corresponding
    to the ICD-10 nodes that should be kept).
    Note: The original given ICD10_tree is not changed.
    '''
    
    if verbose:
        log('Filtering the ICD-10 tree to keep only nodes descending from %d specific codes...' % len(desired_ICD10_codes))
    
    tree_mask = pd.Series(False, index = ICD10_tree.index)
    coding_to_node_id = ICD10_tree.reset_index().set_index('coding')['node_id']
    node_ids_to_visit = {coding_to_node_id[coding] for coding in desired_ICD10_codes}
    visited_node_ids = set()
        
    while len(node_ids_to_visit) > 0:
        node_id = node_ids_to_visit.pop()
        visited_node_ids.add(node_id)
        node_ids_to_visit.update(ICD10_tree.loc[node_id, 'children_ids'] - visited_node_ids)
        tree_mask[node_id] = True
    
    if verbose: 
        log('Remained with %d of %d nodes in the ICD-10 tree.' % (tree_mask.sum(), len(tree_mask)))
    
    return ICD10_tree[tree_mask]

def construct_ICD10_tree(desired_ICD10_codes = None, verbose = True):

    '''
    Constructs an ICD-10 tree (not associated with samples at this point).
    @param desired_ICD10_codes (set of strings or None): If not None, will construct the ICD-10 tree using only nodes descending from these
    ICD-10 codes.
    @param verbose (bool): Whether to log details of the operation of this function.
    @return: A pd.DataFrame representing the constructed ICD-10 tree, where each row represents an ICD-10 node and pointing to its direct
    parent and all descending children.
    '''

    ICD10_tree = pd.read_csv(os.path.join(ukbb_paths.CODINGS_DIR, 'coding19.tsv'), sep = '\t').set_index('node_id')
    _add_children_ids(ICD10_tree)

    if desired_ICD10_codes is not None:
        ICD10_tree = filter_ICD10_tree(ICD10_tree, desired_ICD10_codes, verbose = verbose)
            
    return ICD10_tree
    
def parse_raw_ICD10_dataset_into_tree(raw_ICD10_dataset, desired_ICD10_codes = None, filter_empty_nodes = True, verbose = True):
    
    '''
    Parses a raw dataset with ICD-10 codings (e.g. from the fields specified in ICD10_FIELDS) into a structured ICD-10 tree which is
    indexed by ICD-10 nodes.
    @param raw_ICD10_dataset (pd.DataFrame): The raw dataset with the ICD-10 codings. Must not contain any other values (which are not
    ICD-10 codings).
    @param desired_ICD10_codes (set of strings or None): If not None, will construct the ICD-10 tree using only nodes descending from these
    ICD-10 codes (making the construction go faster).
    @param filter_empty_nodes (bool): Whether to trim the ICD-10 tree of nodes without any samples.
    @param verbose (bool): Whether to log details of the operation of this function.
    @return: A pd.DataFrame describing the created ICD-10 tree, where each row represents an ICD-10 node. The samples associated with each
    node are present in the 'samples' column, which contains, for each ICD-10 node, a set of strings indicating the sample IDs. These sample
    IDs correspond to the index of each sample in the provided raw_ICD10_dataset. 
    '''
    
    ICD10_tree = construct_ICD10_tree(desired_ICD10_codes = desired_ICD10_codes, verbose = verbose)
    relevant_codings = set(ICD10_tree['coding'])
    coding_to_ancestor_node_indices = _get_ICD10_tree_coding_to_ancestor_node_indices(ICD10_tree)
    tree_sample_sets = [set() for _ in range(len(ICD10_tree))]
    
    for sample, raw_sample_ICD10_codings in raw_ICD10_dataset.iterrows():
        for coding in set(raw_sample_ICD10_codings.dropna().unique()) & relevant_codings:
            for node_index in coding_to_ancestor_node_indices[coding]:
                tree_sample_sets[node_index].add(sample)
    
    ICD10_tree['samples'] = tree_sample_sets
    ICD10_tree['n_samples'] = ICD10_tree['samples'].apply(len)
    
    if filter_empty_nodes:
        
        ICD10_tree_mask = (ICD10_tree['n_samples'] > 0)
        ICD10_tree = ICD10_tree[ICD10_tree_mask]
        
        if verbose:
            log('Keeping only %d ICD-10 nodes that contain any samples (of %d nodes).' % (ICD10_tree_mask.sum(), len(ICD10_tree_mask)))
    
    return ICD10_tree

def read_raw_dataset(field_ids, nrows = None, verbose = True):

    '''
    Read a raw dataset from UKBB's primary CSV file for a set of requested fields.
    @param field_ids (set (or something convertable to a set) of ints): Will read only columns relevant to these fields (and to
    covariate fields). 
    @param nrows (int or None): An upper limit to the number of read samples (or None to read all rows).
    @param verbose (bool): Whether to log details of the operation of this function.
    @return: A pd.DataFrame of the raw dataset that was read.
    '''
    
    covariate_field_ids = set(list(zip(*COVARIATE_FIELDS))[1])
    all_field_ids = covariate_field_ids | set(field_ids)

    headers = pd.read_csv(ukbb_paths.PHENOTYPE_DATASET_CSV_FILE_PATH, encoding = ukbb_paths.PHENOTYPE_DATASET_CSV_FILE_ENCODING, nrows = 1).columns
    relevant_headers = [i for i, header in enumerate(headers) if header == 'eid' or \
            any([header.startswith('%d-' % field_id) for field_id in all_field_ids])]
    ICD10_dtypes = {i: str for i in relevant_headers if any(headers[i].startswith('%d-' % field_id) for _, field_id, _ in ICD10_FIELDS)}

    with create_time_measure_if_verbose('Reading %s dataset rows of %d columns (for %d fields)...' % ('all' if nrows is None else nrows, \
            len(relevant_headers), len(all_field_ids)), verbose):
        
        return pd.read_csv(ukbb_paths.PHENOTYPE_DATASET_CSV_FILE_PATH, encoding = ukbb_paths.PHENOTYPE_DATASET_CSV_FILE_ENCODING, usecols = relevant_headers, \
                nrows = nrows, dtype = ICD10_dtypes)

def parse_dataset_covariates(raw_dataset, use_genotyping_metadata = True, include_const = True, include_assessment_centers = True, \
        include_batch = True, verbose = True):

    '''
    Parses a collection of numerical covariates from a raw dataset read from the UKBB's primary CSV file (and may extract additional data from
    the metadata of the genotyping data, see get_sample_genotyping_metadata). The included covariates are:
    - const: a constant 1
    - sex: male (1) or female (0)
    - year_of_birth
    - 40 PCs of the genetic data
    - assessment center (one hot encoding)
    - genotyping batch (one hot encoding)
    @param raw_dataset (pd.DataFrame): The raw data extracted from the phenotypic CSV file of the UKBB, with all the columns associated with
    the required fields provided by COVARIATE_FIELDS (other columns will be ignored).
    @param use_genotyping_metadata (bool): Whether to use genotyping metadata for additional covariates (40 PCs and batch) and to further filter 
    samples lacking genotyping data or with mismatching genetic vs. self-reported sex. Note that include_batch is ignored when this is set to
    False.
    @param include_const (bool): Whether to include the constant (1) in the covariates.
    @param include_assessment_centers (bool): Whether to include the assessment centers in the covariates (using one-hot-encoding).
    @param include_batch (bool): Whether to include the batch in the covariates (using one-hot-encoding). This parameter is ignored when 
    use_genotyping_metadata = False.
    @param verbose (bool): Whether to log details of the operation of this function.
    @return: A pd.DataFrame with the parsed covariates.     
    '''

    basic_covariates = parse_fields(raw_dataset, COVARIATE_FIELDS)
    basic_covariates['sex'] = basic_covariates['sex'].astype(float)
    
    # Each set of covariates will be added to this list as a tuple of (ordering, DataFrame/Series), where ordering is only used to
    # determine the final order of the final dataframe.
    parsed_covariates = [(1, basic_covariates)]
    
    if include_const:
        const = pd.Series(np.ones(len(basic_covariates)), index = basic_covariates.index).rename('const')
        parsed_covariates.append((0, const))
    
    if include_assessment_centers:
        parsed_covariates.append((3, _parse_assessment_center_covariates_as_one_hot_encoding(raw_dataset)))

    if use_genotyping_metadata:
        
        genotyping_metadata, genetic_mask = _parse_genotyping_metadata_for_covariates(raw_dataset['eid'], basic_covariates['sex'], verbose = verbose)
        parsed_covariates = [(ordering, covariates[genetic_mask]) for ordering, covariates in parsed_covariates]
        parsed_covariates.append((2, genotyping_metadata.loc[:, 'PC1':'PC40']))
        
        if include_batch:
            parsed_covariates.append((4, _parse_batch_covariates_as_one_hot_encoding(genotyping_metadata['batch'])))
            
    return pd.concat([covariates for ordering, covariates in sorted(parsed_covariates)], axis = 1)

def parse_fields(raw_dataset, fields_specs, verbose = True):
    
    '''
    Parses a collection of phenotypic fields from a raw dataset read from the UKBB's primary CSV file.
    @param raw_dataset (pd.DataFrame): The raw data extracted from the phenotypic CSV file of the UKBB, with all the columns associated with
    the requested fields (other columns will be ignored). Note that each field typically has more than one column, e.g. to account for multiple
    visits of the sample
    @param fields_specs (list of tuples): The specification of the requested fields. Each field is a tuple of the form (field_name, field_id,
    field_type) where:
    - field_name (str): The name of the column to assign to the parsed values of the field (could be any user-specified name).
    - field_id (int): The ID of the field according to the UKBB (see http://biobank.ctsu.ox.ac.uk/crystal/).
    - field_tpye (str or function): The type of the field, which could take any of the values provided to the field_type argument of the function
    parse_field_data (i.e. 'binary', 'continuous', 'set' or a function; see the documentation of parse_field_data), or 'raw' to include in the
    results all the columns in the raw_dataset associated with the field (potentially more than one column per field), or 'ignore' to ignore the
    field altogether and not include it in the results.
    @param verbose (bool): Whether to log details of the operation of this function.
    @return: A pd.DataFrame with the parsed field values (where the raw values of most fields are aggregated into a single column). 
    '''
        
    parsed_fields = []
    
    for field_name, field_id, field_type in fields_specs:
        
        if field_type == 'ignore':
            continue
            
        raw_field_dataset = raw_dataset[_get_field_columns(raw_dataset.columns, field_id)]
            
        if field_type == 'raw':
            parsed_fields.append(raw_field_dataset)
        else:
        
            if verbose:
                log('Parsing field %s...' % field_name)
            
            parsed_fields.append(parse_field_data(raw_field_dataset, field_type).rename(field_name))
        
    return pd.concat(parsed_fields, axis = 1)

def get_sample_genotyping_metadata():

    '''
    Extracts the following sample metadata fields from from the raw-marker genetic data: batch, genetic_sex, submitted_sex, inferred_sex, PCs.
    @return: A pd.DataFrame of the extracted data, indexed by sample_id.
    '''
    
    PC_COLUMNS = ['PC%d' % (i + 1) for i in range(40)]
    
    result = parse_fam()[['batch', 'genetic_sex']].copy()
    result['genetic_sex'] = result['genetic_sex'].map({2: 0, 1: 1, 0: np.nan})
    
    sample_genotyping_metadata = parse_sample_genotyping_metadata()
    sample_genotyping_metadata.index = result.index

    result['submitted_sex'] = sample_genotyping_metadata['Submitted.Gender'].map({'F': 0, 'M': 1}).values
    result['inferred_sex'] = sample_genotyping_metadata['Inferred.Gender'].map({'F': 0, 'M': 1}).values
    result[PC_COLUMNS] = sample_genotyping_metadata[PC_COLUMNS]
    
    return result

def parse_sample_genotyping_metadata():
    '''
    Parses the metadata of the samples, according to the ukb_sqc_v2.txt file provided by the UKBB as part of the raw-marker genetic data.
    @return: A pd.DataFrame of the parsed data.
    '''
    return pd.read_csv(ukbb_paths.SAMPLE_GENOTYPING_METADATA_FILE_PATH, sep = ' ', names = _get_sample_genotyping_metadata_headers())

def parse_fam():
    '''
    Parses the UKBB's PLINK .fam file, which contains the metadata of the samples in the raw-marker genetic data.
    @return: A pd.DataFrame of the parsed FAM fields, indexed by sample_id.
    '''
    HEADERS = ['family_id', 'sample_id', 'father_id', 'mother_id', 'genetic_sex', 'batch']
    return pd.read_csv(ukbb_paths.FAM_FILE_PATH, sep = ' ', names = HEADERS).set_index('sample_id', drop = True)

def parse_field_data(raw_field_dataset, field_type, null_values = {-1, -3}):

    '''
    Parses the raw columns associated with a specific field (typically each field has more than one column, e.g. to account for multiple visits
    of the sample).
    @param raw_field_dataset (pd.DataFrame): The raw data extracted from the phenotypic CSV file of the UKBB of all the columns associated with
    the given field (it's important not to accidentally provide columns of other fields).
    @param field_type ('binary', 'continuous', 'set' or a function): If a function, will just apply it on the received dataframe (after determining
    null values). If 'binary', will be 1 if there if at least one of the columns is 1, 0 if all the columns are 0, or null if all the columns are null.
    If 'continuous', will take the last not-null value of each row. If 'set', will return a series of sets, where each set contains all the not-null
    values collected from the row.
    @param null_values (set of numbers): The values -1 (Do not know) and -3 (Prefer not to answer) are interpreted as null (np.nan) by default.
    @return: A pd.Series of the parsed values, with an index corresponding to the index of the original dataframe.
    '''
    
    raw_field_dataset = raw_field_dataset.copy()
    raw_field_dataset[raw_field_dataset.isin(null_values)] = np.nan
    
    if callable(field_type):
        return field_type(raw_field_dataset)
    elif field_type == 'binary':
        return raw_field_dataset.max(axis = 1)
    elif field_type == 'continuous':
        return get_row_last_values(raw_field_dataset)
    elif field_type == 'set':
        return raw_field_dataset.apply(lambda raw_values: set(raw_values.dropna().unique()), axis = 1)
    else:
        raise Exception('Unexpected field type: %s' % field_type)
    
def get_withdrawn_eids():
    
    '''
    According to UKBB rules, participants who have wished to withdraw from the database must be removed from any future analysis.
    The IDs (eid) of these participants are recieved through periodic emails containing CSV files. Whenever receiving such an email, just
    put the attached CSV file in the directory specified by the PARTICIPANT_WITHDRAWALS_DIR configuration. If you do not have such directory,
    just set the value PARTICIPANT_WITHDRAWALS_DIR to None, and no samples will be removed (not recommended).
    @return: A set of strings (the strings are the sample IDs, i.e. eid, of the withdrawn samples).
    '''
    
    withdrawn_eids = set()

    if ukbb_paths.PARTICIPANT_WITHDRAWALS_DIR is not None:
        for file_name in os.listdir(ukbb_paths.PARTICIPANT_WITHDRAWALS_DIR):
            withdrawn_eids |= set(pd.read_csv(os.path.join(ukbb_paths.PARTICIPANT_WITHDRAWALS_DIR, file_name), names = ['eid'])['eid'])
        
    return withdrawn_eids
            
def remove_kinships(eid, verbose = True):

    '''
    Determines which samples need to be removed such that the remaining samples will have no kinship connections whatsoever (according to the
    kinship table provided by the UKBB). In order to determine that, kinship groups will first be determined (@see get_kinship_groups), and 
    only one sample will remain within each of the groups. For the sake of determinism, the sample with the lowest eid will be selected within
    each kinship group, and the rest will be discarded.
    @param eid (pd.Series): A series whose values are UKBB sample IDs, from which kinships should be removed.
    @param verbose (bool): Whether to log details of the operation of this function.
    @return: A mask of samples to keep (pd.Series with index corresponding to the eid input, and boolean values).
    '''
    
    all_eids = set(eid)
    kinship_groups = get_kinship_groups()
    
    relevant_kinship_groups = [kinship_group & all_eids for kinship_group in kinship_groups]
    relevant_kinship_groups = [kinship_group for kinship_group in relevant_kinship_groups if len(kinship_group) >= 2]
    unchosen_kinship_representatives = set.union(*[set(sorted(kinship_group)[1:]) for kinship_group in relevant_kinship_groups])
    no_kinship_mask = ~eid.isin(unchosen_kinship_representatives)
    
    if verbose:
        log(('Constructed %d kinship groups (%d samples), of which %d (%d samples) are relevant for the dataset (i.e. containing at least 2 ' + \
                'samples in the dataset). Picking only one representative of each group and removing the %d other samples in those groups ' + \
                'has reduced the dataset from %d to %d samples.') % (len(kinship_groups), len(set.union(*kinship_groups)), \
                len(relevant_kinship_groups), len(set.union(*relevant_kinship_groups)), len(unchosen_kinship_representatives), len(no_kinship_mask), \
                no_kinship_mask.sum()))
    
    return no_kinship_mask
    
def get_kinship_groups():

    '''
    Uses the kinship table provided by the UKBB (as specified by the KINSHIP_TABLE_FILE_PATH configuration) in order to determine kinship groups.
    Each kinship group is a connected component of samples in the graph of kinships (where each node is a UKBB sample, and an edge exists between
    each pair of samples reported in the kinship table).
    @return: A list of sets of strings (the strings are the sample IDs, i.e. eid). Each set of samples is a kinship group.
    '''
    
    kinship_table = pd.read_csv(ukbb_paths.KINSHIP_TABLE_FILE_PATH, sep = ' ')
    kinship_ids = np.array(sorted(set(kinship_table['ID1']) | set(kinship_table['ID2'])))
    n_kinship_ids = len(kinship_ids)
    kinship_id_to_index = pd.Series(np.arange(n_kinship_ids), index = kinship_ids)

    kinship_index1 = kinship_table['ID1'].map(kinship_id_to_index).values
    kinship_index2 = kinship_table['ID2'].map(kinship_id_to_index).values

    symmetric_kinship_index1 = np.concatenate([kinship_index1, kinship_index2])
    symmetric_kinship_index2 = np.concatenate([kinship_index2, kinship_index1])

    kinship_matrix = csr_matrix((np.ones(len(symmetric_kinship_index1), dtype = bool), (symmetric_kinship_index1, \
            symmetric_kinship_index2)), shape = (n_kinship_ids, n_kinship_ids), dtype = bool)

    _, kinship_labels = connected_components(kinship_matrix, directed = False)
    kinship_labels = pd.Series(kinship_labels, index = kinship_ids)
    return [set(group_kinship_labels.index) for _, group_kinship_labels in kinship_labels.groupby(kinship_labels)]
    
def _determine_required_field_ids(fields_specs, only_caucasians):

    field_ids = set(list(zip(*fields_specs))[1])
    
    if only_caucasians:
        field_ids.update({ETHNIC_BACKGROUND_FIELD_ID, GENETIC_ETHNIC_GROUPING_FIELD_ID})
            
    return field_ids
    
def _get_field_columns(columns, field_id):
    return [column for column in columns if column.startswith('%d-' % field_id)]
    
def _filter_caucasians(raw_dataset, verbose):

    # Look at self-reported ethnicity
    raw_ethnic_background = raw_dataset[_get_field_columns(raw_dataset.columns, ETHNIC_BACKGROUND_FIELD_ID)]
    first_white_code = next(iter(WHITE_ETHNIC_BACKGROUND_CODES))
    self_reported_white_mask = pd.notnull(raw_ethnic_background).any(axis = 1) & raw_ethnic_background.fillna(first_white_code)\
            .astype(int).isin(WHITE_ETHNIC_BACKGROUND_CODES).all(axis = 1)
    
    # Look at genetic ethnicity
    raw_genetic_ethnic_grouping = raw_dataset[_get_field_columns(raw_dataset.columns, GENETIC_ETHNIC_GROUPING_FIELD_ID)]
    genetic_caucasian_mask = (raw_genetic_ethnic_grouping == CAUCASIAN_GENETIC_ETHNIC_GROUPING_CODE).any(axis = 1)
    
    # Combine the two masks
    caucasian_mask = self_reported_white_mask & genetic_caucasian_mask
    
    if verbose:
        log(('Of %d samples, %d are self-reported whites and %d are Caucasians according to their genetics. Keeping only ' + \
                'the %d samples that are both.') % (len(raw_dataset), self_reported_white_mask.sum(), genetic_caucasian_mask.sum(), \
                caucasian_mask.sum()))
    
    return raw_dataset[caucasian_mask]
    
def _parse_genotyping_metadata_for_covariates(eid, sex, verbose):

    genotyping_metadata = get_sample_genotyping_metadata().reindex(eid.values).set_index(eid.index)
    mask = pd.notnull(genotyping_metadata).all(axis = 1)
    
    if verbose:
        log('Filtering out %d samples without genotyping metadata.' % (~mask).sum())
                
    sex = sex[mask]
    genotyping_metadata = genotyping_metadata[mask]
    
    assert (genotyping_metadata['genetic_sex'] == genotyping_metadata['inferred_sex']).all()
    sex_mask = (sex == genotyping_metadata['genetic_sex']) & (sex == genotyping_metadata['submitted_sex'])
    
    if verbose:
        log('Filtering out %d samples with mismatching genetic and self-reported sex.' % (~sex_mask).sum())
    
    genotyping_metadata = genotyping_metadata[sex_mask]
    mask[~sex_mask] = False
    return genotyping_metadata, mask
    
def _get_sample_genotyping_metadata_headers():

    with open(ukbb_paths.GENETICS_DESCRIPTION_FILE_PATH, 'r') as f:
        specfile = f.read()
        
    headers = [line.split(' ')[0].strip() for line in specfile.splitlines()[43:70]]
    pc_index = headers.index('PC1-PC40')
    return ['Spurious.Affymetrix%d' % (i + 1) for i in range(2)] + headers[:pc_index] + \
            ['PC%d' % (i + 1) for i in range(40)] + headers[(pc_index + 1):]
    
def _parse_batch_covariates_as_one_hot_encoding(raw_batch):
    all_batches = np.sort(raw_batch.unique())
    batch_to_index = pd.Series(np.arange(len(all_batches)), index = all_batches)
    batch_covariates = pd.DataFrame(0, index = raw_batch.index, columns = all_batches)
    batch_covariates.values[np.arange(len(raw_batch)), batch_to_index.loc[raw_batch].values] = 1
    batch_covariates.rename(columns = lambda colname: 'batch_%s' % colname, inplace = True)
    resolve_dummy_variable_trap(batch_covariates, validate_completeness = False, inplace = True)
    return batch_covariates
    
def _parse_assessment_center_covariates_as_one_hot_encoding(raw_dataset):
    
    # Parse the assessment center covariates as One Hot Encoding (technically, there can be more than one).
    assessment_center_coding_to_meaning = pd.read_csv(os.path.join(ukbb_paths.CODINGS_DIR, 'coding10.tsv'), sep = '\t', \
            index_col = 'coding')['meaning'].apply(lambda meaning: 'AC_' + meaning.lower().replace(' ', '_').replace('(', '').replace(')', ''))
    assessment_center_field_id, = [field_id for field_name, field_id, _ in COVARIATE_FIELDS if field_name == 'assessment_centers']
    raw_assessment_centers = raw_dataset[_get_field_columns(raw_dataset.columns, assessment_center_field_id)]
    assessment_center_covariates = value_df_to_hot_encoding_df(raw_assessment_centers, value_headers = assessment_center_coding_to_meaning)
    
    # Split the assessment center covariates into the standard assessment centers (of which there should be exactly one per sample), and the
    # additional assessment centers (ending with "_imaging" or "_revisit"), which can be added on top of that, in order to resolve the "dummy
    # variable trap" for the standard ones but leave the additional ones.
    additional_assessment_center_filter = lambda assessment_center_name: assessment_center_name.endswith('_imaging') or \
            assessment_center_name.endswith('_revisit')
    additional_assessment_center_covariates = assessment_center_covariates[[column_name for column_name in assessment_center_covariates.columns if \
            additional_assessment_center_filter(column_name)]]
    standard_assessment_center_covariates = assessment_center_covariates[[column_name for column_name in assessment_center_covariates.columns if \
            not additional_assessment_center_filter(column_name)]]
    standard_assessment_center_covariates = resolve_dummy_variable_trap(standard_assessment_center_covariates, validate_completeness = False)
    return pd.concat([standard_assessment_center_covariates, additional_assessment_center_covariates], axis = 1)
    
def _add_children_ids(tree):
    
    tree['children_ids'] = [set() for _ in range(len(tree))]
    
    for node_id, parent_id in tree['parent_id'].iteritems():
        if parent_id != 0:
            tree.loc[parent_id, 'children_ids'].add(node_id)
            
def _get_ICD10_tree_coding_to_ancestor_node_indices(ICD10_tree):

    coding_to_ancestor_node_indices = {}
    coding_to_node_index = swap_series_index_and_value(ICD10_tree.reset_index(drop = True)['coding'])

    for _, ICD10_node in ICD10_tree.iterrows():

        original_node_coding = ICD10_node['coding']
        ancestor_node_indices = []

        while True:

            ancestor_node_indices.append(coding_to_node_index[ICD10_node['coding']])

            if ICD10_node['parent_id'] in ICD10_tree.index:
                ICD10_node = ICD10_tree.loc[ICD10_node['parent_id']]
            else:
                break

        coding_to_ancestor_node_indices[original_node_coding] = ancestor_node_indices
        
    return coding_to_ancestor_node_indices

def _create_chrom_raw_marker_links(chrom):
    chrom_raw_marker_links_path_prefix = _get_chrom_raw_marker_links_path_prefix(chrom)
    safe_symlink(os.path.join(ukbb_paths.CALL_DIR, 'ukb_snp_chr%s_v2.bim' % chrom), '%s.bim' % chrom_raw_marker_links_path_prefix)
    safe_symlink(os.path.join(ukbb_paths.CALL_DIR, 'ukb_cal_chr%s_v2.bed' % chrom), '%s.bed' % chrom_raw_marker_links_path_prefix)
    _create_plink_fixed_fam_file_if_needed('%s.fam' % chrom_raw_marker_links_path_prefix)
    
def _create_plink_fixed_fam_file_if_needed(fixed_fam_file_path):
        
    '''
    Apparently UKBB's FAM file deviates from PLINK's FAM format, by setting the 6th column (of costum values) to be non-numeric (in UKBB, they
    store the batch). It also appears that the pandas-plink module is not flexible enough to handle with such a deviation. As a result, this
    horrible patch is required: removing the 6th column from the FAM file before saving it in a temporary path, so that it could be processed
    by pandas-plink.
    '''
    
    if os.path.exists(fixed_fam_file_path):
        log('%s: already exists.' % fixed_fam_file_path)
    else:
        
        fam_df = pd.read_csv(ukbb_paths.FAM_FILE_PATH, sep = ' ', header = None)
        fam_df.iloc[:, 5] = np.nan
        
        try:
            fam_df.to_csv(fixed_fam_file_path, sep = ' ', index = False, header = False)
        except OSError as e:
            if e.errno == 17:
                log('%s: already exists after all.' % fixed_fam_file_path)
            else:
                raise e
    
def _get_chrom_raw_marker_links_path_prefix(chrom):
    return os.path.join(ukbb_paths.TMP_DIR_FOR_RAW_MARKER_LINKS, 'ukb_%s_chr%s' % (getuser(), chrom))
    
def _load_ukbb_paths():
    return imp.load_source('ukbb_paths', UKBB_PATHS_SETTINGS_FILE_PATH)
    
ukbb_paths = _load_ukbb_paths()