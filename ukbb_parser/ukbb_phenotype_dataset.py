import inspect

import numpy as np
import pandas as pd

from .shared_utils.util import log, set_series_to_hot_encoding_df
from .ukbb_parser import create_dataset, create_ICD10_dataset, pivot_ICD10_tree_to_dataset

USER_SPECIFIED_SEXES = {
    # user_specified_sex_code: (sex_name, ukbb_sex_value)
    'F': ('female', 0),
    'M': ('male', 1),
}
        
def create_phenotype_dataset(full_phenotype_specs, nrows = None, only_caucasians = True, no_kinship = True, \
        parse_dataset_covariates_kwargs = dict(), ignore_samples_without_any_ICD10 = True, verbose = True):
    
    '''
    Create a dataset of phenotypes (and covariates) according to a list of phenotype specifications.
    @param full_phenotype_specs (a list, or anything covertable to a list, of FieldPhenotypeSpec objects, or dictionaries specifying them): The full
    phenotype specifications used to define the phenotypes for the dataset. See the documentation of FullPhenotypeSpec, FieldPhenotypeSpec,
    ICD10PhenotypeSpec and AggregationPhenotypeSpec for further explanation on how to define phenotypes. In short, FieldPhenotypeSpec can be used to define
    a phenotype based on a UKBB field (resulting in binary, continuous or one-hot-encoded phenotype values), ICD10PhenotypeSpec can be used to define 
    a binary phenotype based on a set of ICD-10 codes, and AggregationPhenotypeSpec can be used to recursively define more complex phenotypes based on any
    aggregation function of underlying phenotypes. FieldPhenotypeSpec is then used to wrap the final phenotype spec and provide additional, optional filters
    (currently ony sex_filter is supported). This function can parse a (potentialy hierarchical) dictionary defining a FullPhenotypeSpec. See
    examples/phenotype_specs.py for an example of a specification list defined by dictionaries.
    @param nrows (int or None): An upper limit to the number of read samples (or None to read all rows).
    @param only_caucasians (bool): Whether to filter out non-Caucasian samples.
    @param no_kinship (bool): Whether to filter out related individuals, leaving only one sample from each kinship group (see the function
    remove_kinships).
    @param parse_dataset_covariates_kwargs (dictionary): Additional kwargs to affect the construction of the covariates. These parameters are
    passed to the function parse_dataset_covariates.
    @param ignore_samples_without_any_ICD10 (bool): If True, all phenotypes defined by ICD-10 codes (with ICD10PhenotypeSpec) will be parsed as np.nan
    if they don't have any ICD-10 codes associated with them.
    @param verbose (bool): Whether to log details of the operation of this function.
    @return:
    1. eid (pd.Series): The eid is the primary index used by the UKBB to uniquely identify samples (also used across the genetic datasets).
    2. phenotype_values (pd.DataFrame): The parsed phenotypes, as requested by the full_phenotype_specs argument.
    3. covariates (pd.DataFrame): The parsed covariates (see parse_dataset_covariates).
    Importantly, all three returned dataframes/series have a matching index, which indicates the row number of each sample within the primary
    phenotypic CSV file (PHENOTYPE_DATASET_CSV_FILE_PATH).
    '''
    
    full_phenotype_specs = list(map(_resolve_full_phenotype_spec, full_phenotype_specs))
    required_ukbb_field_specs = [field_spec for full_phenotype_spec in full_phenotype_specs for field_spec in \
            full_phenotype_spec.phenotype_spec.get_required_ukbb_field_specs()]
    required_ICD10_codings = set.union(*(full_phenotype_spec.phenotype_spec.get_required_ICD10_codings() for full_phenotype_spec in \
            full_phenotype_specs))
    
    if len(required_ICD10_codings) == 0:
        eid, fields, covariates = create_dataset(required_ukbb_field_specs, nrows = nrows, only_caucasians = only_caucasians, \
                no_kinship = no_kinship, parse_dataset_covariates_kwargs = parse_dataset_covariates_kwargs, verbose = verbose)
        ICD10_tree = None
        any_ICD10_mask = None
    else:
        eid, ICD10_tree, fields, covariates, any_ICD10_mask = create_ICD10_dataset(required_ukbb_field_specs, nrows = nrows, \
                only_caucasians = only_caucasians, no_kinship = no_kinship, parse_dataset_covariates_kwargs = \
                parse_dataset_covariates_kwargs, filter_samples_without_ICD10 = False, desired_ICD10_codes = required_ICD10_codings, \
                verbose = verbose)
                
    if not ignore_samples_without_any_ICD10:
        any_ICD10_mask = pd.Series(True, index = any_ICD10_mask.index)
        
    phenotype_values = pd.concat([full_phenotype_spec.parse(fields, ICD10_tree, covariates, any_ICD10_mask, verbose = verbose) for \
            full_phenotype_spec in full_phenotype_specs], axis = 1)
    return eid, phenotype_values, covariates
    
class FullPhenotypeSpec:

    '''
    A full specification of a phenotype, determined by providing any of the phenotype specification objects below, and additional potential filters.
    '''

    def __init__(self, phenotype_spec, sex_filter = None):
        '''
        @param phenotype_spec (any spec object): The underlying phenotype specification. 
        @param sex_filter (None, 'M' or 'F'): Whether to filter and include only samples of a given sex, setting the phenotype value of other samples
        to np.nan.
        '''
        self.phenotype_spec = phenotype_spec
        self.sex_filter = sex_filter
        
    def parse(self, fields, ICD10_tree, covariates, any_ICD10_mask, verbose = True):
        
        phenotype_values = self.phenotype_spec.parse(fields, ICD10_tree, any_ICD10_mask)
        
        if self.sex_filter is not None:
            _filter_sex(self.phenotype_spec, phenotype_values, covariates, self.sex_filter, verbose)
            
        return phenotype_values

class FieldPhenotypeSpec:

    '''
    Specification of a phenotype by a UKBB field.
    '''

    def __init__(self, name, field_id, field_type, one_hot_encoding = False):
        '''
        @param name (str): The name of the phenotype.
        @param field_id (int): The ID of the field according to the UKBB (see http://biobank.ctsu.ox.ac.uk/crystal/).
        @param field_type (str or function): The type of the field, which could take any of the values provided to the field_type argument of the function
        parse_field_data (i.e. 'binary', 'continuous', 'set' or a function; see the documentation of parse_field_data). Unlike in the create_dataset function,
        'raw' or 'ignore' shouldn't be provided here as a field type.
        @param one_hot_encoding (bool): If field_type is 'set', one_hot_encoding = True will result in multiple columns for the field, one for each of the
        values it takes (using One Hot Encoding).
        '''
        self.name = name
        self.field_id = field_id
        self.field_type = field_type
        self.one_hot_encoding = one_hot_encoding
        assert not self.one_hot_encoding or self.field_type == 'set', 'one_hot_encoding = True is supported only when field_type = \'set\'.'
        
    def get_required_ukbb_field_specs(self):
        return [(self.name, self.field_id, self.field_type)]
        
    def get_required_ICD10_codings(self):
        return set()
        
    def parse(self, fields, ICD10_tree, any_ICD10_mask):

        phenotype_values = fields[self.name]
        
        if self.one_hot_encoding:
            value_headers = {value: '%s (%s)' % (self.name, value) for value in set.union(*phenotype_values)}
            return set_series_to_hot_encoding_df(phenotype_values, value_headers = value_headers)
        else:
            return phenotype_values
    
class ICD10PhenotypeSpec:

    '''
    Specification of a binary phenotype by ICD-10 codes.
    '''
    
    def __init__(self, name, codings):
        '''
        @param name (str): The name of the phenotype.
        @param codings (set of strings, or anything convertable to it): The ICD-10 codes defining the binary phenotype, which should take the value 1 for
        samples having any of the specified codes, or 0 if they have none of the specified codes (but still have other ICD-10 codes), or np.nan if they
        don't have any ICD-10 code associated with them and ignore_samples_without_any_ICD10 (within create_phenotype_dataset) is set to True (if it's
        set to False, then those samples will be parsed as 0 as well).  
        '''
        self.name = name
        self.codings = set(codings)
        
    def get_required_ukbb_field_specs(self):
        return []
        
    def get_required_ICD10_codings(self):
        return self.codings
        
    def parse(self, fields, ICD10_tree, any_ICD10_mask):
        relevant_subtree = ICD10_tree[ICD10_tree['coding'].isin(self.codings)]
        relevant_dataset = pivot_ICD10_tree_to_dataset(relevant_subtree, fields.index)
        phenotype_values = relevant_dataset.any(axis = 1).rename(self.name)
        phenotype_values.loc[~any_ICD10_mask] = np.nan
        return phenotype_values
        
class AggregationPhenotypeSpec:

    '''
    A phenotype specified by some aggregation of any number of sub-specifications.
    '''

    def __init__(self, name, subspecs, aggregation_function):
        '''
        @param name (str): The name of the phenotype.
        @param subspecs (a list of spec objects or dictionaries): The sub-specifications for the aggregation.
        @param aggregation_function (function): The aggregation function to use. Should recieve the parsed series/dataframes resulted from the
        sub-specification (i.e. the number of arguments the function is expected to get is determined by the number of sub-specifications).
        '''
        self.name = name
        self.subspecs = list(map(_resolve_phenotype_spec, subspecs))
        self.aggregation_function = aggregation_function
        
    def get_required_ukbb_field_specs(self):
        return [field_spec for subspec in self.subspecs for field_spec in subspec.get_required_ukbb_field_specs()]
        
    def get_required_ICD10_codings(self):
        return set.union(*(subspec.get_required_ICD10_codings() for subspec in self.subspecs))
        
    def parse(self, *args, **kwargs):
        return self.aggregation_function(*(subspec.parse(*args, **kwargs) for subspec in self.subspecs)).rename(self.name)
        
def _resolve_full_phenotype_spec(full_phenotype_spec):
    if isinstance(full_phenotype_spec, dict):
        return _parse_full_phenotype_spec_from_dict(full_phenotype_spec)
    elif isinstance(full_phenotype_spec, FullPhenotypeSpec):
        return full_phenotype_spec
    else:
        raise ValueError('Cannot convert type %s into FullPhenotypeSpec.' % type(full_phenotype_spec))
        
def _resolve_phenotype_spec(phenotype_spec):
    if isinstance(phenotype_spec, dict):
        phenotype_spec, extra_kwargs = _parse_phenotype_spec_from_dict(phenotype_spec)
        assert len(extra_kwargs) == 0
        return phenotype_spec
    else:
        return phenotype_spec
        
def _parse_full_phenotype_spec_from_dict(full_phenotype_spec_dict):    
    phenotype_spec, extra_kwargs = _parse_phenotype_spec_from_dict(full_phenotype_spec_dict)
    return FullPhenotypeSpec(phenotype_spec, **extra_kwargs)
    
def _parse_phenotype_spec_from_dict(phenotype_spec_dict):

    SOURCE_TO_PHENOTYPE_SPEC_TYPE = {
        'field': FieldPhenotypeSpec,
        'ICD-10': ICD10PhenotypeSpec,
        'aggregation': AggregationPhenotypeSpec,
    }
    
    phenotype_spec_type = SOURCE_TO_PHENOTYPE_SPEC_TYPE[phenotype_spec_dict['source']]
    phenotype_spec_type_init_args = set(inspect.getfullargspec(phenotype_spec_type.__init__).args[1:])
    phenotype_spec_kwargs = {arg: value for arg, value in phenotype_spec_dict.items() if arg in phenotype_spec_type_init_args}
    extra_kwargs = {arg: value for arg, value in phenotype_spec_dict.items() if arg != 'source' and arg not in \
            phenotype_spec_type_init_args}
    
    return phenotype_spec_type(**phenotype_spec_kwargs), extra_kwargs
        
def _filter_sex(phenotype_spec, phenotype_values, covariates, sex_filter, verbose):
    
    sex_name, required_sex_value = USER_SPECIFIED_SEXES[sex_filter]
    sex_mask = (covariates['sex'] == required_sex_value)
    phenotype_values.loc[~sex_mask] = np.nan
    
    if verbose:
        log('%s: Filtered %d %ss of %d total samples.' % (phenotype_spec.name, sex_mask.sum(), sex_name, len(sex_mask)))
