import os
import numpy as np
import pandas as pd


data_folder = '/data'


def reformat_integer_columns(table):
    for c in table.columns:
        if table[c].dtype == float:
            if (table[c] % 1.0).sum() == 0.0:
                table[c] = table[c].fillna(0).astype(int)
    return table


def load_ena_samples():
    return reformat_integer_columns(pd.read_csv(os.path.join(data_folder, 'sample_ena.csv.gz')))


def load_gisaid_samples():
    return reformat_integer_columns(pd.read_csv(os.path.join(data_folder, 'sample_gisaid.csv.gz')))

# ENA mutations

def load_ena_variants():
    return reformat_integer_columns(pd.read_csv(os.path.join(data_folder, 'variant.csv.gz')))
    
    
def load_ena_variant_observations():
    return reformat_integer_columns(pd.read_csv(os.path.join(data_folder, 'variant_observation.csv.gz')))


def load_ena_variant_cooccurrence():
    return reformat_integer_columns(pd.read_csv(os.path.join(data_folder, 'variant_cooccurrence.csv.gz')))


def load_ena_subclonal_variants():
    return reformat_integer_columns(pd.read_csv(os.path.join(data_folder, 'subclonal_variant.csv.gz')))
    
    
def load_ena_subclonal_variant_observations():
    return reformat_integer_columns(pd.read_csv(os.path.join(data_folder, 'subclonal_variant_observation.csv.gz')))

# GISAID mutations

def load_gisaid_variants():
    return reformat_integer_columns(pd.read_csv(os.path.join(data_folder, 'gisaid_variant.csv.gz')))
    
    
def load_gisaid_variant_observations():
    return reformat_integer_columns(pd.read_csv(os.path.join(data_folder, 'gisaid_variant_observation.csv.gz')))