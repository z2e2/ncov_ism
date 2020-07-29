import logging
import pandas as pd

logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)

def read_alignment(filename):
    """
    Take in a multiple sequence alignment(msa) output file and output a list
    Parameters
    ----------
    filename: str
        filepath to msa file
    Returns
    -------
    seq_list: list
        list containing msa
    """
    seq_list = []
    seq = ''
    with open(filename) as f:
        for line in f:
            if line[0] == '>':
                if len(seq) > 0:
                    seq_list.append((header, seq.upper()))
                    seq = ''
                    header = line[1:].strip('\n')
                else:
                    seq = ''
                    header = line[1:].strip('\n')
            else:
                seq += line.strip('\n')
        if len(seq) > 0:
            seq_list.append((header, seq.upper()))
    return seq_list

def load_msa_data(filename):
    """
    Take in a multiple sequence alignment(msa) output file and output a dataframe
    Parameters
    ----------
    filename: str
        filepath to msa file
    Returns
    -------
    seq_df: pandas.DataFrame
        Pandas dataframe containing msa
    """
    seq_list = read_alignment(filename)
    seq_dict = {'gisaid_epi_isl': [], 'sequence': []}
    for header, seq in seq_list:
        header = header.split('|')[1]
        seq_dict['gisaid_epi_isl'].append(header)
        seq_dict['sequence'].append(seq)
    seq_df = pd.DataFrame.from_dict(seq_dict)
    return seq_df

def merge_metadata(seq_df, meta_df):
    """
    Takes msa dataframe and metadata dataframe as input and output a merged dataframe.
    The two dataframes are merged by gisaid_epi_isl column
    Parameters
    ----------
    seq_df: pandas.DataFrame
        msa dataframe
    meta_df: pandas.DataFrame
        metadata dataframe
    Returns
    -------
    data_df: pandas.DataFrame
        merged Pandas dataframe
    """
    # join sequence with metadate
    data_df = seq_df.join(meta_df.set_index(['gisaid_epi_isl']), on = ['gisaid_epi_isl'],how = 'left')
    # filter by Human, valid date time and convert date column to 'datetime' dtype
    data_df = data_df[data_df.apply(lambda x: (x['host'] == 'Human') and ('X' not in x['date']) and len(x['date'].split('-')) == 3, axis=1)]
    data_df = data_df[data_df.apply(lambda x: (x['host'] == 'Human') and ('X' not in x['date']) and len(x['date'].split('-')) == 3, axis=1)]
    data_df['country/region'] = data_df.apply(lambda x: 'Mainland China' if x['country'] == 'China' else x['country'], axis=1)
    data_df['country/region_exposure'] = data_df.apply(lambda x: 'Mainland China' if x['country_exposure'] == 'China' else x['country_exposure'], axis=1)
    data_df['date'] = pd.to_datetime(data_df['date'])
    data_df = data_df.rename(columns={'region': 'continent', 'region_exposure': 'continent_exposure'})
    data_df = data_df.drop(['virus', 'strain', 'genbank_accession', 'country', 'title', 'country_exposure'], axis = 1)
    return data_df

def load_data(MSA_FILE_NAME, META_FILE_NAME):
    """
    Takes a multiple sequence alignment(msa) output file and metadata tsv file as input and output a merged dataframe.
    The two dataframes are merged by gisaid_epi_isl column 
    Parameters
    ----------
    MSA_FILE_NAME: str
        path to a multiple sequence alignment(msa) output file
    meta_df: pandas.DataFrame
        path to a metadata tsv file
    Returns
    -------
    data_df: pandas.DataFrame
        merged Pandas dataframe
    """
    seq_df = load_msa_data(MSA_FILE_NAME)
    meta_df = pd.read_csv(META_FILE_NAME, sep='\t')
    data_df = merge_metadata(seq_df, meta_df)
    return data_df