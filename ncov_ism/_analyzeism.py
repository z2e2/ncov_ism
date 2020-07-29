import logging
import json
import datetime
import pandas as pd
import numpy as np
from collections import Counter, OrderedDict
from math import log2, ceil

logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)

def regional_analysis(df, region):
    """
    regional analysis of ISM table
    Parameters
    ----------
    df: pandas.DataFrame
        Pandas dataframe containing ISM column
    region: str
        country/region of interest, corresponding to the 'country/region' column in ISM table
    Returns
    -------
    dict_freq: dictionary
        ISM frequency of a region of interest
    """
    df_tmp = df[df['country/region'] == region]
    total = df_tmp.sum()['count']
    ISM_list = list(df_tmp['ISM'].values)
    dict_freq = OrderedDict() 
    for ISM in ISM_list:
        freq = str(df_tmp[df_tmp['ISM'] == ISM]['count'].values[0])
        date = np.datetime_as_string(df_tmp[df_tmp['ISM'] == ISM]['date'].values[0], unit='D')
        dict_freq[ISM] = (date, freq)
    return dict_freq

def time_subset(ISM_df, start, end):
    """
    get a subset of records given start and end date
    Parameters
    ----------
    ISM_df: pandas.DataFrame
        Pandas dataframe containing ISM column
    start: str
        start date in '%Y-%m-%d' format
    end: str
        end date in '%Y-%m-%d' format
    Returns
    -------
    subset_df: pandas.DataFrame
        a subset of Pandas dataframe
    """
    # start: exclusive, end: inclusive
    filter_date_1 = pd.to_datetime(start, format='%Y-%m-%d')
    filter_date_2 = pd.to_datetime(end, format='%Y-%m-%d')
    return ISM_df[(filter_date_1 < ISM_df['date']) & (ISM_df['date'] <= filter_date_2)]

def frequency_count(df):
    """
    count the frequency of each ISM per country/region
    Parameters
    ----------
    df: pandas.DataFrame
        Pandas dataframe containing ISM column
    Returns
    -------
    df_ISM_date: pandas.DataFrame
        Pandas dataframe containing the count of each ISM per country/region
    """
    df_date = df.groupby(['country/region','ISM']).agg({'date': 'min'}).reset_index()
    df_ISM = df.groupby('country/region')['ISM'].value_counts().to_frame()
    df_ISM = df_ISM.rename(columns={'ISM': 'count'}).reset_index()
    df_ISM_date = df_ISM.join(df_date.set_index(['country/region','ISM']), on = ['country/region','ISM'],how = 'left')
    return df_ISM_date

def statewise_analysis(df, state):
    """
    regional analysis of ISM table
    Parameters
    ----------
    df: pandas.DataFrame
        Pandas dataframe containing ISM column
    state: str
        US state of interest, corresponding to the 'division' column in ISM table
    Returns
    -------
    dict_freq: dictionary
        ISM frequency of a state of interest
    """
    df_tmp = df[df['division'] == state]
    total = df_tmp.sum()['count']
    ISM_list = list(df_tmp['ISM'].values)
    dict_freq = OrderedDict() 
    for ISM in ISM_list:
        freq = str(df_tmp[df_tmp['ISM'] == ISM]['count'].values[0])
        date = np.datetime_as_string(df_tmp[df_tmp['ISM'] == ISM]['date'].values[0], unit='D')
        dict_freq[ISM] = (date, freq)   
    return dict_freq

def regional_timeseries_analysis(df, region):
    """
    regional time series analysis of ISM table
    Parameters
    ----------
    df: pandas.DataFrame
        Pandas dataframe containing ISM column
    region: str
        region of interest, corresponding to the 'country/region' column in ISM table
    Returns
    -------
    dict_freq: dictionary
        ISM frequency of a region of interest over time
    """
    df_tmp = df[df['country/region'] == region]
    total = df_tmp.sum()['count']
    ISM_list = list(df_tmp['ISM'].values)
    dict_freq = OrderedDict() 
    for ISM in ISM_list:
        freq = df_tmp[df_tmp['ISM'] == ISM]['count'].values[0]
        dict_freq[ISM] = str(freq)
    return dict_freq

def ISM_filter(dict_freq, threshold):
    """
    collapse low frequency ISMs into "OTHER" per location
    Parameters
    ----------
    dict_freq: dictionary
        ISM frequency of a location of interest
    threshold: float
        ISMs lower than this threshold will be collapsed into "OTHER"
    Returns
    -------
    res_dict: dictionary
        filtered ISM frequency of a location of interest
    """
    res_dict = {'OTHER': [0, 0]}
    total = sum([int(dict_freq[ISM][1]) for ISM in dict_freq])
    for ISM in dict_freq:
        if int(dict_freq[ISM][1])/total < threshold:
            res_dict['OTHER'] = [0, res_dict['OTHER'][1] + int(dict_freq[ISM][1])]
        else:
            res_dict[ISM] = [dict_freq[ISM][0], int(dict_freq[ISM][1]) + res_dict.get(ISM, [0, 0])[1]]
    if res_dict['OTHER'][1] == 0:
        del res_dict['OTHER']
    return res_dict

def ISM_time_series_filter(dict_freq, threshold):
    """
    collapse low frequency ISMs into "OTHER" per location
    Parameters
    ----------
    dict_freq: dictionary
        ISM frequency of a location of interest
    threshold: float
        ISMs lower than this threshold will be collapsed into "OTHER"
    Returns
    -------
    res_dict: dictionary
        filtered ISM frequency of a location of interest
    """
    res_dict = {'OTHER': [0, 0]}
    total = sum([int(dict_freq[ISM]) for ISM in dict_freq])
    for ISM in dict_freq:
        if int(dict_freq[ISM])/total < threshold:
            res_dict['OTHER'] = [0, res_dict['OTHER'][1] + int(dict_freq[ISM])]
        else:
            res_dict[ISM] = [dict_freq[ISM], int(dict_freq[ISM]) + res_dict.get(ISM, [0, 0])[1]]
    if res_dict['OTHER'][1] == 0:
        del res_dict['OTHER']
    return res_dict

def ISM_analysis(ISM_df, output_folder):
    """
    Informative Subtype Marker analysis of ISM table
    Parameters
    ----------
    ISM_df: pandas.DataFrame
        Pandas dataframe containing ISM column
    output_folder: str
        path to the output folder
    Returns
    -------
    region_raw_count: dictionary
        ISM frequency per region
    state_raw_count: dictionary
        ISM frequency per state
    count_dict: dictionary
        ISM frequency time series per region
    """
    logging.info('ISM Analysis in progress: ...')
    
    region_first_date = ISM_df.groupby(['country/region','ISM']).agg({'date': 'min'}).reset_index()
    region_ISM_count = ISM_df.groupby('country/region')['ISM'].value_counts().to_frame()
    region_ISM_count = region_ISM_count.rename(columns={'ISM': 'count'}).reset_index()
    region_ISM_count_date = region_ISM_count.join(region_first_date.set_index(['country/region','ISM']), on = ['country/region','ISM'],how = 'left')
    region_list = ISM_df['country/region'].unique().tolist()

    region_raw_count = OrderedDict() 
    for idx, region in enumerate(region_list):
        dict_freq = regional_analysis(region_ISM_count_date, region)
        region_raw_count[region] = dict_freq

    with open('{}/region_pie_chart.json'.format(output_folder), 'w') as fp:
        json.dump(region_raw_count, fp)
    
    logging.info('ISM Analysis in progress: regional analysis completed.')
    
    intra_use_first_date = ISM_df[ISM_df['country/region'] == 'USA'].groupby(['division','ISM']).agg({'date': 'min'}).reset_index()
    intra_usa_ISM_count = ISM_df[ISM_df['country/region'] == 'USA'].groupby('division')['ISM'].value_counts().to_frame()
    intra_usa_ISM_count = intra_usa_ISM_count.rename(columns={'ISM': 'count'}).reset_index()
    intra_usa_ISM_count_date = intra_usa_ISM_count.join(intra_use_first_date.set_index(['division','ISM']), on = ['division','ISM'],how = 'left')

    state_count_thres = 0
    state_list = [state for state, count in ISM_df[ISM_df['country/region'] == 'USA']['division'].value_counts().items() if count > state_count_thres and state != 'USA']
    
    state_raw_count = OrderedDict() 
    for idx, state in enumerate(state_list):
        dict_freq = statewise_analysis(intra_usa_ISM_count_date, state)
        state_raw_count[state] = dict_freq

    with open('{}/state_pie_chart.json'.format(output_folder), 'w') as fp:
        json.dump(state_raw_count, fp)
        
    logging.info('ISM Analysis in progress: state analysis completed.')
    
    start_date = datetime.date(2019, 12, 1)
    end_date = ISM_df['date'].max().date()
#     end_date = datetime.date(int(GISAID_DATA_DATE[:4]), int(GISAID_DATA_DATE[4:6]), int(GISAID_DATA_DATE[6:]))
    delta = datetime.timedelta(days=1)
    region_list = ISM_df['country/region'].unique().tolist()
    count_dict = OrderedDict() 

    while start_date <= end_date:
        df_tmp = time_subset(ISM_df, '2019-11-01', str(start_date))
        if df_tmp.shape[0] == 0:
            start_date += delta
            continue
        df_tmp_tmp = frequency_count(df_tmp)
        dict_freq = OrderedDict() 
        for region in region_list:
            regional_dict_freq = regional_timeseries_analysis(df_tmp_tmp, region)
            dict_freq[region] = regional_dict_freq
        count_dict[str(start_date)] = dict_freq
        start_date += delta

    with open('{}/region_time_series.json'.format(output_folder), 'w') as fp:
        json.dump(count_dict, fp)
    
    logging.info('ISM Analysis in progress: time series analysis completed.')
    
    logging.info('ISM Analysis in progress: DONE.')
    
    return region_raw_count, state_raw_count, count_dict
