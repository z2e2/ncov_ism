import pandas as pd
import pickle
import numpy as np
import datetime
from ._analyzeism import time_subset
from ._pickism import entropy_analysis

class CustomizedISM:
    def __init__(self, MSA_FILE_NAME, META_FILE_NAME, id_col, loc_col, OUTPUT_FOLDER, 
                 reference_genbank_name, REFERENCE_ID, en_thres, null_thres, region_list):
        self.MSA_FILE_NAME = MSA_FILE_NAME
        self.META_FILE_NAME = META_FILE_NAME
        self.id_col = id_col
        self.loc_col = loc_col
        self.OUTPUT_FOLDER = OUTPUT_FOLDER
        self.reference_genbank_name = reference_genbank_name
        self.REFERENCE_ID = REFERENCE_ID
        self.en_thres = en_thres
        self.null_thres = null_thres
        self.region_list = region_list
        
    def read_alignment(self):
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

        seq_dict = {self.id_col: [], 'sequence': []}
        seq = ''
        with open(self.MSA_FILE_NAME) as f:
            for line in f:
                if line[0] == '>':
                    if len(seq) > 0:
                        seq_dict[self.id_col].append(header.split('|')[1])
                        seq_dict['sequence'].append(seq.upper())
                        seq = ''
                        header = line[1:].strip('\n')
                    else:
                        seq = ''
                        header = line[1:].strip('\n')
                else:
                    seq += line.strip('\n')
            if len(seq) > 0:
                seq_dict[self.id_col].append(header.split('|')[1])
                seq_dict['sequence'].append(seq.upper())
        return pd.DataFrame.from_dict(seq_dict)
    
    def merge_metadata(self, seq_df):
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
        meta_df = pd.read_csv(self.META_FILE_NAME, sep='\t')
        # join sequence with metadate
        data_df = seq_df.join(meta_df.set_index([self.id_col]), on = [self.id_col],how = 'inner')
        
        data_df['date'] = pd.to_datetime(data_df['date'])
        
        return data_df
    
    def entropy_time_series_analysis(self, data_df_all):
        '''
        perform entropy analysis on sequences collected in different time period to identify covarying positions.
        '''

        REFERENCE = (self.REFERENCE_ID, data_df_all[data_df_all[self.id_col] == self.REFERENCE_ID]['sequence'].iloc[0])

        seq_index = []
        index = 0
        for base in REFERENCE[1]:
            if base == '-':
                seq_index.append(index)
            else:
                index += 1
                seq_index.append(index)

        reference_local_index_map = np.array(seq_index)
        MAX_DATE = data_df_all['date'].max().date()
        MIN_DATE = data_df_all['date'].min().date()

        FLAG = True
        start = datetime.datetime.strptime(str(MIN_DATE), "%Y-%m-%d")
        time_list = []
        while FLAG:
            end = start + datetime.timedelta(days=7)
            time_list.append((str(start.date()), str(end.date())))
            start = end
            if end.date() > MAX_DATE:
                FLAG = False

        for start, end in time_list:
            data_df = time_subset(data_df_all, start, end)
            if data_df.shape[0] == 0:
                continue
            print('Entropy Time Series Analysis: processing data between {} and {} ({} sequences).'.format(start, end, data_df.shape[0]))
            H_list, null_freq_list = entropy_analysis(data_df)
            H_list = np.array(H_list)
            null_freq_list = np.array(null_freq_list)
            pickle.dump([H_list, null_freq_list], open('{}/ENTS_{}_{}.pkl'.format(self.OUTPUT_FOLDER, start, end), 'wb'))
