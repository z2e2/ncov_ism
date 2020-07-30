import pickle
import pandas as pd
import numpy as np
from random import sample 
from collections import Counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

def sub_sampling(table, sampling_depth = 50):
    '''
    subsamping a table to get same sampling depth for all samples.
    '''
    sub_sampled_table = np.zeros(table.shape)
    for idx in range(table.shape[0]):
        pool = []
        for ISM_idx, count in enumerate(table[idx].tolist()):
            pool.extend([ISM_idx for i in range(int(count))])
        sub_sampled_dict = Counter(sample(pool,sampling_depth))
        for ISM_idx in sub_sampled_dict:
            sub_sampled_table[idx, ISM_idx] = sub_sampled_dict[ISM_idx]
    return sub_sampled_table

def bray_curtis_distance(table, sample1_id, sample2_id):
    '''
    compute Bray Curtis distance between two samples.
    '''
    numerator = 0
    denominator = 0
    sample1_counts = table[sample1_id]
    sample2_counts = table[sample2_id]
    for sample1_count, sample2_count in zip(sample1_counts, sample2_counts):
        numerator += abs(sample1_count - sample2_count)
        denominator += sample1_count + sample2_count
    return numerator / denominator

def build_ISM_abundance_table(data_df, sampling_depth = 150):
    '''
    convert ISM dataframe to ISM abundance table.
    '''
    NT_SET = set(['A', 'C', 'G', 'T', '-'])
    
    tmp = data_df[['continent', 'country/region']]
    country_to_continent = {}
    for i in range(tmp.shape[0]):
        continent, country = tmp.iloc[i]['continent'], tmp.iloc[i]['country/region']
        if country in country_to_continent:
            pass
        else:
            country_to_continent[country] = continent

    ISM_to_idx, idx_to_ISM, region_to_idx, idx_to_region = {}, {}, {}, {}
    region_list = []
    ISM_set = set([])
    for region, count in data_df['country/region'].value_counts().items():
        if count >= sampling_depth:
            region_list.append(region)
            tmp = data_df[data_df['country/region'] == region]
            ISM_set.update(tmp['ISM'].unique())

    country_to_continent = {country: country_to_continent[country] for country in country_to_continent if country in region_list}
    continent_list = sorted(set(country_to_continent.values()))
    color_list = ['red', 'orange', 'green', 'blue', 'brown', 'black', 'pink', 'mediumpurple','gold', 'lightskyblue', ]
    continent_to_color = {}
    for idx, item in enumerate(continent_list):
        continent_to_color[item] = color_list[idx]

    n_row, n_ISM = len(region_list), len(ISM_set)
    raw_table = np.zeros((n_row, n_ISM))
    for idx, ISM in enumerate(ISM_set):
        ISM_to_idx[ISM] = idx
        idx_to_ISM[idx] = ISM
    for idx, region in enumerate(region_list):
        region_to_idx[region] = idx
        idx_to_region[idx] = region

    for region_idx in range(n_row):
        region = idx_to_region[region_idx]
        tmp = data_df[data_df['country/region'] == region]
        for ISM, count in tmp['ISM'].value_counts().items():
            ISM_idx = ISM_to_idx[ISM]
            raw_table[region_idx, ISM_idx] += count
    return raw_table, n_row, n_ISM, country_to_continent, continent_to_color, idx_to_region

def region_pca_plot(INPUT_FOLDER, OUTPUT_FOLDER, sampling_depth=150):
    '''
    pca plot of ISM abundance table
    '''
    ISM_df = pd.read_csv('{}/ISM_df_with_correction.csv'.format(INPUT_FOLDER))
    ISM_df['date'] = pd.to_datetime(ISM_df['date'])
    
    raw_table, n_row, n_ISM, country_to_continent, continent_to_color, idx_to_region = build_ISM_abundance_table(ISM_df, sampling_depth)

    table = sub_sampling(raw_table, sampling_depth = sampling_depth)
    pickle.dump([raw_table, table, idx_to_region], open('{}/abundance_table.pkl'.format(OUTPUT_FOLDER), 'wb'))
    d_bray_curtis = np.zeros((n_row, n_row))
    for i in range(n_row):
        for j in range(n_row):
            d_bray_curtis[i, j] = bray_curtis_distance(table, i, j)
    region_bray_curtis_dist = []
    for i in range(n_row):
        for j in range(i + 1, n_row):
            region_bray_curtis_dist.append(bray_curtis_distance(table, i, j))

    from sklearn.decomposition import PCA
    pca = PCA(n_components=2)
    X_2d = pca.fit_transform(d_bray_curtis)
    
    DPI = 300
    
    fig = plt.figure(figsize=(2400/DPI, 2100/DPI), dpi=DPI)   
    n_row = X_2d.shape[0]
    seen = set([])
    for i in range(n_row):
        continent = country_to_continent[idx_to_region[i]]
        if continent in seen:
            plt.plot(X_2d[i, 0], X_2d[i, 1], 'o', color=continent_to_color[continent])
        else:
            plt.plot(X_2d[i, 0], X_2d[i, 1], 'o', color=continent_to_color[continent], label=continent)
            seen.add(continent)

    for i in range(n_row):
        plt.text(X_2d[i, 0] + 0.015, X_2d[i, 1] + 0.015, idx_to_region[i], size=8)
        
    plt.xlim([-1.01, 2.5])
    plt.legend(loc='lower right')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.savefig('{}/country_2d.png'.format(OUTPUT_FOLDER), bbox_inches='tight')
    plt.show()

