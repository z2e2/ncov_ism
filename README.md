# *ncov_ism*: Characterizing Geographical and Temporal Dynamics of Novel Coronavirus SARS-CoV-2 using Informative Subtype Markers

*ncov_ism* (Novel Coronavirus Informative Subtype Marker) is an efficient framework for genetic subtyping of a pandemic virus and implement it for SARS-CoV-2, the novel coronavirus that causes COVID-19.        
Drexel University EESI Lab, 2020        
Maintainer: Zhengqiao Zhao, zz374 at drexel dot edu        
Owner: Gail Rosen, gailr at ece dot drexel dot edu        
###### Latest report is available at [here](ISM_report.md)
###### A web map interface of ISM distribution for different regions are available at [here](https://covid19-ism.coe.drexel.edu/)

## Dependencies
The following are required:    
- python=3.5.4
- matplotlib=2.1.1
- pandas=0.20.3
- biopython=1.71
- scipy=1.0.0
- scikit-learn=0.19.1

Conda installation:
```
conda create -n ISM python=3.5.4 biopython=1.71 pandas=0.20.3 matplotlib=2.1.1
source activate ISM
git clone https://github.com/z2e2/ncov_ism.git
```

## Tutorial
The jupyter notebook demo can be found at [source code](demo.ipynb)

## Running ncov_ism
#### Input Data
SARS-CoV-2 (hCoV-19) sequences and metadata are made available at [GISAID](www.gisaid.org). Our pipeline takes the multiple sequence alignment results and metadata as input which are both available at [GISAID](www.gisaid.org).
#### Reference genomes
The default reference sequences is [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2/) in NCBI. It also corresponds to *EPI_ISL_402125* in [GISAID](www.gisaid.org) hCoV-19 database.
#### Quick start
Once *ncov_ism* is installed, you can run the default setting using the following command with minimum number of parameters:
```
ncov_ism default -i <multiple sequence alignment file> -m <metadata file> -o <outpu directory>
```

## About
The novel coronavirus responsible for COVID-19, SARS-CoV-2, expanded to reportedly 8.7 million confirmed cases worldwide by June 21, 2020. The global SARS-CoV-2 pandemic highlights the importance of tracking viral transmission dynamics in real-time. Through June 2020, researchers have obtained genetic sequences of SARS-CoV-2 from over 50 thousand samples from infected individuals worldwide. Since the virus readily mutates, each sequence of an infected individual contains useful information linked to the individual's exposure location and sample date. But, there are over 30,000 bases in the full SARS-CoV-2 genome, so tracking genetic variants on a whole-sequence basis becomes unwieldy. *ncov_ism* is a method to instead efficiently identify and label genetic variants, or "subtypes" of SARS-CoV-2. This method defines a compact set of nucleotide sites that characterize the most variable (and thus most informative) positions in the viral genomes sequenced from different individuals, called an Informative Subtype Marker or *ISM*. This tool defines viral subtypes for each ISM, and analyze the regional distribution of subtypes to track the progress of the pandemic.

## Reference
Preprint in BioRxiv [Genetic Grouping of SARS-CoV-2 Coronavirus Sequences using Informative Subtype Markers for Pandemic Spread Visualization](https://www.biorxiv.org/content/10.1101/2020.04.07.030759v5)
If you find our work helpful, please cite:
```
@article {Zhao2020.04.07.030759,
	author = {Zhao, Zhengqiao and Sokhansanj, Bahrad A. and Malhotra, Charvi and Zheng, Kitty and Rosen, Gail L.},
	title = {Genetic Grouping of SARS-CoV-2 Coronavirus Sequences using Informative Subtype Markers for Pandemic Spread Visualization},
	elocation-id = {2020.04.07.030759},
	year = {2020},
	doi = {10.1101/2020.04.07.030759},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2020/07/10/2020.04.07.030759},
	eprint = {https://www.biorxiv.org/content/early/2020/07/10/2020.04.07.030759.full.pdf},
	journal = {bioRxiv}
}

```
## Acknowledgement
We would like to thank [GISAID](www.gisaid.org) for sharing the sequence data and metadata. We also gratefully acknowledge the authors, originating and submitting laboratories of the sequences from GISAID’s EpiFlu Database on which this research is based. The list is detailed in [here](acknowledgement_table.csv). All submitters of data may be contacted directly via the [GISAID](www.gisaid.org) website.
