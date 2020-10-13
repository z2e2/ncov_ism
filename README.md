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
conda create -n ISM python=3.5.4 biopython=1.71 pandas=0.20.3 matplotlib=2.1.1 scipy=1.0.0 scikit-learn=0.19.1
source activate ISM
git clone https://github.com/z2e2/ncov_ism.git
cd ncov_ism
python setup.py install
```
## Running default ncov_ism
#### Input Data
SARS-CoV-2 (hCoV-19) sequences and metadata are made available at [GISAID](www.gisaid.org). Our pipeline takes the multiple sequence alignment results and metadata as input which are both available at [GISAID](www.gisaid.org). The multiple sequence alignment results of sequences are required. 
#### Reference genomes
The default reference sequences is [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2/) in NCBI. It also corresponds to *EPI_ISL_402125* in [GISAID](www.gisaid.org) hCoV-19 database.
#### Quick start
Once *ncov_ism* is installed, you can run the default setting using the following command with minimum number of parameters:
```
ncov_ism default -i <multiple sequence alignment file> -m <metadata file> -o <output directory>
```
## Running customized ncov_ism
#### Input Data
Viral sequences and metadata that are similar to the format posted in [GISAID](www.gisaid.org) should be supplied. Our pipeline takes the multiple sequence alignment results and metadata as primary input. The multiple sequence alignment results of sequences are required. The metadata should be a Tab Separated Values (TSV) file with a column designated for sequence identification number (e.g., *gisaid_epi_isl* in GISAID metadata), a column designated for geographical location that the sequence was collected from (e.g., *country* and *region* in GISAID metadata) and a column called *date* designated for collection date of the virus genomes (with format "%Y-%m-%d", e.g., 2020-09-24). For the sequence file, the header should be formated the same as the GISAID style. For example:
```
>GISAID_Strain_Name|GISAID_EPI_ISL|DATE
ACGT.....
```
Note *GISAID_Strain_Name* and *DATE* fields are not important but rather a placeholder since we can get the same information from the metadata file. But *GISAID_EPI_ISL* is the sequence Identification numbedr used to cross-reference between differnet files. The pattern above is mandatory so that the tool can parse the file correctly.
#### Other required arguments
A file named *regions_to_visualize.txt* is needed for the tool to determine regions of interest to visualize. For example, if the user is interested in four coutries, namely, USA, Canada, United Kingdom, Iceland. They should put the following content in *regions_to_visualize.txt*. 
```
USA, Canada, United Kingdom, Iceland
```
#### Reference genomes
The default reference sequence is [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2/) in NCBI. It also corresponds to *EPI_ISL_402125* in [GISAID](www.gisaid.org) hCoV-19 database. You can supply with your own reference genome. The corresponding genbank file is also required (e.g., [covid-19-genbank.gb](data/covid-19-genbank.gb) is the genbank file for our default reference gennome).
#### Command
Once *ncov_ism* is installed and required data are ready, you can run the customized ISM analysis with your data using the following command with minimum number of parameters:
```
ncov_ism customized -i <multiple sequence alignment file> -m <metadata file> -o <output directory> -id <column name for accession numbers> -loc <column name for locations> -e <entropy threshold>
```
Here is an simple example, assuming the multiple sequence alignment file is *msa.fasta*, the metadata file is *metadata.tsv*, the output folder is *results*, the accession number column is *gisaid_epi_isl*, the location column is *country* and the entropy threshold is *0.2*:
```
ncov_ism customized -i msa.fasta -m metadata.tsv -o results -id gisaid_epi_isl -loc country -e 0.2
```
## About
The novel coronavirus responsible for COVID-19, SARS-CoV-2, expanded to reportedly 8.7 million confirmed cases worldwide by June 21, 2020. The global SARS-CoV-2 pandemic highlights the importance of tracking viral transmission dynamics in real-time. Through June 2020, researchers have obtained genetic sequences of SARS-CoV-2 from over 50 thousand samples from infected individuals worldwide. Since the virus readily mutates, each sequence of an infected individual contains useful information linked to the individual's exposure location and sample date. But, there are over 30,000 bases in the full SARS-CoV-2 genome, so tracking genetic variants on a whole-sequence basis becomes unwieldy. *ncov_ism* is a method to instead efficiently identify and label genetic variants, or "subtypes" of SARS-CoV-2. This method defines a compact set of nucleotide sites that characterize the most variable (and thus most informative) positions in the viral genomes sequenced from different individuals, called an Informative Subtype Marker or *ISM*. This tool defines viral subtypes for each ISM, and analyze the regional distribution of subtypes to track the progress of the pandemic.

## Reference
Research article in PLOS COMPUTATIONAL BIOLOGY [Genetic Grouping of SARS-CoV-2 Coronavirus Sequences using Informative Subtype Markers for Pandemic Spread Visualization](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008269).
If you find our work helpful, please cite:
```
@article{10.1371/journal.pcbi.1008269,
    author = {Zhao, Zhengqiao AND Sokhansanj, Bahrad A. AND Malhotra, Charvi AND Zheng, Kitty AND Rosen, Gail L.},
    journal = {PLOS Computational Biology},
    publisher = {Public Library of Science},
    title = {Genetic grouping of SARS-CoV-2 coronavirus sequences using informative subtype markers for pandemic spread visualization},
    year = {2020},
    month = {09},
    volume = {16},
    url = {https://doi.org/10.1371/journal.pcbi.1008269},
    pages = {1-32},
    number = {9},
    doi = {10.1371/journal.pcbi.1008269}
}
```
## Acknowledgement
We would like to thank [GISAID](www.gisaid.org) for sharing the sequence data and metadata. We also gratefully acknowledge the authors, originating and submitting laboratories of the sequences from GISAID’s EpiFlu Database on which this research is based. The list is detailed in [here](results/acknowledgement_table.txt). All submitters of data may be contacted directly via the [GISAID](www.gisaid.org) website.
