# ISM
Informative Subtype Markers (ISM) is an efficient framework for genetic subtyping of a pandemic virus and implement it for SARS-CoV-2, the novel coronavirus that causes COVID-19.        
Drexel University EESI Lab, 2020        
Maintainer: Zhengqiao Zhao, zz374 at drexel dot edu  
Owner: Gail Rosen, gailr at ece dot drexel dot edu  

## Abstract
We utilize an entropy-based analysis to identify mutational signatures of SARS-CoV-2 sequences in the GISAID database available as of April 5, 2020. Our  subtyping method identifies nucleotide sites within the viral genome which are highly informative of variation between the viral genomes sequenced in different individuals. These sites are used to characterize individual virus sequence with a characteristic Informative Subtype Marker (ISM). The ISMs provide signatures that can be efficiently and rapidly utilized to quantitatively trace viral dynamics through geography and time. We show that by analyzing the ISM of currently available SARS-CoV-2 sequences, we are able to profile international and interregional differences in viral subtype, and visualize the emergence of viral subtypes in different countries over time. To validate and demonstrate the utility of ISM-based subtyping: (1) We show the distinct genetic subtypes of European infections, in which early on infections are related to the viral subtypes that has become dominant in Italy followed by the development of local subtypes, (2) We distinguish subtypes associated with outbreaks in distinct parts of the United States, identify the development of a local subtype potentially due to community to transmission and distinguish it from the predominant subtype in New York, suggesting that the outbreak in New York is linked to imported cases from Europe. (3) We present results that quantitatively show the temporal behavior of the emergence of SARS-CoV-2 from localization in China to a pattern of distinct regional subtypes as the virus spreads throughout the world over time.
## Entropy of viral sequences
![Fig 0](results/1_overall_entropy.pdf "Entropy analysis")
## ISM distribution worldwide
![Fig 1](results/6_regional_ISM.pdf "Subtype composition in different locations worldwide")
## ISM distribution in US
![Fig 2](results/7_intra-US_ISM.pdf "Subtype composition in different locations in US")
## The dynamic of ISM in different locations
![Fig 3](results/8_ISM_growth_USA.pdf "the dynamic subtype composition in different locations over time")
![Fig 4](results/8_ISM_growth_Japan.pdf "the dynamic subtype composition in different locations over time")
![Fig 5](results/8_ISM_growth_Australia.pdf "the dynamic subtype composition in different locations over time")