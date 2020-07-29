import argparse
import logging

from ._loaddata import load_data
from ._pickism import entropy_analysis, pick_ISM_spots, annotate_ISM, ISM_disambiguation
from ._analyzeism import ISM_analysis
from ._visualization import ISM_visualization, ISM_plot
from .build_ism import build_ISM

logging.basicConfig(format='%(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)

def main():
    
    parser = argparse.ArgumentParser(description='Informative Subtype Marker (ISM) is an efficient framework for genetic subtyping of a pandemic virus and implement it for SARS-CoV-2, the novel coronavirus that causes COVID-19.', prog='main')
    subparsers = parser.add_subparsers(title='subcommands',
                                       description='the following subcommands \
                                    are available: build, analyze, visualize, default', 
                                    dest='subparser_name')
    build_parser = subparsers.add_parser("build")
    build_parser.add_argument("-i", metavar="input",
                        help="Multiple sequence alignment results in fasta format",
                        required=True, type=str)
    build_parser.add_argument("-m", metavar="metadata",
                        help="Metadata in tsv format",
                        required=True, type=str)
    build_parser.add_argument("-gb", metavar="genbank",
                        help="genbank file for reference sequence",
                        required=True, type=str)
    build_parser.add_argument("-o", metavar="output",
                        help="Path to the output directory",
                        required=True, type=str)
    build_parser.add_argument("-r", metavar="referenceId", help="accession number of the reference sequence", 
                        default='EPI_ISL_402125', required=False, type=str)
    build_parser.add_argument("-e", metavar="entropy", help="entropy threshold", 
                        default=0.2, required=False, type=float)
    build_parser.add_argument("-n", metavar="nullFreq", help="null frequency threshold", 
                        default=0.25, required=False, type=float)
    default_parser = subparsers.add_parser("default")
    default_parser.add_argument("-i", metavar="input",
                        help="Multiple sequence alignment results in fasta format",
                        required=True, type=str)
    default_parser.add_argument("-m", metavar="metadata",
                        help="Metadata in tsv format",
                        required=True, type=str)
    default_parser.add_argument("-o", metavar="output",
                        help="Path to the output directory",
                        required=True, type=str)
    args = parser.parse_args()
    if args.subparser_name == "build":
        MSA_FILE_NAME = args.i
        META_FILE_NAME = args.m
        reference_genbank_name = args.gb
        OUTPUT_FOLDER = args.o
        REFERENCE_ID = args.r
        en_thres = args.e
        null_thres = args.n
        ISM_df = build_ISM(MSA_FILE_NAME, META_FILE_NAME, reference_genbank_name, OUTPUT_FOLDER, 
                   REFERENCE_ID, en_thres, null_thres)
    elif args.subparser_name == "default":
        MSA_FILE_NAME = args.i
        META_FILE_NAME = args.m
        OUTPUT_FOLDER = args.o
        reference_genbank_name = 'data/covid-19-genbank.gb'
        REFERENCE_ID = 'EPI_ISL_402125'
        en_thres = 0.2
        null_thres = 0.25
        ISM_df = build_ISM(MSA_FILE_NAME, META_FILE_NAME, reference_genbank_name, OUTPUT_FOLDER, 
                   REFERENCE_ID, en_thres, null_thres)
        region_raw_count, state_raw_count, count_dict = ISM_analysis(ISM_df, OUTPUT_FOLDER)
        
        region_list = ['Mainland China', 'Japan', 'Singapore', 'Hong Kong', 'India',
               'Australia', 'New Zealand', 'Brazil', 'USA', 'Canada', 
               'United Kingdom', 'Iceland', 'Belgium', 'Netherlands', 'Denmark',
               'France', 'Italy', 'Spain', 'Germany', 'Russia',
               ]

        state_list = [
                        'Washington','Oregon','California','Alaska','Idaho',
                        'Michigan','Wisconsin','Minnesota','Illinois','New Mexico', 
                        'Nebraska','Wyoming','Utah','Arizona','Texas',
                        'Massachusetts','Connecticut','New York','New Jersey','Pennsylvania',
                        'Maryland','Washington DC','Virginia','Florida','Louisiana',
                     ]

        time_series_region_list = ['Mainland China', 'Japan', 'USA', 'France', 'Denmark', 
                                   'United Kingdom', 'Netherlands', 'Australia', 'Canada', 'Spain']

        ISM_set, region_pie_chart, state_pie_chart, count_list, date_list = ISM_visualization(region_raw_count, state_raw_count, 
                                                                                              count_dict, region_list, 
                                                                                              state_list,
                                                                                              time_series_region_list,
                                                                                              OUTPUT_FOLDER,
                                                                                              ISM_FILTER_THRESHOLD=0.05,
                                                                                              ISM_TIME_SERIES_FILTER_THRESHOLD=0.025)
        REFERENCE_date = ISM_df[ISM_df['gisaid_epi_isl'] == REFERENCE_ID]['date'].min().date()
        ISM_plot(ISM_df, region_list, region_pie_chart, state_list, state_pie_chart, REFERENCE_date, time_series_region_list, count_list, date_list, OUTPUT_FOLDER)
    
    
if __name__ == "__main__":
    main()
