#!/usr/bin/python3

'''This main script uses the expression data class and the
statistical analysis class to extract data from the liver cancer 
data csv file and reports the highest expressed genes in 
hepatocellular carcinoma and normal tissue 
'''

import argparse
from expression_data_class import GeneExpressionData
from statistical_analysis_class import StatisticalAnalysis
from barplot import barplot_creator


def main():
    '''The main function is used to pass in arguments from
    the command line to extract data 
    '''
    # initialize the argument parser and adding a usefull description to it
    parser = argparse.ArgumentParser(
        description="Process data from microarray expresison csv files")
    # adding arguments to the parser that can be used to contain information via the commandline
    # its uses are specified in the help
    parser.add_argument("--file-path", required=True,
                        help="Please enter the path to the csv data file")
    parser.add_argument("--gene-name", required=True,
                        help="Please enter the gene name of interest, this returns \
                        all the gene expression values associated with it.")
    parser.add_argument("--gene", required=True,
                        help="Please enter a gene name to calculate its mean \
                        expression across all samples")
    parser.add_argument("--expression-values1", required=True,
                        help="Please enter the first tissue type you want to compare the second to")
    parser.add_argument("--expression-values2", required=True,
                        help="Please enter the second tissue type you want to compare the first to")
    parser.add_argument("--top-n", type=int, required=True,
                        help="Please enter the top number of genes you want to review")
    args = parser.parse_args()

    # initialize the expression data
    expression_data = GeneExpressionData(args.file_path)
    # initialize the statistical analysis class on the expression data that is
    # obtained from the gene expression data
    # class
    analysis = StatisticalAnalysis(expression_data, args.expression_values1,
                                   args.expression_values2, args.top_n)

    # get the expression values of a specific gene name that is parsed through the commandline
    exp_dat = expression_data.get_expression(args.gene_name)

    # calculate the mean expression for a certain gene that is parsed through the commandline
    mean_dat = analysis.calculate_mean_expression(args.gene)

    # performing differential expression by calculating the log2 values
    # of certain tissue types and checking the
    # top n genes, that can be parsed through the commandline


    # performing confidence interval calculations on the parsed through tissue types
    analysis.confidence_interval()

    # creating barplots with the specific tissue types
    # and differential expression of the tissue types
    # that have been parsed through in the commandline
    barplot_creator(analysis)

    print(exp_dat)
    print(mean_dat)


if __name__ == "__main__":
    main()
