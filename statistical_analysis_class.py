#!/usr/bin/python3

'''This module processes from the GeneExpressionData class.
It calculates the mean expression of a specific gene across all samples
and it calculates the differential expression between two sets of gene expression values
(normal vs HCC) and returns the result
'''
import math
import scipy.stats as st
from expression_data_class import GeneExpressionData


class StatisticalAnalysis(GeneExpressionData):
    '''With this statistical analysis class the mean expression 
    of a specific gene can be calculated, 
    the differential expression between two sets of gene expression values can be calculated and 
    returns the top 10 genes expressed between two sets in a written report (pdf)
    '''

    def __init__(self, gene_expression_data):
        '''This function initializes gene expression data and stores
        it in a variable
        '''
        # inherit the gene expression data class by using the super() method
        super().__init__(gene_expression_data.file_path)
        # initialize all necessary lists, dictionaries and values that are being used
        # in this script
        self.get_expression(None)
        self.top_n_upregulated_genes = []
        self.top_n_downregulated_genes = []
        self.gene_dict_type1_avg = {}
        self.gene_dict_type2_avg = {}
        self.gene_dict_type1 = {}
        self.gene_dict_type2 = {}
        self.confidence_interval_dict = {}

    def calculate_mean_expression(self, gene):
        '''This function calculates the mean expression of a specific gene
        across all the samples
        '''
        try:
            return f"The gene name with name: '{gene}' is found and has the following \
                mean expression value: {sum(self.gene_dict[gene])/len(self.gene_dict[gene])}"
        except KeyError:
            return f"The gene name with name: '{gene}' is not found in this file.\
                Please enter a different gene name."

    def calculate_differential_expression(self, expression_values1, expression_values2, top_n):
        '''This function calculates the differential expression between two sets of gene expression
        values (normal vs HCC) and returns the result
        '''
        # initialize an empty dictionary to store the gene names
        # and their associated expression values in
        # I chose for dictionaries because that is more clear then a 
        # list or a tuple for further analysis
        self.gene_dict_type1 = {gene_name: [] for gene_name in self.gene_list}
        self.gene_dict_type2 = {gene_name: [] for gene_name in self.gene_list}
        # initialize an empty dictionary to store the differential values in
        differential_dict = {}
        # initialize the parsed through expression values,
        # type 1 and 2, type 1 for the one condition parsed through
        # and type 2 for the other condition that is parsed through
        # to for example compare cancerous tissue to normal tissue
        self.expression_values1 = expression_values1
        self.expression_values2 = expression_values2
        # initialize the top_n variable, which is parsed through the commandline
        # to determine the number of genes that are being shown
        self.top_n = top_n
        # loop through the tissue type and expressoin values
        for _, (tissue_type, expression_data) in self.expression_dict.items():
            # loop through the gene list and expression data
            for gene_name, expression in zip(self.gene_list, expression_data):
                # determine whether the tissue type in the expression dictionary
                # corresponds to the tissue type that has been parsed through the
                # commandline
                if tissue_type == expression_values1:
                    # seperate the one tissue type from the other by using an if
                    # statement followed by appending the expression value
                    # that belongs to the tissue type to the gene dict type(n) list
                    self.gene_dict_type1[gene_name].append(expression)
                elif tissue_type == expression_values2:
                    self.gene_dict_type2[gene_name].append(expression)

        # loop through the gene list
        for gene_name in self.gene_list:
            # calculate the mean values of the genes by dividing the sum of the values
            # by the length of the values (using the values of the gene list as keys)
            self.gene_dict_type1_avg[gene_name] = sum(
                self.gene_dict_type1[gene_name])/len(self.gene_dict_type1[gene_name])
            self.gene_dict_type2_avg[gene_name] = sum(
                self.gene_dict_type2[gene_name])/len(self.gene_dict_type2[gene_name])
            # calculate the log2 fold change by using the math.log2 function of the module math
            # that has been imported previously. The log2 fold change is calculated by
            # taking the difference of the log2(condition1)-log2(condition2)
            differential_dict[gene_name] = (math.log2(
                self.gene_dict_type1_avg[gene_name])-(
                    math.log2(self.gene_dict_type2_avg[gene_name])))
        # for further analysis of the top n upregulated genes, sort the genes by fold change from
        # high to low by using the sorted method and using reverse=True to sort from high to low
        self.top_n_upregulated_genes = sorted(
            differential_dict.items(), key=lambda item: item[1], reverse=True)[:top_n]
        # using the same method to determine the downregulated genes, but now use the reverse=False
        # because the genes need to be sorted from low to high for this
        self.top_n_downregulated_genes = sorted(
            differential_dict.items(), key=lambda item: item[1], reverse=False)[:top_n]

    def confidence_interval(self):
        '''This function calculates the confidence intervals 
        of the log2 expressions of given tissue types
        input: mean of expression values and expression values
        output: confidence intervals of gene expressions of tissue types for error bars in bar plot
        '''
        mean_log2_dict1 = {gene_name: math.log2(
            expression) for gene_name, expression in self.gene_dict_type1_avg.items()}
        mean_log2_dict2 = {gene_name: math.log2(
            expression) for gene_name, expression in self.gene_dict_type2_avg.items()}

        # loop through the gene list
        for gene_name in self.gene_list:
            log2_mean1 = mean_log2_dict1[gene_name]
            log2_mean2 = mean_log2_dict2[gene_name]
            # calculate the confidence interval of the different gene types by using the
            # scipy module, using an alpha of 0.95,
            # and using the numpy module to calculate the standard error
            # mean value
            confidence95_mean1 = st.t.interval(alpha=0.95, df=len(
                self.gene_dict_type1[gene_name]) - 1, loc=log2_mean1, scale=st.sem(
                    self.gene_dict_type1[gene_name]))
            confidence95_mean2 = st.t.interval(alpha=0.95, df=len(
                self.gene_dict_type2[gene_name]) - 1, loc=log2_mean2, scale=st.sem(
                    self.gene_dict_type2[gene_name]))
            # calculate the interval of the fold change
            fold_change_interval = (
                confidence95_mean1[0]-confidence95_mean2[0],
                confidence95_mean1[1] - confidence95_mean2[1])
            # fill the confidence interval dictionary with
            # the gene name as the keys and the confidence intervals
            # of the fold changes
            self.confidence_interval_dict[gene_name] = fold_change_interval
