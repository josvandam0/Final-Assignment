#!/usr/bin/python3

'''This script handles fold change results of one gene vs another gene, 
with its corresponding confidence intervals and writes a pdf file with barplots
in which genes that are potentially interesting to look at are displayed (up- and downregulated)
'''
import math
import scipy.stats as st
import matplotlib.pyplot as plt

class BarPlotCreator:
    '''This class processes data obtained from microarray, which are passed through
    via analysis in the class. With this class, confidence intervals can be calculated
    by using processed data in which the mean and log2 values are passed through.
    It is able to create barplots with its create_barplot function.
    Processed data has to be a list of tuples with as first element the gene name and as
    a second element the corresponding log2 fold change value. 
    '''
    def __init__(self, analysis, expression_values1, expression_values2, top_n):
        self.analysis = analysis
        self.expression_values1 = expression_values1
        self.expression_values2 = expression_values2
        self.top_n = top_n
        self.label_size = 20
        self.confidence_interval_dict = {}

    def confidence_interval(self):
        '''This function calculates the confidence intervals 
        of the log2 expressions of given tissue types
        input: mean of expression values and expression values
        output: confidence intervals of gene expressions of tissue types for error bars in bar plot
        '''
        mean_log2_dict1 = {gene_name: math.log2(
            expression) for gene_name, expression in self.analysis.gene_dict_type1_avg.items()}
        mean_log2_dict2 = {gene_name: math.log2(
            expression) for gene_name, expression in self.analysis.gene_dict_type2_avg.items()}

        # loop through the gene list
        for gene_name in self.analysis.gene_list:
            log2_mean1 = mean_log2_dict1[gene_name]
            log2_mean2 = mean_log2_dict2[gene_name]
            # calculate the confidence interval of the different gene types by using the
            # scipy module, using an alpha of 0.95,
            # and using the numpy module to calculate the standard error
            # mean value
            confidence95_mean1 = st.t.interval(alpha=0.95, df=len(
                self.analysis.gene_dict_type1[gene_name]) - 1, loc=log2_mean1, scale=st.sem(
                    self.analysis.gene_dict_type1[gene_name]))
            confidence95_mean2 = st.t.interval(alpha=0.95, df=len(
                self.analysis.gene_dict_type2[gene_name]) - 1, loc=log2_mean2, scale=st.sem(
                    self.analysis.gene_dict_type2[gene_name]))
            # calculate the interval of the fold change
            fold_change_interval = (
                confidence95_mean1[0]-confidence95_mean2[0],
                confidence95_mean1[1] - confidence95_mean2[1])
            # fill the confidence interval dictionary with
            # the gene name as the keys and the confidence intervals
            # of the fold changes
            self.confidence_interval_dict[gene_name] = fold_change_interval


    def create_barplot(self):
        '''This function creates histograms with gene expression data
        input: gene expression data with fold changes
        output: histogram with fold changes plotted

        '''
        # extract the gene names and the log2 fold changes in lists for further analysis by plotting
        gene_name_upreg, log2_fold_change_upreg = zip(*self.analysis.top_n_upregulated_genes)
        gene_name_downreg, log2_fold_change_downreg = zip(*self.analysis.top_n_downregulated_genes)
        # initialize empty lists to store the confidence intervals
        # of the upregulated and downregulated genes in
        upreg_interval = []
        downreg_interval = []
        # loop through the gene name upregulated list and fill with the upper confidence intervals
        # from the confidence interval dictionary
        # this is done for the upregulated and downregulated intervals seperately by
        # making use of the fact that in both dictionaries the genes are differently ordered
        for gene in gene_name_upreg:
            _, upper = self.confidence_interval_dict[gene]
            upreg_interval.append(upper)
        for gene in gene_name_downreg:
            _, upper = self.confidence_interval_dict[gene]
            downreg_interval.append(upper)

        upreg_ci95 = [abs(ci95 - value) for ci95,
                        value in zip(upreg_interval, log2_fold_change_upreg)]
        downreg_ci95 = [abs(ci95 - value) for ci95,
                        value in zip(downreg_interval, log2_fold_change_downreg)]

        # creating a plot with a certain figuresize by using the imported module
        # called matplotlib for the upregulated genes
        plot1 = plt.figure(figsize=(50, 20))

        # create bars for the plot by passing through the gene names of the upregulated
        # genes on the x axis and the log2 fold change of the upregulated genes on the y
        # axis. Error bars are added by yerr= and as the values of the error bars
        # the confidence intervals of the upregulated genes are being passed through
        plt.bar(gene_name_upreg, log2_fold_change_upreg,
                color='lightgray', width=0.6, yerr=upreg_ci95, capsize=10)
        # add an x label with the name genes names
        plt.xlabel("Gene names", size=self.label_size)
        # add a y label with the name log2 fold changes
        plt.ylabel("Log2 Fold Changes", size=self.label_size)
        # add a title to the plot with the f string formatted
        # passed through arguments from the command line
        plt.title(f"Fold changes of top {self.top_n} upregulated"
                f"{self.expression_values1} vs {self.expression_values2}", size=self.label_size)
        # save the plot as pdf with the f string formatted passed through arguments as the names
        plot1.savefig(
            f"top-{self.top_n}-upregulated"
            f"-{self.expression_values1}-vs-{self.expression_values2}.pdf")
        # creating a plot with a certain figuresize by using the imported module
        # called matplotlib for the downregulated genes
        plot2 = plt.figure(figsize=(50, 20))
        # create bars for the plot by passing through the gene names of the downregulated
        # genes on the x axis and the log2 fold change of the downregulated genes on the y
        # axis. Error bars are added by yerr= and as the values of the error bars
        # the confidence intervals of the downregulated genes are being passed through
        plt.bar(gene_name_downreg, log2_fold_change_downreg,
                color='lightgray', width=0.6, yerr=downreg_ci95, capsize=10)
        # add a y label with the name log2 fold changes
        plt.xlabel("Gene names", size=self.label_size)
        # add a y label with the name log2 fold changes
        plt.ylabel("Log2 Fold Changes", size=self.label_size)
        # add a title to the plot with the f string formatted
        # passed through arguments from the command line
        plt.title(f"Fold changes top {self.top_n} downregulated {self.expression_values1}"
                f"vs {self.expression_values2}", size=self.label_size)
        # save the plot as pdf with the f string formatted passed through arguments as the names
        plot2.savefig(
            f"top-{self.top_n}-downregulated"
            f"-{self.expression_values1}-vs-{self.expression_values2}.pdf")
