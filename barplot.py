import matplotlib.pyplot as plt

def barplot_creator(analysis):
    '''This function creates histograms with gene expression data
    input: gene expression data with fold changes
    output: histogram with fold changes plotted

    '''
    label_size = 20
    # extract the gene names and the log2 fold changes in lists for further analysis by plotting
    gene_name_upreg = [gene[0] for gene in analysis.top_n_upregulated_genes]
    log2_fold_change_upreg = [value[1]
                                for value in analysis.top_n_upregulated_genes]
    gene_name_downreg = [gene[0]
                            for gene in analysis.top_n_downregulated_genes]
    log2_fold_change_downreg = [gene[1]
                                for gene in analysis.top_n_downregulated_genes]
    # initialize empty lists to store the confidence intervals
    # of the upregulated and downregulated genes in
    upreg_interval = []
    downreg_interval = []
    # loop through the gene name upregulated list and fill with the upper confidence intervals
    # from the confidence interval dictionary
    # this is done for the upregulated and downregulated intervals seperately by
    # making use of the fact that in both dictionaries the genes are differently ordered
    for gene in gene_name_upreg:
        _, upper = analysis.confidence_interval_dict[gene]
        upreg_interval.append(upper)
    for gene in gene_name_downreg:
        _, upper = analysis.confidence_interval_dict[gene]
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
    plt.xlabel("Gene names", size=label_size)
    # add a y label with the name log2 fold changes
    plt.ylabel("Log2 Fold Changes", size=label_size)
    # add a title to the plot with the f string formatted
    # passed through arguments from the command line
    plt.title(f"Fold changes of top {analysis.top_n} upregulated genes {analysis.expression_values1} vs {analysis.expression_values2}", size=label_size)
    # save the plot as pdf with the f string formatted passed through arguments as the names
    plot1.savefig(
        f"top-{analysis.top_n}-upregulated_genes-{analysis.expression_values1}-vs-{analysis.expression_values2}.png")
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
    plt.xlabel("Gene names", size=label_size)
    # add a y label with the name log2 fold changes
    plt.ylabel("Log2 Fold Changes", size=label_size)
    # add a title to the plot with the f string formatted
    # passed through arguments from the command line
    plt.title(
        f"Fold changes of top {analysis.top_n} downregulated genes {analysis.expression_values1} vs {analysis.expression_values2}", size=label_size)
    # save the plot as pdf with the f string formatted passed through arguments as the names
    plot2.savefig(
        f"top-{analysis.top_n}-downregulated_genes-{analysis.expression_values1}-vs-{analysis.expression_values2}.png")