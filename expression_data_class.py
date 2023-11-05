#!/usr/bin/python3

'''This script processes data that is obatined from
a csv file that contains liver cancer data with different
genes and expression values from different samples
'''


class GeneExpressionData:
    '''This class can be used to process expression data from liver
    cancer gene expression files It has a load data module where data from
    the file path that is initialized in the initialize method can be
    read. It has a get gene names method, where gene names are return in 
    a list, from the data that has been processed in the load data method.
    And it has a get expression method where given a gene name the expression
    values across different samples is returned
    '''

    def __init__(self, file_path):
        '''This method initializes the data that is needed to
        initialize the file path in which the liver cancer gene expression 
        data is contained
        '''
        # initialize all the necessary information, like file path, lists
        # dictionaries and initialize the load_data method
        self.file_path = file_path
        self.expression_dict = {}
        self.gene_list = []
        self.gene_dict = {}
        self.load_data()

    def load_data(self):
        '''This method loads the data from liver cancer expression files
        and puts all the neccessary information in a dictionary. 
        input: liver cancer data from a csv file
        output: a dictionary that contains as the key the sample names and as value
        a tuple containing tissue_type and expression data with associated gene names
        '''
        # open the file with the with open function
        with open(self.file_path, encoding = "UTF-8") as file:
            # store the header
            header = next(file).split(',')
            # store the gene names in a variable
            self.gene_list = header[2:]
            # initialize empty dictionary
            self.expression_dict = {}
            # iterate through the lines of the file
            for line in file:
                # initialize a variable where the lines of the files are stored in
                read_lines = line.split(',')
                # initialize a variable where the sample names are stored in
                sample_name = read_lines[0]
                # initialize a variable that contains the types of tissue
                tissue_type = read_lines[1]
                # initialize a variable that uses a list comprehension to store
                # the expression data as floats
                expression_data = [float(expr) for expr in read_lines[2:]]
                # fill the expression dictionary with the appropriate values, I chose a
                # dictionary because I found that more clear then a list or a seperate tuple
                self.expression_dict[sample_name] = (
                    tissue_type, expression_data)

    def get_gene_names(self):
        '''This function returns the gene names of the csv file
        that are initiatlized in the load data function
        '''
        return self.gene_list

    def get_expression(self, gene_name):
        '''This function takes in the gene name list and the
        expression values and processes them to return the
        gene name and its associated expression values
        input: gene name and expression data
        output: per gene name the associated expression values
        '''
        # initialize an empty dictionary to store the gene names
        # and their associated expression values in
        # here I chose a dictionary because it would be
        # more handy later on to evaulate the values per gene name
        self.gene_dict = {gene: [] for gene in self.gene_list}
        # loop through the gene names, loop through
        # the expression data and append the expression values to
        # the dictionary with gene names as the keys
        for _, (_, expression_data) in self.expression_dict.items():
            for gene, expression in zip(self.gene_list, expression_data):
                self.gene_dict[gene].append(expression)
        # check whether a gene name has been parsed into the commandline
        if gene_name is not None:
            # try and except KeyError, to see whether the gene name
            # that has been parsed through is correct, or needs to be corrected
            try:
                return f"The gene name with name: '{gene_name}' is found and \
                has the following expression values: {self.gene_dict[gene_name]}"
            except KeyError:
                return f"The gene name with name: '{gene_name}' is not found in this file. \
                Please enter a different gene name."
        else:
            return 0
