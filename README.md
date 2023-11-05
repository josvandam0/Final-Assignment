# Final assignment
This is a project which processes microarray gene expression files and writes a barplot in which the top n expressed genes are shown.

## Table of contents
- [Packages](#packages)
- [What files to use](#what-files-to-use)
- [How to use](#how-to-use)
- [Screenshots](#screenshots)

## Packages
 To use this script you need to have some python libraries installed:

 **argparse**
 
 **matplotlib**
 
 **math**
 
 **scipy** 

## What files to use
To use this script you will need a CSV file containing microarray expression data. This script processes the data of these microarary expressions. The gene names have to be in the first row after the second column, the sample names have to be in the first column of each row, the tissue types in the  second columns of each row and the expression data from the 3rd column to the n'th column of each row. 

## How to use
The main script uses a parser where the path of the csv file and the desired top n (up- or down)expressed of one tissue vs another can be specified. The expression_data_class is imported in the main function, and reads the data of the file, and puts the necessary data in a dictionary. 

The following arguments need to be filled with the desired commands:
- **--file-path** *In here the path of the required csv file needs to be provided.*

- **--gene-name** *To check whether the script works and provides the correct expression data, a specific gene name within the csv file needs to be called on in this argument.*

- **--gene** *To visualize the mean expression of a certain gene, this argument needs to be provided with the desired gene name*

- **--expression-values1** *The final function of this script is to provide a barplot in which the desired tissue types can be compared. The barplot shows the difference in expression of top n genes, that are over- or underexpressed. In this argument the first tissue type needs to be added.*

- **--expression-values2** *In this argument the second tissue type needs to be added (see --expression-values1 for a detailed description of its function)*

- **--top-n** *To determine the number of the top genes that are over- and underexpressed, in this argument the number of top genes can be provided. The created barplot displays the top n genes that are over- and underexpressed. It created two different pdf files in which the barplot with corresponding errorbars are displayed.*

## Screenshots
Examples of output files of down- and upregulated genes in a microarray experiment with normal vs hepatocellular cancer are:

**Downregulated genes**
![Downreg](https://github.com/josvandam0/Final-Assignment/assets/131524850/30c65f6d-0ea6-49f5-9b62-de93e6f309f4)

**Upregulated genes**
![Upreg](https://github.com/josvandam0/Final-Assignment/assets/131524850/4470c399-9806-45fa-9bfe-cbd09a26f8c3)
