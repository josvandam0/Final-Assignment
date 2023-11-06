# Final assignment
This is a project which processes microarray gene expression files and writes a barplot in which the top n expressed genes are shown.
There are 3 classes within this project, that complement each other. 

## Table of contents
- [Packages](#packages)
- [What files to use](#what-files-to-use)
- [Classes](#classes)
- [How to use](#how-to-use)
- [Screenshots](#screenshots)

## Packages
 To use this script you need to have some python libraries installed:

 **```argparse```**
 
 **```matplotlib```**
 
 **```math```**
 
 **```scipy```** 

## What files to use
To use this script you will need a CSV file containing microarray expression data. This script processes the data of these microarary expressions. The gene names have to be in the first row after the second column, the sample names have to be in the first column of each row, the tissue types in the  second columns of each row and the expression data from the 3rd column to the n'th column of each row. 

## Classes
- **```GeneExpressionData```** This class processes microarray expression csv files and puts all necessary information in the so called ```expression_dict``` (dictionary) and the corresponding gene names are stored in a list. With the ```get_expression``` function, the desired gene name needs to be passed through, which is realised via the commandline through parsing. It returns the corresponding gene name with its expression values in f-string format.
 
-  **```StatisticalAnalysis```** This class processes the obtained expression of the csv file further into mean values and subsequently into mean log2 fold change values. It has a couple of functions which require the gene names in a list and the expression values in a ```expression_dict```. With the ```calculate_mean_expression``` function, the desired ```gene``` needs to be passed through via the commandline. With this gene its mean expression values are return in f-string format. In the ```calculate_differential_expression``` function the log2 fold-change values are calculated of the desired ```expression_values1``` and ```expression_values2```, which can be for example ```normal``` or ```HCC``` which are normal tissue type or hepatocellular carcinoma tissue types, respectively. The desired top n of the most expressed up- or downregulated genes, also needs to be passed through as the 3rd argument of this function as ```top_n```.
  
-  **```BarPlotCreator```** With this class the 'Statistical relevance' of samples from microarray expression csv files can be determined. The class uses processed data from the **```StatisticalAnalysis```** class. In the function ```confidence_interval``` the confidence intervals are determined and put in this format ```confidence_interval_dict[gene_name] = fold_change_interval``` for further analysis. The class also contains another function called ```create_barplot``` in which all the obtained and processed data is analysed and pdf files are created which contain barplots and errorbars or up- and downregulated genes (seperate pdf files of up- and downregulated genes are generated).

## How to use
The main script uses a parser where the path of the csv file and the desired top n (up- or down)expressed of one tissue vs another can be specified. The expression_data_class is imported in the main function, and reads the data of the file, and puts the necessary data in a dictionary. 

The following arguments need to be filled with the desired commands:
- **```--file-path```** *In here the path of the required csv file needs to be provided.*

- **```--gene-name```** *To check whether the script works and provides the correct expression data, a specific gene name within the csv file needs to be called on in this argument.*

- **```--gene```** *To visualize the mean expression of a certain gene, this argument needs to be provided with the desired gene name*

- **```--expression-values1```** *The final function of this script is to provide a barplot in which the desired tissue types can be compared. The barplot shows the difference in expression of top n genes, that are over- or underexpressed. In this argument the first tissue type needs to be added.*

- **```--expression-values2```** *In this argument the second tissue type needs to be added (see --expression-values1 for a detailed description of its function)*

- **```--top-n```** *To determine the number of the top genes that are over- and underexpressed, in this argument the number of top genes can be provided. The created barplot displays the top n genes that are over- and underexpressed. It created two different pdf files in which the barplot with corresponding errorbars are displayed.*

## Screenshots
Example output of f-string formatted gene expression values is:

PS \Final assignment> python main_script.py ```--file-path Liver_GSE14520_U133A.csv, --gene-name 209220_at, --gene 209220_at, --exp-val1 normal --exp-val2 HCC, --top-n 10```

The gene name with name: '209220_at' is found and has the following expression values: [10.4940660132875, 10.6636412511969, 11.7427472834902, 11.2649087811095, 9.93381700180943, 10.7557465502368, 11.2494300186965,....]

The gene name with name: '209220_at' is found and has the following mean expression value: 7.427340320636849


Examples of output files of down- and upregulated genes in a microarray experiment with normal vs hepatocellular cancer are:

**Downregulated genes**
![Downreg](https://github.com/josvandam0/Final-Assignment/assets/131524850/30c65f6d-0ea6-49f5-9b62-de93e6f309f4)

**Upregulated genes**
![Upreg](https://github.com/josvandam0/Final-Assignment/assets/131524850/4470c399-9806-45fa-9bfe-cbd09a26f8c3)
