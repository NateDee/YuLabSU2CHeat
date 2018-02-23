# YuLabSU2CHeat
Generating heatmaps for gene sets of interest in Prostate Cancer datasets (SU2C) in the Yu Lab at Northwestern University

This app is built on the Shiny framework in R, it will open in a web browser for user interaction.


## Purpose of Heatmap App:
In RNA-seq and microarray experiments we often get a list of differentially regulated genes.  You can use this heatmap generator to identify if your gene set associates with prostate cancer progression in SU2C and TCGA datasets.  All you need is a list of genes you’re interested in.

### Setup gene list file:

You must prepare your gene list to be accepted by the app.  
•	Gene symbols in the first column, with a header row (i.e. GENE in the first row)
  o	Currently only accepts gene symbol as gene identifiers (i.e. AR, FOXA1; will not work with refseq IDs, i.e. NM_12345)
•	The file must have at least two columns, if you only have a list of genes, just duplicate them in column 2 and name the header row GENE2.
•	Using excel, save your file as a tab-delimited .txt file

### Heatmap app, how to use:

1.    If you do not have R installed, download and install from:
    a.	https://cran.r-project.org/bin/windows/base/rpatched.html
    b.	R install doesn’t require any special permissions, so you should be able to install
2.    Open R
3.    Open the folder:
    a.	R:\Medicine\Hematology-Oncology\Yu_Lab\Nate\scripts_and_tools\YuLab-SU2C-TCGA
4.    Drag the BootSU2C-TCGA-HeatmapApp.R file into your R session.
    a.	The Boot script will test to make sure you have everything installed that is needed, if this is your first time running it may take a couple minutes to set up.
5.    The app will then open in a web browser:
6.    Select the Database you would like to analyse (Either SU2C or TCGA)
    a.	SU2C contains Benign, Primary PCa, CRPC, and CRPC-NE tumors
    b.	TCGA contains Benign, and a large dataset of Primary PCa
7.    Upload your geneset by clicking “Browse” to generate a heatmap.
8.    Test set: sharma-neal_etal2013_cancercell.txt
    a.	R:\Medicine\Hematology-Oncology\Yu_Lab\Nate\scripts_and_tools\ YuLab-SU2C-TCGA \ sharma-neal_etal20113_cancercell.txt
9.    If you wish to sort by a single gene of interest type it in.
    a.	This will also provide a single heatmap for the gene of interest
10.    Click Download Heatmap to save a high resolution version.
    a.	You can also download a heatmap for your gene of interest if you selected one
11.    The “Genes” tab will tell you the order of your genes in the heatmap
    a.	If you select the checkbox for “Sort gene rows by correlation with GOI?” – The genes tab will display Pearson R correlation values for your gene of interest vs. each gene
12.    If any genes in your gene list are not found in the database, they will be included in the “Excluded Genes” tab.

