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

### Heatmap apps, how to use:

1.	If you do not have R installed, download and install from:
    a.	https://cran.r-project.org/bin/windows/base/rpatched.html
    b.	R install doesn’t require any special permissions, so you should be able to install
2.	Open R
3.	Open the folder:
    a.	R:\Medicine\Hematology-Oncology\Yu_Lab\Nate\scripts_and_tools\su2c_script
4.	There are two boot files in the folder, one for creating heatmaps with SU2C data and one with TCGA data
    a.	BootSU2C.R
    b.	BootTCGA.R
5.	Drag the .R file for the dataset you want to load
6.	The Boot script will test to make sure you have everything installed that is needed, if this is your first time running it may take a couple minutes to set up.
7.	It will then read in the data file needed to make heatmaps, it will show you the progress.
8.	The app will then open in a web browser:
9.	Upload your geneset to generate a heatmap.
    a.	Test set: sharma-neal_etal2013_cancercell.txt
        i.	R:\Medicine\Hematology-Oncology\Yu_Lab\Nate\scripts_and_tools\su2c_script\ sharma-neal_etal20113_cancercell.txt
10.	If you wish to sort by a single gene of interest type it in.
    a.	This will also provide a single heatmap for the gene of interest
11.	Click Download Heatmap to save a high resolution version.
