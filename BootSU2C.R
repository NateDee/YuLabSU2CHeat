#Set working directory to SU2C_script
setwd("R:/Medicine/Hematology-Oncology/Yu_Lab/Nate/scripts_and_tools/su2c_script")

#Create function to test for packages that are needed to run app
installPkges <- function(pkg){
	if(!pkg %in% installed.packages()) install.packages(pkg, repos = "https://mirror.las.iastate.edu/CRAN/")
	}

#Required packages, "data.table", "gplots", "RColorBrewer", "shiny"
print("Checking for required packages, installing if needed")

installPkges("data.table")
installPkges("gplots")
installPkges("RColorBrewer")
installPkges("shiny")

#Load shiny
library(shiny)

#Load SU2C app
runApp("su2cHeatmapApp")

#or

#runApp("tcgaHeatmapApp")