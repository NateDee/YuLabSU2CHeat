library(shiny)
library(RColorBrewer)
library(gplots)
library(data.table)
# See above for the definitions of ui and server
pageWithSidebar(
	#App Title
	titlePanel("Yu Lab SU2C Heatmap Generator"),
	
	sidebarPanel(
		fileInput("file1", h3("Gene List File")),
		sliderInput("range",
			label = "Heatmap Scale",
			min = -10, max = 10, value = c(-2,2)),
		textInput("GOI", h4("Gene of interest to sort groups"), value = NULL, placeholder ="Input gene, i.e. AR"),
		downloadButton("downloadHeatmap", "Download Heatmap"),
		downloadButton("downloadGeneHeatmap", "Download Gene Heatmap")
	),
	mainPanel(
		tabsetPanel(
		tabPanel("Heatmap", h4("Heatmap (Gene Z-score across all samples)"),plotOutput("plot1"), plotOutput("plot2")),
		tabPanel("Genes", h4("Included Gene List (Ordered)"), tableOutput("table1")),
		tabPanel("Excluded Genes", h4("Included Gene List (Ordered)"), tableOutput("table2"))
		)
		)
	)



