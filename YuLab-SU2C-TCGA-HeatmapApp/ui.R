library(shiny)
library(RColorBrewer)
library(gplots)
library(data.table)
# See above for the definitions of ui and server
pageWithSidebar(
	#App Title
	titlePanel("Yu Lab SU2C/TCGA Heatmap Generator"),
	
	sidebarPanel(
		selectInput("database", h3("Select Database"),
			choices = list("Select Database...", "TCGA-Prostate", "SU2C-Prostate", "SU2C-Benign-vs-PCa", "SU2C PCa-vs-CRPC", "SU2C-CRPC-vs-NEPC"),
			selected = NULL,
			multiple = FALSE
			),
		fileInput("file1", h3("Gene List File")),
		textInput("GOI", h3("Gene of interest to sort groups"), value = NULL, placeholder ="Input gene, i.e. AR"),
		checkboxInput("correlationSort", h5("Sort gene rows by correlation with GOI?"), value = FALSE),
		h3("Heatmap Options:"),
		div(style="display: inline-block;vertical-align:top; width: 150px;",
			selectInput("lowcolor", "Low Color", choices = list("green", "red", "darkred", "blue", "darkblue"), selected = "green", multiple = FALSE)),
		div(style="display: inline-block;vertical-align:top; width: 150px;",
			selectInput("midcolor", "Middle Color", choices = list("black", "white", "yellow"), selected = "black", multiple = FALSE)),
		div(style="display: inline-block;vertical-align:top; width: 150px;",
			selectInput("highcolor", "High Color", choices = list("red", "darkred", "green", "blue", "darkblue"), selected = "red", multiple =FALSE)),
		sliderInput("xfontsize", "Choose Gene Font Size:", min = 0.3, max = 2, value = 0.5),
		sliderInput("range",
			label = "Heatmap Scale",
			min = -10, max = 10, value = c(-2,2)
		),
		h3("Download Data:"),
		downloadButton("downloadHeatmap", "Download Heatmap"),
		downloadButton("downloadGeneHeatmap", "Download Gene Heatmap"),
		downloadButton("downloadDataMatrix", "Download heatmap data matrix for gene set"),
		downloadButton("downloadGOIMatrix", "Download heatmap data matrix for GOI(s)")
	),
	mainPanel(
		tabsetPanel(
		tabPanel("Heatmap", h4("Heatmap (Gene Z-score across all samples)"),plotOutput("plot1"), plotOutput("plot2")),
		tabPanel("Genes", h4("Included Gene List (Ordered)"), tableOutput("table1")),
		tabPanel("Excluded Genes", h4("Genes not found in selected database"), tableOutput("table2"))
		)
		)
	)



