library(shiny)
library(RColorBrewer)
library(gplots)
library(data.table)
library(DT)
# See above for the definitions of ui and server
pageWithSidebar(
	#App Title
	titlePanel("Yu Lab Human Prostate Heatmap Generator"),
	
	sidebarPanel(
		selectInput("database", h3("Select Database"),
			choices = list("Select Database...", "TCGA-Prostate", "DBGAP-Prostate", "DBGAP-Benign-vs-PCa", "DBGAP-PCa-vs-CRPC", "DBGAP-CRPC-vs-NEPC", "GLINSKY_U133A", "VARAMBALLY_GSE3325", "LAPOINTE_GSE3933", "TOMLINS_GSE6099", "YU_GSE6919", "SBONER_GSE16560", "TAYLOR_GSE21034", "LONG_GSE26367", "CAI_GSE32269", "GRASSO_GSE35988", "ROSS_GSE70770", "VANAJA_U133A"),
			selected = NULL,
			multiple = FALSE
			),
		fileInput("file1", h3("Gene List File")),
		h4("MyData Heatmap Options:"),
		div(style="display: inline-block;vertical-align:top; width: 30%;",
			selectInput("doUserHeatmap", "Heatmap from input?", choices = list("Yes", "No"), selected = "No", multiple = FALSE)),
		div(style="display: inline-block;vertical-align:top; width: 30%;",
			selectInput("normalize", "Normalization", choices = list("Z-Score", "Raw", "Log2"), selected = "Z-Score", multiple = FALSE)),
		div(style="display: inline-block;vertical-align:top; width: 30%;",
			selectInput("userCluster", "Cluster Rows", choices = list("No", "Yes"), selected = "No", multiple = FALSE)),
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
		downloadButton("downloadUserHeatmap", "Download My Data Heatmap"),
		downloadButton("downloadHeatmap", "Download Heatmap"),
		downloadButton("downloadGeneHeatmap", "Download Gene Heatmap"),
		downloadButton("downloadDataMatrix", "Download heatmap data matrix for gene set"),
		downloadButton("downloadGOIMatrix", "Download heatmap data matrix for GOI(s)")
	),
	mainPanel(
		tabsetPanel(
			tabPanel(
				"Instructions", 
				h3(strong("Heatmap App Instructions:"), align="center"),
				htmlOutput("text1"),
				div(style="display: inline-block;vertical-align:top; width: 5%;", ""),
				div(style="display: inline-block;vertical-align:top; width: 20%;",	
					img(src='gene_list.png', height = '100%', width='100%')),
				div(style="display: inline-block;vertical-align:top; width: 70%;",
					htmlOutput("text2")),
				htmlOutput("text3"),
				
				h3(strong("Create a heatmap of your data!"), align="center"),
				htmlOutput("text4"),
				div(style="display: inline-block;vertical-align:middle; width: 4%;", ""),
				div(style="display: inline-block;vertical-align:top; width: 35%;",	
					img(src='ExampleDataset.png', height = '100%', width='100%')),
				div(style="display: inline-block;vertical-align:top; width: 10%;",	
					img(src='ColorKey.png', height = '100%', width='100%')),
				div(style="display: inline-block;vertical-align:top; width: 45%;",
					htmlOutput("text5")),
				htmlOutput("text6")
			), 
			tabPanel("Heatmap", h4("Heatmap (Gene Z-score across all samples)", align="center"),plotOutput("plot1"), plotOutput("plot2")),
			tabPanel("Genes", h4("Included Gene List (Ordered)"), tableOutput("table1")),
			tabPanel("Excluded Genes", h4("Genes not found in selected database"), tableOutput("table2")),
			tabPanel("MyData Heatmap", h4("Heatmap of Input Data (Gene Z-score across all samples)", align="center"),plotOutput("plot3")),
			tabPanel("Datasets Summaries", h4(strong("Summaries of Database Datasets for Analysis:"), align="center"), DT::dataTableOutput("tableSummary"))
		)
	)
)



