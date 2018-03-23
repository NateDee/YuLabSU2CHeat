#Define server logic



server <- function(input,output) {
	
	#Create "Instructions" panel:
	output$text1 <- renderText(
		{ paste( 
		'<b>',
			'<ol>',
				'<li style="margin:10px 0";>To visualize how your gene set is regulated in prostate disease progression, start by selecting a database to query.</li>',
				'<li style="margin:10px 0";>To see what each database contains, navigate to the "Datasets Summaries." tab</li>',
				'<li style="margin:10px 0";>Upload your gene list or choose a list from the "test_gene_sets" folder. Your uploaded gene list must be in the following format:</li>',
			'</ol>',	
		'</b>'
		)}
	)
	
	output$text2 <- renderText(
		{ paste( 
		'<p>',
			'<ul>',
				'<li style="margin:10px 0";>The first row must contain headers (See Row 1 to the left).</li>',
				'<li style="margin:10px 0";>The first column (Column A) must contain gene names (gene symbols).</li>',
				'<li style="margin:10px 0";>If you only have a list of genes, duplicate row one (as shown in Column B).</li>',
				'<li style="margin:10px 0";>If you have other data like gene expression measurements, these can be left as they are and you don\'t need to duplicate column 1.</li>',
				'<li style="margin:10px 0";>To make a heatmap of your own data, see the instructions below under "Create a heatmap of your data!"</li>',
				'<b><li style="margin:10px 0";>Save your data in the "Text (Tab delimited) (*.txt)" format.</li></b>',
			'</ul>',	
		'</p>'
		)}
	)
	output$text3 <- renderText(
		{ paste( 
		'<b>',
			'<ol start="4">',
				'<li style="margin:10px 0";>The "Heatmap" tab will now display the expression of your genes in the tissues of the selected database.</li>',
				'<li style="margin:10px 0";>By default, these genes will be in the same order as your input list.</li>',
				'<li style="margin:10px 0";>If you would like to see how these genes correlate with a gene of interest, type in your gene (or genes) in the side bar.</li>',
				'<ul><li>For example type AR, if you would like to see multiple genes, separate them by a comma: AR,FOXA1,MYCN</li> <li>Typing in a gene will also create a heatmap for that gene that maches up with the main heatmap.</li> <li>Selecting the "Sort gene rows by correlation with GOI?" option will run a pearson R correlation for all genes in the gene set versus your selected gene of interest.  The heatmap rows will be resorted by correlation with your GOI from low to high Pearson R</li></ul>',
				'<li style="margin:10px 0";>The "Genes" panel provides a table of the gene order in the heatmap and the Pearson R if correlation sort is selected.</li>',
				'<li style="margin:10px 0";>The "Excluded Genes" panel provides a table of the genes from your list that were not found in the selected database.</li>',
				'<li style="margin:10px 0";>The generated heatmaps and data can be downloaded by navigating to the bottom of the left sidebar and selecting which data to download.</li>',
			'</ol>',	
		'</b>'
		)}
	)
	
	output$text4 <- renderText(
		{ paste( 
		'<b>',
			'<ol>',
				'<li style="margin:20px 0";>Format your expression data according to specifications below.</li>',
			'</ol>',	
		'</b>'
		)}
	)
	
	output$text5 <- renderText(
		{ paste( 
		'<p>',
			'<ul>',
				'<li style="margin:10px 0";>The first row must contain headers (See Row 1 to the left).</li>',
				'<li style="margin:10px 0";>The first column (Column A) must contain gene names (gene symbols).</li>',
				'<b><li style="margin:10px 0";>The second column must contain a group number, see the color key to choose a group color label.</li></b>',
				'<li style="margin:10px 0";>Your expression measurements can be in multiple formats (RPKM, Log2, Z-Score, Normalized Microarray).  Make sure to adjust the heatmap scale accordingly!</li>',
			'</ul>',	
		'</p>'
		)}
	)
	
	output$text6 <- renderText(
		{ paste( 
		'<b>',
			'<ol start="2">',
				'<li style="margin:10px 0";>Once your data is formatted, upload the .txt file in the same way you would a gene list.</li>',
				'<li style="margin:10px 0";>Under the "MyData Heatmap Options", choose to "Yes" under "Heatmap from Input?"</li>',
				'<li style="margin:10px 0";>Navigate to the "MyData Heatmap" Tab, by default the genes will be in the order you provided.</li>',
				'<li style="margin:10px 0";>You can cluster the genes by choosing "Cluster Rows" = "Yes".</li>',
				'<li style="margin:10px 0";>Choosinge a normalization method:</li>',
				'<ul><li>If your data is in RPKM/FPKM or raw normalized values from a microarray experiment you will probably want to choose "Z-Score" or "Log2" as a normalization method.</li> <li>If your data is already in Log2 or Z-Score format, you should choose "Raw".</li> <b><li>IMPORTANT: If you choose Log2, you should adjust the heatmap scale to have a low value of 0!!! The method is a log2(data + 1) adjustment, log2(0+1) = 0! </li></b></ul>',
				'<li style="margin:10px 0";>The generated heatmaps and data can be downloaded by navigating to the bottom of the left sidebar and selecting "MyData Download Heatmap".</li>',
			'</ol>',	
		'</b>'
		)}
	)
	
	# Import Database, based on user choice
	databaseInput <- reactive({
		validate(
    		need(input$database != "Select Database...", "Select which database you would like to analyze") 
    		)
		inFile <- input$database
		if (inFile == "DBGAP-Prostate") {
			withProgress(message = "Reading DBGAP", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 			
				tbl <- as.data.frame(fread("dbgap_su2c_zscore.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
			})
		} else if (inFile == "TCGA-Prostate") {
			withProgress(message = "Reading TCGA", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("TCGA_normal_primary_zscore.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
			})
		} else if (inFile == "DBGAP-Benign-vs-PCa") {
			withProgress(message = "Reading DBGAP-Benign-vs-PCa", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("dbgap_su2c_zscore_benign-v-pca.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
			})
		} else if (inFile == "DBGAP-PCa-vs-CRPC") {
			withProgress(message = "Reading DBGAP-PCa-vs-CRPC", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("dbgap_su2c_zscore_pca-v-crpc.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
			})
		} else if (inFile == "DBGAP-CRPC-vs-NEPC") {
			withProgress(message = "Reading DBGAP-CRPC-vs-NEPC", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("dbgap_su2c_zscore_crpc-v-nepc.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
			})
		} else if (inFile == "GLINSKY_U133A") {
			withProgress(message = "Reading GLINSKY_U133A", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("./datasets/LAB_GLINSKY_U133A_z-score.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
			})
		} else if (inFile == "VARAMBALLY_GSE3325") {
			withProgress(message = "Reading VARAMBALLY_GSE3325", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("./datasets/LAB_GSE3325_VARAM_z-score.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
			})
		} else if (inFile == "LAPOINTE_GSE3933") {
			withProgress(message = "Reading LAPOINTE_GSE3933", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("./datasets/LAB_GSE3933_LAPOINTE_z-score.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
			})
		} else if (inFile == "TOMLINS_GSE6099") {
			withProgress(message = "Reading TOMLINS_GSE6099", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("./datasets/LAB_GSE6099_TOMLINS_z-score.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
			})
		} else if (inFile == "YU_GSE6919") {
			withProgress(message = "Reading YU_GSE6919", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("./datasets/LAB_GSE6919_YU_z-score.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
			})
		} else if (inFile == "SBONER_GSE16560") {
			withProgress(message = "Reading SBONER_GSE16560", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("./datasets/LAB_GSE16560_SBONER_z-score.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
			})
		} else if (inFile == "TAYLOR_GSE21034") {
			withProgress(message = "Reading TAYLOR_GSE21034", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("./datasets/LAB_GSE21034_TAYLOR_z-score.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
			})
		} else if (inFile == "LONG_GSE26367") {
			withProgress(message = "Reading LONG_GSE26367", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("./datasets/LAB_GSE26367_LONG_z-score.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
			})
		} else if (inFile == "CAI_GSE32269") {
			withProgress(message = "Reading CAI_GSE32269", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("./datasets/LAB_GSE32269_CAI_z-score.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
			})
		} else if (inFile == "GRASSO_GSE35988") {
			withProgress(message = "Reading GRASSO_GSE35988", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("./datasets/LAB_GSE35988_GRASSO_z-score.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
			})
		} else if (inFile == "ROSS_GSE70770") {
			withProgress(message = "Reading ROSS_GSE70770", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("./datasets/LAB_GSE70770_ROSS_z-score.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
			})
		} else if (inFile == "VANAJA_U133A") {
			withProgress(message = "Reading VANAJA_U133A", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("./datasets/LAB_VANAJA_U133A_z-score.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
			})
		}
		return(tbl)
	})
	
	#Upload user file
	datasetInput <- reactive({
		validate(
    		need(input$file1 != 0, "To begin drawing a heatmap, please select a file for input") 
    		)
		inFile <- input$file1
		if (is.null(inFile))
			return(NULL)
		tbl <- as.data.frame(fread(inFile$datapath, sep="\t", header=TRUE, stringsAsFactors=FALSE))
		return(tbl)
	})
	
	#User input Gene of Interest
	geneInput <- reactive({
		validate(
    		need(input$GOI != "", "Enter a gene of interest, or genes seperated by commas (i.e. FOXA1,AR,CXCR7)") 
    		)
		inGOI <- input$GOI
		if (is.null(inGOI))
			return(NULL)
		return(inGOI)
	})
	
	#Does user want to resort heatmap by correlation with gene of interest?
	correlInput <- reactive({
		inCorrel <- input$correlationSort
		return(inCorrel)
	})

	#Make heatmap of user uploaded data
	doInputHeat <- reactive({
		doHeat <- input$doUserHeatmap
		return(doHeat)
	})
	
	#Draw the main heatmap using database selected
	drawHeatmap <- function() {
		selectDB <- databaseInput()
		query <- datasetInput()
		correlate <- correlInput()
				
			data.column = which(colnames(selectDB) == "Gene_name") + 1 #Data starts after Gene_name column, so add 1, for index position
			group = c(rep("green", sum(selectDB[1,] == -2)), rep("orange", sum(selectDB[1,] == -1)), rep("black", sum(selectDB[1,] == 0)), rep("red", sum(selectDB[1,] == 1)), rep("darkred", sum(selectDB[1,] == 2)))
		#If user put in a gene to sort by, use it to resort dbgap:
			if (unlist(strsplit(input$GOI,","))[1] %in% selectDB[,"Gene_name"]) {
				rowOfGene = which(selectDB[,"Gene_name"] == unlist(strsplit(input$GOI,","))[1])
				selectDB.tmp = selectDB[,(data.column):(ncol(selectDB))] #Get just the data for reordering
				#Row 1 contains groupings!!!
				selectDB.tmp = selectDB.tmp[,order(selectDB.tmp[1,], selectDB.tmp[rowOfGene,])]
				selectDB = cbind("Gene_name"=selectDB$Gene_name, selectDB.tmp)
				data.column = which(colnames(selectDB) == "Gene_name") + 1 #Data starts after Gene_name column, so add 1, for index position
			}
			#Filter dbgap with query
			selectDB.query = selectDB[which(selectDB[,"Gene_name"] %in% query[,1]),]
			#Reorder to match query order
			selectDB.query = selectDB.query[match(query[,1],selectDB.query[,"Gene_name"]), ]
			#If correlation plot is wanted, find correlations and sort!
			if (correlate == TRUE & (unlist(strsplit(input$GOI,","))[1] %in% selectDB[,"Gene_name"])) {
				genedf = selectDB[which(selectDB[,"Gene_name"] == unlist(strsplit(input$GOI,","))[1]),]
				t.selectDB.query = t(as.matrix(selectDB.query[,(data.column):(ncol(selectDB.query))]))
				colnames(t.selectDB.query) = selectDB.query[,"Gene_name"]
				t.genedf = t(as.matrix(genedf[,(data.column):(ncol(genedf))]))
				colnames(t.genedf) = genedf[,"Gene_name"]
				correl.mat = t(cor(t.genedf,t.selectDB.query))
				correl.mat = as.data.frame(cbind(correl.mat, correl.mat))
				correl.mat = correl.mat[order(correl.mat[,2]),]
				selectDB.query = selectDB.query[match(rownames(correl.mat),selectDB.query[,"Gene_name"]),]
			}
			
			#Set heatmap arguments
			palette = colorRampPalette(c(input$lowcolor, input$midcolor, input$highcolor))(n=100)
			scalelow = input$range[1]
			scalehigh = input$range[2]
			#Make matrix
			mat.dbq = as.matrix(selectDB.query[,(data.column):(ncol(selectDB))])
			rownames(mat.dbq) = selectDB.query[,"Gene_name"]
			#Get rid of NA values
			mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later
			cexRow = as.numeric(as.character(input$xfontsize))
			heatmap.2(mat.dbq, scale="none", trace="none", dendrogram="none", breaks=seq(scalelow, scalehigh, length.out=101), Rowv=FALSE, Colv=FALSE, col=palette, ColSideColors=group, keysize=0.75, key.par=list(cex=0.5), labCol=FALSE, cexRow=cexRow)
	}

	#Draw heatmap from gene of interest
	drawGeneHeatmap <- function() {
		query <- datasetInput()
		GOI <- geneInput()
		GOI <- unlist(strsplit(GOI, ","))
		selectDB <- databaseInput()
		data.column = which(colnames(selectDB) == "Gene_name") + 1 #Data starts after Gene_name column, so add 1, for index position
		group = c(rep("green", sum(selectDB[1,] == -2)), rep("orange", sum(selectDB[1,] == -1)), rep("black", sum(selectDB[1,] == 0)), rep("red", sum(selectDB[1,] == 1)), rep("darkred", sum(selectDB[1,] == 2)))

		#If user put in a gene to sort by, use it to resort selectDB:
			if (GOI[1] %in% selectDB[,"Gene_name"]) {
				rowOfGene = which(selectDB[,"Gene_name"] == GOI[1])
				selectDB.tmp = selectDB[,(data.column):(ncol(selectDB))]
				selectDB.tmp = selectDB.tmp[,order(selectDB.tmp[1,], selectDB.tmp[rowOfGene,])]
				selectDB = cbind("Gene_name"=selectDB$Gene_name, selectDB.tmp)
				data.column = which(colnames(selectDB) == "Gene_name") + 1 #Data starts after Gene_name column, so add 1, for index position
			}
			#Filter selectDB with query
			selectDB.query = selectDB[which(selectDB[,"Gene_name"] %in% GOI),]
			if (length(GOI) > 1) { 
				selectDB.query = selectDB.query[order(factor(selectDB.query$Gene_name, levels=GOI)),]
			}
			palette = colorRampPalette(c(input$lowcolor, input$midcolor, input$highcolor))(n=100)
			scalelow = input$range[1]
			scalehigh = input$range[2]
			#Make matrix double if only one input gene
			if (length(GOI) == 1) {
				mat.dbq = as.matrix(rbind(selectDB.query[,(data.column):(ncol(selectDB.query))],selectDB.query[,(data.column):(ncol(selectDB.query))]), )
				rownames(mat.dbq) = c(GOI, GOI)
			} else {
				mat.dbq = as.matrix(selectDB.query[,(data.column):(ncol(selectDB.query))])
				rownames(mat.dbq) = selectDB.query[,"Gene_name"]
			}
			#Get rid of NA values
			#mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later,****added, see table output
			cexRow = as.numeric(as.character(input$xfontsize))
			heatmap.2(mat.dbq, scale="none", trace="none", dendrogram="none", breaks=seq(scalelow, scalehigh, length.out=101), Rowv=FALSE, Colv=FALSE, col=palette, ColSideColors=group, keysize=0.75, key.par=list(cex=0.5), labCol=FALSE, cexRow=cexRow)		
	}
	
	#Draw user heatmap from uploaded data
	drawUserHeatmap <- function() {
		#Run only if user checks box to ask:
		if (doInputHeat() == "Yes") {
			query <- datasetInput()
			#GOI <- geneInput()
			#GOI <- unlist(strsplit(GOI, ","))
			group = c(rep("green", sum(query[1,] == -2)), rep("orange", sum(query[1,] == -1)), rep("black", sum(query[1,] == 0)), rep("red", sum(query[1,] == 1)), rep("darkred", sum(query[1,] == 2)))
			#Process from raw values to z-score
			rawToZscore <- function(raw) {
				#Calculate Z-Score: For each row - (Cell - mean)/SD = Z-score units
				row.means = rowMeans(raw[,2:ncol(raw)], na.rm=TRUE)
				row.sd = apply(raw[,2:ncol(raw)], 1, function(x) sd(x, na.rm=TRUE))
				raw[,2:ncol(raw)] = (raw[,2:ncol(raw)] - row.means)/row.sd
				zscore = rbind(raw[2:nrow(raw),])
				return(zscore)
			}
			
			rawToRaw <- function(raw) {
				#Just get rid of group row
				raw <- raw[2:nrow(raw), ]
				return(raw)
			}
			
			rawToLog2 <- function(raw) {
				#Just transform data to log2 + 1
				raw[2:nrow(raw), 2:ncol(raw)] <- log2(raw[2:nrow(raw), 2:ncol(raw)] + 1)
				#Already got groups so drop groups
				log2 <- raw[2:nrow(raw),]
				return(log2)
			}
			
			#Create matrix and options for heatmap - run as function
			returnHeat <- function(df) {
				#Get matrix
				mat = as.matrix(df[,2:ncol(query)])
				rownames(mat) = df[,1]
				colnames(mat) = colnames(df)[2:ncol(query)]
				#Heatmap options:
				palette = colorRampPalette(c(input$lowcolor, input$midcolor, input$highcolor))(n=100)
				scalelow = input$range[1]
				scalehigh = input$range[2]
				cexRow = as.numeric(as.character(input$xfontsize))
				#Cluster rows?
				if (input$userCluster == "Yes") {
					Rowv=TRUE
				} else {
					Rowv=FALSE
				}
				#Cluster columns?
				if (input$userClusterCol == "Yes") {
					Colv=TRUE
					
				} else {
					Colv=FALSE
				}
				#Set dendrogram
				if (input$userCluster == "Yes" & input$userClusterCol == "Yes") {
					dendrogram = "both"
				} else if (input$userCluster == "Yes" & input$userClusterCol == "No"){
					dendrogram = "row"
				} else if (input$userCluster == "No" & input$userClusterCol == "Yes") {
					dendrogram = "column"
				} else {
					dendrogram = "none"
				}
				
				#Run Heatmap:
				userheat = heatmap.2(mat, scale = "none", trace="none", dendrogram=dendrogram, breaks=seq(scalelow, scalehigh, length.out=101), Rowv=Rowv, Colv=Colv, col=palette, ColSideColors=group, keysize=0.75, key.par = list(cex=0.5), labCol=colnames(mat), cexRow=cexRow)
				# Get table of genes (ordered to match heatmap, set global variable to be used by Table 3
				Genes_Order <- rownames(mat)[userheat$rowInd]
				output$table3 <- renderTable({
					return(cbind(Genes_Order))
				})
			}
			
			processData <- function(df) {
				#Normalize data based on choices to user: "Z-Score", "Raw", "Log2"
				if (input$normalize == "Z-Score") {
					df <- rawToZscore(df)
				} else if (input$normalize == "Raw") {
					df <- rawToRaw(df)
				} else if (input$normalize == "Log2") {
					df <- rawToLog2(df)
				}
				return(df)			
			}
			query <- processData(query)	
			# Run heatmap and set as var (will still output heatmap), then grab gene list set as global var to plot
			returnHeat(query)
		}
	}
	
	#Output plot1, database heatmap
	output$plot1 <- renderPlot({
		drawHeatmap()
	})
	
	#Output plot2 Gene(s) of interest heatmap
	output$plot2 <- renderPlot({
		drawGeneHeatmap()
	})
	
	#Output plot3, user heatmap of uploaded data
	output$plot3 <- renderPlot({
		drawUserHeatmap()
	})

	#Table 1 = "Genes included" tab, show genes that matched from user upload list to selected database
	# Can move this into the Heatmap Function**************
	output$table1 <- renderTable({
		query <- datasetInput()
		selectDB <- databaseInput()
		correlate = correlInput()
		data.column = which(colnames(selectDB) == "Gene_name") + 1 #Data starts after Gene_name column, so add 1, for index position
		#If SU2C, use SU2C selectDB:
			#Filter selectDB with query
			selectDB.query = selectDB[which(selectDB[,"Gene_name"] %in% query[,1]),]
			#Reorder to match query order
			selectDB.query = selectDB.query[match(query[,1],selectDB.query[,"Gene_name"]), ]
			#If correlation plot is wanted, find correlations and sort!
			if ((input$GOI != "") & (correlate == TRUE) & (unlist(strsplit(input$GOI,","))[1] %in% selectDB[,"Gene_name"])) {
				genedf = selectDB[which(selectDB[,"Gene_name"] == unlist(strsplit(input$GOI,","))[1]),]
				t.selectDB.query = t(as.matrix(selectDB.query[,(data.column):(ncol(selectDB.query))]))
				colnames(t.selectDB.query) = selectDB.query[,"Gene_name"]
				t.genedf = t(as.matrix(genedf[,(data.column):(ncol(genedf))]))
				colnames(t.genedf) = genedf[,"Gene_name"]
				correl.mat = t(cor(t.genedf,t.selectDB.query))
				correl.mat = as.data.frame(cbind(correl.mat, correl.mat)) # Need two columns to be dataframe, duplicate
				correl.mat = correl.mat[order(correl.mat[,2]),] # Order by low to high correlation
				selectDB.query = selectDB.query[match(rownames(correl.mat),selectDB.query[,"Gene_name"]),]
				correl.mat = correl.mat[complete.cases(correl.mat),] # Get rid of NA values
				return(cbind(Genes=rownames(correl.mat),Pearson.R=correl.mat[,1]))
			}
			#Set heatmap arguments
			genes = selectDB.query[,"Gene_name"]
			#Make matrix
			mat.dbq = as.matrix(selectDB.query[,(data.column):(ncol(selectDB.query))])
			rownames(mat.dbq) = selectDB.query[,"Gene_name"]
			#Get rid of NA values
			mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later
			if (correlate == FALSE) {return(cbind(Genes=rownames(mat.dbq)))}
	})

	#Table 2 - Genes not found from user list in selected database
	output$table2 <- renderTable({
		query <- datasetInput()
		selectDB <- databaseInput()
		data.column = which(colnames(selectDB) == "Gene_name") + 1 #Data starts after Gene_name column, so add 1, for index position
		return(cbind(Genes= subset(query[,1], !(query[,1] %in% selectDB[,"Gene_name"]))))
	})
	
	#Download user heatmap - "MyData")
	output$downloadUserHeatmap <- downloadHandler(
		filename <- function() {
			paste0(
				basename(unlist(strsplit(input$file1$name, split='.txt', fixed=TRUE))[1]), #Get name of file minus .txt
				"_",
				input$normalize,
				'_heatmap',
				'.png', 
				sep=''
				)
			},
		content <- function(file) {
			png(file, width = 6, height = 10, units ="in", res = 300)
			drawUserHeatmap()
			dev.off()
		}
	)
	
	#Download main heatmap
	output$downloadHeatmap <- downloadHandler(
		filename <- function() {
			paste0(
				basename(unlist(strsplit(input$file1$name, split='.txt', fixed=TRUE))[1]), #Get name of file minus .txt
				input$database,
				'-heatmap',
				'.png', 
				sep=''
				)
			},
		content <- function(file) {
			png(file, width = 6, height = 10, units = "in", res = 300)
			drawHeatmap()
			dev.off()
			}
	)
	
	#Download gene of interest heatmap
	output$downloadGeneHeatmap <- downloadHandler(
		filename <- function() {
			paste0(
				input$GOI,
				"_",
				basename(unlist(strsplit(input$file1$name, split='.txt', fixed=TRUE))[1]), #Get name of file minus .txt
				"_",
				input$database,
				'-heatmap', 
				'.png', 
				sep=''
			)
		},
		content <- function(file) {
			png(file, width = 6, height = 10, units = "in", res = 300)
			drawGeneHeatmap()
			dev.off()
		}
	)
	
	#Download data matrix, which is zscore data, for selected DB (matched genes from user list)
	output$downloadDataMatrix <- downloadHandler(
		filename <- function() {
			paste0(
				input$GOI,
				"_",
				basename(unlist(strsplit(input$file1$name, split='.txt', fixed=TRUE))[1]), #Get name of file minus .txt
				"_",
				input$database,
				'-heatmap-matrix', 
				'.txt', 
				sep=''
			)
		},
		
		content <- function(file) {
			query <- datasetInput()
			selectDB <- databaseInput()
			correlate = correlInput()
			data.column = which(colnames(selectDB) == "Gene_name") + 1 #Data starts after Gene_name column, so add 1, for index position
			group = c(rep("Group", 2), rep("Benign", sum(selectDB[1,] == -2)), rep("Primary PCa", sum(selectDB[1,] == -1)), rep("Metastasis", sum(selectDB[1,] == 0)), rep("CRPC", sum(selectDB[1,] == 1)), rep("NE-CRPC", sum(selectDB[1,] == 2)))
			#Unified handling for any database choosen:
				#If user put in a gene to sort by, use it to resort selectDB:
				if (unlist(strsplit(input$GOI,","))[1] %in% selectDB[,"Gene_name"]) {
					rowOfGene = which(selectDB[,"Gene_name"] == unlist(strsplit(input$GOI,","))[1])
					selectDB.tmp = selectDB[,(data.column):(ncol(selectDB))]
					selectDB.tmp = selectDB.tmp[,order(selectDB.tmp[1,], selectDB.tmp[rowOfGene,])]
					selectDB = cbind("Gene_name"=selectDB$Gene_name, selectDB.tmp)
				}
				#Filter selectDB with query
				selectDB.query = selectDB[which(selectDB[,"Gene_name"] %in% query[,1]),]
				#Reorder to match query order
				selectDB.query = selectDB.query[match(query[,1],selectDB.query[,"Gene_name"]), ]
				#If correlation plot is wanted, find correlations and sort!
				if (correlate == TRUE & (unlist(strsplit(input$GOI,","))[1] %in% selectDB[,"Gene_name"])) {
					genedf = selectDB[which(selectDB[,"Gene_name"] == unlist(strsplit(input$GOI,","))[1]),]
					t.selectDB.query = t(as.matrix(selectDB.query[,(data.column):(ncol(selectDB.query))]))
					colnames(t.selectDB.query) = selectDB.query[,"Gene_name"]
					t.genedf = t(as.matrix(genedf[ ,(data.column):(ncol(genedf))]))
					colnames(t.genedf) = genedf[,"Gene_name"]
					correl.mat = t(cor(t.genedf,t.selectDB.query))
					correl.mat = as.data.frame(cbind(correl.mat, correl.mat))
					correl.mat = correl.mat[order(correl.mat[,2]),]
					selectDB.query = selectDB.query[match(rownames(correl.mat),selectDB.query[,"Gene_name"]),] #Match correlation sort order
				}
				#Add sample types to output
				names = colnames(selectDB.query) #Have to reset names after rbind
				selectDB.query = rbind(group, selectDB.query)
				colnames(selectDB.query) = names #Reset colnames 
				#Write table to output
				write.table(selectDB.query, file, sep="\t", row.names=FALSE)
		}
	)
	
	#Download zscore matrix for gene(s) of interest
	output$downloadGOIMatrix <- downloadHandler(
		filename <- function() {
			paste0(
				input$GOI,
				"_",
				basename(unlist(strsplit(input$file1$name, split='.txt', fixed=TRUE))[1]), #Get name of file minus .txt
				"_",
				input$database,
				'-GOIs-matrix', 
				'.txt', 
				sep=''
				)
		},
			
		content <- function(file) {
			query <- datasetInput()
			GOI <- geneInput()
			GOI <- unlist(strsplit(GOI, ","))
			selectDB <- databaseInput()
			data.column = which(colnames(selectDB) == "Gene_name") + 1 #Data starts after Gene_name column, so add 1, for index position
			group = c(rep("Group", 2), rep("Benign", sum(selectDB[1,] == -2)), rep("Primary PCa", sum(selectDB[1,] == -1)), rep("Metastasis", sum(selectDB[1,] == 0)), rep("CRPC", sum(selectDB[1,] == 1)), rep("NE-CRPC", sum(selectDB[1,] == 2)))
			#If user put in a gene to sort by, use it to resort selectDB:
			if (GOI[1] %in% selectDB[,"Gene_name"]) {
				rowOfGene = which(selectDB[,"Gene_name"] == GOI[1])
				selectDB.tmp = selectDB[,(data.column):(ncol(selectDB))]
				selectDB.tmp = selectDB.tmp[,order(selectDB.tmp[1,], selectDB.tmp[rowOfGene,])]
				selectDB = cbind("Gene_name"=selectDB$Gene_name, selectDB.tmp) 
				data.column = which(colnames(selectDB) == "Gene_name") + 1 #Data starts after Gene_name column, so add 1, for index position
				}
			#Filter selectDB with query
			selectDB.query = selectDB[which(selectDB[,"Gene_name"] %in% GOI),]
			if (length(GOI) > 1) { 
				selectDB.query = selectDB.query[order(factor(selectDB.query$Gene_name, levels=GOI)),]
				}
			#Add sample types to output matrix, save colnames then reset after rbind to be sure headers match selectDB
			names = colnames(selectDB.query)
			selectDB.query = rbind(group, selectDB.query)
			colnames(selectDB.query) = names
			write.table(selectDB.query, file, sep="\t", row.names=FALSE)
		}
	)
	
	#Tab of table with database summaries - rendered with "datatable" to jazz it up and provide user interaction
	#Do Datasets summary table
	databaseSummary <- as.data.frame(fread("Datasets-Summary.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE))
	
	output$tableSummary <- DT::renderDataTable(
		datatable(databaseSummary, 
		options = list(
			pageLength = 20, 
			columnDefs = list(list(className = 'dt-center', targets=1:9))
		),	
		rownames=FALSE) %>%
			formatStyle('Benign', color = 'black', backgroundColor = 'green') %>%
			formatStyle('PCa', color = 'black', backgroundColor = 'orange') %>%
			formatStyle('Met', color = 'white', backgroundColor = 'black') %>%
			formatStyle('CRPC', color = 'black', backgroundColor = 'red') %>%
			formatStyle('NEPC', color = 'black', backgroundColor = 'darkred')
	)	
}