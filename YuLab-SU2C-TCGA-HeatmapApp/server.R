#Define server logic



server <- function(input,output) {
	
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
		
	geneInput <- reactive({
		validate(
    		need(input$GOI != "", "Enter a gene of interest, or genes seperated by commas (i.e. FOXA1,AR,CXCR7)") 
    		)
		inGOI <- input$GOI
		if (is.null(inGOI))
			return(NULL)
		return(inGOI)
	})
		
	correlInput <- reactive({
		inCorrel <- input$correlationSort
		return(inCorrel)
	})	
	
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
	
	output$plot1 <- renderPlot({
		drawHeatmap()
	})
	
	output$plot2 <- renderPlot({
		drawGeneHeatmap()
	})

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

	output$table2 <- renderTable({
		query <- datasetInput()
		selectDB <- databaseInput()
		data.column = which(colnames(selectDB) == "Gene_name") + 1 #Data starts after Gene_name column, so add 1, for index position
		return(cbind(Genes= subset(query[,1], !(query[,1] %in% selectDB[,"Gene_name"]))))
	})
	
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