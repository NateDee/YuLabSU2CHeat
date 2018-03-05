#Define server logic



server <- function(input,output) {
	databaseInput <- reactive({
		validate(
    		need(input$database != "Select Database...", "Select which database you would like to analyze") 
    		)
		inFile <- input$database
		if (inFile == "SU2C-Prostate") {
			withProgress(message = "Reading SU2C", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 			
			tbl <- as.data.frame(fread("dbgap_su2c_zscore.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
			incProgress(1)
				})
			} else if (inFile == "TCGA-Prostate") {
				withProgress(message = "Reading TCGA", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("TCGA_normal_primary_zscore.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
				})
			} else if (inFile == "SU2C-Benign-vs-PCa") {
				withProgress(message = "Reading SU2C-Benign-vs-PCa", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("dbgap_su2c_zscore_benign-v-pca.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
				})
			} else if (inFile == "SU2C-PCa-vs-CRPC") {
				withProgress(message = "Reading SU2C-PCa-vs-CRPC", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("dbgap_su2c_zscore_pca-v-crpc.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
				incProgress(1)
				})
			} else if (inFile == "SU2C-CRPC-vs-NEPC") {
				withProgress(message = "Reading SU2C-PCa-vs-CRPC", detail = "Should take < 1 min! Watch your R session for progress.", value = 0, { 
				tbl <- as.data.frame(fread("dbgap_su2c_zscore_crpc-v-nepc.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE))
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
			group = c(rep("green", sum(selectDB[1,] == -2)), rep("orange", sum(selectDB[1,] == -1)), rep("red", sum(selectDB[1,] == 1)), rep("darkred", sum(selectDB[1,] == 2)))
		#If user put in a gene to sort by, use it to resort dbgap:
			if (unlist(strsplit(input$GOI,","))[1] %in% selectDB[,"Gene_name"]) {
				rowOfGene = which(selectDB[,"Gene_name"] == unlist(strsplit(input$GOI,","))[1])
				selectDB.tmp = selectDB[,(data.column):(ncol(selectDB))] #Get just the data for reordering
				#Row 1 contains groupings!!!
				selectDB.tmp = selectDB.tmp[,order(selectDB.tmp[1,], selectDB.tmp[rowOfGene,])]
				selectDB = cbind(selectDB[,1:(data.column -1)], selectDB.tmp)
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
		group = c(rep("green", sum(selectDB[1,] == -2)), rep("orange", sum(selectDB[1,] == -1)), rep("red", sum(selectDB[1,] == 1)), rep("darkred", sum(selectDB[1,] == 2)))

		#If user put in a gene to sort by, use it to resort selectDB:
			if (GOI[1] %in% selectDB[,"Gene_name"]) {
				rowOfGene = which(selectDB[,"Gene_name"] == GOI[1])
				selectDB.tmp = selectDB[,(data.column):(ncol(selectDB))]
				selectDB.tmp = selectDB.tmp[,order(selectDB.tmp[1,], selectDB.tmp[rowOfGene,])]
				selectDB = cbind(selectDB[,1:(data.column -1)], selectDB.tmp)
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
			group = c(rep("Group", 2), rep("Benign", sum(selectDB[1,] == -2)), rep("Primary PCa", sum(selectDB[1,] == -1)), rep("CRPC", sum(selectDB[1,] == 1)), rep("NE-CRPC", sum(selectDB[1,] == 2)))
			#Unified handling for any database choosen:
				#If user put in a gene to sort by, use it to resort selectDB:
				if (unlist(strsplit(input$GOI,","))[1] %in% selectDB[,"Gene_name"]) {
					rowOfGene = which(selectDB[,"Gene_name"] == unlist(strsplit(input$GOI,","))[1])
					selectDB.tmp = selectDB[,(data.column):(ncol(selectDB))]
					selectDB.tmp = selectDB.tmp[,order(selectDB.tmp[1,], selectDB.tmp[rowOfGene,])]
					selectDB = cbind(selectDB[,1:(data.column - 1)], selectDB.tmp)
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
					colnames(t.genedf) = genedf[,2]
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
				group = c(rep("Group", 2), rep("Benign", sum(selectDB[1,] == -2)), rep("Primary PCa", sum(selectDB[1,] == -1)), rep("CRPC", sum(selectDB[1,] == 1)), rep("NE-CRPC", sum(selectDB[1,] == 2)))
				#If user put in a gene to sort by, use it to resort selectDB:
				if (GOI[1] %in% selectDB[,"Gene_name"]) {
					rowOfGene = which(selectDB[,"Gene_name"] == GOI[1])
					selectDB.tmp = selectDB[,(data.column):(ncol(selectDB))]
					selectDB.tmp = selectDB.tmp[,order(selectDB.tmp[1,], selectDB.tmp[rowOfGene,])]
					selectDB = cbind(selectDB[,1:(data.column - 1)], selectDB.tmp) #Grab anything in database before data.column starts and cbind (may have multiple gene db rows)
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
}