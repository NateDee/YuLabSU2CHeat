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
		#If SU2C selected, run SU2C heatmap
		if (input$database == "SU2C-Prostate") {
			dbgap <- selectDB
		#If user put in a gene to sort by, use it to resort dbgap:
			if (unlist(strsplit(input$GOI,","))[1] %in% dbgap[,2]) {
				rowOfGene = which(dbgap[,2] == unlist(strsplit(input$GOI,","))[1])
				dbgap.tmp = dbgap[,3:262]
				#Row 1 contains groupings!!!
				dbgap.tmp = dbgap.tmp[,order(dbgap.tmp[1,], dbgap.tmp[rowOfGene,])]
				dbgap = cbind(dbgap[,1:2], dbgap.tmp)
				}
			#Filter dbgap with query
			dbgap.query = dbgap[which(dbgap[,2] %in% query[,1]),]
			#Reorder to match query order
			dbgap.query = dbgap.query[match(query[,1],dbgap.query[,2]), ]
			#If correlation plot is wanted, find correlations and sort!
			if (correlate == TRUE & (unlist(strsplit(input$GOI,","))[1] %in% dbgap[,2])) {
				genedf = dbgap[which(dbgap[,2] == unlist(strsplit(input$GOI,","))[1]),]
				t.dbgap.query = t(as.matrix(dbgap.query[,3:262]))
				colnames(t.dbgap.query) = dbgap.query[,2]
				t.genedf = t(as.matrix(genedf[,3:262]))
				colnames(t.genedf) = genedf[,2]
				correl.mat = t(cor(t.genedf,t.dbgap.query))
				correl.mat = as.data.frame(cbind(correl.mat, correl.mat))
				correl.mat = correl.mat[order(correl.mat[,2]),]
				dbgap.query = dbgap.query[match(rownames(correl.mat),dbgap.query[,2]),]
				}
			
			group = c(rep("green", 35), rep("orange", 78), rep("red", 132), rep("darkred", 15))
			#Set heatmap arguments
			genes = dbgap.query[,2]
			palette = colorRampPalette(c(input$lowcolor, input$midcolor, input$highcolor))(n=100)
			scalelow = input$range[1]
			scalehigh = input$range[2]
			#Make matrix
			mat.dbq = as.matrix(dbgap.query[,3:262])
			rownames(mat.dbq) = dbgap.query[,2]
			#Get rid of NA values
			mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later
			cexRow = as.numeric(as.character(input$xfontsize))
			heatmap.2(mat.dbq, scale="none", trace="none", dendrogram="none", breaks=seq(scalelow, scalehigh, length.out=101), Rowv=FALSE, Colv=FALSE, col=palette, ColSideColors=group, keysize=0.75, key.par=list(cex=0.5), labCol=FALSE, cexRow=cexRow)
			}
		#If TCGA selected, run TCGA heatmap	
		if (input$database == "TCGA-Prostate") {
			tcga <- selectDB
			if (unlist(strsplit(input$GOI,","))[1] %in% tcga[,2]) {
				rowOfGene = which(tcga[,2] == unlist(strsplit(input$GOI,","))[1])
				tcga.tmp = tcga[,3:553]
				tcga.tmp = tcga.tmp[,order(tcga.tmp[1,], tcga.tmp[rowOfGene,])]
				tcga = cbind(tcga[,1:2], tcga.tmp)
				#tcga[,3:553] <- tcga[order(tcga[1, 3:553], tcga[rowOfGene,3:553]), 3:553]
				}
			#Filter tcga with query
			tcga.query = tcga[which(tcga[,2] %in% query[,1]),]
			#Reorder to match query order
			tcga.query = tcga.query[match(query[,1],tcga.query[,2]), ]
			if (correlate == TRUE & (unlist(strsplit(input$GOI,","))[1] %in% tcga[,2])) {
				genedf = tcga[which(tcga[,2] == unlist(strsplit(input$GOI,","))[1]),]
				t.tcga.query = t(as.matrix(tcga.query[,3:553]))
				colnames(t.tcga.query) = tcga.query[,2]
				t.genedf = t(as.matrix(genedf[,3:553]))
				colnames(t.genedf) = genedf[,2]
				correl.mat = t(cor(t.genedf,t.tcga.query))
				correl.mat = as.data.frame(cbind(correl.mat, correl.mat))
				correl.mat = correl.mat[order(correl.mat[,2]),]
				tcga.query = tcga.query[match(rownames(correl.mat),tcga.query[,2]),]
				}
			
			group = c(rep("green", 52), rep("orange", 499))
			#Set heatmap arguments
			genes = tcga.query[,2]
			palette = colorRampPalette(c(input$lowcolor, input$midcolor, input$highcolor))(n=100)
			scalelow = input$range[1]
			scalehigh = input$range[2]
			#Make matrix
			mat.dbq = as.matrix(tcga.query[,3:553])
			rownames(mat.dbq) = tcga.query[,2]
			#Get rid of NA values
			mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later
			cexRow = as.numeric(as.character(input$xfontsize))
			heatmap.2(mat.dbq, scale="none", trace="none", dendrogram="none", breaks=seq(scalelow, scalehigh, length.out=101), Rowv=FALSE, Colv=FALSE, col=palette, ColSideColors=group, keysize=0.75, key.par=list(cex=0.5), labCol=FALSE, cexRow=cexRow)
			}
	}

	drawGeneHeatmap <- function() {
		query <- datasetInput()
		GOI <- geneInput()
		GOI <- unlist(strsplit(GOI, ","))
		selectDB <- databaseInput()
		#If SU2C selected, use SU2C:
		if (input$database == "SU2C-Prostate") {
			dbgap <- selectDB
			#If user put in a gene to sort by, use it to resort dbgap:
			if (GOI[1] %in% dbgap[,2]) {
				rowOfGene = which(dbgap[,2] == GOI[1])
				dbgap.tmp = dbgap[,3:262]
				dbgap.tmp = dbgap.tmp[,order(dbgap.tmp[1,], dbgap.tmp[rowOfGene,])]
				dbgap = cbind(dbgap[,1:2], dbgap.tmp)
				#dbgap[,3:262] <- dbgap[order(dbgap[1, 3:262], dbgap[rowOfGene,3:262]), 3:262]
				}
			#Filter dbgap with query
			dbgap.query = dbgap[which(dbgap[,2] %in% GOI),]
			if (length(GOI) > 1) { 
				dbgap.query = dbgap.query[order(factor(dbgap.query$Gene_name, levels=GOI)),]
				}
			group = c(rep("green", 35), rep("orange", 78), rep("red", 132), rep("darkred", 15))
			palette = colorRampPalette(c(input$lowcolor, input$midcolor, input$highcolor))(n=100)
			scalelow = input$range[1]
			scalehigh = input$range[2]
			#Make matrix double if only one input gene
			if (length(GOI) == 1) {
				mat.dbq = as.matrix(rbind(dbgap.query[,3:262],dbgap.query[,3:262]), )
				rownames(mat.dbq) = c(GOI, GOI)
				} else {
				mat.dbq = as.matrix(dbgap.query[,3:262])
				rownames(mat.dbq) = dbgap.query[,2]
				}
			#Get rid of NA values
			#mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later,****added, see table output
			cexRow = as.numeric(as.character(input$xfontsize))
			heatmap.2(mat.dbq, scale="none", trace="none", dendrogram="none", breaks=seq(scalelow, scalehigh, length.out=101), Rowv=FALSE, Colv=FALSE, col=palette, ColSideColors=group, keysize=0.75, key.par=list(cex=0.5), labCol=FALSE, cexRow=cexRow)
			}
		#If TCGA selected, use TCGA
		if (input$database == "TCGA-Prostate") {
			tcga <- selectDB
			if (GOI[1] %in% tcga[,2]) {
				rowOfGene = which(tcga[,2] == GOI[1])
				tcga.tmp = tcga[,3:553]
				tcga.tmp = tcga.tmp[,order(tcga.tmp[1,], tcga.tmp[rowOfGene,])]
				tcga = cbind(tcga[,1:2], tcga.tmp)
				#tcga[,3:553] <- tcga[order(tcga[1, 3:553], tcga[rowOfGene,3:553]), 3:553]
				}
			#Filter tcga with query
			tcga.query = tcga[which(tcga[,2] %in% GOI),]
			if (length(GOI) > 1) { 
				tcga.query = tcga.query[order(factor(tcga.query$Gene_name, levels=GOI)),]
				}
			group = c(rep("green", 52), rep("orange", 499))
			palette = colorRampPalette(c(input$lowcolor, input$midcolor, input$highcolor))(n=100)
			scalelow = input$range[1]
			scalehigh = input$range[2]
			#Make matrix
			if (length(GOI) == 1) {
				mat.dbq = as.matrix(rbind(tcga.query[,3:553],tcga.query[,3:553]), )
				rownames(mat.dbq) = c(GOI, GOI)
				} else {
				mat.dbq = as.matrix(tcga.query[,3:553])
				rownames(mat.dbq) = tcga.query[,2]
				}
			#Get rid of NA values
			#mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later,****added, see table output
			cexRow = as.numeric(as.character(input$xfontsize))
			heatmap.2(mat.dbq, scale="none", trace="none", dendrogram="none", breaks=seq(scalelow, scalehigh, length.out=101), Rowv=FALSE, Colv=FALSE, col=palette, ColSideColors=group, keysize=0.75, key.par=list(cex=0.5), labCol=FALSE, cexRow=cexRow)
			}		
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
		#If SU2C, use SU2C dbgap:
		if (input$database == "SU2C-Prostate") {
			dbgap = selectDB
			#Filter dbgap with query
			dbgap.query = dbgap[which(dbgap[,2] %in% query[,1]),]
			#Reorder to match query order
			dbgap.query = dbgap.query[match(query[,1],dbgap.query[,2]), ]
			#If correlation plot is wanted, find correlations and sort!
			if ((input$GOI != "") & (correlate == TRUE) & (unlist(strsplit(input$GOI,","))[1] %in% dbgap[,2])) {
				genedf = dbgap[which(dbgap[,2] == unlist(strsplit(input$GOI,","))[1]),]
				t.dbgap.query = t(as.matrix(dbgap.query[,3:262]))
				colnames(t.dbgap.query) = dbgap.query[,2]
				t.genedf = t(as.matrix(genedf[,3:262]))
				colnames(t.genedf) = genedf[,2]
				correl.mat = t(cor(t.genedf,t.dbgap.query))
				correl.mat = as.data.frame(cbind(correl.mat, correl.mat))
				correl.mat = correl.mat[order(correl.mat[,2]),]
				dbgap.query = dbgap.query[match(rownames(correl.mat),dbgap.query[,2]),]
				return(cbind(genes=rownames(correl.mat),Pearson.R=correl.mat[,1]))
				}
			#Set heatmap arguments
			genes = dbgap.query[,2]
			#Make matrix
			mat.dbq = as.matrix(dbgap.query[,3:262])
			rownames(mat.dbq) = dbgap.query[,2]
			#Get rid of NA values
			mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later
			if (correlate == FALSE) {return(cbind(Genes=rownames(mat.dbq)))}
			}
		#If TCGA, use TCGA
		if (input$database == "TCGA-Prostate") {
			tcga = selectDB
			tcga.query = tcga[which(tcga[,2] %in% query[,1]),]
			#Reorder to match query order
			tcga.query = tcga.query[match(query[,1],tcga.query[,2]), ]
			#Reorder to match query order
			tcga.query = tcga.query[match(query[,1],tcga.query[,2]), ]
			if ((correlate == TRUE) & (unlist(strsplit(input$GOI,","))[1] %in% tcga[,2])) {
				genedf = tcga[which(tcga[,2] == unlist(strsplit(input$GOI,","))[1]),]
				t.tcga.query = t(as.matrix(tcga.query[,3:553]))
				colnames(t.tcga.query) = tcga.query[,2]
				t.genedf = t(as.matrix(genedf[,3:553]))
				colnames(t.genedf) = genedf[,2]
				correl.mat = t(cor(t.genedf,t.tcga.query))
				correl.mat = as.data.frame(cbind(correl.mat, correl.mat))
				correl.mat = correl.mat[order(correl.mat[,2]),]
				tcga.query = tcga.query[match(rownames(correl.mat),tcga.query[,2]),]
				return(cbind(Genes=rownames(correl.mat),Pearson.R=correl.mat[,1]))
				}
			#Set heatmap arguments
			genes = tcga.query[,2]
			#Make matrix
			mat.dbq = as.matrix(tcga.query[,3:553])
			rownames(mat.dbq) = tcga.query[,2]
			#Get rid of NA values
			mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later
			if (correlate == FALSE) {return(cbind(Genes=rownames(mat.dbq)))}
			}
	})

	output$table2 <- renderTable({
		query <- datasetInput()
		selectDB <- databaseInput()
		#If SU2C, use SU2C dbgap:
		if (input$database == "SU2C-Prostate") {
			dbgap = selectDB
			#Filter dbgap with query
			dbgap.query = dbgap[which(dbgap[,2] %in% query[,1]),]
			#Reorder to match query order
			dbgap.query = dbgap.query[match(query[,1],dbgap.query[,2]), ]
			#Set heatmap arguments
			genes = dbgap.query[,2]
			#Make matrix
			mat.dbq = as.matrix(dbgap.query[,3:262])
			rownames(mat.dbq) = dbgap.query[,2]
			#Get rid of NA values
			mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later
			table = subset(query[,1], !(query[,1] %in% rownames(mat.dbq)))
			}
		#If TCGA, use TCGA data:
		if (input$database == "TCGA-Prostate") {
			tcga = selectDB
			#Filter tcga with query
			tcga.query = tcga[which(tcga[,2] %in% query[,1]),]
			#Reorder to match query order
			tcga.query = tcga.query[match(query[,1],tcga.query[,2]), ]
			#Set heatmap arguments
			genes = tcga.query[,2]
			#Make matrix
			mat.dbq = as.matrix(tcga.query[,3:553])
			rownames(mat.dbq) = tcga.query[,2]
			#Get rid of NA values
			mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later
			table = subset(query[,1], !(query[,1] %in% rownames(mat.dbq)))
			}
		return(table)
	})
	
	output$downloadHeatmap <- downloadHandler(
	filename <- function() {
		if (input$database == "SU2C-Prostate") {
			paste0(
				basename(unlist(strsplit(input$file1$name, split='.txt', fixed=TRUE))[1]), #Get name of file minus .txt
				'_SU2C-heatmap', 
				'.png', 
				sep=''
				)
			} else if (input$database == "TCGA-Prostate") {
				paste0(
					basename(unlist(strsplit(input$file1$name, split='.txt', fixed=TRUE))[1]), #Get name of file minus .txt
					'_TCGA-heatmap', 
					'.png', 
					sep=''
					)
				}
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
			#If SU2C, use SU2C data:
			if (input$database == "SU2C-Prostate") {
				dbgap = selectDB
				#If user put in a gene to sort by, use it to resort dbgap:
				if (unlist(strsplit(input$GOI,","))[1] %in% dbgap[,2]) {
					rowOfGene = which(dbgap[,2] == unlist(strsplit(input$GOI,","))[1])
					dbgap.tmp = dbgap[,3:262]
					dbgap.tmp = dbgap.tmp[,order(dbgap.tmp[1,], dbgap.tmp[rowOfGene,])]
					dbgap = cbind(dbgap[,1:2], dbgap.tmp)
					#dbgap[,3:262] <- dbgap[order(dbgap[1, 3:262], dbgap[rowOfGene,3:262]), 3:262]
					}
				#Filter dbgap with query
				dbgap.query = dbgap[which(dbgap[,2] %in% query[,1]),]
				#Reorder to match query order
				dbgap.query = dbgap.query[match(query[,1],dbgap.query[,2]), ]
				#If correlation plot is wanted, find correlations and sort!
				if (correlate == TRUE & (unlist(strsplit(input$GOI,","))[1] %in% dbgap[,2])) {
					genedf = dbgap[which(dbgap[,2] == unlist(strsplit(input$GOI,","))[1]),]
					t.dbgap.query = t(as.matrix(dbgap.query[,3:262]))
					colnames(t.dbgap.query) = dbgap.query[,2]
					t.genedf = t(as.matrix(genedf[,3:262]))
					colnames(t.genedf) = genedf[,2]
					correl.mat = t(cor(t.genedf,t.dbgap.query))
					correl.mat = as.data.frame(cbind(correl.mat, correl.mat))
					correl.mat = correl.mat[order(correl.mat[,2]),]
					dbgap.query = dbgap.query[match(rownames(correl.mat),dbgap.query[,2]),]
					}
				write.table(dbgap.query, file, sep="\t", row.names=FALSE)
				}
			#If TCGA, use TCGA data	
			if (input$database == "TCGA-Prostate") {
				tcga <- selectDB
				if (unlist(strsplit(input$GOI,","))[1] %in% tcga[,2]) {
					rowOfGene = which(tcga[,2] == unlist(strsplit(input$GOI,","))[1])
					tcga.tmp = tcga[,3:553]
					tcga.tmp = tcga.tmp[,order(tcga.tmp[1,], tcga.tmp[rowOfGene,])]
					tcga = cbind(tcga[,1:2], tcga.tmp)
					#tcga[,3:553] <- tcga[order(tcga[1, 3:553], tcga[rowOfGene,3:553]), 3:553]
					}
				#Filter tcga with query
				tcga.query = tcga[which(tcga[,2] %in% query[,1]),]
				#Reorder to match query order
				tcga.query = tcga.query[match(query[,1],tcga.query[,2]), ]
				#If correlation plot is wanted, find correlations and sort!
				if (correlate == TRUE & (unlist(strsplit(input$GOI,","))[1] %in% tcga[,2])) {
					genedf = tcga[which(tcga[,2] == unlist(strsplit(input$GOI,","))[1]),]
					t.tcga.query = t(as.matrix(tcga.query[,3:553]))
					colnames(t.tcga.query) = tcga.query[,2]
					t.genedf = t(as.matrix(genedf[,3:553]))
					colnames(t.genedf) = genedf[,2]
					correl.mat = t(cor(t.genedf,t.tcga.query))
					correl.mat = as.data.frame(cbind(correl.mat, correl.mat))
					correl.mat = correl.mat[order(correl.mat[,2]),]
					tcga.query = tcga.query[match(rownames(correl.mat),tcga.query[,2]),]
					}
				write.table(tcga.query, file, sep="\t", row.names=FALSE)
				}
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
				#If SU2C, use SU2C data:
				if (input$database == "SU2C-Prostate") {
					dbgap <- selectDB
					#If user put in a gene to sort by, use it to resort dbgap:
					if (GOI[1] %in% dbgap[,2]) {
						rowOfGene = which(dbgap[,2] == GOI[1])
						dbgap.tmp = dbgap[,3:262]
						dbgap.tmp = dbgap.tmp[,order(dbgap.tmp[1,], dbgap.tmp[rowOfGene,])]
						dbgap = cbind(dbgap[,1:2], dbgap.tmp)
						}
					#Filter dbgap with query
					dbgap.query = dbgap[which(dbgap[,2] %in% GOI),]
					if (length(GOI) > 1) { 
						dbgap.query = dbgap.query[order(factor(dbgap.query$Gene_name, levels=GOI)),]
						}
					write.table(dbgap.query, file, sep="\t", row.names=FALSE)
					}
				#If TCGA, use TCGA data	
				if (input$database == "TCGA-Prostate") {
					tcga <- selectDB
					if (GOI[1] %in% tcga[,2]) {
						rowOfGene = which(tcga[,2] == GOI[1])
						tcga.tmp = tcga[,3:553]
						tcga.tmp = tcga.tmp[,order(tcga.tmp[1,], tcga.tmp[rowOfGene,])]
						tcga = cbind(tcga[,1:2], tcga.tmp)
						}
					#Filter tcga with query
					tcga.query = tcga[which(tcga[,2] %in% GOI),]
					if (length(GOI) > 1) { 
						tcga.query = tcga.query[order(factor(tcga.query$Gene_name, levels=GOI)),]
						}
					write.table(tcga.query, file, sep="\t", row.names=FALSE)
					}
				}
		)
}