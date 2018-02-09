#Define server logic



server <- function(input,output) {
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
    		need(input$GOI %in% dbgap[,2], "To enter a valid gene of interest for a matched heatmap") 
    		)
		inGOI <- input$GOI
		if (is.null(inGOI))
			return(NULL)
		return(inGOI)
		})
	
	drawHeatmap <- function() {
		query <- datasetInput()
		#If user put in a gene to sort by, use it to resort dbgap:
			if (input$GOI %in% dbgap[,2]) {
				rowOfGene = which(dbgap[,2] == input$GOI)
				dbgap.tmp = dbgap[,3:262]
				dbgap.tmp = dbgap.tmp[,order(dbgap.tmp[1,], dbgap.tmp[rowOfGene,])]
				dbgap = cbind(dbgap[,1:2], dbgap.tmp)
				#dbgap[,3:262] <- dbgap[order(dbgap[1, 3:262], dbgap[rowOfGene,3:262]), 3:262]
				}
			#Filter dbgap with query
			dbgap.query = dbgap[which(dbgap[,2] %in% query[,1]),]
			#Reorder to match query order
			dbgap.query = dbgap.query[match(query[,1],dbgap.query[,2]), ]
			group = c(rep("green", 35), rep("orange", 78), rep("red", 132), rep("darkred", 15))
			#Set heatmap arguments
			genes = dbgap.query[,2]
			palette = colorRampPalette(c("green", "black", "red"))(n=100)
			scalelow = input$range[1]
			scalehigh = input$range[2]
			#Make matrix
			mat.dbq = as.matrix(dbgap.query[,3:262])
			rownames(mat.dbq) = dbgap.query[,2]
			#Get rid of NA values
			mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later

			heatmap.2(mat.dbq, scale="none", trace="none", dendrogram="none", breaks=seq(scalelow, scalehigh, length.out=101), Rowv=FALSE, Colv=FALSE, col=palette, ColSideColors=group, keysize=0.75, key.par=list(cex=0.5), labCol=FALSE, cexRow=0.5)
	}

	drawGeneHeatmap <- function() {
		query <- datasetInput()
		GOI <- geneInput()
		#If user put in a gene to sort by, use it to resort dbgap:
			if (GOI %in% dbgap[,2]) {
				rowOfGene = which(dbgap[,2] == GOI)
				dbgap.tmp = dbgap[,3:262]
				dbgap.tmp = dbgap.tmp[,order(dbgap.tmp[1,], dbgap.tmp[rowOfGene,])]
				dbgap = cbind(dbgap[,1:2], dbgap.tmp)
				#dbgap[,3:262] <- dbgap[order(dbgap[1, 3:262], dbgap[rowOfGene,3:262]), 3:262]
				}
			#Filter dbgap with query
			dbgap.query = dbgap[which(dbgap[,2] == input$GOI),]
			group = c(rep("green", 35), rep("orange", 78), rep("red", 132), rep("darkred", 15))
			palette = colorRampPalette(c("green", "black", "red"))(n=100)
			scalelow = input$range[1]
			scalehigh = input$range[2]
			#Make matrix
			mat.dbq = as.matrix(rbind(dbgap.query[,3:262],dbgap.query[,3:262]), )
			assign('mat', mat.dbq, envir = globalenv())
			#Get rid of NA values
			#mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later,****added, see table output
			heatmap.2(mat.dbq, scale="none", trace="none", dendrogram="none", breaks=seq(scalelow, scalehigh, length.out=101), Rowv=FALSE, Colv=FALSE, col=palette, ColSideColors=group, keysize=0.75, key.par=list(cex=0.5), labCol=FALSE, cexRow=0.5)
	}
	
	output$plot1 <- renderPlot({
		drawHeatmap()
	})
	
	output$plot2 <- renderPlot({
		drawGeneHeatmap()
	})

	output$table1 <- renderTable({
		query <- datasetInput()
		#Filter dbgap with query
		dbgap.query = dbgap[which(dbgap[,2] %in% query[,1]),]
		#Reorder to match query order
		dbgap.query = dbgap.query[match(query[,1],dbgap.query[,2]), ]
		group = c(rep("green", 35), rep("orange", 78), rep("red", 132), rep("darkred", 15))
		#Set heatmap arguments
		genes = dbgap.query[,2]
		palette = colorRampPalette(c("green", "black", "red"))(n=100)
		scalelow = input$range[1]
		scalehigh = input$range[2]
		#Make matrix
		mat.dbq = as.matrix(dbgap.query[,3:262])
		rownames(mat.dbq) = dbgap.query[,2]
		#Get rid of NA values
		mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later
		rownames(mat.dbq)
	})

	output$table2 <- renderTable({
		query <- datasetInput()
		#Filter dbgap with query
		dbgap.query = dbgap[which(dbgap[,2] %in% query[,1]),]
		#Reorder to match query order
		dbgap.query = dbgap.query[match(query[,1],dbgap.query[,2]), ]
		group = c(rep("green", 35), rep("orange", 78), rep("red", 132), rep("darkred", 15))
		#Set heatmap arguments
		genes = dbgap.query[,2]
		palette = colorRampPalette(c("green", "black", "red"))(n=100)
		scalelow = input$range[1]
		scalehigh = input$range[2]
		#Make matrix
		mat.dbq = as.matrix(dbgap.query[,3:262])
		rownames(mat.dbq) = dbgap.query[,2]
		#Get rid of NA values
		mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later
		subset(query[,1], !(query[,1] %in% rownames(mat.dbq)))
	})
	
	output$downloadHeatmap <- downloadHandler(
	filename <- function() {
		paste0(
			basename(unlist(strsplit(input$file1$name, split='.txt', fixed=TRUE))[1]), #Get name of file minus .txt
			'_SU2C-heatmap', 
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
			'_SU2C-heatmap', 
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
}