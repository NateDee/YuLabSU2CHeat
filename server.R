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
		#___________________________________________________________________________________
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
		#___________________________________________________________________________________
		#SU2C-Benign-vs-PCa
		if (input$database == "SU2C-Benign-vs-PCa") {
			dbgap.bp <- selectDB
		#If user put in a gene to sort by, use it to resort dbgap.bp:
			if (unlist(strsplit(input$GOI,","))[1] %in% dbgap.bp[,2]) {
				rowOfGene = which(dbgap.bp[,2] == unlist(strsplit(input$GOI,","))[1])
				dbgap.bp.tmp = dbgap.bp[,3:115]
				#Row 1 contains groupings!!!
				dbgap.bp.tmp = dbgap.bp.tmp[,order(dbgap.bp.tmp[1,], dbgap.bp.tmp[rowOfGene,])]
				dbgap.bp = cbind(dbgap.bp[,1:2], dbgap.bp.tmp)
				}
			#Filter dbgap.bp with query
			dbgap.bp.query = dbgap.bp[which(dbgap.bp[,2] %in% query[,1]),]
			#Reorder to match query order
			dbgap.bp.query = dbgap.bp.query[match(query[,1],dbgap.bp.query[,2]), ]
			#If correlation plot is wanted, find correlations and sort!
			if (correlate == TRUE & (unlist(strsplit(input$GOI,","))[1] %in% dbgap.bp[,2])) {
				genedf = dbgap.bp[which(dbgap.bp[,2] == unlist(strsplit(input$GOI,","))[1]),]
				t.dbgap.bp.query = t(as.matrix(dbgap.bp.query[,3:115]))
				colnames(t.dbgap.bp.query) = dbgap.bp.query[,2]
				t.genedf = t(as.matrix(genedf[,3:115]))
				colnames(t.genedf) = genedf[,2]
				correl.mat = t(cor(t.genedf,t.dbgap.bp.query))
				correl.mat = as.data.frame(cbind(correl.mat, correl.mat))
				correl.mat = correl.mat[order(correl.mat[,2]),]
				dbgap.bp.query = dbgap.bp.query[match(rownames(correl.mat),dbgap.bp.query[,2]),]
				}
			
			group = c(rep("green", 35), rep("orange", 78))
			#Set heatmap arguments
			genes = dbgap.bp.query[,2]
			palette = colorRampPalette(c(input$lowcolor, input$midcolor, input$highcolor))(n=100)
			scalelow = input$range[1]
			scalehigh = input$range[2]
			#Make matrix
			mat.dbq = as.matrix(dbgap.bp.query[,3:115])
			rownames(mat.dbq) = dbgap.bp.query[,2]
			#Get rid of NA values
			mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later
			cexRow = as.numeric(as.character(input$xfontsize))
			heatmap.2(mat.dbq, scale="none", trace="none", dendrogram="none", breaks=seq(scalelow, scalehigh, length.out=101), Rowv=FALSE, Colv=FALSE, col=palette, ColSideColors=group, keysize=0.75, key.par=list(cex=0.5), labCol=FALSE, cexRow=cexRow)
			}
		#___________________________________________________________________________________
		#SU2C-PCa-vs-CRPC
		if (input$database == "SU2C-PCa-vs-CRPC") {
			dbgap.pc <- selectDB
		#If user put in a gene to sort by, use it to resort dbgap.pc:
			if (unlist(strsplit(input$GOI,","))[1] %in% dbgap.pc[,2]) {
				rowOfGene = which(dbgap.pc[,2] == unlist(strsplit(input$GOI,","))[1])
				dbgap.pc.tmp = dbgap.pc[,3:212]
				#Row 1 contains groupings!!!
				dbgap.pc.tmp = dbgap.pc.tmp[,order(dbgap.pc.tmp[1,], dbgap.pc.tmp[rowOfGene,])]
				dbgap.pc = cbind(dbgap.pc[,1:2], dbgap.pc.tmp)
				}
			#Filter dbgap.pc with query
			dbgap.pc.query = dbgap.pc[which(dbgap.pc[,2] %in% query[,1]),]
			#Reorder to match query order
			dbgap.pc.query = dbgap.pc.query[match(query[,1],dbgap.pc.query[,2]), ]
			#If correlation plot is wanted, find correlations and sort!
			if (correlate == TRUE & (unlist(strsplit(input$GOI,","))[1] %in% dbgap.pc[,2])) {
				genedf = dbgap.pc[which(dbgap.pc[,2] == unlist(strsplit(input$GOI,","))[1]),]
				t.dbgap.pc.query = t(as.matrix(dbgap.pc.query[,3:212]))
				colnames(t.dbgap.pc.query) = dbgap.pc.query[,2]
				t.genedf = t(as.matrix(genedf[,3:212]))
				colnames(t.genedf) = genedf[,2]
				correl.mat = t(cor(t.genedf,t.dbgap.pc.query))
				correl.mat = as.data.frame(cbind(correl.mat, correl.mat))
				correl.mat = correl.mat[order(correl.mat[,2]),]
				dbgap.pc.query = dbgap.pc.query[match(rownames(correl.mat),dbgap.pc.query[,2]),]
				}
			
			group = c(rep("orange", 78), rep("red", 132))
			#Set heatmap arguments
			genes = dbgap.pc.query[,2]
			palette = colorRampPalette(c(input$lowcolor, input$midcolor, input$highcolor))(n=100)
			scalelow = input$range[1]
			scalehigh = input$range[2]
			#Make matrix
			mat.dbq = as.matrix(dbgap.pc.query[,3:212])
			rownames(mat.dbq) = dbgap.pc.query[,2]
			#Get rid of NA values
			mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later
			cexRow = as.numeric(as.character(input$xfontsize))
			heatmap.2(mat.dbq, scale="none", trace="none", dendrogram="none", breaks=seq(scalelow, scalehigh, length.out=101), Rowv=FALSE, Colv=FALSE, col=palette, ColSideColors=group, keysize=0.75, key.par=list(cex=0.5), labCol=FALSE, cexRow=cexRow)
			}
		#___________________________________________________________________________________
		#SU2C-CRPC-vs-NEPC
		if (input$database == "SU2C-CRPC-vs-NEPC") {
			dbgap.cn <- selectDB
		#If user put in a gene to sort by, use it to resort dbgap.cn:
			if (unlist(strsplit(input$GOI,","))[1] %in% dbgap.cn[,2]) {
				rowOfGene = which(dbgap.cn[,2] == unlist(strsplit(input$GOI,","))[1])
				dbgap.cn.tmp = dbgap.cn[,3:149]
				#Row 1 contains groupings!!!
				dbgap.cn.tmp = dbgap.cn.tmp[,order(dbgap.cn.tmp[1,], dbgap.cn.tmp[rowOfGene,])]
				dbgap.cn = cbind(dbgap.cn[,1:2], dbgap.cn.tmp)
				}
			#Filter dbgap.cn with query
			dbgap.cn.query = dbgap.cn[which(dbgap.cn[,2] %in% query[,1]),]
			#Reorder to match query order
			dbgap.cn.query = dbgap.cn.query[match(query[,1],dbgap.cn.query[,2]), ]
			#If correlation plot is wanted, find correlations and sort!
			if (correlate == TRUE & (unlist(strsplit(input$GOI,","))[1] %in% dbgap.cn[,2])) {
				genedf = dbgap.cn[which(dbgap.cn[,2] == unlist(strsplit(input$GOI,","))[1]),]
				t.dbgap.cn.query = t(as.matrix(dbgap.cn.query[,3:149]))
				colnames(t.dbgap.cn.query) = dbgap.cn.query[,2]
				t.genedf = t(as.matrix(genedf[,3:149]))
				colnames(t.genedf) = genedf[,2]
				correl.mat = t(cor(t.genedf,t.dbgap.cn.query))
				correl.mat = as.data.frame(cbind(correl.mat, correl.mat))
				correl.mat = correl.mat[order(correl.mat[,2]),]
				dbgap.cn.query = dbgap.cn.query[match(rownames(correl.mat),dbgap.cn.query[,2]),]
				}
			
			group = c(rep("red", 132), rep("darkred", 15))
			#Set heatmap arguments
			genes = dbgap.cn.query[,2]
			palette = colorRampPalette(c(input$lowcolor, input$midcolor, input$highcolor))(n=100)
			scalelow = input$range[1]
			scalehigh = input$range[2]
			#Make matrix
			mat.dbq = as.matrix(dbgap.cn.query[,3:149])
			rownames(mat.dbq) = dbgap.cn.query[,2]
			#Get rid of NA values
			mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later
			cexRow = as.numeric(as.character(input$xfontsize))
			heatmap.2(mat.dbq, scale="none", trace="none", dendrogram="none", breaks=seq(scalelow, scalehigh, length.out=101), Rowv=FALSE, Colv=FALSE, col=palette, ColSideColors=group, keysize=0.75, key.par=list(cex=0.5), labCol=FALSE, cexRow=cexRow)
			}
		#___________________________________________________________________________________
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
		
		#___________________________________________________________________________________
		#SU2C-Benign-vs-PCa
		if (input$database == "SU2C-Benign-vs-PCa") {
			dbgap.bp <- selectDB
		#If user put in a gene to sort by, use it to resort dbgap.bp:
			if (GOI[1] %in% dbgap.bp[,2]) {
				rowOfGene = which(dbgap.bp[,2] == GOI[1])
				dbgap.bp.tmp = dbgap.bp[,3:115]
				#Row 1 contains groupings!!!
				dbgap.bp.tmp = dbgap.bp.tmp[,order(dbgap.bp.tmp[1,], dbgap.bp.tmp[rowOfGene,])]
				dbgap.bp = cbind(dbgap.bp[,1:2], dbgap.bp.tmp)
				}
			#Filter dbgap.bp with query
			dbgap.bp.query = dbgap.bp[which(dbgap.bp[,2] %in% GOI),]
			if (length(GOI) > 1) { 
				dbgap.bp.query = dbgap.bp.query[order(factor(dbgap.bp.query$Gene_name, levels=GOI)),]
				}
			group = c(rep("green", 35), rep("orange", 78))
			palette = colorRampPalette(c(input$lowcolor, input$midcolor, input$highcolor))(n=100)
			scalelow = input$range[1]
			scalehigh = input$range[2]
			#Make matrix double if only one input gene
			if (length(GOI) == 1) {
				mat.dbq = as.matrix(rbind(dbgap.bp.query[,3:115],dbgap.bp.query[,3:115]), )
				rownames(mat.dbq) = c(GOI, GOI)
				} else {
				mat.dbq = as.matrix(dbgap.bp.query[,3:115])
				rownames(mat.dbq) = dbgap.bp.query[,2]
				}
			#Get rid of NA values
			#mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later,****added, see table output
			cexRow = as.numeric(as.character(input$xfontsize))
			heatmap.2(mat.dbq, scale="none", trace="none", dendrogram="none", breaks=seq(scalelow, scalehigh, length.out=101), Rowv=FALSE, Colv=FALSE, col=palette, ColSideColors=group, keysize=0.75, key.par=list(cex=0.5), labCol=FALSE, cexRow=cexRow)
			}
		#___________________________________________________________________________________
		#SU2C-PCa-vs-CRPC
		if (input$database == "SU2C-PCa-vs-CRPC") {
			dbgap.pc <- selectDB
			#If user put in a gene to sort by, use it to resort dbgap.pc:
			if (GOI[1] %in% dbgap.pc[,2]) {
				rowOfGene = which(dbgap.pc[,2] == GOI[1])
				dbgap.pc.tmp = dbgap.pc[,3:212]
				dbgap.pc.tmp = dbgap.pc.tmp[,order(dbgap.pc.tmp[1,], dbgap.pc.tmp[rowOfGene,])]
				dbgap.pc = cbind(dbgap.pc[,1:2], dbgap.pc.tmp)
				#dbgap.pc[,3:212] <- dbgap.pc[order(dbgap.pc[1, 3:212], dbgap.pc[rowOfGene,3:212]), 3:212]
				}
			#Filter dbgap.pc with query
			dbgap.pc.query = dbgap.pc[which(dbgap.pc[,2] %in% GOI),]
			if (length(GOI) > 1) { 
				dbgap.pc.query = dbgap.pc.query[order(factor(dbgap.pc.query$Gene_name, levels=GOI)),]
				}
			group = c(rep("orange", 78), rep("red", 132))
			palette = colorRampPalette(c(input$lowcolor, input$midcolor, input$highcolor))(n=100)
			scalelow = input$range[1]
			scalehigh = input$range[2]
			#Make matrix double if only one input gene
			if (length(GOI) == 1) {
				mat.dbq = as.matrix(rbind(dbgap.pc.query[,3:212],dbgap.pc.query[,3:212]), )
				rownames(mat.dbq) = c(GOI, GOI)
				} else {
				mat.dbq = as.matrix(dbgap.pc.query[,3:212])
				rownames(mat.dbq) = dbgap.pc.query[,2]
				}
			#Get rid of NA values
			#mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later,****added, see table output
			cexRow = as.numeric(as.character(input$xfontsize))
			heatmap.2(mat.dbq, scale="none", trace="none", dendrogram="none", breaks=seq(scalelow, scalehigh, length.out=101), Rowv=FALSE, Colv=FALSE, col=palette, ColSideColors=group, keysize=0.75, key.par=list(cex=0.5), labCol=FALSE, cexRow=cexRow)
			}
		#___________________________________________________________________________________
		#SU2C-CRPC-vs-NEPC
		if (input$database == "SU2C-CRPC-vs-NEPC") {
			dbgap.cn <- selectDB
			#If user put in a gene to sort by, use it to resort dbgap.cn:
			if (GOI[1] %in% dbgap.cn[,2]) {
				rowOfGene = which(dbgap.cn[,2] == GOI[1])
				dbgap.cn.tmp = dbgap.cn[,3:149]
				dbgap.cn.tmp = dbgap.cn.tmp[,order(dbgap.cn.tmp[1,], dbgap.cn.tmp[rowOfGene,])]
				dbgap.cn = cbind(dbgap.cn[,1:2], dbgap.cn.tmp)
				#dbgap.cn[,3:149] <- dbgap.cn[order(dbgap.cn[1, 3:149], dbgap.cn[rowOfGene,3:149]), 3:149]
				}
			#Filter dbgap.cn with query
			dbgap.cn.query = dbgap.cn[which(dbgap.cn[,2] %in% GOI),]
			if (length(GOI) > 1) { 
				dbgap.cn.query = dbgap.cn.query[order(factor(dbgap.cn.query$Gene_name, levels=GOI)),]
				}
			group = c(rep("red", 132), rep("darkred", 15))
			palette = colorRampPalette(c(input$lowcolor, input$midcolor, input$highcolor))(n=100)
			scalelow = input$range[1]
			scalehigh = input$range[2]
			#Make matrix double if only one input gene
			if (length(GOI) == 1) {
				mat.dbq = as.matrix(rbind(dbgap.cn.query[,3:149],dbgap.cn.query[,3:149]), )
				rownames(mat.dbq) = c(GOI, GOI)
				} else {
				mat.dbq = as.matrix(dbgap.cn.query[,3:149])
				rownames(mat.dbq) = dbgap.cn.query[,2]
				}
			#Get rid of NA values
			#mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later,****added, see table output
			cexRow = as.numeric(as.character(input$xfontsize))
			heatmap.2(mat.dbq, scale="none", trace="none", dendrogram="none", breaks=seq(scalelow, scalehigh, length.out=101), Rowv=FALSE, Colv=FALSE, col=palette, ColSideColors=group, keysize=0.75, key.par=list(cex=0.5), labCol=FALSE, cexRow=cexRow)
			}
		#___________________________________________________________________________________
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
		#___________________________________________________________________________________
		#SU2C-Benign-vs-PCa
		if (input$database == "SU2C-Benign-vs-PCa") {
			dbgap.bp = selectDB
			#Filter dbgap.bp with query
			dbgap.bp.query = dbgap.bp[which(dbgap.bp[,2] %in% query[,1]),]
			#Reorder to match query order
			dbgap.bp.query = dbgap.bp.query[match(query[,1],dbgap.bp.query[,2]), ]
			#If correlation plot is wanted, find correlations and sort!
			if ((input$GOI != "") & (correlate == TRUE) & (unlist(strsplit(input$GOI,","))[1] %in% dbgap.bp[,2])) {
				genedf = dbgap.bp[which(dbgap.bp[,2] == unlist(strsplit(input$GOI,","))[1]),]
				t.dbgap.bp.query = t(as.matrix(dbgap.bp.query[,3:115]))
				colnames(t.dbgap.bp.query) = dbgap.bp.query[,2]
				t.genedf = t(as.matrix(genedf[,3:115]))
				colnames(t.genedf) = genedf[,2]
				correl.mat = t(cor(t.genedf,t.dbgap.bp.query))
				correl.mat = as.data.frame(cbind(correl.mat, correl.mat))
				correl.mat = correl.mat[order(correl.mat[,2]),]
				dbgap.bp.query = dbgap.bp.query[match(rownames(correl.mat),dbgap.bp.query[,2]),]
				return(cbind(genes=rownames(correl.mat),Pearson.R=correl.mat[,1]))
				}
			#Set heatmap arguments
			genes = dbgap.bp.query[,2]
			#Make matrix
			mat.dbq = as.matrix(dbgap.bp.query[,3:115])
			rownames(mat.dbq) = dbgap.bp.query[,2]
			#Get rid of NA values
			mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later
			if (correlate == FALSE) {return(cbind(Genes=rownames(mat.dbq)))}
			}
		#___________________________________________________________________________________
		#SU2C-PCa-vs-CRPC
		if (input$database == "SU2C-PCa-vs-CRPC") {
			dbgap.pc = selectDB
			#Filter dbgap.pc with query
			dbgap.pc.query = dbgap.pc[which(dbgap.pc[,2] %in% query[,1]),]
			#Reorder to match query order
			dbgap.pc.query = dbgap.pc.query[match(query[,1],dbgap.pc.query[,2]), ]
			#If correlation plot is wanted, find correlations and sort!
			if ((input$GOI != "") & (correlate == TRUE) & (unlist(strsplit(input$GOI,","))[1] %in% dbgap.pc[,2])) {
				genedf = dbgap.pc[which(dbgap.pc[,2] == unlist(strsplit(input$GOI,","))[1]),]
				t.dbgap.pc.query = t(as.matrix(dbgap.pc.query[,3:212]))
				colnames(t.dbgap.pc.query) = dbgap.pc.query[,2]
				t.genedf = t(as.matrix(genedf[,3:212]))
				colnames(t.genedf) = genedf[,2]
				correl.mat = t(cor(t.genedf,t.dbgap.pc.query))
				correl.mat = as.data.frame(cbind(correl.mat, correl.mat))
				correl.mat = correl.mat[order(correl.mat[,2]),]
				dbgap.pc.query = dbgap.pc.query[match(rownames(correl.mat),dbgap.pc.query[,2]),]
				return(cbind(genes=rownames(correl.mat),Pearson.R=correl.mat[,1]))
				}
			#Set heatmap arguments
			genes = dbgap.pc.query[,2]
			#Make matrix
			mat.dbq = as.matrix(dbgap.pc.query[,3:212])
			rownames(mat.dbq) = dbgap.pc.query[,2]
			#Get rid of NA values
			mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later
			if (correlate == FALSE) {return(cbind(Genes=rownames(mat.dbq)))}
			}
		#___________________________________________________________________________________
		#SU2C-CRPC-vs-NEPC
		if (input$database == "SU2C-CRPC-vs-NEPC") {
			dbgap.cn = selectDB
			#Filter dbgap.cn with query
			dbgap.cn.query = dbgap.cn[which(dbgap.cn[,2] %in% query[,1]),]
			#Reorder to match query order
			dbgap.cn.query = dbgap.cn.query[match(query[,1],dbgap.cn.query[,2]), ]
			#If correlation plot is wanted, find correlations and sort!
			if ((input$GOI != "") & (correlate == TRUE) & (unlist(strsplit(input$GOI,","))[1] %in% dbgap.cn[,2])) {
				genedf = dbgap.cn[which(dbgap.cn[,2] == unlist(strsplit(input$GOI,","))[1]),]
				t.dbgap.cn.query = t(as.matrix(dbgap.cn.query[,3:149]))
				colnames(t.dbgap.cn.query) = dbgap.cn.query[,2]
				t.genedf = t(as.matrix(genedf[,3:149]))
				colnames(t.genedf) = genedf[,2]
				correl.mat = t(cor(t.genedf,t.dbgap.cn.query))
				correl.mat = as.data.frame(cbind(correl.mat, correl.mat))
				correl.mat = correl.mat[order(correl.mat[,2]),]
				dbgap.cn.query = dbgap.cn.query[match(rownames(correl.mat),dbgap.cn.query[,2]),]
				return(cbind(genes=rownames(correl.mat),Pearson.R=correl.mat[,1]))
				}
			#Set heatmap arguments
			genes = dbgap.cn.query[,2]
			#Make matrix
			mat.dbq = as.matrix(dbgap.cn.query[,3:149])
			rownames(mat.dbq) = dbgap.cn.query[,2]
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
		#If SU2C, use SU2C-Benign-vs-PCa dbgap.bp:
		if (input$database == "SU2C-Benign-vs-PCa") {
			dbgap.bp = selectDB
			#Filter dbgap.bp with query
			dbgap.bp.query = dbgap.bp[which(dbgap.bp[,2] %in% query[,1]),]
			#Reorder to match query order
			dbgap.bp.query = dbgap.bp.query[match(query[,1],dbgap.bp.query[,2]), ]
			#Set heatmap arguments
			genes = dbgap.bp.query[,2]
			#Make matrix
			mat.dbq = as.matrix(dbgap.bp.query[,3:115])
			rownames(mat.dbq) = dbgap.bp.query[,2]
			#Get rid of NA values
			mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later
			table = subset(query[,1], !(query[,1] %in% rownames(mat.dbq)))
			}
		#If SU2C, use SU2C-PCa-vs-CRPC dbgap.pc:
		if (input$database == "SU2C-PCa-vs-CRPC") {
			dbgap.pc = selectDB
			#Filter dbgap.pc with query
			dbgap.pc.query = dbgap.pc[which(dbgap.pc[,2] %in% query[,1]),]
			#Reorder to match query order
			dbgap.pc.query = dbgap.pc.query[match(query[,1],dbgap.pc.query[,2]), ]
			#Set heatmap arguments
			genes = dbgap.pc.query[,2]
			#Make matrix
			mat.dbq = as.matrix(dbgap.pc.query[,3:212])
			rownames(mat.dbq) = dbgap.pc.query[,2]
			#Get rid of NA values
			mat.dbq = mat.dbq[complete.cases(mat.dbq), ]  ###Add something to show NA values later
			table = subset(query[,1], !(query[,1] %in% rownames(mat.dbq)))
			}
		#If SU2C, use SU2C-CRPC-vs-NEPC dbgap.cn:
		if (input$database == "SU2C-CRPC-vs-NEPC") {
			dbgap.cn = selectDB
			#Filter dbgap.cn with query
			dbgap.cn.query = dbgap.cn[which(dbgap.cn[,2] %in% query[,1]),]
			#Reorder to match query order
			dbgap.cn.query = dbgap.cn.query[match(query[,1],dbgap.cn.query[,2]), ]
			#Set heatmap arguments
			genes = dbgap.cn.query[,2]
			#Make matrix
			mat.dbq = as.matrix(dbgap.cn.query[,3:149])
			rownames(mat.dbq) = dbgap.cn.query[,2]
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
				#Add sample types to output
				group = c(rep("Group", 2), rep("Benign", 35), rep("Primary-PCa", 78), rep("CRPC", 132), rep("NE-CRPC", 15))
				names = colnames(dbgap.query)
				dbgap.query = rbind(group, dbgap.query)
				colnames(dbgap.query) = names
				#Write table to output
				write.table(dbgap.query, file, sep="\t", row.names=FALSE)
				}
			#___________________________________________________________________________________
			#If SU2C-Benign-vs-PCa use Benign vs. PCa:
			if (input$database == "SU2C-Benign-vs-PCa") {
				dbgap.bp = selectDB
				#If user put in a gene to sort by, use it to resort dbgap.bp:
				if (unlist(strsplit(input$GOI,","))[1] %in% dbgap.bp[,2]) {
					rowOfGene = which(dbgap.bp[,2] == unlist(strsplit(input$GOI,","))[1])
					dbgap.bp.tmp = dbgap.bp[,3:115]
					dbgap.bp.tmp = dbgap.bp.tmp[,order(dbgap.bp.tmp[1,], dbgap.bp.tmp[rowOfGene,])]
					dbgap.bp = cbind(dbgap.bp[,1:2], dbgap.bp.tmp)
					#dbgap.bp[,3:115] <- dbgap.bp[order(dbgap.bp[1, 3:115], dbgap.bp[rowOfGene,3:115]), 3:115]
					}
				#Filter dbgap.bp with query
				dbgap.bp.query = dbgap.bp[which(dbgap.bp[,2] %in% query[,1]),]
				#Reorder to match query order
				dbgap.bp.query = dbgap.bp.query[match(query[,1],dbgap.bp.query[,2]), ]
				#If correlation plot is wanted, find correlations and sort!
				if (correlate == TRUE & (unlist(strsplit(input$GOI,","))[1] %in% dbgap.bp[,2])) {
					genedf = dbgap.bp[which(dbgap.bp[,2] == unlist(strsplit(input$GOI,","))[1]),]
					t.dbgap.bp.query = t(as.matrix(dbgap.bp.query[,3:115]))
					colnames(t.dbgap.bp.query) = dbgap.bp.query[,2]
					t.genedf = t(as.matrix(genedf[,3:115]))
					colnames(t.genedf) = genedf[,2]
					correl.mat = t(cor(t.genedf,t.dbgap.bp.query))
					correl.mat = as.data.frame(cbind(correl.mat, correl.mat))
					correl.mat = correl.mat[order(correl.mat[,2]),]
					dbgap.bp.query = dbgap.bp.query[match(rownames(correl.mat),dbgap.bp.query[,2]),]
					}
				#Add sample types to output
				group = c(rep("Group", 2), rep("Benign", 35), rep("Primary-PCa", 78))
				names = colnames(dbgap.bp.query)
				dbgap.bp.query = rbind(group, dbgap.bp.query)
				colnames(dbgap.bp.query) = names
				#Write table to output
				write.table(dbgap.bp.query, file, sep="\t", row.names=FALSE)
				}	
			#___________________________________________________________________________________
			#If SU2C-PCa-vs-CRPC:
			if (input$database == "SU2C-PCa-vs-CRPC") {
				dbgap.pc = selectDB
				#If user put in a gene to sort by, use it to resort dbgap.pc:
				if (unlist(strsplit(input$GOI,","))[1] %in% dbgap.pc[,2]) {
					rowOfGene = which(dbgap.pc[,2] == unlist(strsplit(input$GOI,","))[1])
					dbgap.pc.tmp = dbgap.pc[,3:212]
					dbgap.pc.tmp = dbgap.pc.tmp[,order(dbgap.pc.tmp[1,], dbgap.pc.tmp[rowOfGene,])]
					dbgap.pc = cbind(dbgap.pc[,1:2], dbgap.pc.tmp)
					#dbgap.pc[,3:212] <- dbgap.pc[order(dbgap.pc[1, 3:212], dbgap.pc[rowOfGene,3:212]), 3:212]
					}
				#Filter dbgap.pc with query
				dbgap.pc.query = dbgap.pc[which(dbgap.pc[,2] %in% query[,1]),]
				#Reorder to match query order
				dbgap.pc.query = dbgap.pc.query[match(query[,1],dbgap.pc.query[,2]), ]
				#If correlation plot is wanted, find correlations and sort!
				if (correlate == TRUE & (unlist(strsplit(input$GOI,","))[1] %in% dbgap.pc[,2])) {
					genedf = dbgap.pc[which(dbgap.pc[,2] == unlist(strsplit(input$GOI,","))[1]),]
					t.dbgap.pc.query = t(as.matrix(dbgap.pc.query[,3:212]))
					colnames(t.dbgap.pc.query) = dbgap.pc.query[,2]
					t.genedf = t(as.matrix(genedf[,3:212]))
					colnames(t.genedf) = genedf[,2]
					correl.mat = t(cor(t.genedf,t.dbgap.pc.query))
					correl.mat = as.data.frame(cbind(correl.mat, correl.mat))
					correl.mat = correl.mat[order(correl.mat[,2]),]
					dbgap.pc.query = dbgap.pc.query[match(rownames(correl.mat),dbgap.pc.query[,2]),]
					}
				#Add sample types to output
				group = c(rep("Group", 2),rep("Primary-PCa", 78), rep("CRPC", 132))
				names = colnames(dbgap.pc.query)
				dbgap.pc.query = rbind(group, dbgap.pc.query)
				colnames(dbgap.pc.query) = names
				#Write table to output
				write.table(dbgap.pc.query, file, sep="\t", row.names=FALSE)
				}
			#___________________________________________________________________________________
			#If SU2C-CRPC-vs-NEPC:
			if (input$database == "SU2C-CRPC-vs-NEPC") {
				dbgap.cn = selectDB
				#If user put in a gene to sort by, use it to resort dbgap.cn:
				if (unlist(strsplit(input$GOI,","))[1] %in% dbgap.cn[,2]) {
					rowOfGene = which(dbgap.cn[,2] == unlist(strsplit(input$GOI,","))[1])
					dbgap.cn.tmp = dbgap.cn[,3:149]
					dbgap.cn.tmp = dbgap.cn.tmp[,order(dbgap.cn.tmp[1,], dbgap.cn.tmp[rowOfGene,])]
					dbgap.cn = cbind(dbgap.cn[,1:2], dbgap.cn.tmp)
					#dbgap.cn[,3:149] <- dbgap.cn[order(dbgap.cn[1, 3:149], dbgap.cn[rowOfGene,3:149]), 3:149]
					}
				#Filter dbgap.cn with query
				dbgap.cn.query = dbgap.cn[which(dbgap.cn[,2] %in% query[,1]),]
				#Reorder to match query order
				dbgap.cn.query = dbgap.cn.query[match(query[,1],dbgap.cn.query[,2]), ]
				#If correlation plot is wanted, find correlations and sort!
				if (correlate == TRUE & (unlist(strsplit(input$GOI,","))[1] %in% dbgap.cn[,2])) {
					genedf = dbgap.cn[which(dbgap.cn[,2] == unlist(strsplit(input$GOI,","))[1]),]
					t.dbgap.cn.query = t(as.matrix(dbgap.cn.query[,3:149]))
					colnames(t.dbgap.cn.query) = dbgap.cn.query[,2]
					t.genedf = t(as.matrix(genedf[,3:149]))
					colnames(t.genedf) = genedf[,2]
					correl.mat = t(cor(t.genedf,t.dbgap.cn.query))
					correl.mat = as.data.frame(cbind(correl.mat, correl.mat))
					correl.mat = correl.mat[order(correl.mat[,2]),]
					dbgap.cn.query = dbgap.cn.query[match(rownames(correl.mat),dbgap.cn.query[,2]),]
					}
				#Add sample types to output
				group = c(rep("Group", 2), rep("CRPC", 132), rep("NE-CRPC", 15))
				names = colnames(dbgap.cn.query)
				dbgap.cn.query = rbind(group, dbgap.cn.query)
				colnames(dbgap.cn.query) = names
				#Write table to output
				write.table(dbgap.cn.query, file, sep="\t", row.names=FALSE)
				}
			#___________________________________________________________________________________
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
				group = c(rep("Group", 2), rep("Benign", 52), rep("Primary-PCa", 499))
				names = colnames(tcga.query)
				tcga.query = rbind(group, tcga.query)
				colnames(tcga.query) = names
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
					#Add sample types to output matrix
					group = c(rep("Group", 2), rep("Benign", 35), rep("Primary-PCa", 78), rep("CRPC", 132), rep("NE-CRPC", 15))
					names = colnames(dbgap.query)
					dbgap.query = rbind(group, dbgap.query)
					colnames(dbgap.query) = names
					write.table(dbgap.query, file, sep="\t", row.names=FALSE)
					}
				#If SU2C Benign vs. Cancer:
				if (input$database == "SU2C-Benign-vs-PCa") {
					dbgap.pc <- selectDB
					#If user put in a gene to sort by, use it to resort dbgap.pc:
					if (GOI[1] %in% dbgap.pc[,2]) {
						rowOfGene = which(dbgap.pc[,2] == GOI[1])
						dbgap.pc.tmp = dbgap.pc[,3:115]
						dbgap.pc.tmp = dbgap.pc.tmp[,order(dbgap.pc.tmp[1,], dbgap.pc.tmp[rowOfGene,])]
						dbgap.pc = cbind(dbgap.pc[,1:2], dbgap.pc.tmp)
						}
					#Filter dbgap.pc with query
					dbgap.pc.query = dbgap.pc[which(dbgap.pc[,2] %in% GOI),]
					if (length(GOI) > 1) { 
						dbgap.pc.query = dbgap.pc.query[order(factor(dbgap.pc.query$Gene_name, levels=GOI)),]
						}
					#Add sample types to output matrix
					group = c(rep("Group", 2), rep("Benign", 35), rep("Primary-PCa", 78))
					names = colnames(dbgap.pc.query)
					dbgap.pc.query = rbind(group, dbgap.pc.query)
					colnames(dbgap.pc.query) = names
					write.table(dbgap.pc.query, file, sep="\t", row.names=FALSE)
					}
				#If SU2C Cancer vs. CRPC:
				if (input$database == "SU2C-PCa-vs-CRPC") {
					dbgap.pc <- selectDB
					#If user put in a gene to sort by, use it to resort dbgap.pc:
					if (GOI[1] %in% dbgap.pc[,2]) {
						rowOfGene = which(dbgap.pc[,2] == GOI[1])
						dbgap.pc.tmp = dbgap.pc[,3:212]
						dbgap.pc.tmp = dbgap.pc.tmp[,order(dbgap.pc.tmp[1,], dbgap.pc.tmp[rowOfGene,])]
						dbgap.pc = cbind(dbgap.pc[,1:2], dbgap.pc.tmp)
						}
					#Filter dbgap.pc with query
					dbgap.pc.query = dbgap.pc[which(dbgap.pc[,2] %in% GOI),]
					if (length(GOI) > 1) { 
						dbgap.pc.query = dbgap.pc.query[order(factor(dbgap.pc.query$Gene_name, levels=GOI)),]
						}
					#Add sample types to output matrix
					group = c(rep("Group", 2), rep("Primary-PCa", 78), rep("CRPC", 132))
					names = colnames(dbgap.pc.query)
					dbgap.pc.query = rbind(group, dbgap.pc.query)
					colnames(dbgap.pc.query) = names
					write.table(dbgap.pc.query, file, sep="\t", row.names=FALSE)
					}
				#If SU2C-CRPC-vs-NEPC:
				if (input$database == "SU2C-CRPC-vs-NEPC") {
					dbgap.cn <- selectDB
					#If user put in a gene to sort by, use it to resort dbgap.cn:
					if (GOI[1] %in% dbgap.cn[,2]) {
						rowOfGene = which(dbgap.cn[,2] == GOI[1])
						dbgap.cn.tmp = dbgap.cn[,3:149]
						dbgap.cn.tmp = dbgap.cn.tmp[,order(dbgap.cn.tmp[1,], dbgap.cn.tmp[rowOfGene,])]
						dbgap.cn = cbind(dbgap.cn[,1:2], dbgap.cn.tmp)
						}
					#Filter dbgap.cn with query
					dbgap.cn.query = dbgap.cn[which(dbgap.cn[,2] %in% GOI),]
					if (length(GOI) > 1) { 
						dbgap.cn.query = dbgap.cn.query[order(factor(dbgap.cn.query$Gene_name, levels=GOI)),]
						}
					#Add sample types to output matrix
					group = c(rep("Group", 2), rep("CRPC", 132), rep("NE-CRPC", 15))
					names = colnames(dbgap.cn.query)
					dbgap.cn.query = rbind(group, dbgap.cn.query)
					colnames(dbgap.cn.query) = names
					write.table(dbgap.cn.query, file, sep="\t", row.names=FALSE)
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
					#Add sample types to output matrix
					group = c(rep("Group", 2), rep("Benign", 52), rep("Primary-PCa", 499))
					names = colnames(tcga.query)
					tcga.query = rbind(group, tcga.query)
					colnames(tcga.query) = names
					write.table(tcga.query, file, sep="\t", row.names=FALSE)
					}
				}
		)
}