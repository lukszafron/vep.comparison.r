#! /usr/bin/env Rscript
cat("Program path:", unlist(strsplit(grep(commandArgs(), pattern = "file=", value = T), split = "="))[2], "\n")

arguments2 <- commandArgs(trailingOnly = T)

if(length(arguments2) != 8) {stop("This program requires eight arguments. The first one is the destination folder, the second is the number of threads to use, the third should point to the CSV file used in the vep.r analysis, the fourth is the grouping variable, the fifth is the independent factor used, the sixth determines if the samples are paired, the seventh indicates the pair indicator column, while the eighth determines whether the FDR correction should be performed.")}
arguments2

library(doMC)
library(foreach)
library(dplyr)
library(matrixStats)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(pals)
library(openxlsx)
library(data.table)
library(pdftools)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(limma)
library(topGO)
library(ReactomePA)
library(clusterProfiler)

f.vep.comparison <- function(Robject1, Robject2) {
  f.topGO <- function(geneList, sigGenes) {
    topDiffGenes <- function(x) {x %in% sigGenes}
    colMap <- function(x) {
      .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
      return(.col[match(1:length(x), order(x))])
    }
    suffix <- paste(runid, group, groupid, "impact", sub(impact, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE"), sep = ".")
    wb <- createWorkbook()
    for(GOtype in c("BP", "CC", "MF")) { 
      GOdata.GOtype <- new("topGOdata",
                           ontology = GOtype,
                           allGenes = geneList,
                           geneSel = topDiffGenes,
                           annotationFun=annFUN.org, mapping = "org.Hs.eg.db", ID = "symbol",
                           nodeSize = 10)
      assign(paste("GOdata", GOtype, sep = "."), value = GOdata.GOtype)
      
      if(length(GOdata.GOtype@graph@nodes) >= 50) {nodesn.GOtype <- 50} else {nodesn.GOtype <- length(GOdata.GOtype@graph@nodes)}
      
      resultFisher.GOtype <- runTest(GOdata.GOtype, algorithm = "classic", statistic = "fisher")
      resultFisher.GOtype
      assign(paste("resultFisher", GOtype, sep = "."), value = resultFisher.GOtype)
      resultKS.GOtype <- runTest(GOdata.GOtype, algorithm = "classic", statistic = "ks")
      resultKS.GOtype
      assign(paste("resultKS", GOtype, sep = "."), value = resultKS.GOtype)
      resultKS.elim.GOtype <- tryCatch(runTest(GOdata.GOtype, algorithm = "elim", statistic = "ks"), error = function(e){NULL})
      resultKS.elim.GOtype
      assign(paste("resultKS.elim", GOtype, sep = "."), value = resultKS.elim.GOtype)
      if(!is.null(resultKS.elim.GOtype)) {
        allRes.GOtype.classicKS <- GenTable(GOdata.GOtype, classicFisher = resultFisher.GOtype,
                                            classicKS = resultKS.GOtype, elimKS = resultKS.elim.GOtype,
                                            orderBy = "classicKS", ranksOf = "classicFisher", topNodes = nodesn.GOtype)
        allRes.GOtype.elimKS <- GenTable(GOdata.GOtype, classicFisher = resultFisher.GOtype,
                                         classicKS = resultKS.GOtype, elimKS = resultKS.elim.GOtype,
                                         orderBy = "elimKS", ranksOf = "classicFisher", topNodes = nodesn.GOtype)
        allRes.GOtype.Fisher <- GenTable(GOdata.GOtype, classicFisher = resultFisher.GOtype,
                                         classicKS = resultKS.GOtype, elimKS = resultKS.elim.GOtype,
                                         orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = nodesn.GOtype)} else
                                         {
                                           allRes.GOtype.classicKS <- GenTable(GOdata.GOtype, classicFisher = resultFisher.GOtype,
                                                                               classicKS = resultKS.GOtype,
                                                                               orderBy = "classicKS", ranksOf = "classicFisher", topNodes = nodesn.GOtype)
                                           allRes.GOtype.Fisher <- GenTable(GOdata.GOtype, classicFisher = resultFisher.GOtype,
                                                                            classicKS = resultKS.GOtype,
                                                                            orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = nodesn.GOtype)
                                         }
      geneSig <- geneList[topDiffGenes(geneList)]
      
      AnnotatedGenes.GOtype.classicKS <- sapply(allRes.GOtype.classicKS$GO.ID, function(x) { as.character(unlist(genesInTerm(object = GOdata.GOtype, whichGO = x)))})
      genes_with_GO_term <- function(x) {AnnotatedGenes.GOtype.classicKS[[x]][AnnotatedGenes.GOtype.classicKS[[x]] %in% names(geneSig)]}
      GeneList.GOtype.classicKS <- sapply(allRes.GOtype.classicKS$GO.ID,genes_with_GO_term)
      
      sink(paste("GO_top50-significant_genes.", GOtype, ".classicKS",".(",Indfactor,").",suffix,".txt", sep = ""))

      print(GeneList.GOtype.classicKS)
      sink()
      if(exists("allRes.GOtype.elimKS")) {  
        AnnotatedGenes.GOtype.elimKS <- sapply(allRes.GOtype.elimKS$GO.ID, function(x) { as.character(unlist(genesInTerm(object = GOdata.GOtype, whichGO = x)))})
        genes_with_GO_term <- function(x) {AnnotatedGenes.GOtype.elimKS[[x]][AnnotatedGenes.GOtype.elimKS[[x]] %in% names(geneSig)]}
        GeneList.GOtype.elimKS <- sapply(allRes.GOtype.elimKS$GO.ID,genes_with_GO_term)

        sink(paste("GO_top50-significant_genes.", GOtype, ".elimKS",".(",Indfactor,").",suffix,".txt", sep = ""))

        print(GeneList.GOtype.elimKS)
        sink()
      }
      AnnotatedGenes.GOtype.Fisher <- sapply(allRes.GOtype.Fisher$GO.ID, function(x) { as.character(unlist(genesInTerm(object = GOdata.GOtype, whichGO = x)))})
      genes_with_GO_term <- function(x) {AnnotatedGenes.GOtype.Fisher[[x]][AnnotatedGenes.GOtype.Fisher[[x]] %in% names(geneSig)]}
      GeneList.GOtype.Fisher <- sapply(allRes.GOtype.Fisher$GO.ID,genes_with_GO_term)

      sink(paste("GO_top50-significant_genes.", GOtype, ".Fisher",".(",Indfactor,").",suffix,".txt", sep = ""))

      print(GeneList.GOtype.Fisher)
      sink()
      
      sheet.name <- paste("GO-top_50", GOtype, "classicKS", sep = ".")
      addWorksheet(wb = wb, sheetName = sheet.name)
      writeData(allRes.GOtype.classicKS, wb = wb, sheet = sheet.name, rowNames = F)
      
      sheet.name <- paste("GO-top_50", GOtype, "elimKS", sep = ".")
      addWorksheet(wb = wb, sheetName = sheet.name)
      
      if(exists("allRes.GOtype.elimKS")) {
        writeData(allRes.GOtype.elimKS, wb = wb, sheet = sheet.name, rowNames = F)
      }    
      sheet.name <- paste("GO-top_50", GOtype, "Fisher", sep = ".")
      addWorksheet(wb = wb, sheetName = sheet.name)
      writeData(allRes.GOtype.Fisher, wb = wb, sheet = sheet.name, rowNames = F)
      suppressWarnings(rm(resultKS.elim.GOtype, allRes.GOtype.elimKS))
    }

      saveWorkbook(wb, file = paste("GO-top_50",".(",Indfactor,").",suffix,".xlsx", sep = ""), overwrite = T)

      pdf(paste("GO-p-values_comparison",".(",Indfactor,").",suffix,".pdf", sep = ""))

    for(GOtype in c("BP", "CC", "MF")) {
      if(!is.null(paste("resultKS.elim", GOtype, sep = "."))) {
        pValue.classic <- score(get(paste("resultKS", GOtype, sep = ".")))
        pValue.elim <- score(get(paste("resultKS.elim", GOtype, sep = ".")))[names(pValue.classic)]
        gstat <- termStat(get(paste("GOdata", GOtype, sep = ".")), names(pValue.classic))
        gSize <- gstat$Annotated / max(gstat$Annotated) * 4
        gCol <- colMap(gstat$Significant)
        plot(x=pValue.classic, y=pValue.elim, xlab = "p-value - KS test (classic)", ylab = "p-value - KS test (elim)", pch = 19, cex = gSize, col = gCol, title(main=paste("Gene ontology analysis", GOtype, sep = " - ")))
      } else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene ontology", GOtype, "resultKS.elim ERROR", sep = " - "))}
    }
    dev.off()

    pdf(paste("GO-diagrams_top_10_test_KS_classic",".(",Indfactor,").",suffix,".pdf", sep = ""))

    for(GOtype in c("BP", "CC", "MF")) {
      tryCatch(expr={showSigOfNodes(get(paste("GOdata", GOtype, sep = ".")), score(get(paste("resultKS",GOtype, sep = "."))), firstSigNodes = 10, useInfo ='all'); title(main=paste("Gene ontology", GOtype, sep = " - "))},
               error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene ontology", GOtype, "ERROR", sep = " - "))}})
    }
    dev.off()
    
    pdf(paste("GO-diagrams_top_10_test_KS_elim",".(",Indfactor,").",suffix,".pdf", sep = ""))

    for(GOtype in c("BP", "CC", "MF")) {
      tryCatch(expr={showSigOfNodes(get(paste("GOdata", GOtype, sep = ".")), score(get(paste("resultKS.elim", GOtype, sep = "."))), firstSigNodes = 10, useInfo ='all'); title(main=paste("Gene ontology", GOtype, sep = " - "))},
               error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene ontology", GOtype, "ERROR", sep = " - "))}})
    }
    dev.off()
    
    pdf(paste("GO-diagrams_top_10_test_Fisher",".(",Indfactor,").",suffix,".pdf", sep = ""))

    for(GOtype in c("BP", "CC", "MF")) {
      tryCatch(expr={showSigOfNodes(get(paste("GOdata", GOtype, sep = ".")), score(get(paste("resultFisher", GOtype, sep = "."))), firstSigNodes = 10, useInfo ='all'); title(main=paste("Gene ontology", GOtype, sep = " - "))},
               error={function(e) {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene ontology", GOtype, "ERROR", sep = " - "))}})
    }
    dev.off()
  }
  f.Reactome <- function(group) {
    subgroups <- unlist(strsplit(group, split = "_vs_"))
    for(impact in impacts) {
      all.genes <- get(paste("allGenes.list", groupid, impact, sep = "."))[[group]]
      all.genes_g1 <- all.genes[grepl(names(all.genes), pattern = ".*#1$")]
      all.genes_g2 <- all.genes[grepl(names(all.genes), pattern = ".*#2$")]
      sig.genes.name <- paste("sigGenes.list", groupid, impact, sep = ".")
      if(! exists(sig.genes.name)) {next}
      sig.genes <- get(sig.genes.name)[[group]]
      sig.genes_g1 <- sig.genes[grepl(names(sig.genes), pattern = ".*#1$")]
      sig.genes_g2 <- sig.genes[grepl(names(sig.genes), pattern = ".*#2$")]
      all.genes.list <- c('all.genes_g1', 'all.genes_g2')
      sig.genes.list <- c('sig.genes_g1', 'sig.genes_g2')
      if(impact == "X.HIGH.MODERATE.") {impact <- "HIGH_or_MODERATE"}
      
      for(index in c(1,2)) {
        var.subgroup <- subgroups[index]
        all.genes <- get(all.genes.list[index])
        sig.genes <- get(sig.genes.list[index])
        if(length(all.genes) > 0 & length(sig.genes) > 0) {
          names(all.genes) <- sub(names(get(all.genes.list[index])), pattern = "#(1|2)$", replacement = "")
          names(sig.genes) <- sub(names(get(sig.genes.list[index])), pattern = "#(1|2)$", replacement = "")
          sig.genes[sig.genes == 0] <- 2.2e-16
          minus.log10.sig.genes <- -log10(sig.genes)
          names(minus.log10.sig.genes) <- tryCatch(expr = {names(minus.log10.sig.genes) <- foreach(gname = names(minus.log10.sig.genes), .combine = c) %do% {as.vector(tryCatch(expr=mapIds(org.Hs.eg.db, keys=gname, column = "ENTREZID", keytype="SYMBOL", multiVals="first"), error = function(e) {if(length(alias2Symbol(sub(gname, pattern = "-", replacement = ""))) > 0) {mapIds(org.Hs.eg.db, keys = alias2Symbol(sub(gname, pattern = "-", replacement = ""), species = "Hs"), keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")} else{NA}}))}},
                   error = {function(e) {ambiguous.gene.names.list <- foreach(gname = names(minus.log10.sig.genes)) %do% {as.vector(tryCatch(expr=mapIds(org.Hs.eg.db, keys=gname, column = "ENTREZID", keytype="SYMBOL", multiVals="first"), error = function(e) {if(length(alias2Symbol(sub(gname, pattern = "-", replacement = ""))) > 0) {mapIds(org.Hs.eg.db, keys = alias2Symbol(sub(gname, pattern = "-", replacement = ""), species = "Hs"), keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")} else{NA}}))}
                            names(ambiguous.gene.names.list) <- names(minus.log10.sig.genes)
                            ambiguous.gene.names.list <- ambiguous.gene.names.list[sapply(ambiguous.gene.names.list, length) > 1]
                            sink(file = "WARNING:Ambiguous.gene.names.list")
                            cat("INFO: For each gene symbol, the first ENTREZID was used in further analyses.\n\n")
                            print(ambiguous.gene.names.list)
                            sink()
                            foreach(gname = names(minus.log10.sig.genes), .combine = c) %do% {as.vector(tryCatch(expr=mapIds(org.Hs.eg.db, keys=gname, column = "ENTREZID", keytype="SYMBOL", multiVals="first"), error = function(e) {if(length(alias2Symbol(sub(gname, pattern = "-", replacement = ""))) > 0) {mapIds(org.Hs.eg.db, keys = alias2Symbol(sub(gname, pattern = "-", replacement = ""), species = "Hs")[1], keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")} else{NA}}))}}})
          minus.log10.sig.genes <- minus.log10.sig.genes[!is.na(names(minus.log10.sig.genes))]
          geneList <- minus.log10.sig.genes
          geneList <- geneList[!is.na(names(geneList))]
          geneList <- sort(geneList, decreasing = TRUE)
          geneList <- geneList[!duplicated(names(geneList))]
          geneList.df <- as.data.frame(geneList)
          colnames(geneList.df) <- "-log10(p-value)"
          if (nrow(geneList.df) > 0) {
            rownames(geneList.df) <- as.vector(mapIds(org.Hs.eg.db,
                                                      keys=rownames(geneList.df),
                                                      column="SYMBOL",
                                                      keytype="ENTREZID",
                                                      multiVals="first"))}
          geneList.df.name <- paste("geneList", var.subgroup, groupid, paste( "impact", sub(impact, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE"), sep = ":"), sep = ".")
          assign(geneList.df.name, value = geneList.df)
          
          de <- names(geneList)
          if (length(de) > 0) {
            de.name <- paste("de", var.subgroup, groupid, paste( "impact", sub(impact, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE"), sep = ":"), sep = ".")
            assign(de.name, value = de)
            
            temp.x <- enrichPathway(gene=get(de.name), organism = "human", pAdjustMethod = "BH", pvalueCutoff=0.05, readable=TRUE)
            x.name <- paste("x", var.subgroup, groupid, paste( "impact", sub(impact, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE"), sep = ":"), sep = ".")
            x.name.df <- paste("x", var.subgroup, groupid, paste( "impact", sub(impact, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE"), sep = ":"), "df", sep = ".")
            assign(x.name, value = temp.x)
            rm(temp.x)
            try(p1 <- dotplot(get(x.name), showCategory=50) + labs(title = paste("Pathway enrichment analysis", paste(var.subgroup, "up-altered"), sep = "-")) + theme(plot.title = element_text(hjust = 0.5, face = "bold")))
            p1.name <- paste('p1', var.subgroup, groupid, paste( "impact", sub(impact, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE"), sep = ":"), sep = ".")
            
            if(exists('p1')) {if(nrow(p1$data) > 0) {assign(p1.name, value = p1)
              p1.pathways <- as.character(p1$data[order(p1$data$GeneRatio, decreasing = T),][["Description"]][1:if(nrow(p1$data)<3){nrow(p1$data)} else {3}])
              if(index == 1) {
                pdf(paste("Reactome analysis", var.subgroup, groupid, paste("impact", impact, sep = ":"), "01#1", "pdf", sep = "."), height = 10, width = 13)} else if(index == 2) {
                  pdf(paste("Reactome analysis", var.subgroup, groupid, paste("impact", impact, sep = ":"), "07#1", "pdf", sep = "."), height = 10, width = 13)}
              
              for(i in p1.pathways) {try(expr = {pp1 <- viewPathway(i)
              pp1 <- pp1 + labs(title = paste(paste(var.subgroup, "up-altered"), i, sep = ": ")) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
              print(pp1)})}
              dev.off()
            }}
            try(p2 <- heatplot(get(x.name), showCategory = 50, foldChange = geneList) + labs(title = paste("Pathway enrichment heatmap", paste(var.subgroup, "up-altered"), sep = "-"), subtitle = "(Info: Fold change values correspond to the -log10(p-values)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid.major.x  = element_line(color = "grey", linetype = "solid", size = 0.5), panel.grid.major.y  = element_line(color = "grey", linetype = "solid", size = 0.5)))
            p2.name <- paste('p2', var.subgroup, groupid, paste("impact", sub(impact, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE"), sep = ":"), sep = ".")
            if(exists('p2')) {assign(p2.name, value = p2)}
            try(p3 <- emapplot(get(x.name), showCategory = 50, color = "p.adjust") + labs(title = paste("Pathway enrichment map", paste(var.subgroup, "up-altered"), sep = "-")) + theme(plot.title = element_text(hjust = 0.5, face = "bold")))
            p3.name <- paste('p3', var.subgroup, groupid, paste("impact", sub(impact, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE"), sep = ":"), sep = ".")
            if(exists('p3')) {assign(p3.name, value = p3)}
            try(p4 <- cnetplot(get(x.name), categorySize="qvalue", foldChange=geneList) + labs(title = paste("Complex associations map", paste(var.subgroup, "up-altered"), sep = "-"), subtitle = "(Info: Fold change values correspond to the -log10(p-values)") + theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5)))
            p4.name <- paste('p4', var.subgroup, groupid, paste( "impact", sub(impact, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE"), sep = ":"), sep = ".")
            if(exists('p4')) {assign(p4.name, value = p4)}
            
            tmp.df <- as.data.frame(get(x.name))
            assign(x.name.df, value = tmp.df)
            
            try({y <- gsePathway(geneList, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "BH", by = "fgsea")
            res.fgsea <- as.data.frame(y)
            p5 <- emapplot(y, showCategory = 50, color = "p.adjust") + labs(title = paste("Gene set enrichment analysis", paste(var.subgroup, "up-altered"), sep = "-")) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))})
            if(exists("y")) { if(nrow(y) > 0) {
              for(i in seq(1,length(y@result$core_enrichment))) {y@result$core_enrichment[i] <- paste(as.character(mapIds(org.Hs.eg.db, keys = unlist(strsplit(y@result$core_enrichment[i], split = "/")), keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")),collapse = "/")}}}
            p5.name <- paste('p5', var.subgroup, groupid, paste( "impact", sub(impact, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE"), sep = ":"), sep = ".")
            if(exists('p5')) {assign(p5.name, value = p5)}
            
            y.name <- paste("y", var.subgroup, groupid, paste( "impact", sub(impact, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE"), sep = ":"), sep = ".")
            y.name.df <- paste("y", var.subgroup, groupid, paste( "impact", sub(impact, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE"), sep = ":"), "df", sep = ".")
            if(exists("y")) {assign(y.name, value = y)}
            if(exists("y")) {assign(y.name.df, value = as.data.frame(get(y.name)))}
            geneList.Symbols <- geneList
            if (length(geneList.Symbols) > 0) {
              names(geneList.Symbols) <- mapIds(org.Hs.eg.db, keys = names(geneList.Symbols), keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")}
            if(exists('y')) {if(nrow(get(y.name.df)) > 0) { 
              p6 <- heatplot(get(y.name), showCategory = 50, foldChange = geneList.Symbols) + labs(title = paste("Gene set enrichment heatmap", paste(var.subgroup, "up-altered"), sep = "-")) + theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.x = element_text(angle = 90, vjust = 0.5), panel.grid.major.x  = element_line(color = "grey", linetype = "solid", size = 0.5), panel.grid.major.y  = element_line(color = "grey", linetype = "solid", size = 0.5))}}
            p6.name <- paste('p6', var.subgroup, groupid, paste( "impact", sub(impact, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE"), sep = ":"), sep = ".")
            if(exists('p6')) {assign(p6.name, value = p6)}
            rm(p1, p2, p3, p4, p5, p6, y)
          }}}
      de.1 <- paste("de", subgroups[1], groupid, paste("impact", impact, sep = ":"), sep = ".")
      de.2 <- paste("de", subgroups[2], groupid, paste("impact", impact, sep = ":"), sep = ".")
      
      if(!exists(de.1)) {assign(de.1, value = NULL)}
      if(!exists(de.2)) {assign(de.2, value = NULL)}
      
      de.list.full <- do.call(mget, list(c(de.1, de.2)))
      names(de.list.full) <- sub(subgroups, pattern = "(.*)", replacement = "\\1-up_altered")
      
      try(compareClustersRes <- compareCluster(de.list.full, fun="enrichPathway", organism = "human", pAdjustMethod = "BH", pvalueCutoff=0.05, readable=TRUE))
      try(p7 <- dotplot(compareClustersRes, showCategory=50) + labs(title = paste("Pathway enrichment analysis - group comparison", sep = "-")) + theme(plot.title = element_text(hjust = 0.5, face = "bold")))
      p7.name <- paste('p7', "ZZZ", groupid, paste( "impact", sub(impact, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE"), sep = ":"), sep = ".")
      if(exists('p7')) {assign(p7.name, value = p7)}
      rm(p7)
      
      try(compareClustersRes.df <- as.data.frame(compareClustersRes))
      compareClustersRes.df.name <- paste("compareClustersRes.df", groupid, paste("impact", impact, sep = ":"), sep = ".")
      if(exists('compareClustersRes.df')) {
        assign(compareClustersRes.df.name, value = compareClustersRes.df)}
      rm(compareClustersRes, compareClustersRes.df)
      
      suffix <- paste(subgroups[1], groupid, paste("impact", impact, sep = ":"), sep = ".")
      if (exists(paste("p1", suffix, sep = "."))){
        pdf(paste("Reactome analysis", suffix, "01", "pdf", sep = "."), height = 10, width = 13)
        print(get(paste("p1", suffix, sep = ".")))
        dev.off()}
      if(exists(paste("p2", suffix, sep = "."))) { if(nrow(get(paste("p2", suffix, sep = "."))$data) > 260){
        pdf(paste("Reactome analysis", suffix, "02", "pdf", sep = "."), height = 10, width = nrow(get(paste("p2", suffix, sep = "."))$data)/20)} else {
          pdf(paste("Reactome analysis", suffix, "02", "pdf", sep = "."), height = 10, width = 13)
        }} else {
          pdf(paste("Reactome analysis", suffix, "02", "pdf", sep = "."), height = 10, width = 10)}
      if(exists(paste("p2", suffix, sep = "."))){
        print(get(paste("p2", suffix, sep = ".")))} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Pathway enrichment heatmap: not enough terms to draw a map", paste("subgroup", suffix, sep = ":"), sep = "-"))}
      dev.off()
      
      pdf(paste("Reactome analysis", suffix, "03", "pdf", sep = "."), height = 13, width = 13)
      if(exists(paste("p3", suffix, sep = "."))) {
        print(get(paste("p3", suffix, sep = ".")))} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Pathway enrichment map: not enough terms to draw a map", paste("subgroup", suffix, sep = ":"), sep = "-"))}
      dev.off()
      
      pdf(paste("Reactome analysis", suffix, "04", "pdf", sep = "."), height = 13, width = 13)
      if(exists(paste("p4", suffix, sep = "."))) {
        print(get(paste("p4", suffix, sep = ".")))} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Complex assocations map: not enough terms to draw a map", paste("subgroup", suffix, sep = ":"), sep = "-"))}
      dev.off()
      
      pdf(paste("Reactome analysis", suffix, "05", "pdf", sep = "."), height = 13, width = 13)
      if(exists(paste("p5", suffix, sep = "."))) {
        print(get(paste("p5", suffix, sep = ".")))} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene set enrichment analysis: no term enriched under specific pvalueCutoff", paste("subgroup", suffix, sep = ":"), sep = "-"))}
      dev.off()
      
      if(exists(paste("p6", suffix, sep = "."))) { if (nrow(get(paste("p6", suffix, sep = "."))$data) > 260){
        pdf(paste("Reactome analysis", suffix, "06", "pdf", sep = "."), height = 10, width = nrow(get(paste("p6", suffix, sep = "."))$data)/20)} else {
          pdf(paste("Reactome analysis", suffix, "06", "pdf", sep = "."), height = 10, width = 13)
        }} else {
          pdf(paste("Reactome analysis", suffix, "06", "pdf", sep = "."), height = 10, width = 10)}
      if(exists(paste("p6", suffix, sep = "."))) {
        print(get(paste("p6", suffix, sep = ".")))} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene set enrichment heatmap: not enough terms to draw a map", paste("subgroup", suffix, sep = ":"), sep = "-"))}
      dev.off()
      
      suffix <- paste(subgroups[2], groupid, paste("impact", impact, sep = ":"), sep = ".")
      if (exists(paste("p1", suffix, sep = "."))){
        pdf(paste("Reactome analysis", suffix, "07", "pdf", sep = "."), height = 10, width = 13)
        print(get(paste("p1", suffix, sep = ".")))
        dev.off()}
      if(exists(paste("p2", suffix, sep = "."))) { if(nrow(get(paste("p2", suffix, sep = "."))$data) > 260){
        pdf(paste("Reactome analysis", suffix, "08", "pdf", sep = "."), height = 10, width = nrow(get(paste("p2", suffix, sep = "."))$data)/20)} else {
          pdf(paste("Reactome analysis", suffix, "08", "pdf", sep = "."), height = 10, width = 13)
        }} else {
          pdf(paste("Reactome analysis", suffix, "08", "pdf", sep = "."), height = 10, width = 10)}
      if(exists(paste("p2", suffix, sep = "."))){
        print(get(paste("p2", suffix, sep = ".")))} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Pathway enrichment heatmap: not enough terms to draw a map", paste("subgroup", suffix, sep = ":"), sep = "-"))}
      dev.off()
      
      pdf(paste("Reactome analysis", suffix, "09", "pdf", sep = "."), height = 13, width = 13)
      if(exists(paste("p3", suffix, sep = "."))) {
        print(get(paste("p3", suffix, sep = ".")))} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Pathway enrichment map: not enough terms to draw a map", paste("subgroup", suffix, sep = ":"), sep = "-"))}
      dev.off()
      
      pdf(paste("Reactome analysis", suffix, "10", "pdf", sep = "."), height = 13, width = 13)
      if(exists(paste("p4", suffix, sep = "."))) {
        print(get(paste("p4", suffix, sep = ".")))} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Complex assocations map: not enough terms to draw a map", paste("subgroup", suffix, sep = ":"), sep = "-"))}
      dev.off()
      
      pdf(paste("Reactome analysis", suffix, "11", "pdf", sep = "."), height = 13, width = 13)
      if(exists(paste("p5", suffix, sep = "."))) {
        print(get(paste("p5", suffix, sep = ".")))} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene set enrichment analysis: no term enriched under specific pvalueCutoff", paste("subgroup", suffix, sep = ":"), sep = "-"))}
      dev.off()
      
      if(exists(paste("p6", suffix, sep = "."))) { if (nrow(get(paste("p6", suffix, sep = "."))$data) > 260){
        pdf(paste("Reactome analysis", suffix, "12", "pdf", sep = "."), height = 10, width = nrow(get(paste("p6", suffix, sep = "."))$data)/20)} else {
          pdf(paste("Reactome analysis", suffix, "12", "pdf", sep = "."), height = 10, width = 13)
        }} else {
          pdf(paste("Reactome analysis", suffix, "12", "pdf", sep = "."), height = 10, width = 10)}
      if(exists(paste("p6", suffix, sep = "."))) {
        print(get(paste("p6", suffix, sep = ".")))} else {plot.new() + plot.window(xlim=c(-5,5), ylim=c(-5,5)); title(main=paste("Gene set enrichment heatmap: not enough terms to draw a map", paste("subgroup", suffix, sep = ":"), sep = "-"))}
      dev.off()      
      
      suffix <- paste("ZZZ", groupid, paste("impact", impact, sep = ":"), sep = ".")
      if(exists(paste("p7", suffix, sep = "."))) {
        if(nrow(get(paste("p7", suffix, sep = "."))$data)>60) {p7.height = floor(nrow(get(paste("p7", suffix, sep = "."))$data)/6)} else {p7.height = 10}
        pdf(paste("Reactome analysis", suffix, "13", "pdf", sep = "."), height = p7.height, width = 13)
        print(get(paste("p7", suffix, sep = ".")))
        dev.off()}
      
      pdffiles <- sort(list.files(pattern = paste("Reactome analysis", paste0("(", paste(c(subgroups, "ZZZ"), collapse = "|"),")"), make.names(groupid), paste("impact", impact, sep = ":"), "[0-9#]+", "pdf", sep = ".")))

      reactome.pdf.name = paste0("Reactome_analysis.", Indfactor, ".", group, ".", groupid, ".", paste("impact", impact, sep = ":"), ".", runid,".pdf", sep = "")

      pdf_combine(pdffiles, output = reactome.pdf.name)
      unlink(pdffiles)
      
      reactome.xlsx.name <- paste0("Reactome_analysis.", Indfactor, ".", group, ".", groupid, ".", paste("impact", impact, sep = ":"), ".", runid,".xlsx", sep = "")

      unlink(reactome.xlsx.name)
      
      wb <- createWorkbook()
      
      geneList.final <- paste("geneList", subgroups[1], groupid, paste("impact", impact, sep = ":"), sep = ".")
      if(exists(geneList.final)) {
        sheet.name <- paste(subgroups[1], "up-altered", sep = ".")
        sheet.name <- substr(sheet.name, 1, 31)
        addWorksheet(wb = wb, sheetName = sheet.name)
        writeData(x = get(geneList.final), wb = wb, sheet = sheet.name, rowNames = T)
      }
      x.df.final <- paste("x", subgroups[1], groupid, paste("impact", impact, sep = ":"), "df", sep = ".")
      if(exists(x.df.final)) {
        sheet.name <- paste("Path.enrich", subgroups[1], "up-altered", sep = ".")
        sheet.name <- substr(sheet.name, 1, 31)
        addWorksheet(wb = wb, sheetName = sheet.name)
        writeData(x = get(x.df.final), wb = wb, sheet = sheet.name, rowNames = F)
      }
      
      y.df.final <- paste("y", subgroups[1], groupid, paste("impact", impact, sep = ":"), "df", sep = ".")
      if(exists(y.df.final)) {
        sheet.name <- paste("Gene.set.enrich", subgroups[1], "up-altered", sep = ".")
        sheet.name <- substr(sheet.name, 1, 31)
        addWorksheet(wb = wb, sheetName = sheet.name)
        writeData(x = get(y.df.final), wb = wb, sheet = sheet.name, rowNames = F)
      }
      
      geneList.final <- paste("geneList", subgroups[2], groupid, paste("impact", impact, sep = ":"), sep = ".")
      if(exists(geneList.final)) {
        sheet.name <- paste(subgroups[2], "up-altered", sep = ".")
        sheet.name <- substr(sheet.name, 1, 31)
        addWorksheet(wb = wb, sheetName = sheet.name)
        writeData(x = get(geneList.final), wb = wb, sheet = sheet.name, rowNames = T)
      }
      x.df.final <- paste("x", subgroups[2], groupid, paste("impact", impact, sep = ":"), "df", sep = ".")
      if(exists(x.df.final)) {
        sheet.name <- paste("Path.enrich", subgroups[2], "up-altered", sep = ".")
        sheet.name <- substr(sheet.name, 1, 31)
        addWorksheet(wb = wb, sheetName = sheet.name)
        writeData(x = get(x.df.final), wb = wb, sheet = sheet.name, rowNames = F)
      }
      
      y.df.final <- paste("y", subgroups[2], groupid, paste("impact", impact, sep = ":"), "df", sep = ".")
      if(exists(y.df.final)) {
        sheet.name <- paste("Gene.set.enrich", subgroups[2], "up-altered", sep = ".")
        sheet.name <- substr(sheet.name, 1, 31)
        addWorksheet(wb = wb, sheetName = sheet.name)
        writeData(x = get(y.df.final), wb = wb, sheet = sheet.name, rowNames = F)
      }  
      
      compareClustersRes.df.final <- paste("compareClustersRes.df", groupid, paste("impact", impact, sep = ":"), sep = ".")
      if(exists(compareClustersRes.df.final)) {
        sheet.name <- "Path.enrich.group.comparison"
        addWorksheet(wb = wb, sheetName = sheet.name)
        writeData(x = get(compareClustersRes.df.final), wb = wb, sheet = sheet.name, rowNames = F)
      }
      
      saveWorkbook(wb = wb, file = reactome.xlsx.name, overwrite = T)
    }}

  if(all(sapply(c(Robject1, Robject2), file.exists))) {
  Robject1.list <<- load(Robject1, envir = .GlobalEnv)
  Robject2.list <<- load(Robject2, envir = .GlobalEnv)
  
  if(txt.file != "NA") {
  con <- file(txt.file)
  gene.list1 <- unique(gsub(readLines(con), pattern = "\\s", replacement = ""))
  gene.signature.name <- gsub(txt.file, pattern = "^.*\\/", replacement = "")
  close(con)

  gene.signature <- NULL
  for(i in gene.list1) {if(length(alias2Symbol(i, species = "Hs")) == 0) {gene.signature <- append(gene.signature, values = i)} else {gene.signature <- append(gene.signature, values = alias2Symbol(i, species = "Hs"))}}
  gene.signature <- unique(gene.signature)
  } else {
  gene.signature.name <- NULL
  }
  
  anno1 <- paste("anno", "SNP", groupid, "HIGH", sep = ".")
  anno2 <- paste("anno", "SNP", groupid, "X.HIGH.MODERATE.", sep = ".")
  anno3 <- paste("anno", "NON_SNP", groupid, "HIGH", sep = ".")
  anno4 <- paste("anno", "NON_SNP", groupid, "X.HIGH.MODERATE.", sep = ".")
  
  anno.checker <- function(annos) {
    anno.rownames <- foreach(anno = annos, .combine = cbind) %dopar% {
      rownames(get(anno))}
    row.values <- foreach(i = seq(1, nrow(anno.rownames)), .combine = c) %dopar% {
      length(unique(anno.rownames[i,1:ncol(anno.rownames)]))
    }
    all(row.values == 1)}
  if(! anno.checker(c(anno1, anno2, anno3, anno4))) {stop("Not all annotation files match.")}
  anno <- get(anno1)
  
  impacts <- c("HIGH", "X.HIGH.MODERATE.")
  
  needed.objects <- c(paste("res1.bool.SNP", groupid, "X.HIGH.MODERATE.", sep = "."), 
                      paste("res1.bool.SNP", groupid, "HIGH", sep = "."),
                      paste("res1.bool.NON_SNP", groupid, "X.HIGH.MODERATE.", sep = "."), 
                      paste("res1.bool.NON_SNP", groupid, "HIGH", sep = "."),
                      
                      paste("res1.SNP", groupid, "X.HIGH.MODERATE.", sep = "."), 
                      paste("res1.SNP", groupid, "HIGH", sep = "."),
                      paste("res1.NON_SNP", groupid, "X.HIGH.MODERATE.", sep = "."), 
                      paste("res1.NON_SNP", groupid, "HIGH", sep = "."),
                      
                      paste("new.vars.SNP", groupid, "X.HIGH.MODERATE.", sep = "."), 
                      paste("new.vars.SNP", groupid, "HIGH", sep = "."),
                      paste("new.vars.NON_SNP", groupid, "X.HIGH.MODERATE.", sep = "."), 
                      paste("new.vars.NON_SNP", groupid, "HIGH", sep = ".")
                      )

  if(! all(sapply(needed.objects, exists))) {stop(paste("At least one needed object does not exist,", "groupid:", groupid, ", objects:", paste(needed.objects[!sapply(needed.objects, exists)], collapse = ", ")))}
  
  wb <- createWorkbook()
  if(txt.file != "NA") {
    wb2 <- createWorkbook()
  }

  for(impact in impacts) {
    impact.name <- sub(impact, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE")  
    
    df.snp.name <- paste("res1", "SNP", groupid, impact, sep = ".")
    colSums.snp.df <- as.data.frame(colSums(get(df.snp.name)))
    colSums.snp.df[["Genes"]] <- rownames(colSums.snp.df)
    
    df.snp.new.vars.name <- paste("new.vars", "SNP", groupid, impact, sep = ".")
    unique.var.genes.snp <- unique(get(df.snp.new.vars.name)[!duplicated(get(df.snp.new.vars.name)[,"HGVSg"]), , drop = F][,"SYMBOL"])
    unique.var.genes.snp.numbers <- foreach(gene = unique.var.genes.snp, .combine = c) %dopar% {sum(get(df.snp.new.vars.name)[!duplicated(get(df.snp.new.vars.name)[,"HGVSg"]), , drop = F][,"SYMBOL"] == gene)}
    if(length(unique.var.genes.snp.numbers) >0) {names(unique.var.genes.snp.numbers) <- unique.var.genes.snp
    unique.var.genes.snp.numbers.df <- as.data.frame(unique.var.genes.snp.numbers)
    unique.var.genes.snp.numbers.df[["Genes"]] <- rownames(unique.var.genes.snp.numbers.df)} else {
    unique.var.genes.snp.numbers.df <- data.frame()
    }
    
    df.non_snp.name <- paste("res1", "NON_SNP", groupid, impact, sep = ".")
    colSums.non_snp.df <- as.data.frame(colSums(get(df.non_snp.name)))
    colSums.non_snp.df[["Genes"]] <- rownames(colSums.non_snp.df)
    
    df.non_snp.new.vars.name <- paste("new.vars", "NON_SNP", groupid, impact, sep = ".")
    unique.var.genes.non_snp <- unique(get(df.non_snp.new.vars.name)[!duplicated(get(df.non_snp.new.vars.name)[,"HGVSg"]), , drop = F][,"SYMBOL"])
    unique.var.genes.non_snp.numbers <- foreach(gene = unique.var.genes.non_snp, .combine = c) %dopar% {sum(get(df.non_snp.new.vars.name)[!duplicated(get(df.non_snp.new.vars.name)[,"HGVSg"]), , drop = F][,"SYMBOL"] == gene)}
    if(length(unique.var.genes.non_snp.numbers) >0) {names(unique.var.genes.non_snp.numbers) <- unique.var.genes.non_snp
    unique.var.genes.non_snp.numbers.df <- as.data.frame(unique.var.genes.non_snp.numbers)
    unique.var.genes.non_snp.numbers.df[["Genes"]] <- rownames(unique.var.genes.non_snp.numbers.df)} else {
    unique.var.genes.non_snp.numbers.df <- data.frame()
    }

    df.snp.name <- paste("res1.bool", "SNP", groupid, impact, sep = ".")
    df.non.snp.name <- paste("res1.bool", "NON_SNP", groupid, impact, sep = ".")
    res1.bool.snp <- get(df.snp.name)
    res1.bool.non.snp <- get(df.non.snp.name)
    res1.bool.snp.t <- t(res1.bool.snp)
    res1.bool.non.snp.t <- t(res1.bool.non.snp)
    if(!all(colnames(res1.bool.non.snp.t) == colnames(res1.bool.snp.t))) {stop(paste("Sample names in SNP and NON_SNP data frames do not match, groupid:", groupid, ", SNP:", paste(colnames(res1.bool.snp.t), collapse = ", "), ", NON_SNP:", paste(colnames(res1.bool.non.snp.t), collapse = ", ")))}
    sampleIDs <- colnames(res1.bool.non.snp.t)
    res1.bool.merged <- merge(x = res1.bool.snp.t, y = res1.bool.non.snp.t, by.x = 0, by.y = 0, all.x = T, all.y = T)
    
    if(nrow(res1.bool.merged) > 0) {
    res1.bool.merged[is.na(res1.bool.merged)] <- 0
    rownames(res1.bool.merged) <- res1.bool.merged$Row.names
    res1.bool.merged <- res1.bool.merged %>% dplyr::select(-c("Row.names"))
    res1.bool.merged.sum <- foreach( i = sampleIDs, .combine="cbind") %dopar% {res1.bool.merged %>% transmute(z = rowSums(.[grep(pattern = paste0("^", i, "\\.(x|y)$"), colnames(.))]))}
    rownames(res1.bool.merged.sum) <- rownames(res1.bool.merged)
    colnames(res1.bool.merged.sum) <- sampleIDs
    indexes <- which(res1.bool.merged.sum == 2, arr.ind=T)
    snp.non.snp.df <- as.data.frame(cbind(rownames(res1.bool.merged.sum)[indexes[,1]],colnames(res1.bool.merged.sum)[indexes[,2]]))
    colnames(snp.non.snp.df) <- c("Gene", "Sample")
    snp.non.snp.df[["Sample"]] <- sub(snp.non.snp.df[["Sample"]], pattern = "#.*$", replacement = "")
    snp.non.snp.df <- merge(x = snp.non.snp.df, y = anno, by.x = "Sample", by.y = 0)
    snp.non.snp.df <- snp.non.snp.df[order(snp.non.snp.df[[Indfactor]]),]
    
    sheet.name <- paste("SNP&NON_SNP", impact.name, sep = ".")
    addWorksheet(wb = wb, sheetName = sheet.name)
    writeData(x = snp.non.snp.df, wb = wb, sheet = sheet.name, rowNames = F)
    
    if(txt.file != "NA") {
      snp.non.snp.df.signature <- snp.non.snp.df[snp.non.snp.df[["Gene"]] %in% gene.signature, , drop = F]
      addWorksheet(wb = wb2, sheetName = sheet.name)
      writeData(x = snp.non.snp.df.signature, wb = wb2, sheet = sheet.name, rowNames = F)
    }
    
    res1.bool.merged.sum[res1.bool.merged.sum == 2] <- 1
    res1.bool.merged.sum.t <- as.data.frame(t(res1.bool.merged.sum))

    old.names <- colnames(res1.bool.merged.sum.t)[grepl(colnames(res1.bool.merged.sum.t), pattern = "^[0-9]+$")]
    if(length(old.names) > 0) {
    new.names <- mapIds(x = org.Hs.eg.db, keys = old.names, keytype = "ENTREZID", column = "ENSEMBL", multiVals = "first")
    new.names[is.na(new.names)] <- names(new.names[is.na(new.names)])
    
    colnames(res1.bool.merged.sum.t)[grepl(colnames(res1.bool.merged.sum.t), pattern = "^[0-9]+$")] <- new.names
    dup.cols <- colnames(res1.bool.merged.sum.t)[duplicated(colnames(res1.bool.merged.sum.t))]
    pattern <- paste0("(", dup.cols, "(\\.[0-9]+$)?)", collapse = "|")
    colnames(res1.bool.merged.sum.t) <- make.unique(colnames(res1.bool.merged.sum.t))
    corr.df <- foreach(col = dup.cols, .combine = cbind) %do% {
      zcol <- res1.bool.merged.sum.t %>% transmute(z = rowSums(.[grep(pattern = paste0("^", col, "(\\.[0-9]+)?"), colnames(.))]))
      colnames(zcol) <- col
      zcol}
    corr.df[corr.df > 1] <- 1
    res1.bool.merged.sum.t <- res1.bool.merged.sum.t %>% dplyr::select(-colnames(res1.bool.merged.sum.t)[grepl(colnames(res1.bool.merged.sum.t), pattern = pattern)])
    res1.bool.merged.sum.t <- merge(x = res1.bool.merged.sum.t, y = corr.df, by.x = 0, by.y = 0, sort = F)
    rownames(res1.bool.merged.sum.t) <- res1.bool.merged.sum.t[["Row.names"]]
    res1.bool.merged.sum.t <- res1.bool.merged.sum.t[,-1]}
    
    no.var.samples.df <- as.data.frame(colSums(res1.bool.merged.sum.t))
    no.var.samples.df[["Genes"]] <- rownames(no.var.samples.df)
    freq.var.samples.df <- as.data.frame(colSums(res1.bool.merged.sum.t)/nrow(res1.bool.merged.sum.t))
    freq.var.samples.df[["Genes"]] <- rownames(freq.var.samples.df)
    
    if(!all(rownames(res1.bool.merged.sum.t) == rownames(anno))) {stop("Samples do not match annotations.")}
    res1.bool.merged.sum.t.list <- as.list(with(res1.bool.merged.sum.t, by(data = res1.bool.merged.sum.t, INDICES = anno[[Indfactor]], FUN = print, simplify = F)))
    res1.bool.merged.sum.t.list <- res1.bool.merged.sum.t.list[sapply(res1.bool.merged.sum.t.list, function(x){!is.null(x)})]
    res1.bool.merged.sum.t.comb <- foreach(subset1 = res1.bool.merged.sum.t.list, .combine = cbind) %do% {data.frame(colSums(subset1), colMeans(subset1), colSds(as.matrix(subset1)))}
    colnames(res1.bool.merged.sum.t.comb) <- as.vector(outer(c("Altered samples", "Fraction of altered samples", "SD"), names(res1.bool.merged.sum.t.list), paste, sep="."))
    res1.bool.merged.sum.t.comb.res <- reshape(res1.bool.merged.sum.t.comb, direction = "long", ids = rownames(res1.bool.merged.sum.t.comb), varying = colnames(res1.bool.merged.sum.t.comb))
    res1.bool.merged.sum.t.comb.res[is.na(res1.bool.merged.sum.t.comb.res)] <- 0
    res1.bool.merged.sum.t.comb.res <- setNames(res1.bool.merged.sum.t.comb.res, nm = c(Indfactor, "Altered samples", "Fraction of altered samples", "SD", "Gene"))
    res1.bool.merged.sum.t.comb.res <- res1.bool.merged.sum.t.comb.res[, c(5, 1:4), drop = F]
    res1.bool.merged.sum.t.comb.res[[Indfactor]] <- as.factor(res1.bool.merged.sum.t.comb.res[[Indfactor]])
    
    sheet.name <- paste("Var.freq.sum", impact.name, sep = ".")
    addWorksheet(wb = wb, sheetName = sheet.name)
    writeData(x = res1.bool.merged.sum.t.comb.res, wb = wb , sheet = sheet.name, rowNames = F)
    
    if(txt.file != "NA") {
      res1.bool.merged.sum.t.comb.res.signature <- res1.bool.merged.sum.t.comb.res[res1.bool.merged.sum.t.comb.res[["Gene"]] %in% gene.signature, , drop = F]
      addWorksheet(wb = wb2, sheetName = sheet.name)
      writeData(x = res1.bool.merged.sum.t.comb.res.signature, wb = wb2 , sheet = sheet.name, rowNames = F)
    }
    
  if(length(unique(res1.bool.merged.sum.t.comb.res[["Gene"]]))>56) {width <- length(unique(res1.bool.merged.sum.t.comb.res[["Gene"]]))/8} else {width <- 7}
  pdf(file = paste(runid, "VEP.analysis.cumulative.res", groupid, impact, "pdf", sep = "."), width = width)
    ggplot1 <- ggplot(res1.bool.merged.sum.t.comb.res, aes(x = Gene, y = `Fraction of altered samples`, fill = get(Indfactor))) + 
      geom_bar(stat = "identity", position = position_dodge()) + 
      labs(title = paste0("Cumulative frequency of variants ", "(", "impact: ", impact.name, ")"), x = "Gene names", y = "Fraction of altered samples", fill = Indfactor) +
      theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      scale_fill_manual(values = as.vector(cols25()))
    print(ggplot1)
  dev.off()
  
  if(txt.file != "NA") {
    res1.bool.merged.sum.t.comb.res.signature <- res1.bool.merged.sum.t.comb.res[res1.bool.merged.sum.t.comb.res[["Gene"]] %in% gene.signature, , drop = F]
    if(nrow(res1.bool.merged.sum.t.comb.res.signature) >0)
    {
      if(length(unique(res1.bool.merged.sum.t.comb.res.signature[["Gene"]]))>56) {width <- length(unique(res1.bool.merged.sum.t.comb.res.signature[["Gene"]]))/8} else {width <- 7}
      pdf(file = paste(runid, "VEP.analysis.cumulative.res", groupid, "gene.signature", gene.signature.name, impact, "pdf", sep = "."), width = width)
      ggplot1 <- ggplot(res1.bool.merged.sum.t.comb.res.signature, aes(x = Gene, y = `Fraction of altered samples`, fill = get(Indfactor))) + 
        geom_bar(stat = "identity", position = position_dodge()) + 
        labs(title = paste0("Cumulative frequency of variants ", "(", "impact: ", impact.name, ")\ngene signature: ", gene.signature.name), x = "Gene names", y = "Fraction of altered samples", fill = Indfactor) +
        theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        scale_fill_manual(values = as.vector(cols25()))
      print(ggplot1)
      dev.off()
    }}
  
    geneList <- colnames(res1.bool.merged.sum.t)
    
    level.names <- sort(as.character(unique(anno[[Indfactor]])))
    if(!all(rownames(res1.bool.merged.sum.t) == rownames(anno))) {stop("Annotations do not match sample names.")}
    res1.bool.merged.sum.t.colSums <- with(res1.bool.merged.sum.t, by(res1.bool.merged.sum.t, INDICES = anno[[Indfactor]], FUN = colSums, simplify = F))
    res1.bool.merged.sum.t.colSums.list <- foreach(level = level.names) %dopar% {
      df.tmp <- as.data.frame(cbind(res1.bool.merged.sum.t.colSums[[level]], nrow(subset(res1.bool.merged.sum.t, subset = anno[[Indfactor]] == level)) - res1.bool.merged.sum.t.colSums[[level]]))
      df.tmp <- setNames(df.tmp, nm = c("Altered", "Not_altered"))}
    names(res1.bool.merged.sum.t.colSums.list) <- level.names
    
    rownames(res1.bool.merged.sum.t) <- sub(rownames(res1.bool.merged.sum.t), pattern = "#.*$", replacement = "")
    write.table(x = res1.bool.merged.sum.t, row.names = T, col.names = NA, quote = F, sep = ";", file = paste(runid, "VEP.analysis.results", groupid, Indfactor, "impact", sub(impact, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE"), "csv", sep = "."))
    
    if(txt.file != "NA") {
      res1.bool.merged.sum.t.signature <- res1.bool.merged.sum.t[, colnames(res1.bool.merged.sum.t) %in% gene.signature, drop = F]
      write.table(x = res1.bool.merged.sum.t.signature, row.names = T, col.names = NA, quote = F, sep = ";", file = paste(runid, "VEP.analysis.results", "gene.signature", gene.signature.name, groupid, Indfactor, "impact", sub(impact, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE"), "csv", sep = "."))
    }
    
    if(length(level.names) > 1) {
      combinations <- combn(level.names, m = 2)
      comp.list <- NULL
      for(i in seq(1, ncol(combinations))) {
        comb.res.list.name <- paste("comb.res.list", groupid, Indfactor, paste(combinations[,i], collapse = "_vs_"), "impact", impact, sep = ".")
        comb.res.list.tmp <- foreach(gr1 = combinations[, i][1], gr2 = combinations[, i][2]) %do% {
          l1.tmp <- foreach(gene = geneList) %dopar% {
            vec1 <- res1.bool.merged.sum.t.colSums.list[[gr1]][rownames(res1.bool.merged.sum.t.colSums.list[[gr1]]) == gene, , drop = F]
            vec2 <- res1.bool.merged.sum.t.colSums.list[[gr2]][rownames(res1.bool.merged.sum.t.colSums.list[[gr2]]) == gene, , drop = F]
            df.tmp <- rbind(vec1, vec2)
            rownames(df.tmp) <- c(paste(gene, gr1, sep = "."), paste(gene, gr2, sep = "."))
            df.tmp
          }
          geneList.final.tmp <- foreach(gene = geneList, .combine = c) %dopar% {
            vec1 <- res1.bool.merged.sum.t.colSums.list[[gr1]][rownames(res1.bool.merged.sum.t.colSums.list[[gr1]]) == gene, , drop = F]
            vec2 <- res1.bool.merged.sum.t.colSums.list[[gr2]][rownames(res1.bool.merged.sum.t.colSums.list[[gr2]]) == gene, , drop = F]
            df.tmp <- rbind(vec1, vec2)
            var.freq <- NULL
            for(i in seq(1, nrow(df.tmp))) {
              var.freq <- append(var.freq, values = df.tmp[i,"Altered"] / sum(df.tmp[i,]))}
            gene.new.name <- paste(gene, which.max(var.freq), sep = "#")
            gene.new.name
          }
          names(l1.tmp) <- geneList.final.tmp
          l1.tmp
        }
        comp1 <- paste(gr1, gr2, sep = "_vs_")
        comp.list <- append(x = comp.list, values = comp1)
        names(comb.res.list.tmp) <- comp1
        assign(comb.res.list.name, value = comb.res.list.tmp)
        sink(file = paste(runid, sub(sub(comb.res.list.name, pattern = "comb.res.list.", replacement = "Two_by_two_tables."), groupid, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE"), "txt", sep = "."))
        print(get(comb.res.list.name))
        sink()
      }
      
      chi.fisher.res <- foreach(comp = comp.list) %do% {
        comp.name <- paste("comb.res.list", groupid, Indfactor, comp, "impact", impact, sep = ".")
        geneListNew <- NULL
        for(i in seq(1, length(geneList))) {
          geneListNew <- append(geneListNew, values = (grep(names(get(comp.name)[[comp]]), pattern = paste(geneList[i], "(1|2)", sep = "#"), value = T)))
          geneListNew <- unique(geneListNew)}
        gene.list <- foreach(gene = geneListNew) %dopar% {
          if(min(get(comp.name)[[comp]][[gene]])<5) {fisher.test(get(comp.name)[[comp]][[gene]])} else {chisq.test(get(comp.name)[[comp]][[gene]])}
        }
        names(gene.list) <- geneListNew
        gene.list
      }
      names(chi.fisher.res) <- comp.list
      chi.fisher.res.name <- paste("chi.fisher.res", groupid, "impact", impact, sep = ".")
      assign(chi.fisher.res.name, value = chi.fisher.res)
      sink(file = paste(runid, sub(sub(chi.fisher.res.name, pattern = "chi.fisher.res.", replacement = "Chi_square_or_Fisher_exact_tests."), Indfactor, groupid, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE"), "txt", sep = "."))
      print(get(chi.fisher.res.name))
      sink()
      
      allGenes.list <- foreach(group = comp.list) %do% {
        geneListNew <- names(get(paste("chi.fisher.res", groupid, "impact", impact, sep = "."))[[group]])
        vec <- foreach(gene = geneListNew, .combine = "c") %dopar% {
          gene.p.val <- get(paste("chi.fisher.res", groupid, "impact", impact, sep = "."))[[group]][[gene]][["p.value"]]}
        names(vec) <- geneListNew
        vec}
      names(allGenes.list) <- comp.list
      allGenes.list.name <- paste("allGenes.list", groupid, impact, sep = ".")
      assign(allGenes.list.name, value = allGenes.list)
      
      allGenes.df.list <- foreach(group1 = comp.list) %dopar% {
        group1_cats <- unlist(strsplit(group1, split = "_vs_"))
        up_altered.group <- sub(names(allGenes.list[[group1]]), pattern = "(^.*#)", replacement = "")
        up_altered_var <- foreach(up_group = up_altered.group, .combine = c) %do% {
          if(up_group == "1") {group1_cats[1]} else 
          if(up_group == "2") {group1_cats[2]} else 
            {stop("The was a problem with the determination of the up-altered group.")}}
        df.new <- data.frame(up_altered_var, allGenes.list[[group1]],
                   p.adjust(allGenes.list[[group1]], method = "fdr"))
        
        rownames(df.new) <- sub(names(allGenes.list[[group1]]), pattern = "#.*$", replacement = "")
        df.new[["Genes"]] <- rownames(df.new)
        colnames(df.new) <- c(as.character(outer(group1, c("up-altered.group", "p-value", "padj"), paste, sep = ".")), "Genes")
        df.new[,c(4,1:3)]
        }
      names(allGenes.df.list) <- comp.list
      
      allGenes.df <- Reduce(function(x, y) merge(x = x, y = y, by.x = "Genes", by.y = "Genes", all = T), allGenes.df.list)
      allGenes.df[["Ensembl.ID"]] <- mapIds(x = org.Hs.eg.db, keys = allGenes.df[["Genes"]], keytype = "SYMBOL", column = "ENSEMBL", multiVals = "first")
      allGenes.df[["Entrez.ID"]] <- mapIds(x = org.Hs.eg.db, keys = allGenes.df[["Genes"]], keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
      allGenes.df <- allGenes.df[, c(1,ncol(allGenes.df)-1, ncol(allGenes.df),2:(ncol(allGenes.df)-2))]
      
      addWorksheet(wb = wb, sheetName = paste("All_stats", impact.name, sep = "."))
      writeData(x = allGenes.df, wb = wb, sheet = paste("All_stats", impact.name, sep = "."), rowNames = F, keepNA = T, na.string = "NA")
      
      if(txt.file != "NA") {
        allGenes.df.signature <- allGenes.df[allGenes.df[["Genes"]] %in% gene.signature, , drop = F]
        addWorksheet(wb = wb2, sheetName = paste("All_stats", impact.name, sep = "."))
        writeData(x = allGenes.df.signature, wb = wb2, sheet = paste("All_stats", impact.name, sep = "."), rowNames = F, keepNA = T, na.string = "NA")
      }
      
      sigGenes.list <- foreach(group1 = comp.list) %dopar% {
      if(fdr) {
        p.values <- p.adjust(allGenes.list[[group1]], method = "fdr")} else {
        p.values <- allGenes.list[[group1]]}
      p.values[p.values < 0.05]
      }
      names(sigGenes.list) <- comp.list
      sigGenes.list.name <- paste("sigGenes.list", groupid, impact, sep = ".")
      assign(sigGenes.list.name, value = sigGenes.list)
      
      sink(file = paste(runid, "Significant_genes", Indfactor, groupid, "impact", sub(impact, pattern = "X.HIGH.MODERATE.", replacement = "HIGH_or_MODERATE"), paste("fdr", fdr, sep = ":"), "txt", sep = "."))
      cat("Genes with at least one variant in at least one sample (the list of genes submitted to the topGO app):\n")
      cat(geneList, sep = ", ")
      cat("\n\nGenes with significant differences in variant frequencies between the displayed groups (p-values of the Chi-square test or the Fisher exact test are included). Values after the hash sign(#) show in which group the frequency of variants is higher (1 and 2 stand for the group before and after 'vs', respectively):\n")
      print(sigGenes.list)
      sink()
      
      for(group in comp.list) {
        group.dir <- paste(Indfactor, group, paste("fdr", fdr, sep = ":"), sep = ".")
        dir.create(group.dir)
        setwd(group.dir)
        all.genes <- allGenes.list[[group]]
        names(all.genes) <- sub(names(all.genes), pattern = "#(1|2)", replacement = "")
        sig.genes <- sigGenes.list[[group]]
        if(!is.null(sig.genes)) {names(sig.genes) <- sub(names(sig.genes), pattern = "#(1|2)", replacement = "")}
        if(length(all.genes) > 0 & length(sig.genes) > 0) {
          try(f.topGO(geneList = all.genes, sigGenes = sig.genes))
        }
        setwd(workdir)
      }}
    summary.table.colnames <- c("Genes", paste0("SNP variants", " (", impact.name, ")"), paste0("New unique SNP variants", " (", impact.name, ")"), paste0("NON_SNP variants", " (", impact.name, ")"), paste0("New unique NON_SNP variants", " (", impact.name, ")"), paste0("Altered samples", " (", impact.name, ")"), paste0("Fraction of altered samples", " (", impact.name, ")"))
    summary.table.list <- list(colSums.snp.df, unique.var.genes.snp.numbers.df, colSums.non_snp.df, unique.var.genes.non_snp.numbers.df, no.var.samples.df, freq.var.samples.df)
    summary.table.colnames.bool <- c(TRUE, sapply(summary.table.list,nrow) != 0)
    summary.table.list <- summary.table.list[sapply(summary.table.list, nrow) != 0]
    summary.table.colnames <- summary.table.colnames[summary.table.colnames.bool]
    summary.table <- Reduce(function(x,y) merge(x = x, y = y, all = T, by.x = "Genes", by.y = "Genes"), summary.table.list)
    summary.table <- setNames(summary.table, nm = summary.table.colnames)
    summary.table[is.na(summary.table)] <- 0
    
    summary.table.name <- paste("summary.table", impact.name, sep = ".")
    assign(summary.table.name, value = summary.table)
    }}
  
  if(length(ls(pattern = "summary\\.table\\.HIGH")) == 2) {
    summary.table.final <- merge(x = summary.table.HIGH, y = summary.table.HIGH_or_MODERATE, by.x = "Genes", by.y = "Genes", all = T)
    summary.table.final[is.na(summary.table.final)] <- 0
    } else
  if(length(ls(pattern = "summary\\.table\\.HIGH")) == 1) {
    summary.table.final <- get(ls(pattern = "summary\\.table\\.HIGH"))}

  if(exists("summary.table.final")) {
    addWorksheet(wb = wb, sheetName = "Summary_table")
    writeData(wb = wb, x = summary.table.final, sheet = "Summary_table")
    
    if(exists("gene.signature")) {
      summary.table.final.gene.signature <- summary.table.final[summary.table.final[["Genes"]] %in% gene.signature, , drop = F]
      addWorksheet(wb = wb2, sheetName = "Summary_table")
      writeData(wb = wb2, x = summary.table.final.gene.signature, sheet = "Summary_table")
    }}

  pdffiles2 <- grep(list.files(pattern = paste0("^", runid, "\\.VEP.analysis.cumulative.res\\.", groupid, ".*", "\\.pdf$")), pattern = paste("\\.gene\\.signature", gene.signature.name, sep = "."), value = T, invert = T)
  pdf_combine(pdffiles2, output = paste(runid, "VEP.analysis.cumulative.results.ALL.GENES", groupid, Indfactor, "pdf", sep = "."))
  unlink(pdffiles2)
  
  if(txt.file != "NA") {
    pdffiles3 <- grep(list.files(pattern = paste0("^", runid, "\\.VEP.analysis.cumulative.res\\.", groupid, ".*", "\\.pdf$")), pattern = paste("\\.gene\\.signature", gene.signature.name, sep = "."), value = T)
    pdf_combine(pdffiles3, output = paste(runid, "VEP.analysis.cumulative.results.gene.signature", gene.signature.name, groupid, Indfactor, "pdf", sep = "."))
    unlink(pdffiles3)
  }
  
  saveWorkbook(wb = wb, file = gsub(paste(runid, "VEP.analysis.cumulative.results", groupid, Indfactor, "xlsx", sep = "."), pattern = "\\.\\.", replacement = "."), overwrite = T)
  
  if(txt.file != "NA") {
  saveWorkbook(wb = wb2, file = gsub(paste(runid, "VEP.analysis.cumulative.results", groupid, Indfactor, paste0("gene_signature:", gene.signature.name), "xlsx", sep = "."), pattern = "\\.\\.", replacement = "."), overwrite = T)
  }
  
  if(length(level.names) > 1) {
    for(group in comp.list) {
      group.dir <- paste(Indfactor, group, paste("fdr", fdr, sep = ":"), sep = ".")
      dir.create(group.dir)
      setwd(group.dir)
      f.Reactome(group = group)
      setwd(workdir)
    }
  }
} else {stop(paste("At least one of the needed R objects was not found:", paste(c(Robject1, Robject2)[!sapply(c(Robject1, Robject2), file.exists)], collapse = ", ")))}
}

workdir <- arguments2[1]
runid <- sub(arguments2[1], pattern = "(.*\\/RUNS\\/)(.*)(\\/(MAPPINGS|MAPPINGS_TRIMMED)\\/.*)", replacement = "\\2")
threads <- as.numeric(arguments2[2])
csvfile <- arguments2[3]
grouping.col <- arguments2[4]
grouping.cols <- unlist(strsplit(arguments2[4], split = "\\,"))
Indfactor <- arguments2[5]
pairedSamples <- arguments2[6]
pair.ident <- arguments2[7]
fdr <- as.logical(arguments2[8])

registerDoMC(threads)

csvtable.full <- fread(file = csvfile, sep = ";", header = T, stringsAsFactors = F)
csvtable.full <- csvtable.full[!grepl(csvtable.full[[Indfactor]], pattern = "^([[:space:]]?)+$"),]
csvtable.full <- csvtable.full[!is.na(csvtable.full[[Indfactor]]),]

if(grouping.col != "ALL_SAMPLES") {
csvtable.full <- csvtable.full[csvtable.full %>% dplyr::select(all_of(grouping.cols)) %>% is.na() %>% rowSums() == 0,]
}

if(pairedSamples == "TRUE") {csvtable.full <- csvtable.full[!is.na(csvtable.full[[pair.ident]]),]}

Robject1.dir <- sub(arguments2[1], pattern = "\\/COMPARISON\\/", replacement = "/SNP/Summary/")
Robject2.dir <- sub(arguments2[1], pattern = "\\/COMPARISON\\/", replacement = "/NON_SNP/Summary/")

dir.create(workdir, recursive = T)
setwd(workdir)

if(grouping.col != "ALL_SAMPLES") {
  if(length(grouping.cols) == 1) {
    groupids <- unique(csvtable.full[[grouping.col]])} else if(length(grouping.cols) == 2) {
      col.pairs <- unique(csvtable.full %>% dplyr::select(grouping.cols))
      groupids <- foreach(i = seq(1,nrow(col.pairs)), .combine = c) %dopar% {paste(as.character(col.pairs[i,]), collapse = "+")}
    } else {stop("This program supports up to two grouping variables only.")}
} else {
  groupids <- "ALL_SAMPLES"
}

save.image(paste(runid, "VEP.analysis.comparison.RData", sep = "."))

for(groupid in groupids) {
  Robject1.file <- paste(runid, "SNP", groupid, Indfactor, paste("paired", pairedSamples, sep = ":"), "VEP_analysis.final.RData", sep = ".")
  Robject2.file <- paste(runid, "NON_SNP", groupid, Indfactor, paste("paired", pairedSamples, sep = ":"), "VEP_analysis.final.RData", sep = ".")
  Robject1.path <- paste(Robject1.dir, Robject1.file, sep = "/")
  Robject2.path <- paste(Robject2.dir, Robject2.file, sep = "/")
  f.vep.comparison(Robject1 = Robject1.path, Robject2 = Robject2.path)
  rm(list = Robject1.list)
  rm(list = Robject2.list)
}

sessionInfo()
proc.time()
date()

cat("All done.\n")
