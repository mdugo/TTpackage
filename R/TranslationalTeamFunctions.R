
#*************************************************************
#* This function performs singscore calculation for gene set
#* collections that include gene sets with only up-regulated
#* genes or both up- and down-regulated genes
#*************************************************************

simpleScoreMod<-function(rankData, mysetlist, knownDir = TRUE, centerScore = TRUE){

  require(Biobase)
  require(singscore)

  score<-matrix(0, nrow=length(mysetlist),ncol=ncol(rankData))
  rownames(score)<-names(mysetlist)
  colnames(score)<-colnames(rankData)
  for(i in 1:length(mysetlist)){
    if(class(mysetlist[[i]])=="character"){
      sl<-mysetlist[[i]]
      sl<-sl[sl%in%rownames(rankData)]
      if(length(sl)>=3){
        scoretemp<-simpleScore(rankData,upSet = sl, knownDirection = knownDir, centerScore = centerScore)
        score[rownames(score)==names(mysetlist)[i],]<-scoretemp$TotalScore
      }
    } else {
      sl.up<-mysetlist[[i]][[grep("UP",names(mysetlist[[i]]),ignore.case=T)]]
      sl.up<-sl.up[sl.up%in%rownames(rankData)]
      sl.dn<-mysetlist[[i]][[grep("DN|DOWN",names(mysetlist[[i]]),ignore.case=T)]]
      sl.dn<-sl.dn[sl.dn%in%rownames(rankData)]
      if(length(sl.up)>=3 & length(sl.dn)>=3){
        scoretemp<-simpleScore(rankData,upSet = sl.up, downSet = sl.dn)
        score[rownames(score)==names(mysetlist)[i],]<-scoretemp$TotalScore
      }
    }
  }
  score<-score[rowSums(score!=0)>0,]
  return(score)
}



#**************************************************************************************
#* Given a list of gene sets and an expression matrix this function performs singscore,
#* evaluates the correlation between gene sets, generate a network and detect gene sets
#* clusters of highly correlated gene sets
#**************************************************************************************



correlateGeneSets <- function(cluster_by = c("gene_content","correlation"),
                               similarity_method = c("kappa", "jaccard", "dice", "overlap"),
                               similarity_th = 0,
                               expMat, enrichment_table, directional=TRUE,
                               GSEA_style=TRUE, gene_set_list=NULL,
                               correl_th=0.9, FDR_th=0.05, FDR_column,
                               gene_set_name_column, NES_column="NES", NES_th=0,
                               direction_column = NULL,
                               gene_content_column=NULL,
                               igraph.vertex.label = NA, igraph.vertex.label.cex = 0.2,
                               igraph.vertex.frame.color="black", igraph.edge.color = "black",
                               igraph.mark.border = "black", legend_cex=0.4){


  require(Biobase)
  require(circlize)
  require(ComplexHeatmap)
  require(dplyr)
  require(ggplot2)
  require(igraph)
  require(lattice)
  require(msigdbr)
  require(randomcoloR)
  require(rcartocolor)
  require(RColorBrewer)
  require(simplifyEnrichment)
  require(singscore)
  require(xlsx)

  enrichment_table<-as.data.frame(enrichment_table)
  sig<-enrichment_table[enrichment_table[,FDR_column] < FDR_th & abs(enrichment_table[,NES_column]) >= NES_th, gene_set_name_column]

  if(length(sig) == 0){
    stop("Error: no significant gene sets found at provided thresholds.")
  }

  if(cluster_by == "gene_content"){
    if(is.null(gene_content_column)){
      print("No leading edge provided. Proceeding with the full gene content of gene sets.")

      if(length(grep("HALLMARK|REACTOME|KEGG|BIOCARTA|PID|WP|GOBP", unique(gsub("_.*", "", sig)))) >= 1){
        print("Retrieving gene sets from MSigDB...")
        gs1 <- msigdbr(species = "Homo sapiens", category = "H")
        gs2 <- msigdbr(species = "Homo sapiens", category = "C2")
        gs3 <- msigdbr(species = "Homo sapiens", category = "C5")
        print("Done.")
        gs<-rbind(gs1, gs2, gs3)
        gs<-gs[gs$gs_name%in%sig,]

        print("Converting to a list...")
        setlist<-vector("list", length=length(unique(gs$gs_name)))
        names(setlist)<-unique(gs$gs_name)
        for(i in 1:length(setlist)){
          setlist[[i]]<-unique(gs$human_gene_symbol[gs$gs_name==names(setlist)[i]])
        }
        if (nrow(gs) == 0){
          stop("Error: provided gene sets are not found in gene set list.")
        }
      }
      if(!all(unique(gsub("_.*", "", sig)) %in% c("HALLMARK", "REACTOME", "KEGG", "BIOCARTA", "PID", "WP", "GOBP"))){
        if(is.null(gene_set_list)){
          stop("Error: additional gene sets detected. Please provide the custom gene set list with gene symbols.")
        }
        gs_custom<-gene_set_list
        toUnlist<-gs_custom[names(gs_custom) %in% names(which(sapply(gs_custom, class) == "list"))]
        if(length(toUnlist)>0){
          gs_custom<-gs_custom[!names(gs_custom) %in% names(which(sapply(gs_custom, class) == "list"))]
          for(i in 1:length(toUnlist)){
            gs_custom<-c(gs_custom, toUnlist[[i]])
          }
        }
        gs_custom<-gs_custom[names(gs_custom) %in% sig]
        gs_custom<-gs_custom[-grep("HALLMARK|REACTOME|KEGG|BIOCARTA|PID|WP|GOBP", names(gs_custom))]
        names(gs_custom)<-paste0("CUSTOM_", names(gs_custom))
        setlist<-c(setlist, gs_custom)
      }
      if(length(grep("HALLMARK|REACTOME|KEGG|BIOCARTA|PID|WP|GOBP", unique(gsub("_.*", "", sig)))) == 0){
        if(is.null(gene_set_list)){
          stop("Error: No gene set list provided. Please provide the custom gene set list with gene symbols.")
        }
        setlist<-gene_set_list
        names(setlist)<-paste0("CUSTOM_", names(setlist))
      }
      print("Done.")

      ts<-term_similarity(setlist, method = similarity_method)

    } else {
      le<-enrichment_table[enrichment_table[,gene_set_name_column] %in% sig,gene_content_column]
      names(le)<-enrichment_table[enrichment_table[,gene_set_name_column] %in% sig,gene_set_name_column]
      ts<-term_similarity(le, method = similarity_method)
    }

    ig<-graph_from_adjacency_matrix(adjmatrix = ts,
                                    mode="upper",
                                    diag=F,
                                    weighted = T)
    ig<-delete_edges(ig, edges=which(E(ig)$weight <= similarity_th))
    if(gsize(ig) == 0){
      stop("Error: no edges left after filtering for similarity threshold.")
    }

    cl = membership(cluster_fast_greedy(ig))
    cl[cl%in%names(which(table(cl)==1))]<-"Unclustered"
    names(cl)<-gsub("CUSTOM_","",names(cl))

    res2<-enrichment_table[enrichment_table[,gene_set_name_column]%in%names(cl),]
    res2<-res2[match(names(cl), res2[,gene_set_name_column]),]
    identical(names(cl), res2[,gene_set_name_column])
    V(ig)$Cluster<-cl
    V(ig)$FDR<--log10(ifelse(res2[,FDR_column]<1e-10, 1e-10, res2[,FDR_column]))

    if(length(table(cl))>1){
      if(directional == FALSE){
        set.seed(333)
        clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = c(distinctColorPalette(k=length(unique(V(ig)$Cluster))-1),"grey80"))
        names(clustercol)<-factor(V(ig)$Cluster)
        clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
        clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
        set.seed(123)
        coords<-layout_nicely(ig)
        plot.igraph(ig,
                    layout=coords,
                    vertex.size=V(ig)$FDR,
                    vertex.label=igraph.vertex.label,
                    vertex.label.cex=igraph.vertex.label.cex,
                    vertex.frame.color=igraph.vertex.frame.color,
                    vertex.color=clustercol,
                    edge.width = E(ig)$weight,
                    edge.color = igraph.edge.color,
                    margin=c(0,0,0,0))
        legend("topleft",
               legend = names(clustercolLegend)[names(clustercolLegend) != "Unclustered"],
               fill = clustercolLegend[names(clustercolLegend) != "Unclustered"],
               border = T,
               title ="Clusters",
               title.adj=0,
               bty = "n",
               cex = legend_cex)
      }
    }
    if(directional == TRUE & GSEA_style == TRUE){
      V(ig)$Direction<-res2[,NES_column]
      set.seed(333)
      clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = c(distinctColorPalette(k=length(unique(V(ig)$Cluster))-1),"grey80"))
      names(clustercol)<-factor(V(ig)$Cluster)
      clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
      clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
      NEScol<-colorRamp2(breaks = c(min(as.numeric(V(ig)$Direction), na.rm = T), 0, as.numeric(max(V(ig)$Direction), na.rm = T)), colors = c("royalblue", "white", "red2"))
      NEScol<-NEScol(as.numeric(V(ig)$Direction))
      markg<-list()
      for(i in 1:length(unique(cl))){
        markg[[i]]<-which(cl%in%unique(cl)[i])
      }
      names(markg)<-unique(cl)
      markg<-markg[names(markg)!="Unclustered"]

      set.seed(123)
      coords<-layout_nicely(ig)
      plot.igraph(ig,
                  layout=coords,
                  vertex.size=V(ig)$FDR,
                  vertex.label=igraph.vertex.label,
                  vertex.label.cex=igraph.vertex.label.cex,
                  vertex.frame.color=igraph.vertex.frame.color,
                  vertex.color=NEScol,
                  edge.width = E(ig)$weight,
                  edge.color = igraph.edge.color,
                  mark.groups = markg,
                  mark.shape=0.9,
                  mark.expand=10,
                  mark.border = igraph.mark.border,
                  mark.col = clustercolLegend,
                  margin=c(0,0,0,0))
      legend("topleft",
             legend = names(clustercolLegend)[names(clustercolLegend) != "Unclustered"],
             fill = clustercolLegend[names(clustercolLegend) != "Unclustered"],
             border = T,
             title ="Clusters",
             title.adj=0,
             bty = "n",
             cex = legend_cex)
    }
    if(directional == TRUE & GSEA_style == FALSE){
      if(sum(colnames(res) == "Direction") == 0){
        stop("Error: No column providing the directionality of the enrichment was found.")
      } else {
        V(ig)$Direction<-res2$Direction
        set.seed(333)
        clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = c(distinctColorPalette(k=length(unique(V(ig)$Cluster))-1),"grey80"))
        names(clustercol)<-factor(V(ig)$Cluster)
        clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
        clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
        nodecol<-V(ig)$Direction
        nodecol[grep("up", V(ig)$Direction, ignore.case = T)]<-"red2"
        nodecol[grep("down|dn", V(ig)$Direction, ignore.case = T)]<-"royalblue"

        markg<-list()
        for(i in 1:length(unique(cl))){
          markg[[i]]<-which(cl%in%unique(cl)[i])
        }
        names(markg)<-unique(cl)
        markg<-markg[names(markg)!="Unclustered"]

        set.seed(123)
        coords<-layout_nicely(ig)
        plot.igraph(ig,
                    layout=coords,
                    vertex.size=V(ig)$FDR,
                    vertex.label=igraph.vertex.label,
                    vertex.label.cex=igraph.vertex.label.cex,
                    vertex.frame.color=igraph.vertex.frame.color,
                    vertex.color=nodecol,
                    edge.width = E(ig)$weight,
                    edge.color = igraph.edge.color,
                    mark.groups = markg,
                    mark.shape=0.9,
                    mark.expand=10,
                    mark.border = igraph.mark.border,
                    mark.col = clustercolLegend,
                    margin=c(0,0,0,0))
        legend("topleft",
               legend = names(clustercolLegend)[names(clustercolLegend) != "Unclustered"],
               fill = clustercolLegend[names(clustercolLegend) != "Unclustered"],
               border = T,
               title ="Clusters",
               title.adj=0,
               bty = "n",
               cex = legend_cex)
      }
    }
    community<-data.frame(geneSet = names(cl), Cluster = as.vector(cl))
    colnames(community)[1]<-gene_set_name_column
    enrichment_table<-merge(enrichment_table, community, by=gene_set_name_column, all.x=T)
    enrichment_table<-enrichment_table[order(enrichment_table$Cluster, enrichment_table[,FDR_column]),]
    return(enrichment_table)

    if(directional == FALSE){
      set.seed(123)
      coords<-layout_nicely(ig)
      plot.igraph(ig,
                  layout=coords,
                  vertex.size=V(ig)$FDR,
                  vertex.label=igraph.vertex.label,
                  vertex.label.cex=igraph.vertex.label.cex,
                  vertex.frame.color=igraph.vertex.frame.color,
                  vertex.color="grey80",
                  edge.width = E(ig)$weight,
                  edge.color = igraph.edge.color,
                  margin=c(0,0,0,0))
    }
    if(directional == TRUE & GSEA_style == TRUE){
      V(ig)$Direction<-res2[,NES_column]
      NEScol<-colorRamp2(breaks = c(min(as.numeric(V(ig)$Direction), na.rm = T), 0, as.numeric(max(V(ig)$Direction), na.rm = T)), colors = c("royalblue", "white", "red2"))
      NEScol<-NEScol(as.numeric(V(ig)$Direction))

      set.seed(123)
      coords<-layout_nicely(ig)
      plot.igraph(ig,
                  layout=coords,
                  vertex.size=V(ig)$FDR,
                  vertex.label=igraph.vertex.label,
                  vertex.label.cex=igraph.vertex.label.cex,
                  vertex.frame.color=igraph.vertex.frame.color,
                  vertex.color=NEScol,
                  edge.width = E(ig)$weight,
                  edge.color = igraph.edge.color,
                  margin=c(0,0,0,0))
    }
    if(directional == TRUE & GSEA_style == FALSE){
      if(sum(colnames(res) == "Direction") == 0){
        stop("Error: No column providing the directionality of the enrichment was found.")
      } else {
        V(ig)$Direction<-res2$Direction
        nodecol<-V(ig)$Direction
        nodecol[grep("up", V(ig)$Direction, ignore.case = T)]<-"red2"
        nodecol[grep("down|dn", V(ig)$Direction, ignore.case = T)]<-"royalblue"

        set.seed(123)
        coords<-layout_nicely(ig)
        plot.igraph(ig,
                    layout=coords,
                    vertex.size=V(ig)$FDR,
                    vertex.label=igraph.vertex.label,
                    vertex.label.cex=igraph.vertex.label.cex,
                    vertex.frame.color=igraph.vertex.frame.color,
                    vertex.color=nodecol,
                    edge.width = E(ig)$weight,
                    edge.color = igraph.edge.color,
                    margin=c(0,0,0,0))
      }
    }
  }

  if(cluster_by == "correlation"){

    if(length(grep("HALLMARK|REACTOME|KEGG|BIOCARTA|PID|WP|GOBP", unique(gsub("_.*", "", sig)))) >= 1){
      print("Retrieving gene sets from MSigDB...")
      gs1 <- msigdbr(species = "Homo sapiens", category = "H")
      gs2 <- msigdbr(species = "Homo sapiens", category = "C2")
      gs3 <- msigdbr(species = "Homo sapiens", category = "C5")
      print("Done.")
      gs<-rbind(gs1, gs2, gs3)
      gs<-gs[gs$gs_name%in%sig,]

      print("Converting to a list...")
      setlist<-vector("list", length=length(unique(gs$gs_name)))
      names(setlist)<-unique(gs$gs_name)
      for(i in 1:length(setlist)){
        setlist[[i]]<-unique(gs$human_gene_symbol[gs$gs_name==names(setlist)[i]])
      }
      if (nrow(gs) == 0){
        stop("Error: provided gene sets are not found in gene set list.")
      }
    }
    if(!all(unique(gsub("_.*", "", sig)) %in% c("HALLMARK", "REACTOME", "KEGG", "BIOCARTA", "PID", "WP", "GOBP"))){
      if(is.null(gene_set_list)){
        stop("Error: additional gene sets detected. Please provide the custom gene set list with gene symbols.")
      }
      gs_custom<-gene_set_list
      toUnlist<-gs_custom[names(gs_custom) %in% names(which(sapply(gs_custom, class) == "list"))]
      if(length(toUnlist)>0){
        gs_custom<-gs_custom[!names(gs_custom) %in% names(which(sapply(gs_custom, class) == "list"))]
        for(i in 1:length(toUnlist)){
          gs_custom<-c(gs_custom, toUnlist[[i]])
        }
      }
      gs_custom<-gs_custom[names(gs_custom) %in% sig]
      gs_custom<-gs_custom[-grep("HALLMARK|REACTOME|KEGG|BIOCARTA|PID|WP|GOBP", names(gs_custom))]
      names(gs_custom)<-paste0("CUSTOM_", names(gs_custom))
      setlist<-c(setlist, gs_custom)
    }
    if(length(grep("HALLMARK|REACTOME|KEGG|BIOCARTA|PID|WP|GOBP", unique(gsub("_.*", "", sig)))) == 0){
      if(is.null(gene_set_list)){
        stop("Error: No gene set list provided. Please provide the custom gene set list with gene symbols.")
      }
      setlist<-gene_set_list
      names(setlist)<-paste0("CUSTOM_", names(setlist))
    }
    print("Done.")

    if(all(colSums(apply(expMat,2,function(x)sapply(x,is.integer))) == nrow(expMat))){
      stop("Error: all data in the expression matrix are integers. Please provide counts normalized for gene length (FPKM/RPKM or TPM)")
    }

    detectSymbols<-grep("^ACTB$|^GAPDH$|^EGFR$", rownames(expMat), value=T)
    if(length(detectSymbols) == 0){
      stop("Error: row names of expression matrix must be gene symbols.")
    }

    print("Performing singscore...")

    rankMat<-rankGenes(expMat)

    if(length(grep("HALLMARK|CUSTOM",names(setlist)))>0){
      scoreMat1<-simpleScoreMod(rankData = rankMat,
                                mysetlist = setlist[grep("HALLMARK|CUSTOM",names(setlist))],
                                knownDir = T)
      if(length(grep("REACTOME|KEGG|BIOCARTA|PID|WP|GOBP",names(setlist)))>0){
        scoreMat2<-simpleScoreMod(rankData = rankMat,
                                  mysetlist = setlist[grep("REACTOME|KEGG|BIOCARTA|PID|WP|GOBP",names(setlist))],
                                  knownDir = F,
                                  centerScore = F)
        scoreMat<-rbind(scoreMat1, scoreMat2)
      } else {
        scoreMat<-scoreMat1
      }
    } else {
      scoreMat<-simpleScoreMod(rankData = rankMat,
                               mysetlist = setlist,
                               knownDir = F,
                               centerScore = F)
    }
    print("Done.")
    rownames(scoreMat)<-gsub("CUSTOM_", "", rownames(scoreMat))

    correl<-cor(t(scoreMat), method="spearman")

    if(directional == FALSE){

      ig<-graph_from_adjacency_matrix(adjmatrix = correl,
                                      mode="upper",
                                      diag=F,
                                      weighted = T)
      ig<-delete_edges(ig, edges=which(E(ig)$weight <= correl_th))
      if(gsize(ig) == 0){
        stop("Error: no edges left after filtering for similarity threshold.")
      }

      cl = membership(cluster_fast_greedy(ig))
      cl[cl%in%names(which(table(cl)==1))]<-"Unclustered"

      res2<-enrichment_table[enrichment_table[,gene_set_name_column]%in%names(cl),]
      res2<-res2[match(names(cl), res2[,gene_set_name_column]),]
      identical(names(cl), res2[,gene_set_name_column])
      V(ig)$Cluster<-cl
      V(ig)$FDR<--log10(ifelse(res2[,FDR_column]<1e-10, 1e-10, res2[,FDR_column]))
      if(length(table(cl))>1){
        set.seed(444)
        clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = c(distinctColorPalette(k=length(unique(V(ig)$Cluster))-1),"grey80"))
        names(clustercol)<-factor(V(ig)$Cluster)
        clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
        clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
        set.seed(123)
        coords<-layout_nicely(ig)
        plot.igraph(ig,
                    layout=coords,
                    vertex.size=V(ig)$FDR,
                    vertex.label=igraph.vertex.label,
                    vertex.label.cex=igraph.vertex.label.cex,
                    vertex.frame.color=igraph.vertex.frame.color,
                    vertex.color=clustercol,
                    edge.width = E(ig)$weight,
                    edge.color = igraph.edge.color,
                    margin=c(0,0,0,0))
        legend("topleft",
               legend = names(clustercolLegend)[names(clustercolLegend) != "Unclustered"],
               fill = clustercolLegend[names(clustercolLegend) != "Unclustered"],
               border = T,
               title ="Clusters",
               title.adj=0,
               bty = "n",
               cex = legend_cex)
      } else {
        set.seed(123)
        coords<-layout_nicely(ig)
        plot.igraph(ig,
                    layout=coords,
                    vertex.size=V(ig)$FDR,
                    vertex.label=igraph.vertex.label,
                    vertex.label.cex=igraph.vertex.label.cex,
                    vertex.frame.color=igraph.vertex.frame.color,
                    vertex.color="grey80",
                    edge.width = E(ig)$weight,
                    edge.color = igraph.edge.color,
                    margin=c(0,0,0,0))
      }
    }
    if(directional == TRUE & GSEA_style == TRUE){
      up<-enrichment_table[enrichment_table[,gene_set_name_column] %in% sig & enrichment_table[,NES_column] > 0,gene_set_name_column]
      dn<-enrichment_table[enrichment_table[,gene_set_name_column] %in% sig & enrichment_table[,NES_column] < 0,gene_set_name_column]

      tempList<-list(up,dn)
      layout(matrix(c(1:2), ncol = sum(sapply(tempList, length) > 0)))

      if(length(up) > 0){

        correl.up<-correl[rownames(correl) %in% up , colnames(correl) %in% up]

        ig<-graph_from_adjacency_matrix(adjmatrix = correl.up,
                                        mode="upper",
                                        diag=F,
                                        weighted = T)
        ig<-delete_edges(ig, edges=which(E(ig)$weight <= correl_th))
        if(gsize(ig) == 0){
          stop("Error: no edges left after filtering for similarity threshold.")
        }

        cl.up <- membership(cluster_fast_greedy(ig))
        cl.up[cl.up%in%names(which(table(cl.up)==1))]<-"Unclustered"
        if(length(table(cl.up))>1){
          res2<-enrichment_table[enrichment_table[,gene_set_name_column]%in%names(cl.up),]
          res2<-res2[match(names(cl.up), res2[,gene_set_name_column]),]
          identical(names(cl.up), res2[,gene_set_name_column])
          V(ig)$Cluster<-cl.up
          V(ig)$FDR<--log10(ifelse(res2[,FDR_column]<1e-10, 1e-10, res2[,FDR_column]))
          V(ig)$Direction<-res2[,NES_column]
          set.seed(333)
          clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = c(distinctColorPalette(k=length(unique(V(ig)$Cluster))-1),"grey80"))
          names(clustercol)<-factor(V(ig)$Cluster)
          clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
          clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
          NEScol<-colorRamp2(breaks = c(min(as.numeric(V(ig)$Direction), na.rm = T), mean(as.numeric(V(ig)$Direction), na.rm = T), max(as.numeric(V(ig)$Direction), na.rm = T)), colors = c("#FB6A4A", "#CB181D", "#67000D"))
          NEScol<-NEScol(as.numeric(V(ig)$Direction))
          markg<-list()
          for(i in 1:length(unique(cl.up))){
            markg[[i]]<-which(cl.up%in%unique(cl.up)[i])
          }
          names(markg)<-unique(cl.up)
          markg<-markg[names(markg)!="Unclustered"]

          set.seed(123)
          coords<-layout_nicely(ig)
          plot.igraph(ig,
                      layout=coords,
                      vertex.size=V(ig)$FDR,
                      vertex.label=igraph.vertex.label,
                      vertex.label.cex=igraph.vertex.label.cex,
                      vertex.frame.color=igraph.vertex.frame.color,
                      vertex.color=NEScol,
                      edge.width = E(ig)$weight,
                      edge.color = igraph.edge.color,
                      mark.groups = markg,
                      mark.shape=0.9,
                      mark.expand=10,
                      mark.border = igraph.mark.border,
                      mark.col = clustercolLegend,
                      margin=c(0,0,0,0))
          title("Up-regulated",cex.main=2,col.main="black")
          legend("topleft",
                 legend = names(clustercolLegend)[names(clustercolLegend) != "Unclustered"],
                 fill = clustercolLegend[names(clustercolLegend) != "Unclustered"],
                 border = T,
                 title ="Clusters",
                 title.adj=0,
                 bty = "n",
                 cex = legend_cex)
        } else {
          res2<-enrichment_table[enrichment_table[,gene_set_name_column]%in%names(cl.up),]
          res2<-res2[match(names(cl.up), res2[,gene_set_name_column]),]
          V(ig)$Cluster<-cl.up
          V(ig)$FDR<--log10(ifelse(res2[,FDR_column]<1e-10, 1e-10, res2[,FDR_column]))
          V(ig)$Direction<-res2[,NES_column]
          set.seed(333)
          NEScol<-colorRamp2(breaks = c(min(as.numeric(V(ig)$Direction), na.rm = T), mean(as.numeric(V(ig)$Direction), na.rm = T), max(as.numeric(V(ig)$Direction), na.rm = T)), colors = c("#FB6A4A", "#CB181D", "#67000D"))
          NEScol<-NEScol(as.numeric(V(ig)$Direction))

          set.seed(123)
          coords<-layout_nicely(ig)
          plot.igraph(ig,
                      layout=coords,
                      vertex.size=V(ig)$FDR,
                      vertex.label=igraph.vertex.label,
                      vertex.label.cex=igraph.vertex.label.cex,
                      vertex.frame.color=igraph.vertex.frame.color,
                      vertex.color=NEScol,
                      edge.width = E(ig)$weight,
                      edge.color = igraph.edge.color,
                      margin=c(0,0,0,0))
          title("Up-regulated",cex.main=2,col.main="black")
        }
      } else {
        cl.up<-NA
      }

      if(length(dn) > 0){
        correl.dn<-correl[rownames(correl) %in% dn , colnames(correl) %in% dn]

        ig<-graph_from_adjacency_matrix(adjmatrix = correl.dn,
                                        mode="upper",
                                        diag=F,
                                        weighted = T)
        ig<-delete_edges(ig, edges=which(E(ig)$weight <= correl_th))
        if(gsize(ig) == 0){
          stop("Error: no edges left after filtering for similarity threshold.")
        }

        cl.dn <- membership(cluster_fast_greedy(ig))
        cl.dn[cl.dn%in%names(which(table(cl.dn)==1))]<-"Unclustered"
        if(length(table(cl.dn))>1){
          res2<-enrichment_table[enrichment_table[,gene_set_name_column]%in%names(cl.dn),]
          res2<-res2[match(names(cl.dn), res2[,gene_set_name_column]),]
          V(ig)$Cluster<-cl.dn
          V(ig)$FDR<--log10(ifelse(res2[,FDR_column]<1e-10, 1e-10, res2[,FDR_column]))
          V(ig)$Direction<-res2[,NES_column]
          set.seed(343)
          clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = c(distinctColorPalette(k=length(unique(V(ig)$Cluster))-1),"grey80"))
          names(clustercol)<-factor(V(ig)$Cluster)
          clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
          clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
          NEScol<-colorRamp2(breaks = c(min(as.numeric(V(ig)$Direction), na.rm = T), mean(as.numeric(V(ig)$Direction), na.rm = T), max(as.numeric(V(ig)$Direction), na.rm = T)), colors = c("#6BAED6", "#2171B5", "#08306B"))
          NEScol<-NEScol(as.numeric(V(ig)$Direction))
          markg<-list()
          for(i in 1:length(unique(cl.dn))){
            markg[[i]]<-which(cl.dn%in%unique(cl.dn)[i])
          }
          names(markg)<-unique(cl.dn)
          markg<-markg[names(markg)!="Unclustered"]

          set.seed(123)
          coords<-layout_nicely(ig)
          plot.igraph(ig,
                      layout=coords,
                      vertex.size=V(ig)$FDR,
                      vertex.label=igraph.vertex.label,
                      vertex.label.cex=igraph.vertex.label.cex,
                      vertex.frame.color=igraph.vertex.frame.color,
                      vertex.color=NEScol,
                      edge.width = E(ig)$weight,
                      edge.color = igraph.edge.color,
                      mark.groups = markg,
                      mark.shape=0.9,
                      mark.expand=10,
                      mark.border = igraph.mark.border,
                      mark.col = clustercolLegend,
                      margin=c(0,0,0,0))
          title("Down-regulated",cex.main=2,col.main="black")
          legend("topleft",
                 legend = names(clustercolLegend)[names(clustercolLegend) != "Unclustered"],
                 fill = clustercolLegend[names(clustercolLegend) != "Unclustered"],
                 border = T,
                 title ="Clusters",
                 title.adj=0,
                 bty = "n",
                 cex = legend_cex)
        } else {
          res2<-enrichment_table[enrichment_table[,gene_set_name_column]%in%names(cl.dn),]
          res2<-res2[match(names(cl.dn), res2[,gene_set_name_column]),]
          V(ig)$Cluster<-cl.dn
          V(ig)$FDR<--log10(ifelse(res2[,FDR_column]<1e-10, 1e-10, res2[,FDR_column]))
          V(ig)$Direction<-res2[,NES_column]
          set.seed(333)
          NEScol<-colorRamp2(breaks = c(min(as.numeric(V(ig)$Direction), na.rm = T), mean(as.numeric(V(ig)$Direction), na.rm = T), max(as.numeric(V(ig)$Direction), na.rm = T)), colors = c("#6BAED6", "#2171B5", "#08306B"))
          NEScol<-NEScol(as.numeric(V(ig)$Direction))

          set.seed(123)
          coords<-layout_nicely(ig)
          plot.igraph(ig,
                      layout=coords,
                      vertex.size=V(ig)$FDR,
                      vertex.label=igraph.vertex.label,
                      vertex.label.cex=igraph.vertex.label.cex,
                      vertex.frame.color=igraph.vertex.frame.color,
                      vertex.color=NEScol,
                      edge.width = E(ig)$weight,
                      edge.color = igraph.edge.color,
                      margin=c(0,0,0,0))
          title("Down-regulated",cex.main=2,col.main="black")
        }
      } else {
        cl.dn<-NA
      }
      cl<-c(cl.up, cl.dn)
    }
    if(directional == TRUE & GSEA_style == FALSE){
      if(is.null(direction_column)){
        stop("Error: No column providing the directionality of the enrichment was found.")
      } else {
        up<-enrichment_table[enrichment_table[,gene_set_name_column] %in% sig,]
        up<-up[grep("up", enrichment_table[,direction_column], ignore.case = T),gene_set_name_column]
        dn<-enrichment_table[enrichment_table[,gene_set_name_column] %in% sig,]
        dn<-dn[grep("dn|down", enrichment_table[,direction_column], ignore.case = T),gene_set_name_column]

        tempList<-list(up, dn)
        layout(matrix(c(1:2), ncol = sum(sapply(tempList, length) > 0)))

        if(length(up) > 0){
          correl.up<-correl[rownames(correl) %in% up , colnames(correl) %in% up]

          ig<-graph_from_adjacency_matrix(adjmatrix = correl.up,
                                          mode="upper",
                                          diag=F,
                                          weighted = T)
          ig<-delete_edges(ig, edges=which(E(ig)$weight <= correl_th))
          if(gsize(ig) == 0){
            stop("Error: no edges left after filtering for similarity threshold.")
          }

          cl.up <- membership(cluster_fast_greedy(ig))
          cl.up[cl.up%in%names(which(table(cl.up)==1))]<-"Unclustered"
          if(length(table(cl.up))>1){
            res2<-enrichment_table[enrichment_table[,gene_set_name_column]%in%names(cl.up),]
            res2<-res2[match(names(cl.up), res2[,gene_set_name_column]),]
            identical(names(cl.up), res2[,gene_set_name_column])
            V(ig)$Cluster<-cl.up
            V(ig)$FDR<--log10(ifelse(res2[,FDR_column]<1e-10, 1e-10, res2[,FDR_column]))
            set.seed(333)
            clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = c(distinctColorPalette(k=length(unique(V(ig)$Cluster))-1),"grey80"))
            names(clustercol)<-factor(V(ig)$Cluster)
            clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
            clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
            markg<-list()
            for(i in 1:length(unique(cl.up))){
              markg[[i]]<-which(cl.up%in%unique(cl.up)[i])
            }
            names(markg)<-unique(cl.up)
            markg<-markg[names(markg)!="Unclustered"]

            set.seed(123)
            coords<-layout_nicely(ig)
            plot.igraph(ig,
                        layout=coords,
                        vertex.size=V(ig)$FDR,
                        vertex.label=igraph.vertex.label,
                        vertex.label.cex=igraph.vertex.label.cex,
                        vertex.frame.color=igraph.vertex.frame.color,
                        vertex.color="red2",
                        edge.width = E(ig)$weight,
                        edge.color = igraph.edge.color,
                        mark.groups = markg,
                        mark.shape=0.9,
                        mark.expand=10,
                        mark.border = igraph.mark.border,
                        mark.col = clustercolLegend,
                        margin=c(0,0,0,0))
            title("Up-regulated",cex.main=2,col.main="black")
            legend("topleft",
                   legend = names(clustercolLegend)[names(clustercolLegend) != "Unclustered"],
                   fill = clustercolLegend[names(clustercolLegend) != "Unclustered"],
                   border = T,
                   title ="Clusters",
                   title.adj=0,
                   bty = "n",
                   cex = legend_cex)
          } else {
            res2<-enrichment_table[enrichment_table[,gene_set_name_column]%in%names(cl.up),]
            res2<-res2[match(names(cl.up), res2[,gene_set_name_column]),]
            identical(names(cl.up), res2[,gene_set_name_column])
            V(ig)$Cluster<-cl.up
            V(ig)$FDR<--log10(ifelse(res2[,FDR_column]<1e-10, 1e-10, res2[,FDR_column]))

            set.seed(123)
            coords<-layout_nicely(ig)
            plot.igraph(ig,
                        layout=coords,
                        vertex.size=V(ig)$FDR,
                        vertex.label=igraph.vertex.label,
                        vertex.label.cex=igraph.vertex.label.cex,
                        vertex.frame.color=igraph.vertex.frame.color,
                        vertex.color="red2",
                        edge.width = E(ig)$weight,
                        edge.color = igraph.edge.color,
                        margin=c(0,0,0,0))
            title("Up-regulated",cex.main=2,col.main="black")
          }
        } else {
          cl.up<-NA
        }

        if(length(dn) > 0){
          correl.dn<-correl[rownames(correl) %in% dn , colnames(correl) %in% dn]

          ig<-graph_from_adjacency_matrix(adjmatrix = correl.dn,
                                          mode="upper",
                                          diag=F,
                                          weighted = T)
          ig<-delete_edges(ig, edges=which(E(ig)$weight <= correl_th))
          if(gsize(ig) == 0){
            stop("Error: no edges left after filtering for similarity threshold.")
          }

          cl.dn <- membership(cluster_fast_greedy(ig))
          cl.dn[cl.dn%in%names(which(table(cl.dn)==1))]<-"Unclustered"
          if(length(table(cl.dn))>1){
            res2<-enrichment_table[enrichment_table[,gene_set_name_column]%in%names(cl.dn),]
            res2<-res2[match(names(cl.dn), res2[,gene_set_name_column]),]
            identical(names(cl.dn), res2[,gene_set_name_column])
            V(ig)$Cluster<-cl.dn
            V(ig)$FDR<--log10(ifelse(res2[,FDR_column]<1e-10, 1e-10, res2[,FDR_column]))
            set.seed(333)
            clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = c(distinctColorPalette(k=length(unique(V(ig)$Cluster))-1),"grey80"))
            names(clustercol)<-factor(V(ig)$Cluster)
            clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
            clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
            markg<-list()
            for(i in 1:length(unique(cl.dn))){
              markg[[i]]<-which(cl.dn%in%unique(cl.dn)[i])
            }
            names(markg)<-unique(cl.dn)
            markg<-markg[names(markg)!="Unclustered"]

            set.seed(123)
            coords<-layout_nicely(ig)
            plot.igraph(ig,
                        layout=coords,
                        vertex.size=V(ig)$FDR,
                        vertex.label=igraph.vertex.label,
                        vertex.label.cex=igraph.vertex.label.cex,
                        vertex.frame.color=igraph.vertex.frame.color,
                        vertex.color="royalblue",
                        edge.width = E(ig)$weight,
                        edge.color = igraph.edge.color,
                        mark.groups = markg,
                        mark.shape=0.9,
                        mark.expand=10,
                        mark.border = igraph.mark.border,
                        mark.col = clustercolLegend,
                        margin=c(0,0,0,0))
            title("Down-regulated",cex.main=2,col.main="black")
            legend("topleft",
                   legend = names(clustercolLegend)[names(clustercolLegend) != "Unclustered"],
                   fill = clustercolLegend[names(clustercolLegend) != "Unclustered"],
                   border = T,
                   title ="Clusters",
                   title.adj=0,
                   bty = "n",
                   cex = legend_cex)
          } else {
            res2<-enrichment_table[enrichment_table[,gene_set_name_column]%in%names(cl.dn),]
            res2<-res2[match(names(cl.dn), res2[,gene_set_name_column]),]
            identical(names(cl.dn), res2[,gene_set_name_column])
            V(ig)$Cluster<-cl.dn
            V(ig)$FDR<--log10(ifelse(res2[,FDR_column]<1e-10, 1e-10, res2[,FDR_column]))

            set.seed(123)
            coords<-layout_nicely(ig)
            plot.igraph(ig,
                        layout=coords,
                        vertex.size=V(ig)$FDR,
                        vertex.label=igraph.vertex.label,
                        vertex.label.cex=igraph.vertex.label.cex,
                        vertex.frame.color=igraph.vertex.frame.color,
                        vertex.color="royalblue",
                        edge.width = E(ig)$weight,
                        edge.color = igraph.edge.color,
                        margin=c(0,0,0,0))
            title("Down-regulated",cex.main=2,col.main="black")
          }
        } else {
          cl.dn<-NA
        }
        cl<-c(cl.up, cl.dn)
      }
    }
    community<-data.frame(geneSet = names(cl), Cluster = as.vector(cl))
    colnames(community)[1]<-gene_set_name_column
    enrichment_table<-merge(enrichment_table, community, by=gene_set_name_column, all.x=T)
    enrichment_table<-enrichment_table[order(enrichment_table$Cluster, enrichment_table[,FDR_column]),]
    return(enrichment_table)
  }
}

# example

# res<-read.table("/Users/dugo.matteo/Desktop/hsr_bioinfo/projects/NAPHER2/data_balazs/rnaseq_expression/results/analysis_8_p21/data_driven/gsea/reactome/fgsea_p21High_vs_p21Low_surgery_allArms_reactome.txt", header=T, sep="\t", as.is=T)
# load("/Users/dugo.matteo/Desktop/hsr_bioinfo/projects/cell_lines_PFHP/data/cell_lines_PFHP_tpm.RData")
# dataset<-cell_lines_PFHP_tpm[!is.na(fData(cell_lines_PFHP_tpm)$gene_name),]
# expmat<-exprs(dataset)
# expmat<-aggregate(expmat, by=list(gene_name=fData(dataset)$gene_name), sum)
# rownames(expmat)<-expmat$gene_name
# expmat<-as.matrix(expmat[,-1])
#
# cc<-correlateGeneSets(expMat = expmat,
#                   enrichment_table = res,
#                   FDR_column = "padj",
#                   gene_set_name_column = "pathway")



#**************************************************************************************
#* generate sankey plots
#* input is a data frame in wide format of paired observations
#**************************************************************************************

sankeyPlot<-function(df, node_width=0.2, node_space=10, node_color="black", node_alpha=0.4,
                     flow_smooth=8, node_label_size=4, node_label_color="black",
                     node_label_fill="white", n_size=4, x_label_size=13, x_label_color="black",
                     legend_name="", fill_colors=c("blue","red"), legend_text_size=12, legend_title_size=14,
                     plot_title="", plot_subtitle_size=13,plot_title_size=18,print_numbers=TRUE, factor.levels=NULL){

  library(dplyr)
  library(ggplot2)
  library(ggsankey)

  df <- df %>%
    make_long(colnames(df))
  if(is.null(factor.levels)){
    factor.levels<-union(levels(factor(df$node)), levels(factor(df$next_node)))
  }
  df$node<-factor(df$node,levels=factor.levels)
  df$next_node<-factor(df$next_node,levels=factor.levels)
  #df$node<-droplevels(df$node)
  #df$next_node<-droplevels(df$next_node)

  # prepare labels and labels positions

  stages<-levels(df$x)
  count.stages<-vector("list", length = length(stages))
  names(count.stages)<-stages
  for(i in 1:length(count.stages)){
    count.stages[[i]]<-table(df$node[df$x%in%stages[i]])
  }
  count.stages<-lapply(count.stages, function(x)x[x>0])

  lab.pos.stages<-lapply(count.stages, cumsum)
  for(i in 1:length(lab.pos.stages)){
    lab.pos.stages[[i]][1]<-lab.pos.stages[[i]][1]/2
    for(j in 2:length(lab.pos.stages[[i]])){
      lab.pos.stages[[i]][j]<-(lab.pos.stages[[i]][j] - (count.stages[[i]][j]/2)) + (node_space*(j-1))
    }
  }

  # calculate p-value

  tnxn<-matrix(0,ncol=length(stages), nrow=length(factor.levels))
  rownames(tnxn)<-factor.levels
  colnames(tnxn)<-stages
  for(i in 1:length(stages)){
    tnxn[,i]<-as.vector(table(df$node[df$x==stages[i]]))
  }
  if(sum(tnxn<=5)>0){
    if(ncol(tnxn)<3){
      p.test<-fisher.test(tnxn)$p.value
    } else {
      p.test<-fisher.test(tnxn, simulate.p.value=TRUE)$p.value
    }
  } else {
    p.test<-chisq.test(tnxn)$p.value
  }

  if(p.test<0.001){
    p.test<-format(p.test, digits=3, scientific=T)
  } else {
    p.test<-round(p.test,3)
  }

  # plot

  p<-ggplot(data=df, aes(x = x,
                         next_x = next_x,
                         node = node,
                         next_node = next_node,
                         fill=node)) +
    geom_sankey(smooth=flow_smooth,type="alluvial",
                node.color = node_color,
                alpha=node_alpha,
                space=node_space,
                width=node_width,
                show.legend = F) +
    theme_sankey() +
    theme(axis.text.x = element_text(size=x_label_size, color=x_label_color),
          legend.text = element_text(size=legend_text_size),
          legend.title = element_text(size=legend_title_size),
          plot.title = element_text(face="bold", size=plot_title_size, hjust=0.5),
          plot.subtitle = element_text(size=plot_subtitle_size, hjust=0.5)) +
    scale_fill_manual(name=legend_name,values=fill_colors) +
    xlab("") +
    ylab("") +
    labs(title=plot_title, subtitle = paste0("p-value = ",p.test))

  # print labels

  if(print_numbers){
    for(i in 1:length(lab.pos.stages)){
      for(j in 1:length(lab.pos.stages[[i]])){
        p <- p + geom_label(x = i, y = lab.pos.stages[[i]][j], label = paste0(names(lab.pos.stages[[i]])[j],"\nN = ",count.stages[[i]][j]), inherit.aes = F, fill = node_label_fill, color = node_label_color, size = node_label_size)
      }
    }
  } else {
    for(i in 1:length(lab.pos.stages)){
      for(j in 1:length(lab.pos.stages[[i]])){
        p <- p + geom_label(x = 1, y = lab.pos.stages[[i]][j], label = names(lab.pos.stages[[i]])[j], inherit.aes = F, fill = node_label_fill, color = node_label_color, size = node_label_size)
      }
    }
  }

  return(p)
  print(p)
}


#*************************************************************
# Univariate and multivariate logistic regression.
# Response must be a numeric vector where 0 indicates failure and 1 success.
# It an also be specified as a factor (when the first level denotes failure
# and all others success).
#*************************************************************

logistic_regression<-function(response, covariate1=NULL, covariate2=NULL, covariate3=NULL, covariate4=NULL, data, rounding_factor = 3, sort = FALSE){
  res.df<-data.frame(Variable = rownames(data),
                     Odds_Ratio_Variable = 0,
                     CI2.5_Variable = 0,
                     CI97.5_Variable = 0,
                     z_value_Variable = 0,
                     P_value_Variable = 0,
                     Odds_Ratio_Covariate1 = 0,
                     CI2.5_Covariate1 = 0,
                     CI97.5_Covariate1 = 0,
                     z_value_Covariate1 = 0,
                     P_value_Covariate1 = 0,
                     Odds_Ratio_Covariate2 = 0,
                     CI2.5_Covariate2 = 0,
                     CI97.5_Covariate2 = 0,
                     z_value_Covariate2 = 0,
                     P_value_Covariate2 = 0,
                     Odds_Ratio_Covariate3 = 0,
                     CI2.5_Covariate3 = 0,
                     CI97.5_Covariate3 = 0,
                     z_value_Covariate3 = 0,
                     P_value_Covariate3 = 0,
                     Odds_Ratio_Covariate4 = 0,
                     CI2.5_Covariate4 = 0,
                     CI97.5_Covariate4 = 0,
                     z_value_Covariate4 = 0,
                     P_value_Covariate4 = 0)
  for(i in 1:nrow(data)){
    if(is.null(covariate1)){
      covariate1<-rep(NA, ncol(data))
    }
    if(is.null(covariate2)){
      covariate2<-rep(NA, ncol(data))
    }
    if(is.null(covariate3)){
      covariate3<-rep(NA, ncol(data))
    }
    if(is.null(covariate4)){
      covariate4<-rep(NA, ncol(data))
    }
    df<-data.frame(Y=response,
                   Gene=data[i,],
                   Covariate1=covariate1,
                   Covariate2=covariate2,
                   Covariate3=covariate3,
                   Covariate4=covariate4,
                   stringsAsFactors = F)
    df<-df[,colSums(is.na(df))<nrow(df)]
    res.df<-res.df[,grep(paste(c("Variable",grep("Covariate",  colnames(df), value=T)), collapse="|"), colnames(res.df))]
    nCov<-grep("Covariate",  colnames(df), value=T)
    if(length(nCov) == 0){
      fit<-glm(Y ~ Gene, family = binomial, data = df)
      res.df[i,grep("_Variable", colnames(res.df))]<-c(round(exp(coef(fit))[2],rounding_factor),
                                                       round(exp(confint.default(fit, level = 0.95))[2,1],rounding_factor),
                                                       round(exp(confint.default(fit, level = 0.95))[2,2],rounding_factor),
                                                       summary(fit)$coefficients[2,3],
                                                       summary(fit)$coefficients[2,4])
    }
    if(length(nCov) > 0){
      f1 <- as.formula(paste("Y ~ Gene +", paste(grep("Covariate", colnames(df), value=T), collapse="+")))
      fit<-glm(f1, family = binomial, data = df)
      res.df[i,grep("_Variable", colnames(res.df))]<-c(round(exp(coef(fit))[2],rounding_factor),
                                                       round(exp(confint.default(fit, level = 0.95))[2,1],rounding_factor),
                                                       round(exp(confint.default(fit, level = 0.95))[2,2],rounding_factor),
                                                       summary(fit)$coefficients[2,3],
                                                       summary(fit)$coefficients[2,4])
      k=3
      for(j in 1:length(nCov)){
        res.df[i,grep(nCov[j], colnames(res.df))]<-c(round(exp(coef(fit))[k],rounding_factor),
                                                     round(exp(confint.default(fit, level = 0.95))[k,1],rounding_factor),
                                                     round(exp(confint.default(fit, level = 0.95))[k,2],rounding_factor),
                                                     summary(fit)$coefficients[k,3],
                                                     summary(fit)$coefficients[k,4])
        k<-k+1
      }
    }
  }
  if(length(nCov) > 0){
    fdr<-matrix(0, ncol=length(nCov), nrow=nrow(res.df))
    colnames(fdr)<-paste("FDR",nCov,sep="_")
    res.df$FDR_Variable<-p.adjust(res.df$P_value_Variable, method="BH")
    for(z in 1:ncol(fdr)){
      fdr[,z]<-p.adjust(res.df[,grep(paste0("P_value_", nCov[z]), colnames(res.df))], method="BH")
    }
    res.df<-as.data.frame(cbind(res.df, fdr))
    colOrder<-grep("Variable", colnames(res.df), value=T)
    for(z in 1:length(nCov)){
      colOrder<-c(colOrder, grep(nCov[z], colnames(res.df), value=T))
    }
    res.df<-res.df[,colOrder]
  } else {
    res.df$FDR_Variable<-p.adjust(res.df$P_value_Variable, method="BH")
  }
  if(sort){
    res.df<-res.df[order(res.df$P_value_Variable),]
  }
  return(res.df)
}

#*************************************************************
# Interaction test for logistic regression.
# Response must be a numeric vector where 0 indicates failure and 1 success.
# It an also be specified as a factor (when the first level denotes failure
# and all others success).
#*************************************************************

interaction_logistic_regression<-function(response, covariate, data, rounding_factor = 3, sort = FALSE){
  res.df<-data.frame(Variable = rownames(data),
                     Odds_Ratio = 0,
                     CI2.5 = 0,
                     CI97.5 = 0,
                     z_value=0,
                     P_value = 0)
  for(i in 1:nrow(data)){
    df<-data.frame(Y=response,
                   Gene=data[i,],
                   Covariate=covariate,
                   stringsAsFactors = F)
    fit<-glm(Y ~ Gene*Covariate, family = binomial, data = df)
    res.df[i,2:ncol(res.df)]<-c(round(exp(coef(fit))[4],rounding_factor),
                                round(exp(confint.default(fit, level = 0.95))[4,1],rounding_factor),
                                round(exp(confint.default(fit, level = 0.95))[4,2],rounding_factor),
                                summary(fit)$coefficients[4,3],
                                summary(fit)$coefficients[4,4])
  }
  res.df$FDR<-p.adjust(res.df$P_value, method="BH")
  if(sort){
    res.df<-res.df[order(res.df$P_value),]
  }
  return(res.df)
}





