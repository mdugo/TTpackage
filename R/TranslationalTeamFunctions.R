
#*************************************************************
#* This function performs singscore calculation for gene set
#* collections that include gene sets with only up-regulated
#* genes or both up- and down-regulated genes
#*************************************************************

simpleScoreMod<-function(rankData, mysetlist, knownDir = TRUE, centerScore = TRUE, minSize = 0){

  require(Biobase)
  require(singscore)

  score<-matrix(0, nrow=length(mysetlist),ncol=ncol(rankData))
  rownames(score)<-names(mysetlist)
  colnames(score)<-colnames(rankData)
  for(i in 1:length(mysetlist)){
    if(class(mysetlist[[i]])=="character"){
      sl<-mysetlist[[i]]
      sl<-sl[sl%in%rownames(rankData)]
      if(length(sl)>=minSize){
        scoretemp<-simpleScore(rankData,upSet = sl, knownDirection = knownDir, centerScore = centerScore)
        score[rownames(score)==names(mysetlist)[i],]<-scoretemp$TotalScore
      }
    } else {
      sl.up<-mysetlist[[i]][[grep("UP",names(mysetlist[[i]]),ignore.case=T)]]
      sl.up<-sl.up[sl.up%in%rownames(rankData)]
      sl.dn<-mysetlist[[i]][[grep("DN|DOWN",names(mysetlist[[i]]),ignore.case=T)]]
      sl.dn<-sl.dn[sl.dn%in%rownames(rankData)]
      if(length(sl.up)>=minSize & length(sl.dn)>=minSize){
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
                              correlation_method = c("pearson", "spearman"),
                              correl_th=0.9, FDR_th=0.05, FDR_column,
                              gene_set_name_column, NES_column="NES", NES_th=0,
                              direction_column = NULL,
                              gene_content_column=NULL,
                              igraph.vertex.label = NA, igraph.vertex.label.cex = 0.2,
                              igraph.vertex.size.cex = 1,
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
    print("No significant gene sets found at provided thresholds.")
  } else {

    if(cluster_by == "gene_content"){
      if(is.null(gene_content_column)){

        ### retrieving gene set content

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
        print("No edges left after filtering for similarity threshold.")
      } else {

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
            clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = c(c(carto_pal(12,"Bold")[-12],carto_pal(12,"Safe")[-12])[1:(length(table(V(ig)$Cluster))-1)],"grey"))
            names(clustercol)<-factor(V(ig)$Cluster)
            clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
            clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
            set.seed(123)
            coords<-layout_nicely(ig)
            plot.igraph(ig,
                        layout=coords,
                        vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                        vertex.label=igraph.vertex.label,
                        vertex.label.cex=igraph.vertex.label.cex,
                        vertex.frame.color=igraph.vertex.frame.color,
                        vertex.color=clustercol,
                        edge.width = E(ig)$weight,
                        edge.color = igraph.edge.color,
                        margin=c(0,0,0,0),
                        main=paste0("Gene content\n(Method = ",similarity_method, "; Threshold = ",similarity_th,")"))
            legend("topleft",
                   legend = names(clustercolLegend)[names(clustercolLegend) != "Unclustered"],
                   fill = clustercolLegend[names(clustercolLegend) != "Unclustered"],
                   border = T,
                   title ="Clusters",
                   title.adj=0,
                   bty = "n",
                   cex = legend_cex)
          }

          if(directional == TRUE & GSEA_style == TRUE){
            V(ig)$Direction<-res2[,NES_column]
            set.seed(333)
            clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = c(c(carto_pal(12,"Bold")[-12],carto_pal(12,"Safe")[-12])[1:(length(table(V(ig)$Cluster))-1)],"grey"))
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
                        vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
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
                        margin=c(0,0,0,0),
                        main=paste0("Gene content\n(Method = ",similarity_method, "; Threshold = ",similarity_th,")"))
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
              clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = c(c(carto_pal(12,"Bold")[-12],carto_pal(12,"Safe")[-12])[1:(length(table(V(ig)$Cluster))-1)],"grey"))
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
                          vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
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
                          margin=c(0,0,0,0),
                          main=paste0("Gene content\n(Method = ",similarity_method, "; Threshold = ",similarity_th,")"))
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
        }

        if(directional == FALSE){
          set.seed(123)
          coords<-layout_nicely(ig)
          plot.igraph(ig,
                      layout=coords,
                      vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                      vertex.label=igraph.vertex.label,
                      vertex.label.cex=igraph.vertex.label.cex,
                      vertex.frame.color=igraph.vertex.frame.color,
                      vertex.color="grey80",
                      edge.width = E(ig)$weight,
                      edge.color = igraph.edge.color,
                      margin=c(0,0,0,0),
                      main=paste0("Gene content\n(Method = ",similarity_method, "; Threshold = ",similarity_th,")"))
        }
        if(directional == TRUE & GSEA_style == TRUE){
          V(ig)$Direction<-res2[,NES_column]
          NEScol<-colorRamp2(breaks = c(min(as.numeric(V(ig)$Direction), na.rm = T), 0, as.numeric(max(V(ig)$Direction), na.rm = T)), colors = c("royalblue", "white", "red2"))
          NEScol<-NEScol(as.numeric(V(ig)$Direction))

          set.seed(123)
          coords<-layout_nicely(ig)
          plot.igraph(ig,
                      layout=coords,
                      vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                      vertex.label=igraph.vertex.label,
                      vertex.label.cex=igraph.vertex.label.cex,
                      vertex.frame.color=igraph.vertex.frame.color,
                      vertex.color=NEScol,
                      edge.width = E(ig)$weight,
                      edge.color = igraph.edge.color,
                      margin=c(0,0,0,0),
                      main=paste0("Gene content\n(Method = ",similarity_method, "; Threshold = ",similarity_th,")"))
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
                        vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                        vertex.label=igraph.vertex.label,
                        vertex.label.cex=igraph.vertex.label.cex,
                        vertex.frame.color=igraph.vertex.frame.color,
                        vertex.color=nodecol,
                        edge.width = E(ig)$weight,
                        edge.color = igraph.edge.color,
                        margin=c(0,0,0,0),
                        main=paste0("Gene content\n(Method = ",similarity_method, "; Threshold = ",similarity_th,")"))
          }
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
                                    knownDir = T,
                                    centerScore = T)
          scoreMat<-rbind(scoreMat1, scoreMat2)
        } else {
          scoreMat<-scoreMat1
        }
      } else {
        scoreMat<-simpleScoreMod(rankData = rankMat,
                                 mysetlist = setlist,
                                 knownDir = T,
                                 centerScore = T)
      }
      print("Done.")
      rownames(scoreMat)<-gsub("CUSTOM_", "", rownames(scoreMat))

      correl<-cor(t(scoreMat), method=correlation_method)

      if(directional == FALSE){
        ig<-graph_from_adjacency_matrix(adjmatrix = correl,
                                        mode="upper",
                                        diag=F,
                                        weighted = T)
        ig<-delete_edges(ig, edges=which(E(ig)$weight <= correl_th))
        if(gsize(ig) == 0){
          print("No edges left after filtering for similarity threshold.")
        } else {

          cl = membership(cluster_fast_greedy(ig))
          cl[cl%in%names(which(table(cl)==1))]<-"Unclustered"

          res2<-enrichment_table[enrichment_table[,gene_set_name_column]%in%names(cl),]
          res2<-res2[match(names(cl), res2[,gene_set_name_column]),]
          identical(names(cl), res2[,gene_set_name_column])
          V(ig)$Cluster<-cl
          V(ig)$FDR<--log10(ifelse(res2[,FDR_column]<1e-10, 1e-10, res2[,FDR_column]))

          if(length(table(cl))>1){
            set.seed(444)
            clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = c(c(carto_pal(12,"Bold")[-12],carto_pal(12,"Safe")[-12])[1:(length(table(V(ig)$Cluster))-1)],"grey"))
            names(clustercol)<-factor(V(ig)$Cluster)
            clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
            clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
            set.seed(123)
            coords<-layout_nicely(ig)
            plot.igraph(ig,
                        layout=coords,
                        vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                        vertex.label=igraph.vertex.label,
                        vertex.label.cex=igraph.vertex.label.cex,
                        vertex.frame.color=igraph.vertex.frame.color,
                        vertex.color=clustercol,
                        edge.width = E(ig)$weight,
                        edge.color = igraph.edge.color,
                        margin=c(0,0,0,0),
                        main=paste0("Correlation\n(Method = ",correlation_method, "; Threshold = ",correl_th,")"))
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
                        vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                        vertex.label=igraph.vertex.label,
                        vertex.label.cex=igraph.vertex.label.cex,
                        vertex.frame.color=igraph.vertex.frame.color,
                        vertex.color="grey80",
                        edge.width = E(ig)$weight,
                        edge.color = igraph.edge.color,
                        margin=c(0,0,0,0),
                        main=paste0("Correlation\n(Method = ",correlation_method, "; Threshold = ",correl_th,")"))
          }
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
            print("No edges left after filtering for similarity threshold.")
            cl.up<-NA
            plot(0:10, 0:10, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", bty = "n")
            text(x = 5, y = 5, labels = "No edges left after\nfiltering for similarity threshold", cex=1.4)
          } else {

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
              clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = c(c(carto_pal(12,"Bold")[-12],carto_pal(12,"Safe")[-12])[1:(length(table(V(ig)$Cluster))-1)],"grey"))
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
                          vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
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
                          margin=c(0,0,0,0),
                          main=paste0("Correlation\n(Method = ",correlation_method, "; Threshold = ",correl_th,")"))
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
                          vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                          vertex.label=igraph.vertex.label,
                          vertex.label.cex=igraph.vertex.label.cex,
                          vertex.frame.color=igraph.vertex.frame.color,
                          vertex.color=NEScol,
                          edge.width = E(ig)$weight,
                          edge.color = igraph.edge.color,
                          margin=c(0,0,0,0),
                          main=paste0("Correlation\n(Method = ",correlation_method, "; Threshold = ",correl_th,")"))
              title("Up-regulated",cex.main=2,col.main="black")
            }
          }
        } else {
          cl.up<-NA
          plot(0:10, 0:10, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", bty = "n")
          text(x = 5, y = 5, labels = "No positively\nenriched gene sets", cex=1.4)
        }

        if(length(dn) > 0){
          correl.dn<-correl[rownames(correl) %in% dn , colnames(correl) %in% dn]
          ig<-graph_from_adjacency_matrix(adjmatrix = correl.dn,
                                          mode="upper",
                                          diag=F,
                                          weighted = T)
          ig<-delete_edges(ig, edges=which(E(ig)$weight <= correl_th))
          if(gsize(ig) == 0){
            print("No edges left after filtering for similarity threshold.")
            cl.dn<-NA
            plot(0:10, 0:10, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", bty = "n")
            text(x = 5, y = 5, labels = "No edges left after\nfiltering for similarity threshold", cex=1.4)
          } else {

            cl.dn <- membership(cluster_fast_greedy(ig))
            cl.dn[cl.dn%in%names(which(table(cl.dn)==1))]<-"Unclustered"

            if(length(table(cl.dn))>1){
              res2<-enrichment_table[enrichment_table[,gene_set_name_column]%in%names(cl.dn),]
              res2<-res2[match(names(cl.dn), res2[,gene_set_name_column]),]
              V(ig)$Cluster<-cl.dn
              V(ig)$FDR<--log10(ifelse(res2[,FDR_column]<1e-10, 1e-10, res2[,FDR_column]))
              V(ig)$Direction<-res2[,NES_column]
              set.seed(343)
              clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = c(c(carto_pal(12,"Bold")[-12],carto_pal(12,"Safe")[-12])[1:(length(table(V(ig)$Cluster))-1)],"grey"))
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
                          vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
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
                          margin=c(0,0,0,0),
                          main=paste0("Correlation\n(Method = ",correlation_method, "; Threshold = ",correl_th,")"))
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
                          vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                          vertex.label=igraph.vertex.label,
                          vertex.label.cex=igraph.vertex.label.cex,
                          vertex.frame.color=igraph.vertex.frame.color,
                          vertex.color=NEScol,
                          edge.width = E(ig)$weight,
                          edge.color = igraph.edge.color,
                          margin=c(0,0,0,0),
                          main=paste0("Correlation\n(Method = ",correlation_method, "; Threshold = ",correl_th,")"))
              title("Down-regulated",cex.main=2,col.main="black")
            }
          }
        } else {
          cl.dn<-NA
          plot(0:10, 0:10, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", bty = "n")
          text(x = 5, y = 5, labels = "No negatively\nenriched gene sets", cex=1.4)
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
              print("No edges left after filtering for similarity threshold.")
              cl.up<-NA
              plot(0:10, 0:10, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", bty = "n")
              text(x = 5, y = 5, labels = "No edges left after\nfiltering for similarity threshold", cex=1.4)
            } else {

              cl.up <- membership(cluster_fast_greedy(ig))
              cl.up[cl.up%in%names(which(table(cl.up)==1))]<-"Unclustered"
              if(length(table(cl.up))>1){
                res2<-enrichment_table[enrichment_table[,gene_set_name_column]%in%names(cl.up),]
                res2<-res2[match(names(cl.up), res2[,gene_set_name_column]),]
                identical(names(cl.up), res2[,gene_set_name_column])
                V(ig)$Cluster<-cl.up
                V(ig)$FDR<--log10(ifelse(res2[,FDR_column]<1e-10, 1e-10, res2[,FDR_column]))
                set.seed(333)
                clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = c(c(carto_pal(12,"Bold")[-12],carto_pal(12,"Safe")[-12])[1:(length(table(V(ig)$Cluster))-1)],"grey"))
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
                            vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
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
                            margin=c(0,0,0,0),
                            main=paste0("Correlation\n(Method = ",correlation_method, "; Threshold = ",correl_th,")"))
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
                            vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                            vertex.label=igraph.vertex.label,
                            vertex.label.cex=igraph.vertex.label.cex,
                            vertex.frame.color=igraph.vertex.frame.color,
                            vertex.color="red2",
                            edge.width = E(ig)$weight,
                            edge.color = igraph.edge.color,
                            margin=c(0,0,0,0),
                            main=paste0("Correlation\n(Method = ",correlation_method, "; Threshold = ",correl_th,")"))
                title("Up-regulated",cex.main=2,col.main="black")
              }
            }
          } else {
            cl.up<-NA
            plot(0:10, 0:10, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", bty = "n")
            text(x = 5, y = 5, labels = "No positively\nenriched gene sets", cex=1.4)
          }

          if(length(dn) > 0){
            correl.dn<-correl[rownames(correl) %in% dn , colnames(correl) %in% dn]
            ig<-graph_from_adjacency_matrix(adjmatrix = correl.dn,
                                            mode="upper",
                                            diag=F,
                                            weighted = T)
            ig<-delete_edges(ig, edges=which(E(ig)$weight <= correl_th))
            if(gsize(ig) == 0){
              print("No edges left after filtering for similarity threshold.")
              cl.dn<-NA
              plot(0:10, 0:10, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", bty = "n")
              text(x = 5, y = 5, labels = "No edges left after\nfiltering for similarity threshold", cex=1.4)
            } else {

              cl.dn <- membership(cluster_fast_greedy(ig))
              cl.dn[cl.dn%in%names(which(table(cl.dn)==1))]<-"Unclustered"
              if(length(table(cl.dn))>1){
                res2<-enrichment_table[enrichment_table[,gene_set_name_column]%in%names(cl.dn),]
                res2<-res2[match(names(cl.dn), res2[,gene_set_name_column]),]
                identical(names(cl.dn), res2[,gene_set_name_column])
                V(ig)$Cluster<-cl.dn
                V(ig)$FDR<--log10(ifelse(res2[,FDR_column]<1e-10, 1e-10, res2[,FDR_column]))
                set.seed(333)
                clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = c(c(carto_pal(12,"Bold")[-12],carto_pal(12,"Safe")[-12])[1:(length(table(V(ig)$Cluster))-1)],"grey"))
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
                            vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
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
                            margin=c(0,0,0,0),
                            main=paste0("Correlation\n(Method = ",correlation_method, "; Threshold = ",correl_th,")"))
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
                            vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                            vertex.label=igraph.vertex.label,
                            vertex.label.cex=igraph.vertex.label.cex,
                            vertex.frame.color=igraph.vertex.frame.color,
                            vertex.color="royalblue",
                            edge.width = E(ig)$weight,
                            edge.color = igraph.edge.color,
                            margin=c(0,0,0,0),
                            main=paste0("Correlation\n(Method = ",correlation_method, "; Threshold = ",correl_th,")"))
                title("Down-regulated",cex.main=2,col.main="black")
              }
            }
          } else {
            cl.dn<-NA
            plot(0:10, 0:10, xaxt = "n", yaxt = "n", xlab = "", ylab = "", col = "white", bty = "n")
            text(x = 5, y = 5, labels = "No negatively\nenriched gene sets", cex=1.4)
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

logistic_regression<-function(response, covariates=NULL, data, rounding_factor = 3, sort = FALSE){
  res.df<-data.frame(Variable = rownames(data), matrix(0, nrow = nrow(data), ncol=(length(covariates)+1)*5))
  colnames(res.df)[-1]<-paste(rep(c("OR","CI2.5","CI97.5","Z","P"),length(covariates)+1),rep(c("Variable",covariates), each=5),sep="_")

  Y<-pData(data)[,response]
  covs<-data.frame(pData(data)[,covariates])
  if(is.null(covariates)){
    resTemp<-apply(exprs(data), 1, function(x){
      z<-glm(Y ~ x, family = binomial)
      ee<-c(round(exp(coef(z))[2],rounding_factor),
            round(exp(confint.default(z, level = 0.95))[2,1],rounding_factor),
            round(exp(confint.default(z, level = 0.95))[2,2],rounding_factor),
            summary(z)$coefficients[2,3],
            summary(z)$coefficients[2,4])
    })
    res.df[,-1]<-t(resTemp)
    res.df$FDR_Variable<-p.adjust(res.df$P_Variable, method="BH")
  }
  if(!is.null(covariates)){
    colnames(covs)<-covariates
    resTemp<-apply(exprs(data), 1, function(x){
      df<-data.frame(Y,
                     Gene=x,
                     covs,
                     stringsAsFactors = F)
      f1 <- as.formula(paste("Y ~ Gene +", paste(colnames(df)[3:ncol(df)], collapse="+")))
      z<-glm(f1, family = binomial, data = df)
      ee<-matrix(as.vector(rbind(round(exp(coef(z))[-1],rounding_factor),
                                 round(exp(confint.default(z, level = 0.95))[-1,1],rounding_factor),
                                 round(exp(confint.default(z, level = 0.95))[-1,2],rounding_factor),
                                 summary(z)$coefficients[-1,3],
                                 summary(z)$coefficients[-1,4])),nrow=1)
    })
    res.df[,-1]<-t(resTemp)
    fdr<-apply(res.df[,grep("^P_",colnames(res.df))], 2, p.adjust, method="BH")
    colnames(fdr)<-gsub("P_","FDR_",colnames(fdr))
    res.df<-as.data.frame(cbind(res.df, fdr))
    index<-unname(unlist(sapply(c("Variable",covariates), function(x){
      grep(x, colnames(res.df))})))
    res.df<-res.df[,index]
  }
  if(sort){
    res.df<-res.df[order(res.df$P_Variable),]
  }
  return(res.df)
}

#*************************************************************
# Interaction test for logistic regression.
# Response must be a numeric vector where 0 indicates failure and 1 success.
# It an also be specified as a factor (when the first level denotes failure
# and all others success).
#*************************************************************


interaction_logistic_regression<-function(response, covariate_interaction, covariates_adjustment = NULL, data, rounding_factor = 3, sort = FALSE){
  res.df<-data.frame(Variable = rownames(data), matrix(0, nrow = nrow(data), ncol=(length(covariates_adjustment)+3)*5))
  colnames(res.df)[-1]<-paste(rep(c("OR","CI2.5","CI97.5","Z","P"),length(covariates_adjustment)+3),rep(c("Variable",covariate_interaction,"Interaction",covariates_adjustment), each=5),sep="_")

  for(i in 1:nrow(data)){
    df<-data.frame(Y=pData(data)[,response],
                   Gene=exprs(data)[i,],
                   Covariate_interaction=pData(data)[,covariate_interaction],
                   pData(data)[,covariates_adjustment],
                   stringsAsFactors = F)
    colnames(df)[grep("pData",colnames(df))]<-covariates_adjustment
    if(is.null(covariates_adjustment)){

      fit<-glm(Y ~ Gene*Covariate_interaction, family = binomial, data = df)
      res.df[i,grep("_Variable", colnames(res.df))]<-c(round(exp(coef(fit))[2],rounding_factor),
                                                       round(exp(confint.default(fit, level = 0.95))[2,1],rounding_factor),
                                                       round(exp(confint.default(fit, level = 0.95))[2,2],rounding_factor),
                                                       summary(fit)$coefficients[2,3],
                                                       summary(fit)$coefficients[2,4])
      res.df[i,grep(covariate_interaction, colnames(res.df))]<-c(round(exp(coef(fit))[3],rounding_factor),
                                                                 round(exp(confint.default(fit, level = 0.95))[3,1],rounding_factor),
                                                                 round(exp(confint.default(fit, level = 0.95))[3,2],rounding_factor),
                                                                 summary(fit)$coefficients[3,3],
                                                                 summary(fit)$coefficients[3,4])
      res.df[i,grep("Interaction", colnames(res.df))]<-c(round(exp(coef(fit))[4],rounding_factor),
                                                         round(exp(confint.default(fit, level = 0.95))[4,1],rounding_factor),
                                                         round(exp(confint.default(fit, level = 0.95))[4,2],rounding_factor),
                                                         summary(fit)$coefficients[4,3],
                                                         summary(fit)$coefficients[4,4])
    }
    if(!is.null(covariates_adjustment)){
      f1 <- as.formula(paste("Y ~ Gene * Covariate_interaction +", paste(colnames(df)[4:ncol(df)], collapse="+")))
      fit<-glm(f1, family = binomial, data = df)
      res.df[i,grep("_Variable", colnames(res.df))]<-c(round(exp(coef(fit))[2],rounding_factor),
                                                       round(exp(confint.default(fit, level = 0.95))[2,1],rounding_factor),
                                                       round(exp(confint.default(fit, level = 0.95))[2,2],rounding_factor),
                                                       summary(fit)$coefficients[2,3],
                                                       summary(fit)$coefficients[2,4])
      res.df[i,grep(covariate_interaction, colnames(res.df))]<-c(round(exp(coef(fit))[3],rounding_factor),
                                                                 round(exp(confint.default(fit, level = 0.95))[3,1],rounding_factor),
                                                                 round(exp(confint.default(fit, level = 0.95))[3,2],rounding_factor),
                                                                 summary(fit)$coefficients[3,3],
                                                                 summary(fit)$coefficients[3,4])
      res.df[i,grep("Interaction", colnames(res.df))]<-c(round(exp(coef(fit))[4+length(covariates_adjustment)],rounding_factor),
                                                         round(exp(confint.default(fit, level = 0.95))[4+length(covariates_adjustment),1],rounding_factor),
                                                         round(exp(confint.default(fit, level = 0.95))[4+length(covariates_adjustment),2],rounding_factor),
                                                         summary(fit)$coefficients[4+length(covariates_adjustment),3],
                                                         summary(fit)$coefficients[4+length(covariates_adjustment),4])
      k=4
      for(j in 1:length(covariates_adjustment)){
        res.df[i,grep(covariates_adjustment[j], colnames(res.df))]<-c(round(exp(coef(fit))[k],rounding_factor),
                                                                      round(exp(confint.default(fit, level = 0.95))[k,1],rounding_factor),
                                                                      round(exp(confint.default(fit, level = 0.95))[k,2],rounding_factor),
                                                                      summary(fit)$coefficients[k,3],
                                                                      summary(fit)$coefficients[k,4])
        k<-k+1
      }
    }
  }
  fdr<-res.df[,grep("^P_",colnames(res.df))]
  colnames(fdr)<-gsub("^P_","FDR_",colnames(fdr))
  fdr<-apply(fdr, 2, p.adjust, "BH")
  res.df<-as.data.frame(cbind(res.df, fdr))
  colOrder<-grep("Variable", colnames(res.df), value=T)
  colOrder<-c(colOrder,grep(covariate_interaction, colnames(res.df), value=T))
  colOrder<-c(colOrder,grep("Interaction", colnames(res.df), value=T))

  if(!is.null(covariates_adjustment)){
    for(z in 1:length(covariates_adjustment)){
      colOrder<-c(colOrder, grep(covariates_adjustment[z], colnames(res.df), value=T))
    }
  }
  res.df<-res.df[,colOrder]
  if(sort){
    res.df<-res.df[order(res.df$P_Interaction),]
  }
  return(res.df)
}



#****************************************************************************
# Multivariate logistic regression for paired gene expression data
# Response must be a numeric vector where 0 indicates failure and 1 success.
# It an also be specified as a factor (when the first level denotes failure
# and all others success).
#****************************************************************************


paired_logistic_regression_additive<-function(response, data, data2, covariates_adjustment=NULL, suffixes, rounding_factor = 3){
  res.df<-data.frame(Variable = rownames(data), matrix(0, nrow = nrow(data), ncol=(length(covariates_adjustment)+2)*5))
  colnames(res.df)[-1]<-paste(rep(c("OR","CI2.5","CI97.5","Z","P"),length(covariates_adjustment)+2),rep(c(suffixes,covariates_adjustment), each=5),sep="_")

  for(i in 1:nrow(data)){
    df<-data.frame(Y=pData(data)[,response],
                   Gene_1=exprs(data)[i,],
                   Gene_2=exprs(data2)[i,],
                   pData(data)[,covariates_adjustment],
                   stringsAsFactors = F)
    colnames(df)[grep("pData",colnames(df))]<-covariates_adjustment
    if(is.null(covariates_adjustment)){
      fit<-glm(Y ~ Gene_1 + Gene_2, family = binomial, data = df)
      res.df[i,grep(suffixes[1],colnames(res.df))]<-c(round(exp(coef(fit))[2],rounding_factor),
                                                      round(exp(confint.default(fit, level = 0.95))[2,1],rounding_factor),
                                                      round(exp(confint.default(fit, level = 0.95))[2,2],rounding_factor),
                                                      summary(fit)$coefficients[2,3],
                                                      summary(fit)$coefficients[2,4])
      res.df[i,grep(suffixes[2],colnames(res.df))]<-c(round(exp(coef(fit))[3],rounding_factor),
                                                      round(exp(confint.default(fit, level = 0.95))[3,1],rounding_factor),
                                                      round(exp(confint.default(fit, level = 0.95))[3,2],rounding_factor),
                                                      summary(fit)$coefficients[3,3],
                                                      summary(fit)$coefficients[3,4])
    }
    if(!is.null(covariates_adjustment)){
      f1 <- as.formula(paste("Y ~ Gene_1 + Gene_2 +", paste(colnames(df)[3:ncol(df)], collapse="+")))
      fit<-glm(f1, family = binomial, data = df)
      res.df[i,grep(suffixes[1], colnames(res.df))]<-c(round(exp(coef(fit))[2],rounding_factor),
                                                       round(exp(confint.default(fit, level = 0.95))[2,1],rounding_factor),
                                                       round(exp(confint.default(fit, level = 0.95))[2,2],rounding_factor),
                                                       summary(fit)$coefficients[2,3],
                                                       summary(fit)$coefficients[2,4])
      res.df[i,grep(suffixes[2], colnames(res.df))]<-c(round(exp(coef(fit))[3],rounding_factor),
                                                       round(exp(confint.default(fit, level = 0.95))[3,1],rounding_factor),
                                                       round(exp(confint.default(fit, level = 0.95))[3,2],rounding_factor),
                                                       summary(fit)$coefficients[3,3],
                                                       summary(fit)$coefficients[3,4])
      k=4
      for(j in 1:length(covariates_adjustment)){
        res.df[i,grep(covariates_adjustment[j], colnames(res.df))]<-c(round(exp(coef(fit))[k],rounding_factor),
                                                                      round(exp(confint.default(fit, level = 0.95))[k,1],rounding_factor),
                                                                      round(exp(confint.default(fit, level = 0.95))[k,2],rounding_factor),
                                                                      summary(fit)$coefficients[k,3],
                                                                      summary(fit)$coefficients[k,4])
        k<-k+1
      }
    }
  }
  fdr<-res.df[,grep("^P_",colnames(res.df))]
  colnames(fdr)<-gsub("^P_","FDR_",colnames(fdr))
  fdr<-apply(fdr, 2, p.adjust, "BH")
  res.df<-as.data.frame(cbind(res.df, fdr))
  colOrder<-c("Variable",grep(suffixes[1], colnames(res.df), value=T))
  colOrder<-c(colOrder,grep(suffixes[2], colnames(res.df), value=T))

  if(!is.null(covariates_adjustment)){
    for(z in 1:length(covariates_adjustment)){
      colOrder<-c(colOrder, grep(covariates_adjustment[z], colnames(res.df), value=T))
    }
  }
  res.df<-res.df[,colOrder]
  return(res.df)
}

#****************************************************************************
# Multivariate logistic regression for paired gene expression data
# Response must be a numeric vector where 0 indicates failure and 1 success.
# It an also be specified as a factor (when the first level denotes failure
# and all others success).
# This function consider interaction between the paired data
#****************************************************************************

paired_logistic_regression_interaction<-function(response, data, data2, covariates_adjustment=NULL, suffixes, rounding_factor = 3){
  res.df<-data.frame(Variable = rownames(data), matrix(0, nrow = nrow(data), ncol=(length(covariates_adjustment)+3)*5))
  colnames(res.df)[-1]<-paste(rep(c("OR","CI2.5","CI97.5","Z","P"),length(covariates_adjustment)+3),rep(c(suffixes,"Interaction",covariates_adjustment), each=5),sep="_")

  for(i in 1:nrow(data)){
    df<-data.frame(Y=pData(data)[,response],
                   Gene_1=exprs(data)[i,],
                   Gene_2=exprs(data2)[i,],
                   pData(data)[,covariates_adjustment],
                   stringsAsFactors = F)
    colnames(df)[grep("pData",colnames(df))]<-covariates_adjustment
    if(is.null(covariates_adjustment)){
      fit<-glm(Y ~ Gene_1 * Gene_2, family = binomial, data = df)
      res.df[i,grep(suffixes[1],colnames(res.df))]<-c(round(exp(coef(fit))[2],rounding_factor),
                                                      round(exp(confint.default(fit, level = 0.95))[2,1],rounding_factor),
                                                      round(exp(confint.default(fit, level = 0.95))[2,2],rounding_factor),
                                                      summary(fit)$coefficients[2,3],
                                                      summary(fit)$coefficients[2,4])
      res.df[i,grep(suffixes[2],colnames(res.df))]<-c(round(exp(coef(fit))[3],rounding_factor),
                                                      round(exp(confint.default(fit, level = 0.95))[3,1],rounding_factor),
                                                      round(exp(confint.default(fit, level = 0.95))[3,2],rounding_factor),
                                                      summary(fit)$coefficients[3,3],
                                                      summary(fit)$coefficients[3,4])
      res.df[i,grep("Interaction",colnames(res.df))]<-c(round(exp(coef(fit))[4],rounding_factor),
                                                        round(exp(confint.default(fit, level = 0.95))[4,1],rounding_factor),
                                                        round(exp(confint.default(fit, level = 0.95))[4,2],rounding_factor),
                                                        summary(fit)$coefficients[4,3],
                                                        summary(fit)$coefficients[4,4])
    }
    if(!is.null(covariates_adjustment)){
      f1 <- as.formula(paste("Y ~ Gene_1 * Gene_2 +", paste(colnames(df)[3:ncol(df)], collapse="+")))
      fit<-glm(f1, family = binomial, data = df)
      res.df[i,grep(suffixes[1], colnames(res.df))]<-c(round(exp(coef(fit))[2],rounding_factor),
                                                       round(exp(confint.default(fit, level = 0.95))[2,1],rounding_factor),
                                                       round(exp(confint.default(fit, level = 0.95))[2,2],rounding_factor),
                                                       summary(fit)$coefficients[2,3],
                                                       summary(fit)$coefficients[2,4])
      res.df[i,grep(suffixes[2], colnames(res.df))]<-c(round(exp(coef(fit))[3],rounding_factor),
                                                       round(exp(confint.default(fit, level = 0.95))[3,1],rounding_factor),
                                                       round(exp(confint.default(fit, level = 0.95))[3,2],rounding_factor),
                                                       summary(fit)$coefficients[3,3],
                                                       summary(fit)$coefficients[3,4])
      res.df[i,grep("Interaction",colnames(res.df))]<-c(round(exp(coef(fit))[4+length(covariates_adjustment)],rounding_factor),
                                                        round(exp(confint.default(fit, level = 0.95))[4+length(covariates_adjustment),1],rounding_factor),
                                                        round(exp(confint.default(fit, level = 0.95))[4+length(covariates_adjustment),2],rounding_factor),
                                                        summary(fit)$coefficients[4+length(covariates_adjustment),3],
                                                        summary(fit)$coefficients[4+length(covariates_adjustment),4])
      k=4
      for(j in 1:length(covariates_adjustment)){
        res.df[i,grep(covariates_adjustment[j], colnames(res.df))]<-c(round(exp(coef(fit))[k],rounding_factor),
                                                                      round(exp(confint.default(fit, level = 0.95))[k,1],rounding_factor),
                                                                      round(exp(confint.default(fit, level = 0.95))[k,2],rounding_factor),
                                                                      summary(fit)$coefficients[k,3],
                                                                      summary(fit)$coefficients[k,4])
        k<-k+1
      }
    }
  }
  fdr<-res.df[,grep("^P_",colnames(res.df))]
  colnames(fdr)<-gsub("^P_","FDR_",colnames(fdr))
  fdr<-apply(fdr, 2, p.adjust, "BH")
  res.df<-as.data.frame(cbind(res.df, fdr))
  colOrder<-c("Variable",grep(suffixes[1], colnames(res.df), value=T))
  colOrder<-c(colOrder,grep(suffixes[2], colnames(res.df), value=T))
  colOrder<-c(colOrder,grep("Interaction", colnames(res.df), value=T))

  if(!is.null(covariates_adjustment)){
    for(z in 1:length(covariates_adjustment)){
      colOrder<-c(colOrder, grep(covariates_adjustment[z], colnames(res.df), value=T))
    }
  }
  res.df<-res.df[,colOrder]
  return(res.df)
}

#**************************************************************************************
#* Comparison of gene set enrichment networks between multiple contrasts
#**************************************************************************************


compareEnrichmentNetworks <- function(cluster_by = c("gene_content", "correlation"),
                                      similarity_method = c("kappa", "jaccard", "dice", "overlap"),
                                      similarity_th = 0,
                                      expMat,
                                      enrichment_list,
                                      directional = TRUE,
                                      GSEA_style = TRUE,
                                      gene_set_list = NULL,
                                      correlation_method = c("spearman", "pearson"),
                                      correl_th = 0.9,
                                      FDR_th = 0.05,
                                      FDR_column,
                                      gene_set_name_column,
                                      NES_column = "NES",
                                      NES_th = 0,
                                      direction_column = NULL,
                                      gene_content_column = NULL,
                                      igraph.vertex.label = NA,
                                      igraph.vertex.size.cex = 1,
                                      igraph.vertex.label.cex = 0.2,
                                      igraph.vertex.frame.color = "black",
                                      igraph.edge.color = "black",
                                      igraph.mark.border = "black",
                                      legend_cex = 0.4,
                                      nColLayout = 3) {



  require(Biobase)
  require(circlize)
  require(ComplexHeatmap)
  require(dplyr)
  require(ggplot2)
  require(igraph)
  require(msigdbr)
  require(rcartocolor)
  require(RColorBrewer)
  require(simplifyEnrichment)
  require(singscore)
  require(xlsx)

  enrichList_sig<-lapply(enrichment_list, function(x)x[x[,FDR_column] < FDR_th & abs(x[,NES_column]) >= NES_th,])
  sig<-unique(unlist(lapply(enrichList_sig, function(x)x[,gene_set_name_column])))
  sig<-sig[!is.na(sig)]

  if(length(sig) == 0){
    print("No significant gene sets found at provided thresholds.")
  } else {

    if(cluster_by == "gene_content"){

      ### retrieving gene set content

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


      ig<-graph_from_adjacency_matrix(adjmatrix = ts,
                                      mode="upper",
                                      diag=F,
                                      weighted = T)
      ig<-delete_edges(ig, edges=which(E(ig)$weight <= similarity_th))
      if(gsize(ig) == 0){
        print("No edges left after filtering for similarity threshold.")
      } else {

        cl = membership(cluster_fast_greedy(ig))
        cl[cl%in%names(which(table(cl)==1))]<-"Unclustered"
        names(cl)<-gsub("CUSTOM_","",names(cl))
        V(ig)$Cluster<-cl

        if(length(table(cl))>1){
          for(i in 1:length(enrichList_sig)){
            cl.temp<-data.frame(x = names(cl), Cluster_gene_content = as.vector(cl))
            colnames(cl.temp)[1]<-gene_set_name_column
            enrichment_list[[i]] <- merge(enrichment_list[[i]], cl.temp, by = gene_set_name_column, all.x = T)
          }
          if(directional == FALSE){
            layout(matrix(1:length(enrichList_sig), ncol = nColLayout))
            for(i in 1:length(enrichList_sig)){
              temp<-enrichment_list[[i]]
              temp<-temp[,colnames(temp) %in% c(gene_set_name_column, FDR_column)]
              foo<-data.frame(Node = names(V(ig)), pathway = gsub("CUSTOM_", "", names(V(ig))), stringsAsFactors = F)
              colnames(foo)[2]<-gene_set_name_column
              foo<-merge(foo, temp, by = gene_set_name_column, all.x = T)
              foo<-foo[match(names(V(ig)), foo$Node),]
              V(ig)$FDR<--log10(ifelse(foo[,FDR_column] < 1e-10, 1e-10, foo[,FDR_column]))
              V(ig)$Sig<-ifelse(foo[,FDR_column] < FDR_th, "gold", "grey80")
              markg<-list()
              for(j in 1:length(unique(cl))){
                markg[[j]]<-which(cl%in%unique(cl)[j])
              }
              names(markg)<-unique(cl)
              markg<-markg[names(markg)!="Unclustered"]
              set.seed(333)
              clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = carto_pal(12,"Bold")[c(1:(length(table(V(ig)$Cluster))-1),12)])
              names(clustercol)<-factor(V(ig)$Cluster)
              clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
              clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
              set.seed(123)
              coords<-layout_nicely(ig)
              plot.igraph(ig,
                          layout=coords,
                          vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                          vertex.label=igraph.vertex.label,
                          vertex.label.cex=igraph.vertex.label.cex,
                          vertex.frame.color=igraph.vertex.frame.color,
                          vertex.color=V(ig)$Sig,
                          edge.width = E(ig)$weight,
                          edge.color = igraph.edge.color,
                          mark.groups = markg,
                          mark.shape=0.9,
                          mark.expand=10,
                          mark.border = igraph.mark.border,
                          mark.col = clustercolLegend,
                          margin=c(0,0,0,0))
              title(names(enrichList_sig)[i],cex.main=1.5,col.main="black")
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
            layout(matrix(1:length(enrichList_sig), ncol = nColLayout))
            for(i in 1:length(enrichList_sig)){
              temp<-enrichment_list[[i]]
              temp<-temp[,colnames(temp) %in% c(gene_set_name_column, NES_column, FDR_column)]
              foo<-data.frame(Node = names(V(ig)), pathway = gsub("CUSTOM_", "", names(V(ig))), stringsAsFactors = F)
              colnames(foo)[2]<-gene_set_name_column
              foo<-merge(foo, temp, by = gene_set_name_column, all.x = T)
              foo<-foo[match(names(V(ig)), foo$Node),]
              V(ig)$FDR<--log10(ifelse(foo[,FDR_column] < 1e-10, 1e-10, foo[,FDR_column]))
              V(ig)$Sig<-ifelse(foo[,FDR_column] < FDR_th & foo[,NES_column] >= NES_th, "red2", ifelse(foo[,FDR_column] < FDR_th & foo[,NES_column] <= -NES_th, "royalblue", "grey80"))
              markg<-list()
              for(j in 1:length(unique(cl))){
                markg[[j]]<-which(cl%in%unique(cl)[j])
              }
              names(markg)<-unique(cl)
              markg<-markg[names(markg)!="Unclustered"]
              set.seed(333)
              clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = carto_pal(12,"Bold")[c(1:(length(table(V(ig)$Cluster))-1),12)])
              names(clustercol)<-factor(V(ig)$Cluster)
              clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
              clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
              set.seed(123)
              coords<-layout_nicely(ig)
              plot.igraph(ig,
                          layout=coords,
                          vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                          vertex.label=igraph.vertex.label,
                          vertex.label.cex=igraph.vertex.label.cex,
                          vertex.frame.color=igraph.vertex.frame.color,
                          vertex.color=V(ig)$Sig,
                          edge.width = E(ig)$weight,
                          edge.color = igraph.edge.color,
                          mark.groups = markg,
                          mark.shape=0.9,
                          mark.expand=10,
                          mark.border = igraph.mark.border,
                          mark.col = clustercolLegend,
                          margin=c(0,0,0,0))
              title(names(enrichList_sig)[i],cex.main=1.5,col.main="black")
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

          if(directional == TRUE & GSEA_style == FALSE){
            layout(matrix(1:length(enrichList_sig), ncol = nColLayout))
            for(i in 1:length(enrichList_sig)){
              temp<-enrichment_list[[i]]
              temp<-temp[,colnames(temp) %in% c(gene_set_name_column, direction_column, FDR_column)]
              foo<-data.frame(Node = names(V(ig)), pathway = gsub("CUSTOM_", "", names(V(ig))), stringsAsFactors = F)
              colnames(foo)[2]<-gene_set_name_column
              foo<-merge(foo, temp, by = gene_set_name_column, all.x = T)
              foo<-foo[match(names(V(ig)), foo$Node),]
              foo[,FDR_column][is.na(foo[,FDR_column])]<-0.5
              foo[,direction_column][is.na(foo[,direction_column])]<-"n.s."
              foo$Color<-"grey80"
              foo$Color[grep("up", foo[,direction_column], ignore.case = T)]<-"red2"
              foo$Color[grep("dn|down", foo[,direction_column], ignore.case = T)]<-"royalblue"
              V(ig)$FDR<--log10(ifelse(foo[,FDR_column] < 1e-10, 1e-10, foo[,FDR_column]))
              V(ig)$Sig<-foo$Color
              markg<-list()
              for(j in 1:length(unique(cl))){
                markg[[j]]<-which(cl%in%unique(cl)[j])
              }
              names(markg)<-unique(cl)
              markg<-markg[names(markg)!="Unclustered"]
              set.seed(333)
              clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = carto_pal(12,"Bold")[c(1:(length(table(V(ig)$Cluster))-1),12)])
              names(clustercol)<-factor(V(ig)$Cluster)
              clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
              clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
              set.seed(123)
              coords<-layout_nicely(ig)
              plot.igraph(ig,
                          layout=coords,
                          vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                          vertex.label=igraph.vertex.label,
                          vertex.label.cex=igraph.vertex.label.cex,
                          vertex.frame.color=igraph.vertex.frame.color,
                          vertex.color=V(ig)$Sig,
                          edge.width = E(ig)$weight,
                          edge.color = igraph.edge.color,
                          mark.groups = markg,
                          mark.shape=0.9,
                          mark.expand=10,
                          mark.border = igraph.mark.border,
                          mark.col = clustercolLegend,
                          margin=c(0,0,0,0))
              title(names(enrichList_sig)[i],cex.main=1.5,col.main="black")
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
        } else {
          if(directional == FALSE){
            layout(matrix(1:length(enrichList_sig), ncol = nColLayout))
            for(i in 1:length(enrichList_sig)){
              temp<-enrichment_list[[i]]
              temp<-temp[,colnames(temp) %in% c(gene_set_name_column, FDR_column)]
              foo<-data.frame(Node = names(V(ig)), pathway = gsub("CUSTOM_", "", names(V(ig))), stringsAsFactors = F)
              colnames(foo)[2]<-gene_set_name_column
              foo<-merge(foo, temp, by = gene_set_name_column, all.x = T)
              foo<-foo[match(names(V(ig)), foo$Node),]
              V(ig)$FDR<--log10(ifelse(foo[,FDR_column] < 1e-10, 1e-10, foo[,FDR_column]))
              V(ig)$Sig<-ifelse(foo[,FDR_column] < FDR_th, "gold", "grey80")
              set.seed(333)
              clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = carto_pal(12,"Bold")[c(1:(length(table(V(ig)$Cluster))-1),12)])
              names(clustercol)<-factor(V(ig)$Cluster)
              clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
              clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
              set.seed(123)
              coords<-layout_nicely(ig)
              plot.igraph(ig,
                          layout=coords,
                          vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                          vertex.label=igraph.vertex.label,
                          vertex.label.cex=igraph.vertex.label.cex,
                          vertex.frame.color=igraph.vertex.frame.color,
                          vertex.color=V(ig)$Sig,
                          edge.width = E(ig)$weight,
                          edge.color = igraph.edge.color,
                          margin=c(0,0,0,0))
              title(names(enrichList_sig)[i],cex.main=1.5,col.main="black")
            }
          }

          if(directional == TRUE & GSEA_style == TRUE){
            layout(matrix(1:length(enrichList_sig), ncol = nColLayout))
            for(i in 1:length(enrichList_sig)){
              temp<-enrichment_list[[i]]
              temp<-temp[,colnames(temp) %in% c(gene_set_name_column, NES_column, FDR_column)]
              foo<-data.frame(Node = names(V(ig)), pathway = gsub("CUSTOM_", "", names(V(ig))), stringsAsFactors = F)
              colnames(foo)[2]<-gene_set_name_column
              foo<-merge(foo, temp, by = gene_set_name_column, all.x = T)
              foo<-foo[match(names(V(ig)), foo$Node),]
              V(ig)$FDR<--log10(ifelse(foo[,FDR_column] < 1e-10, 1e-10, foo[,FDR_column]))
              V(ig)$Sig<-ifelse(foo[,FDR_column] < FDR_th & foo[,NES_column] >= NES_th, "red2", ifelse(foo[,FDR_column] < FDR_th & foo[,NES_column] <= -NES_th, "royalblue", "grey80"))
              set.seed(333)
              clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = carto_pal(12,"Bold")[c(1:(length(table(V(ig)$Cluster))-1),12)])
              names(clustercol)<-factor(V(ig)$Cluster)
              clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
              clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
              set.seed(123)
              coords<-layout_nicely(ig)
              plot.igraph(ig,
                          layout=coords,
                          vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                          vertex.label=igraph.vertex.label,
                          vertex.label.cex=igraph.vertex.label.cex,
                          vertex.frame.color=igraph.vertex.frame.color,
                          vertex.color=V(ig)$Sig,
                          edge.width = E(ig)$weight,
                          edge.color = igraph.edge.color,
                          margin=c(0,0,0,0))
              title(names(enrichList_sig)[i],cex.main=1.5,col.main="black")
            }
          }

          if(directional == TRUE & GSEA_style == FALSE){
            layout(matrix(1:length(enrichList_sig), ncol = nColLayout))
            for(i in 1:length(enrichList_sig)){
              temp<-enrichment_list[[i]]
              temp<-temp[,colnames(temp) %in% c(gene_set_name_column, direction_column, FDR_column)]
              foo<-data.frame(Node = names(V(ig)), pathway = gsub("CUSTOM_", "", names(V(ig))), stringsAsFactors = F)
              colnames(foo)[2]<-gene_set_name_column
              foo<-merge(foo, temp, by = gene_set_name_column, all.x = T)
              foo<-foo[match(names(V(ig)), foo$Node),]
              foo[,FDR_column][is.na(foo[,FDR_column])]<-0.5
              foo[,direction_column][is.na(foo[,direction_column])]<-"n.s."
              foo$Color<-"grey80"
              foo$Color[grep("up", foo[,direction_column], ignore.case = T)]<-"red2"
              foo$Color[grep("dn|down", foo[,direction_column], ignore.case = T)]<-"royalblue"
              V(ig)$FDR<--log10(ifelse(foo[,FDR_column] < 1e-10, 1e-10, foo[,FDR_column]))
              V(ig)$Sig<-foo$Color
              set.seed(333)
              clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = carto_pal(12,"Bold")[c(1:(length(table(V(ig)$Cluster))-1),12)])
              names(clustercol)<-factor(V(ig)$Cluster)
              clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
              clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
              set.seed(123)
              coords<-layout_nicely(ig)
              plot.igraph(ig,
                          layout=coords,
                          vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                          vertex.label=igraph.vertex.label,
                          vertex.label.cex=igraph.vertex.label.cex,
                          vertex.frame.color=igraph.vertex.frame.color,
                          vertex.color=V(ig)$Sig,
                          edge.width = E(ig)$weight,
                          edge.color = igraph.edge.color,
                          margin=c(0,0,0,0))
              title(names(enrichList_sig)[i],cex.main=1.5,col.main="black")
            }
          }
        }
      }
    }


    #################################ù


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
                                    knownDir = T,
                                    centerScore = T)
          scoreMat<-rbind(scoreMat1, scoreMat2)
        } else {
          scoreMat<-scoreMat1
        }
      } else {
        scoreMat<-simpleScoreMod(rankData = rankMat,
                                 mysetlist = setlist,
                                 knownDir = T,
                                 centerScore = T)
      }
      print("Done.")
      rownames(scoreMat)<-gsub("CUSTOM_", "", rownames(scoreMat))

      correl<-cor(t(scoreMat), method=correlation_method)


      ig<-graph_from_adjacency_matrix(adjmatrix = correl,
                                      mode="upper",
                                      diag=F,
                                      weighted = T)
      ig<-delete_edges(ig, edges=which(E(ig)$weight <= correl_th))
      if(gsize(ig) == 0){
        print("No edges left after filtering for similarity threshold.")
      } else {
        cl = membership(cluster_fast_greedy(ig))
        cl[cl%in%names(which(table(cl)==1))]<-"Unclustered"
        #names(cl)<-gsub("CUSTOM_","",names(cl))
        V(ig)$Cluster<-cl

        if(length(table(cl))>1){
          for(i in 1:length(enrichList_sig)){
            cl.temp<-data.frame(x = names(cl), Cluster_correlation = as.vector(cl))
            colnames(cl.temp)[1]<-gene_set_name_column
            enrichment_list[[i]] <- merge(enrichment_list[[i]], cl.temp, by = gene_set_name_column, all.x = T)
          }
          if(directional == FALSE){
            layout(matrix(1:length(enrichList_sig), ncol = nColLayout))
            for(i in 1:length(enrichList_sig)){
              temp<-enrichment_list[[i]]
              temp<-temp[,colnames(temp) %in% c(gene_set_name_column, FDR_column)]
              foo<-data.frame(Node = names(V(ig)), pathway = gsub("CUSTOM_", "", names(V(ig))), stringsAsFactors = F)
              colnames(foo)[2]<-gene_set_name_column
              foo<-merge(foo, temp, by = gene_set_name_column, all.x = T)
              foo<-foo[match(names(V(ig)), foo$Node),]
              V(ig)$FDR<--log10(ifelse(foo[,FDR_column] < 1e-10, 1e-10, foo[,FDR_column]))
              V(ig)$Sig<-ifelse(foo[,FDR_column] < FDR_th, "gold", "grey80")
              markg<-list()
              for(j in 1:length(unique(cl))){
                markg[[j]]<-which(cl%in%unique(cl)[j])
              }
              names(markg)<-unique(cl)
              markg<-markg[names(markg)!="Unclustered"]
              set.seed(333)
              clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = carto_pal(12,"Bold")[c(1:(length(table(V(ig)$Cluster))-1),12)])
              names(clustercol)<-factor(V(ig)$Cluster)
              clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
              clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
              set.seed(123)
              coords<-layout_nicely(ig)
              plot.igraph(ig,
                          layout=coords,
                          vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                          vertex.label=igraph.vertex.label,
                          vertex.label.cex=igraph.vertex.label.cex,
                          vertex.frame.color=igraph.vertex.frame.color,
                          vertex.color=V(ig)$Sig,
                          edge.width = E(ig)$weight,
                          edge.color = igraph.edge.color,
                          mark.groups = markg,
                          mark.shape=0.9,
                          mark.expand=10,
                          mark.border = igraph.mark.border,
                          mark.col = clustercolLegend,
                          margin=c(0,0,0,0))
              title(names(enrichList_sig)[i],cex.main=1.5,col.main="black")
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
            layout(matrix(1:length(enrichList_sig), ncol = nColLayout))
            for(i in 1:length(enrichList_sig)){
              temp<-enrichment_list[[i]]
              temp<-temp[,colnames(temp) %in% c(gene_set_name_column, NES_column, FDR_column)]
              foo<-data.frame(Node = names(V(ig)), pathway = gsub("CUSTOM_", "", names(V(ig))), stringsAsFactors = F)
              colnames(foo)[2]<-gene_set_name_column
              foo<-merge(foo, temp, by = gene_set_name_column, all.x = T)
              foo<-foo[match(names(V(ig)), foo$Node),]
              V(ig)$FDR<--log10(ifelse(foo[,FDR_column] < 1e-10, 1e-10, foo[,FDR_column]))
              V(ig)$Sig<-ifelse(foo[,FDR_column] < FDR_th & foo[,NES_column] >= NES_th, "red2", ifelse(foo[,FDR_column] < FDR_th & foo[,NES_column] <= -NES_th, "royalblue", "grey80"))
              markg<-list()
              for(j in 1:length(unique(cl))){
                markg[[j]]<-which(cl%in%unique(cl)[j])
              }
              names(markg)<-unique(cl)
              markg<-markg[names(markg)!="Unclustered"]
              set.seed(333)
              clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = carto_pal(12,"Bold")[c(1:(length(table(V(ig)$Cluster))-1),12)])
              names(clustercol)<-factor(V(ig)$Cluster)
              clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
              clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
              set.seed(123)
              coords<-layout_nicely(ig)
              plot.igraph(ig,
                          layout=coords,
                          vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                          vertex.label=igraph.vertex.label,
                          vertex.label.cex=igraph.vertex.label.cex,
                          vertex.frame.color=igraph.vertex.frame.color,
                          vertex.color=V(ig)$Sig,
                          edge.width = E(ig)$weight,
                          edge.color = igraph.edge.color,
                          mark.groups = markg,
                          mark.shape=0.9,
                          mark.expand=10,
                          mark.border = igraph.mark.border,
                          mark.col = clustercolLegend,
                          margin=c(0,0,0,0))
              title(names(enrichList_sig)[i],cex.main=1.5,col.main="black")
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

          if(directional == TRUE & GSEA_style == FALSE){
            layout(matrix(1:length(enrichList_sig), ncol = nColLayout))
            for(i in 1:length(enrichList_sig)){
              temp<-enrichment_list[[i]]
              temp<-temp[,colnames(temp) %in% c(gene_set_name_column, direction_column, FDR_column)]
              foo<-data.frame(Node = names(V(ig)), pathway = gsub("CUSTOM_", "", names(V(ig))), stringsAsFactors = F)
              colnames(foo)[2]<-gene_set_name_column
              foo<-merge(foo, temp, by = gene_set_name_column, all.x = T)
              foo<-foo[match(names(V(ig)), foo$Node),]
              foo[,FDR_column][is.na(foo[,FDR_column])]<-0.5
              foo[,direction_column][is.na(foo[,direction_column])]<-"n.s."
              foo$Color<-"grey80"
              foo$Color[grep("up", foo[,direction_column], ignore.case = T)]<-"red2"
              foo$Color[grep("dn|down", foo[,direction_column], ignore.case = T)]<-"royalblue"
              V(ig)$FDR<--log10(ifelse(foo[,FDR_column] < 1e-10, 1e-10, foo[,FDR_column]))
              V(ig)$Sig<-foo$Color
              markg<-list()
              for(j in 1:length(unique(cl))){
                markg[[j]]<-which(cl%in%unique(cl)[j])
              }
              names(markg)<-unique(cl)
              markg<-markg[names(markg)!="Unclustered"]
              set.seed(333)
              clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = carto_pal(12,"Bold")[c(1:(length(table(V(ig)$Cluster))-1),12)])
              names(clustercol)<-factor(V(ig)$Cluster)
              clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
              clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
              set.seed(123)
              coords<-layout_nicely(ig)
              plot.igraph(ig,
                          layout=coords,
                          vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                          vertex.label=igraph.vertex.label,
                          vertex.label.cex=igraph.vertex.label.cex,
                          vertex.frame.color=igraph.vertex.frame.color,
                          vertex.color=V(ig)$Sig,
                          edge.width = E(ig)$weight,
                          edge.color = igraph.edge.color,
                          mark.groups = markg,
                          mark.shape=0.9,
                          mark.expand=10,
                          mark.border = igraph.mark.border,
                          mark.col = clustercolLegend,
                          margin=c(0,0,0,0))
              title(names(enrichList_sig)[i],cex.main=1.5,col.main="black")
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
        } else {
          if(directional == FALSE){
            layout(matrix(1:length(enrichList_sig), ncol = nColLayout))
            for(i in 1:length(enrichList_sig)){
              temp<-enrichment_list[[i]]
              temp<-temp[,colnames(temp) %in% c(gene_set_name_column, FDR_column)]
              foo<-data.frame(Node = names(V(ig)), pathway = gsub("CUSTOM_", "", names(V(ig))), stringsAsFactors = F)
              colnames(foo)[2]<-gene_set_name_column
              foo<-merge(foo, temp, by = gene_set_name_column, all.x = T)
              foo<-foo[match(names(V(ig)), foo$Node),]
              V(ig)$FDR<--log10(ifelse(foo[,FDR_column] < 1e-10, 1e-10, foo[,FDR_column]))
              V(ig)$Sig<-ifelse(foo[,FDR_column] < FDR_th, "gold", "grey80")
              set.seed(333)
              clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = carto_pal(12,"Bold")[c(1:(length(table(V(ig)$Cluster))-1),12)])
              names(clustercol)<-factor(V(ig)$Cluster)
              clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
              clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
              set.seed(123)
              coords<-layout_nicely(ig)
              plot.igraph(ig,
                          layout=coords,
                          vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                          vertex.label=igraph.vertex.label,
                          vertex.label.cex=igraph.vertex.label.cex,
                          vertex.frame.color=igraph.vertex.frame.color,
                          vertex.color=V(ig)$Sig,
                          edge.width = E(ig)$weight,
                          edge.color = igraph.edge.color,
                          margin=c(0,0,0,0))
              title(names(enrichList_sig)[i],cex.main=1.5,col.main="black")
            }
          }

          if(directional == TRUE & GSEA_style == TRUE){
            layout(matrix(1:length(enrichList_sig), ncol = nColLayout))
            for(i in 1:length(enrichList_sig)){
              temp<-enrichment_list[[i]]
              temp<-temp[,colnames(temp) %in% c(gene_set_name_column, NES_column, FDR_column)]
              foo<-data.frame(Node = names(V(ig)), pathway = gsub("CUSTOM_", "", names(V(ig))), stringsAsFactors = F)
              colnames(foo)[2]<-gene_set_name_column
              foo<-merge(foo, temp, by = gene_set_name_column, all.x = T)
              foo<-foo[match(names(V(ig)), foo$Node),]
              V(ig)$FDR<--log10(ifelse(foo[,FDR_column] < 1e-10, 1e-10, foo[,FDR_column]))
              V(ig)$Sig<-ifelse(foo[,FDR_column] < FDR_th & foo[,NES_column] >= NES_th, "red2", ifelse(foo[,FDR_column] < FDR_th & foo[,NES_column] <= -NES_th, "royalblue", "grey80"))
              set.seed(333)
              clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = carto_pal(12,"Bold")[c(1:(length(table(V(ig)$Cluster))-1),12)])
              names(clustercol)<-factor(V(ig)$Cluster)
              clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
              clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
              set.seed(123)
              coords<-layout_nicely(ig)
              plot.igraph(ig,
                          layout=coords,
                          vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                          vertex.label=igraph.vertex.label,
                          vertex.label.cex=igraph.vertex.label.cex,
                          vertex.frame.color=igraph.vertex.frame.color,
                          vertex.color=V(ig)$Sig,
                          edge.width = E(ig)$weight,
                          edge.color = igraph.edge.color,
                          margin=c(0,0,0,0))
              title(names(enrichList_sig)[i],cex.main=1.5,col.main="black")
            }
          }

          if(directional == TRUE & GSEA_style == FALSE){
            layout(matrix(1:length(enrichList_sig), ncol = nColLayout))
            for(i in 1:length(enrichList_sig)){
              temp<-enrichment_list[[i]]
              temp<-temp[,colnames(temp) %in% c(gene_set_name_column, direction_column, FDR_column)]
              foo<-data.frame(Node = names(V(ig)), pathway = gsub("CUSTOM_", "", names(V(ig))), stringsAsFactors = F)
              colnames(foo)[2]<-gene_set_name_column
              foo<-merge(foo, temp, by = gene_set_name_column, all.x = T)
              foo<-foo[match(names(V(ig)), foo$Node),]
              foo[,FDR_column][is.na(foo[,FDR_column])]<-0.5
              foo[,direction_column][is.na(foo[,direction_column])]<-"n.s."
              foo$Color<-"grey80"
              foo$Color[grep("up", foo[,direction_column], ignore.case = T)]<-"red2"
              foo$Color[grep("dn|down", foo[,direction_column], ignore.case = T)]<-"royalblue"
              V(ig)$FDR<--log10(ifelse(foo[,FDR_column] < 1e-10, 1e-10, foo[,FDR_column]))
              V(ig)$Sig<-foo$Color
              set.seed(333)
              clustercol<-level.colors(x=as.numeric(factor(V(ig)$Cluster)), at=c(0:length(unique(V(ig)$Cluster))), col.regions = carto_pal(12,"Bold")[c(1:(length(table(V(ig)$Cluster))-1),12)])
              names(clustercol)<-factor(V(ig)$Cluster)
              clustercolLegend<-structure(unique(clustercol), names = unique(names(clustercol)))
              clustercolLegend<-clustercolLegend[order(names(clustercolLegend))]
              set.seed(123)
              coords<-layout_nicely(ig)
              plot.igraph(ig,
                          layout=coords,
                          vertex.size=V(ig)$FDR*igraph.vertex.size.cex,
                          vertex.label=igraph.vertex.label,
                          vertex.label.cex=igraph.vertex.label.cex,
                          vertex.frame.color=igraph.vertex.frame.color,
                          vertex.color=V(ig)$Sig,
                          edge.width = E(ig)$weight,
                          edge.color = igraph.edge.color,
                          margin=c(0,0,0,0))
              title(names(enrichList_sig)[i],cex.main=1.5,col.main="black")
            }
          }
        }
      }
    }
    return(enrichment_list)
  }
}



#********************************************************
#gsea interaction
#********************************************************
gseaInteraction<-function(data, phenotype, subsetting_variable, nPerm=1000, pheno_factor_levels, pathways, seed=123, simulate_ranking=TRUE, test=c("t-test","lm","logistic"), ranks=NULL){
  if(simulate_ranking){
    require(Rfast)
    require(Rfast2)
    nClasses<-names(table(pData(data)[,subsetting_variable]))
    ESlist<-vector("list",length=length(nClasses))
    names(ESlist)<-nClasses
    rankingList<-ESlist
    for(i in names(ESlist)){
      subData<-data[,pData(data)[,subsetting_variable]==i]
      Y<-as.numeric(factor(pData(subData)[,phenotype],levels=pheno_factor_levels))-1
      dataTransp<-t(exprs(subData))
      set.seed(seed)
      if(test=="t-test"){
        realRanking<-ttests(x=dataTransp,ina=Y)[,1]
        rankingMat<-replicate(nPerm,ttests(x=dataTransp,ina=sample(Y))[,1])
      }
      if(test=="lm"){
        realRanking<-regression(x=dataTransp,y=Y)[,1]
        rankingMat<-replicate(nPerm,regression(x=dataTransp,y=sample(Y))[,1])
      }
      if(test=="logistic"){
        realRanking<-sp.logiregs(x=dataTransp,y=Y)[,1]
        rankingMat<-replicate(nPerm,sp.logiregs(x=dataTransp,y=sample(Y))[,1])
      }
      rankingMat<-as.data.frame(cbind(realRanking, rankingMat))
      colnames(rankingMat)<-c("Real",paste0("Random",1:nPerm))
      rankingMat<-as.list(rankingMat)
      rankingMat<-lapply(rankingMat, function(x)structure(x,names=rownames(subData)))
      rankingMat<-lapply(rankingMat,sort, decreasing=T)
      geneIndexes<-lapply(rankingMat,function(x)lapply(pathways,function(z)which(names(x)%in%z)))
      calcES<-function(geneIndexList, randomRanking){
        unlist(lapply(geneIndexList,function(x)calcGseaStat(stats = randomRanking, selectedStats = x)))
      }
      rankingList[[i]]<-rankingMat
      ESlist[[i]]<-mapply(calcES, geneIndexList=geneIndexes, randomRanking=rankingMat)
    }
  } else {
    ESlist<-vector("list",length=length(ranks))
    names(ESlist)<-names(ranks)
    for(i in names(ESlist)){
      geneIndexes<-lapply(ranks,function(x)lapply(pathways,function(z)which(names(x)%in%z)))
      calcES<-function(geneIndexList, randomRanking){
        unlist(lapply(geneIndexList,function(x)calcGseaStat(stats = randomRanking, selectedStats = x)))
      }
      ESlist[[i]]<-mapply(calcES, geneIndexList=geneIndexes, randomRanking=ranks)
    }
  }
  realDiff<-ESlist[[1]][,1]-ESlist[[2]][,1]
  simDiff<-mapply(function(X,Y){sapply(2:(nPerm+1), function(col) X[,col]-Y[,col])}, X=ESlist[1], Y=ESlist[2], SIMPLIFY=F)[[1]]
  diffs<-cbind(realDiff, simDiff)
  pval<-apply(diffs, 1, function(x){
    if(x[1]>0){
      sum(x[x>0][-1]>=x[1])/sum(x[-1]>0)
    } else {
      sum(x[x<0][-1]<=x[1])/sum(x[-1]<0)
    }
  }
  )
  pvalCalc<-list(ranks=rankingList, ES=ESlist, differences = diffs, pvalues=pval)
  return(pvalCalc)
}


#**********************************************
# random limma
#**********************************************

randomLimma<-function(phenotype, covariates, dgeObject, contrasts){
  f<-factor(sample(phenotype))
  covs<-dgeObject$samples[,covariates]
  covs<-cbind(f,covs)
  if(sum(colSums(is.na(covs)))>0){
    stop("Error: NAs present in the sample covariates")
  }
  if(!is.null(covariates)){
    formula<-as.formula(paste0("~0+f+",paste(covariates, collapse="+")))
  } else {
    formula<-as.formula("~0+f")
  }
  design<-model.matrix(formula, data=covs)
  colnames(design)[1:length(levels(f))]<-levels(f)
  v<-voom(dgeObject, design, plot=F)
  fit <- lmFit(v, design)
  contrast.matrix <- makeContrasts(contrasts=contrasts,
                                   levels=design)
  fit.cont <- contrasts.fit(fit, contrast.matrix)
  fit.eb <- eBayes(fit.cont)
  tList<-vector("list", length=length(contrasts))
  names(tList)<-contrasts
  tStats<-lapply(contrasts, function(x)structure(topTable(fit.eb,coef=x,number=nrow(v), adjust="BH", sort.by="none")$t, names=rownames(dgeObject)))
  names(tStats)<-contrasts
  return(tStats)
}


#********************************************************
# gsea plot
#********************************************************

gseaPlot<-function(rankStat, geneSet, geneSetList, gseaTable, margin.r=0.2, margin.l=0.2, title="", labelPos="", labelNeg="", labelSize=4){
  library(fgsea)
  library(ggplot2)
  library(cowplot)
  statsAdj <- sort(rankStat,decreasing=T)
  mypathway <- unname(as.vector(na.omit(match(geneSetList[[geneSet]], names(statsAdj)))))
  mypathway <- sort(mypathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = mypathway,
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(mypathway - 1, mypathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  sigLabel<-paste0("ES = ",round(gseaTable$ES[gseaTable$pathway==geneSet],2),"\nNES = ",round(gseaTable$NES[gseaTable$pathway==geneSet],2),"\nP = ",round(gseaTable$pval[gseaTable$pathway==geneSet],3),"\nFDR = ",round(gseaTable$padj[gseaTable$pathway==geneSet],3))

  if(gseaTable$ES[gseaTable$pathway==geneSet]<0){
    g1<-ggplot(data=toPlot, aes(x=x, y=y)) +
      geom_line(size=1.2, color="green") +
      xlab(NULL) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            plot.margin = ggplot2::margin(t = 0.2, r = margin.r, b = -0.2, l = margin.l, unit = "cm")) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(breaks=seq(from=min(toPlot$y),to= max(toPlot$y), length.out=5),
                         labels = round(seq(from=min(toPlot$y),to= max(toPlot$y), length.out=5),2)) +
      ylab("Running ES") +
      labs(title=paste0("Geneset: ", geneSet, title)) +
      ggplot2::annotate("text", x= min(toPlot$x)+500, y =min(toPlot$y),label=sigLabel,hjust=0, vjust=0)
  } else {
    g1<-ggplot(data=toPlot, aes(x=x, y=y)) +
      geom_line(size=1.2, color="green") +
      xlab(NULL) +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            plot.margin = ggplot2::margin(t = 0.2, r = margin.r, b = -0.2, l = margin.l, unit = "cm")) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(breaks=seq(from=min(toPlot$y),to= max(toPlot$y), length.out=5),
                         labels = paste0(" ",round(seq(from=min(toPlot$y),to= max(toPlot$y), length.out=5),2))) +
      ylab("Running ES") +
      labs(title=paste0("Geneset: ", geneSet, title)) +
      ggplot2::annotate("text", x= max(toPlot$x)-500, y =max(toPlot$y),label=sigLabel,hjust=1, vjust=1)
  }

  ticks<-data.frame(idx=1:length(statsAdj),ymin=0,ymax=0,ymax2=1)
  ticks$ymax[mypathway]<-1
  g2<-ggplot(data=ticks, aes(x=idx)) +
    geom_linerange(aes(ymin=ymin, ymax=ymax)) +
    #geom_linerange(aes(ymin=ymin, ymax=ymax2,color=statsAdj)) +
    #scale_color_gradient2(low=tickcol[1],mid=tickcol[2],high=tickcol[3], guide="none") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line.x = element_blank(),
          plot.margin = ggplot2::margin(t = 0, r = margin.r, b = 0, l = margin.l, unit = "cm")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0))

  cut1<-cut(statsAdj[statsAdj>0],10, labels=1:10)
  cut2<-cut(statsAdj[statsAdj<0],10, labels=-10:-1)
  cutpoints1<-cumsum(rev(table(cut1)))
  cutpoints2<-cumsum(rev(table(cut2)))+length(statsAdj[statsAdj>0])
  cutpoints<-c(cutpoints1, cutpoints2)
  rects<-data.frame(xmin=c(0,cutpoints[-length(cutpoints)]), xmax=cutpoints, ymin=0, ymax=1, FillVal=seq(10,-9,-1))
  midpoint<-rects[rects$FillVal==0,]
  midpoint$xmax<-midpoint$xmin
  rects2<-rects[rects$FillVal<=0,]
  rects2$FillVal<-rects2$FillVal-1
  rects<-rbind(rects[rects$FillVal>0,],midpoint,rects2)
  g3<-ggplot(aes(x=xmin), data=rects) +
    geom_rect(aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=FillVal), data=rects) +
    scale_fill_gradient2(low="blue",high="red2", guide="none") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line.x = element_blank(),
          plot.margin = ggplot2::margin(t = -0.2, r = margin.r, b = 0, l = margin.l, unit = "cm")) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0)) +
    geom_text(aes(x=500, y=0.5, label=labelPos), size=labelSize, hjust=0, vjust=0.5) +
    geom_text(aes(x=nrow(ticks)-500, y=0.5, label=labelNeg), size=labelSize, hjust=1, vjust=0.5)


  rnkPlot<-data.frame(idx=1:length(statsAdj), rnkStat=statsAdj)
  g4<-ggplot(data=rnkPlot, aes(x=idx, y=rnkStat)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          legend.position = c(0.85, 0.85),
          plot.margin = ggplot2::margin(t = -0.2, r = margin.r, b = 0.2, l = margin.l, unit = "cm")) +
    xlab("Gene rank") +
    ylab("Ranking statitisc") +
    scale_fill_gradient2(low="blue",mid="grey80",high="red2", guide="none") +
    scale_x_continuous(breaks = seq(0,length(statsAdj),5000),expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0))

  gs1<-plot_grid(g1, g2, g3, g4, ncol=1, rel_heights = c(1.5, 0.3, 0.2, 1), align = "v")
  return(gs1)
}


#********************************************************
# NCBI annotation
#********************************************************

ncbiAnnot<-function(fdata, column_name_symbol){
  ncbi<-read.table("/Users/dugo.matteo/Desktop/hsr_bioinfo/databases/annotations/NCBI/Homo_sapiens_gene_info_Jan23", sep="\t", as.is=T, quote="", comment.char = "", header=T)
  ncbi$Synonyms<-paste0(ncbi$Symbol,"|",ncbi$Synonyms)
  geneList<-sapply(ncbi$Synonyms, function(x)strsplit(x,"\\|"))
  names(geneList)<-ncbi$Symbol
  geneList<-lapply(geneList, function(x)data.frame(Official_Symbol=x[1],Alias=x))
  geneDF<-do.call(rbind.data.frame, geneList)
  ncbi2<-unique(ncbi[,c("Symbol","Synonyms")])
  colnames(ncbi2)[1]<-"Official_Symbol"
  geneDF<-merge(geneDF, ncbi2, by="Official_Symbol", all.x=T)
  colnames(fdata)[colnames(fdata)==column_name_symbol]<-"Alias"
  fdata<-merge(fdata, geneDF, by="Alias", all.x=T)
  rownames(fdata)<-fdata$Official_Symbol

}
