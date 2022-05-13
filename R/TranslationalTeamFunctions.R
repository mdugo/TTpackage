
require(Biobase)
require(circlize)
require(ComplexHeatmap)
require(igraph)
require(lattice)
require(msigdbr)
require(randomcoloR)
require(rcartocolor)
require(RColorBrewer)
require(singscore)
require(xlsx)
require(dplyr)
require(ggplot2)
require(ggsankey)

#*************************************************************
#* This function performs singscore calculation for gene set
#* collections that include gene sets with only up-regulated
#* genes or both up- and down-regulated genes
#*************************************************************

simpleScoreMod<-function(rankData, mysetlist, knownDir = TRUE, centerScore = TRUE){
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

correlateGeneSets <- function(expMat, enrichment_table, directional=TRUE, GSEA_style=TRUE, gene_set_list=NULL,
                              correl_th=0.9, FDR_th=0.05, FDR_column,
                              gene_set_name_column, NES_column="NES", NES_th=0,
                              igraph.vertex.label = NA, igraph.vertex.label.cex = 0.2,
                              igraph.vertex.frame.color="black", igraph.edge.color = "black",
                              igraph.mark.border = "black", legend_cex=0.4){

  sig<-enrichment_table[enrichment_table[,FDR_column] < FDR_th & abs(enrichment_table[,NES_column]) >= NES_th, gene_set_name_column]

  if(all(colSums(apply(expMat,2,function(x)sapply(x,is.integer))) == nrow(expMat))){
    stop("Error: all data in the expression matrix are integers. Please provide counts normalized for gene length (FPKM/RPKM or TPM)")
  }

  detectSymbols<-grep("^ACTB$|^GAPDH$|^EGFR$", rownames(expMat), value=T)
  if(length(detectSymbols) == 0){
    stop("Error: row names of expression matrix must be gene symbols.")
  }


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
    names(gs_custom)<-paste0("CUSTOM_", names(gs_custom))
    setlist<-list(setlist, gs_custom)
  }
  if(length(grep("HALLMARK|REACTOME|KEGG|BIOCARTA|PID|WP|GOBP", unique(gsub("_.*", "", sig)))) == 0){
    if(is.null(gene_set_list)){
      stop("Error: No gene set list provided. Please provide the custom gene set list with gene symbols.")
    }
    setlist<-gene_set_list
    names(setlist)<-paste0("CUSTOM_", names(setlist))
  }
  print("Done.")

  print("Performing singscore...")

  rankMat<-rankGenes(expMat)

  if(length(grep("HALLMARK|CUSTOM",names(setlist)))>0){
    scoreMat1<-mySimpleScore(rankData = rankMat,
                             mysetlist = setlist[grep("HALLMARK|CUSTOM",names(setlist))],
                             knownDir = T)
    if(length(grep("REACTOME|KEGG|BIOCARTA|PID|WP|GOBP",names(setlist)))>0){
      scoreMat2<-mySimpleScore(rankData = rankMat,
                               mysetlist = setlist[grep("REACTOME|KEGG|BIOCARTA|PID|WP|GOBP",names(setlist))],
                               knownDir = F,
                               centerScore = F)
      scoreMat<-rbind(scoreMat1, scoreMat2)
    } else {
      scoreMat<-scoreMat1
    }
  } else {
    scoreMat<-mySimpleScore(rankData = rankMat,
                            mysetlist = setlist,
                            knownDir = F,
                            centerScore = F)
  }
  print("Done.")

  correl<-cor(t(scoreMat), method="spearman")

  ig<-graph_from_adjacency_matrix(adjmatrix = correl,
                                  mode="upper",
                                  diag=F,
                                  weighted = T)
  ig<-delete_edges(ig, edges=which(E(ig)$weight <= correl_th))

  cl = membership(cluster_fast_greedy(ig))
  cl[cl%in%names(which(table(cl)==1))]<-"Unclustered"
  res2<-res[res[,gene_set_name_column]%in%names(cl),]
  res2<-res2[match(names(cl), res2[,gene_set_name_column]),]

  V(ig)$Cluster<-cl
  V(ig)$FDR<--log10(ifelse(res2[,FDR_column]<1e-10, 1e-10, res2[,FDR_column]))

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
  if(directional == TRUE & GSEA_style == TRUE){
    V(ig)$Direction<-res2$NES
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
                mark.col = unique(clustercol),
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
                  mark.col = unique(clustercol),
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
  res<-merge(res, community, by=gene_set_name_column, all.x=T)
  res<-res[order(res$Cluster, res[,FDR_column]),]
  return(res)
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
    factor.levels<-levels(factor(df$node))
  }
  df$node<-factor(df$node,levels=factor.levels)
  df$next_node<-factor(df$next_node,levels=factor.levels)

  stages<-levels(df$x)
  count.stages<-matrix(0,nrow=length(stages),ncol=length(factor.levels))
  colnames(count.stages)<-factor.levels
  rownames(count.stages)<-stages
  for(i in 1:nrow(count.stages)){
    count.stages[i,]<-table(df$node[df$x%in%stages[i]])
  }

  lab.pos.stages<-t(apply(count.stages,1,cumsum))
  lab.pos.stages[,1]<-lab.pos.stages[,1]/2
  for(i in 2:ncol(lab.pos.stages)){
    lab.pos.stages[,i]<-(lab.pos.stages[,i]-(count.stages[,i]/2))+(node_space*(i-1))
  }
  if(sum(count.stages<=5)>0){
    if(ncol(count.stages)<3){
      p.test<-fisher.test(count.stages)$p.value
    } else {
      p.test<-fisher.test(count.stages, simulate.p.value=TRUE)$p.value
    }
  } else {
    p.test<-chisq.test(count.stages)$p.value
  }

  if(p.test<0.001){
    p.test<-format(p.test, digits=3, scientific=T)
  } else {
    p.test<-round(p.test,3)
  }



  p<-ggplot(data=df, aes(x = x,
                         next_x = next_x,
                         node = node,
                         next_node = next_node,
                         fill=node)) +
    geom_sankey(smooth=flow_smooth,type="alluvial",node.color = node_color, alpha=node_alpha,space=node_space,width=node_width, show.legend = F) +
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

  if(print_numbers){
    for(i in 1:nrow(lab.pos.stages)){
      for(j in 1:ncol(lab.pos.stages)){
        p <- p + geom_label(x = i, y = lab.pos.stages[i,j], label = paste0(colnames(lab.pos.stages)[j],"\nN = ",count.stages[i,j]), inherit.aes = F, fill = node_label_fill, color = node_label_color, size = node_label_size)
      }
    }
  } else {
    for(i in 1:nrow(lab.pos.stages)){
      for(j in 1:ncol(lab.pos.stages)){
        p <- p + geom_label(x = 1, y = lab.pos.stages[i,j], label = colnames(lab.pos.stages)[j], inherit.aes = F, fill = node_label_fill, color = node_label_color, size = node_label_size)
      }
    }
  }

  return(p)
  print(p)
}


