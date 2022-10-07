\name{correlateGeneSets}
\alias{correlateGeneSets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{correlateGeneSets}
\description{
Given a list of gene sets this function cluster them using a network approach based on similarity in terms of gene content or correlation.
}
\usage{
correlateGeneSets(cluster_by = c("gene_content","correlation"),
                              similarity_method = c("kappa", "jaccard", "dice", "overlap"),
                              similarity_th = 0,
                              expMat,
                              enrichment_table,
                              directional=TRUE,
                              GSEA_style=TRUE,
                              gene_set_list=NULL,
                              correlation_method = c("pearson", "spearman"),
                              correl_th=0.9,
                              FDR_th=0.05,
                              FDR_column,
                              gene_set_name_column,
                              NES_column="NES",
                              NES_th=0,
                              gene_content_column=NULL,
                              igraph.vertex.label = NA,
                              igraph.vertex.label.cex = 0.2,
                              igraph.vertex.size.cex = 1,
                              igraph.vertex.frame.color="black",
                              igraph.edge.color = "black",
                              igraph.mark.border = "black",
                              legend_cex=0.4)
}
\arguments{
  \item{cluster_by}{Character. Whether to cluster gene sets by correlation or by gene content. Values = c("gene_content", "correlation")}
  \item{similarity_method}{Character. Method used to calculate similarity between gene sets when cluster_by is "gene_content". Values = c("kappa", "jaccard", "dice", "overlap"). Ignored if cluster_by = "correlation"}
  \item{similarity_th}{Numerical value. Similarity threshold to filter out network edges. Ignored if cluster_by = "correlation"}
  \item{expMat}{Gene expression matrix with gene symbols as row names. Ignored if cluster_by = "correlation"}
  \item{enrichment_table}{data frame resulting from a GSEA or ORA analysis}
  \item{directional}{If TRUE, a column indicating if gene sets are positively or negatively enriched is required. Default: TRUE}
  \item{GSEA_style}{If TRUE, a column indicating the NES is required. Default: TRUE}
  \item{gene_set_list}{A list of gene sets with gene symbols required when custom gene sets not included in msigdb were used. Default: NULL}
  \item{correl_method}{Character. Correlation method used to calculate correlation between gene sets. Options: pearson or spearman. Ignored if cluster_by = "gene_content"}
  \item{correl_th}{Numerical value. Correlation threshold used to filter edges in the network. Default: 0.9. Ignored if cluster_by = "gene_content"}
  \item{FDR_th}{Numerical value. FDR threshold used to select significant gene sets. Default: 0.05}
  \item{FDR_column}{string indicating the name of the column of the enrichment_table containing FDR. No default}
  \item{gene_set_name_column}{string indicating the name of the column of the enrichment_table containing gene set names. No default}
  \item{NES_column}{string indicating the name of the column of the enrichment_table containing the NES. Required only when GSEA_style = TRUE. No default}
  \item{NES_th}{Numerical value. NES threshold used to filter significant gene sets. Required only when GSEA_style = TRUE. Default: 0}
  \item{gene_content_column}{string indicating the name of the column of the enrichment_table containing leading edge genes or genes of the gene set contained in the input list for ORA analysis. Default = NULL}
  \item{igraph.}{Graphical parameters for igraph network}
  \item{legend_cex}{A numerical value giving the amount by which legend should be magnified relative to the default. A value of 1 means no magnification. Default: 0.4.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An igraph network and a data.frame similar to the input enrichment table with the "Cluster" column added.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Matteo Dugo
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
No example
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
