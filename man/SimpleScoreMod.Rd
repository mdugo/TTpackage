\name{simpleScoreMod}
\alias{simpleScoreMod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{simpleScoreMod}
\description{
This function automatically performs singscore for gene set collections that include both gene sets defined by a single list of genes (e.g.: HALLMARK, REACTOME) or two lists of up- and down-regulated genes.
}
\usage{
simpleScoreMod(rankData,
  mysetlist,
  knownDir = TRUE,
  centerScore = TRUE
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rankData}{A matrix object, ranked gene expression matrix data generated using the rankGenes function from singscore package}
  \item{mysetlist}{A list of gene sets with gene symbols}
  \item{knownDir}{A boolean, determining whether the gene set should be considered to be directional or not. A gene set is directional if the type of genes in it are known i.e. up- or down-regulated. Needs to be specified if the gene set is made of one list only. Set to TRUE by default when the gene set includes two lists of up- and down-regulated genes.}
  \item{centerScore}{A Boolean, specifying whether scores should be centered around 0, default as TRUE. Note: scores never centered if knownDirection = FALSE. Set to TRUE by default when the gene set includes two lists of up- and down-regulated genes.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A matrix of singscores for all samples
}
\references{
Foroutan, M., Bhuva, D. D., Lyu, R., Horan, K., Cursons, J., & Davis, M. J. (2018). Single sample scoring of molecular phenotypes. BMC bioinformatics, 19(1), 1-10.
}
\author{
Matteo Dugo. Adapted from sinscore Bioconductor package.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
No examples
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
