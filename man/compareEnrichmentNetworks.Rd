\name{compareEnrichmentNetworks}
\alias{compareEnrichmentNetworks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{compareEnrichmentNetworks}
\description{Comparison of gene set enrichment networks between multiple contrasts}
\usage{
compareEnrichmentNetworks(x)
}
\arguments{
  \item{enrichment_list}{list of data frames resulting from a GSEA or ORA analysis. For other parameters refer to correlateGeneSets function}
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
