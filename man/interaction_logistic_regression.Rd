\name{interaction_logistic_regression}
\alias{interaction_logistic_regression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
interaction_logistic_regression
}
\description{
This function performs interaction test for logistic regression between genes and a covariate of interest.
}
\usage{
interaction_logistic_regression(response,
  covariate,
  data,
  rounding_factor = 3,
  sort = FALSE
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{response}{Response must be a numeric vector where 0 indicates failure and 1 success. It an also be specified as a factor, where the first level denotes failure and all others success}
  \item{covariate}{Vector. Covariate of interest}
  \item{data}{Matrix. Gene expression matrix with genes in rows and samples in columns}
  \item{rounding_factor}{Integer. Rounding factor of the numerical results}
  \item{sort}{Boolean. Sort results by p-value}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
