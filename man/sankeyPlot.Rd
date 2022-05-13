\name{sankeyPlot}
\alias{sankeyPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{sankeyPlot}
\description{
This function generate sankey plots. The input is a data frame of paired observations in wide format.
}
\usage{
sankeyPlot(x)
}
sankeyPlot(df,
  node_width=0.2,
  node_space=10,
  node_color="black",
  node_alpha=0.4,
  flow_smooth=8,
  node_label_size=4,
  node_label_color="black",
  node_label_fill="white",
  n_size=4,
  x_label_size=13,
  x_label_color="black",
  legend_name="",
  fill_colors=c("blue","red"),
  legend_text_size=12,
  legend_title_size=14,
  plot_title="",
  plot_subtitle_size=13,
  plot_title_size=18,
  print_numbers=TRUE,
  factor.levels=NULL
  )
\arguments{
  \item{df}{Data frame of paired observations in wide format.}
  \item{factor.levels}{Character. Levels of the variables. Useful to specify the order of appearance of the classes in the sankey plot.}
  \item{...}{Graphical parameters taken from ggplot2 and ggsankey.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A sankey plot
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
