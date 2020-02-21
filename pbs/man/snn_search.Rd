\name{snn_search}
\alias{snn_search}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ snn_search
%%  ~~function to do ... ~~
}
\description{ Base on the snn matrix to find the exemplars for each groups.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
snn_search(cl, s_seurat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cl}{ The cell names of target group.
%%     ~~Describe \code{cl} here~~
}
  \item{s_seurat}{ The seurat object of the target group.
%%     ~~Describe \code{s_seurat} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ The exemplars cell names of the target group.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Gerrard Mak
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cl, s_seurat)
{
    snn_table = as.matrix(s_seurat@graphs$RNA_snn[cl, cl])
    colnames(snn_table) = cl
    snn_table[snn_table > 0] = 1
    center_name = names(colSums(snn_table))[which(colSums(snn_table) ==
        max(colSums(snn_table)))]
    if (length(center_name) > 1) {
        snn_colsum = colSums(as.matrix(s_seurat@graphs$RNA_snn[cl,
            center_name]))
        center_name = names(snn_colsum)[which(snn_colsum == max(snn_colsum))][1]
    }
    center_name
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
