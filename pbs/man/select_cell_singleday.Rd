\name{select_cell_singleday}
\alias{select_cell_singleday}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ select_cell_singleday
%%  ~~function to do ... ~~
}
\description{ This function selects cells from large single cell dataset functionally.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
select_cell_singleday(count_data = NULL, count_logdata = NULL, pca_k = NULL, cell_anno = NULL, pca_method = "varianceplot", cluster_method = "kmeans", center_method = "snn_search", k1 = NULL, clustern = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{count_data}{ Raw count matrix. Row is gene and column is cell.
%%     ~~Describe \code{count_data} here~~
}
  \item{count_logdata}{ log count matrix. Row is gene and column is cell.
%%     ~~Describe \code{count_logdata} here~~
}
  \item{pca_k}{ Defaul is null. And pca_k will be calculated by inner function ("varianceplot" or "jackstraw"). pca_k could be set any number greater than 1.
%%     ~~Describe \code{pca_k} here~~
}
  \item{cell_anno}{ Cell annotation matrix. the matrix should have the column name "Day".
%%     ~~Describe \code{cell_anno} here~~
}
  \item{pca_method}{ Defaul is "varianceplot". The other choice is "jackstraw" (not suitable for large dataset). If parameter pca_k is null.
%%     ~~Describe \code{pca_method} here~~
}
  \item{cluster_method}{ Defaul is "kmeans". The other choice are "apcluster", "hierarchy", "pam", "lpc".
%%     ~~Describe \code{cluster_method} here~~
}
  \item{center_method}{ Defaul is "snn_search". The other choice is "nearpoint_search", only for cluster_method is "kmeans".
%%     ~~Describe \code{center_method} here~~
}
  \item{k1}{ Number of groups in second clustering. Defaul is null and it will be calculated by inner function. k1 could be set any number greater than 1.
y
%%     ~~Describe \code{k1} here~~
}
  \item{clustern}{ Number of cluster for parallel computing. Defaul is 4.
%%     ~~Describe \code{clustern} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{ A cell_anno matrix with extra colunms.
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
function (count_data = NULL, count_logdata = NULL, pca_k = NULL,
    cell_anno = NULL, pca_method = "varianceplot", cluster_method = "kmeans",
    center_method = "snn_search", k1 = NULL, clustern = 4)
{
    if (is.null(count_data) & is.null(count_logdata)) {
        stop("empty input")
    }
    if (is.null(cell_anno)) {
        stop("empty cell_anno")
    }
    if (length(unique(cell_anno$Day)) > 1) {
        stop("use other function names select_cell_multiday")
    }
    i <- unique(cell_anno$Day)
    if (!is.null(count_data)) {
        sm <- count_data
        s_seurat <- CreateSeuratObject(counts = sm, project = i)
        s_seurat <- NormalizeData(s_seurat)
    }
    else if (is.null(count_data) & !is.null(count_logdata)) {
        sm_log <- count_logdata
        s_seurat <- CreateSeuratObject(counts = sm_log, project = i)
    }
    s_seurat <- FindVariableFeatures(s_seurat, selection.method = "vst")
    s_seurat <- ScaleData(s_seurat)
    s_seurat <- RunPCA(s_seurat, features = VariableFeatures(object = s_seurat))
    print("finish pca")
    if (is.null(pca_k) & pca_method == "varianceplot") {
        ep <- ElbowPlot(s_seurat, ndims = 50)
        mdim <- sum(abs(diff(ep$data$stdev)) > 0.01) + 1
        if (mdim < 10) {
            mdim <- 10
        }
    }
    else if (is.null(pca_k) & pca_method == "jackstraw") {
        s_seurat <- JackStraw(s_seurat, num.replicate = 100)
        s_seurat <- ScoreJackStraw(s_seurat, dims = 1:20)
        js_p <- s_seurat@reductions$pca@jackstraw$overall.p.values
        mdim <- max(which(js_p[, 2] < 1e-100))
        if (mdim < 10) {
            mdim <- 10
        }
    }
    else {
        mdim = pca_k
    }
    print("finish select pc")
    s_seurat <- FindNeighbors(s_seurat, dims = 1:mdim)
    s_seurat <- FindClusters(s_seurat, resolution = 0.1)
    print("finish step one cluster")
    s_seurat@meta.data$seurat_clusters = paste(i, s_seurat@meta.data$seurat_clusters,
        sep = "_")
    sub_cluster_lab = s_seurat@meta.data$seurat_clusters
    names(sub_cluster_lab) <- rownames(s_seurat@meta.data)
    cn1 <- makeCluster(clustern)
    registerDoParallel(cn1)
    r2 = foreach(j = unique(sub_cluster_lab), .combine = rbind,
        .packages = c("Seurat", "apcluster", "flexclust", "dummies",
            "fastcluster", "cluster", "igraph")) \%dopar\% {
        Label = c()
        cn2 <- names(sub_cluster_lab)[sub_cluster_lab == j]
        sm2 <- s_seurat@reductions$pca@cell.embeddings[cn2, 1:mdim]
        snn2 <- s_seurat@graphs$RNA_snn[cn2, cn2]
        l_snn2 <- snn2
        diag(l_snn2) = 0
        if (is.null(k1)) {
            k1 = round(sqrt(dim(sm2)[1]))
            if (k1 < 50) {
                k1 = 50
            }
        }
        if (k1 < dim(sm2)[1]) {
            if (cluster_method == "lpc") {
                g_r <- graph_from_adjacency_matrix(l_snn2, mode = "undirected",
                  weighted = TRUE)
                lpc_r <- label.propagation.community(g_r)
                cls <- communities(lpc_r)
                cell_num <- sapply(cls, length)
                lpc_l = unlist(lapply(communities(lpc_r), function(cl) length(cl)))
                idx <- rep(1:length(cls), cell_num)
                names(idx) = as.character(unlist(cls))
                idx = idx[rownames(sm2)]
                source("R/snn_search.R")
                gcsample = cls
                exemplars = rep("FALSE", nrow(sm2))
                exemplars[which(rownames(sm2) \%in\% unlist(lapply(gcsample,
                  function(cl) snn_search(cl, s_seurat))))] = "TRUE"
                labels = paste(j, idx, sep = "_")
                labels = cbind(as.character(labels), as.character(exemplars))
                rownames(labels) = rownames(sm2)
                colnames(labels) = c("labels", "exemplars")
                labels = as.data.frame(labels)
                Label = rbind(Label, labels)
            }
            if (cluster_method == "apcluster") {
                apk_r <- apclusterK(negDistMat(r = 2), x = sm2,
                  K = k1)
                cls <- apk_r@clusters
                cell_num <- sapply(cls, length)
                idx <- rep(1:length(apk_r@exemplars), cell_num)
                labels <- data.frame(labels = idx, exemplars = names(unlist(apk_r@clusters)) \%in\%
                  names(unlist(apk_r@exemplars)), row.names = names(unlist(apk_r@clusters)))
                labels$labels = paste(j, labels$labels, sep = "_")
                Label <- rbind(Label, labels)
            }
            if (cluster_method == "kmeans") {
                set.seed(1)
                km_r = kcca(sm2, k = k1, family = kccaFamily("kmeans"),
                  control = list(iter.max = 1000))
                if (center_method == "snn_search") {
                  source("R/snn_search.R")
                  gcsample = lapply(unique(km_r@cluster), function(cl) names(km_r@cluster)[km_r@cluster ==
                    cl])
                  exemplars = rep("FALSE", nrow(sm2))
                  exemplars[which(rownames(sm2) \%in\% unlist(lapply(gcsample,
                    function(cl) snn_search(cl, s_seurat))))] = "TRUE"
                }
                if (center_method == "nearpoint_search") {
                  dsm = rbind(km_r@centers, sm2)
                  dsm_r = as.matrix(dist(dsm))[-c(1:k1), 1:k1]
                  exemplars = rep("FALSE", nrow(sm2))
                  exemplars[apply(dsm_r, 2, function(cl) which(cl ==
                    min(cl)))] = "TRUE"
                }
                labels = paste(j, km_r@cluster, sep = "_")
                labels = cbind(as.character(labels), as.character(exemplars))
                rownames(labels) = rownames(sm2)
                colnames(labels) = c("labels", "exemplars")
                labels = as.data.frame(labels)
                Label = rbind(Label, labels)
            }
            if (cluster_method == "hierarchy") {
                source("R/snn_search.R")
                hr = fastcluster::hclust(as.dist(1 - snn2)^2,
                  method = "ward.D")
                hr_cutree = Cutree(hr, k1)
                names(hr_cutree) = hr$labels
                gcsample = lapply(unique(hr_cutree), function(cl) names(hr_cutree)[hr_cutree ==
                  cl])
                exemplars = rep("FALSE", nrow(sm2))
                exemplars[which(rownames(sm2) \%in\% unlist(lapply(gcsample,
                  function(cl) snn_search(cl, s_seurat))))] = "TRUE"
                labels = paste(j, hr_cutree, sep = "_")
                labels = cbind(as.character(labels), as.character(exemplars))
                rownames(labels) = rownames(sm2)
                colnames(labels) = c("labels", "exemplars")
                labels = as.data.frame(labels)
                Label = rbind(Label, labels)
            }
            if (cluster_method == "pam") {
                pam_r <- pam(sm2, k = k1)
                exemplars = rep("FALSE", nrow(sm2))
                exemplars[which(rownames(sm2) \%in\% rownames(pam_r$medoids))] = "TRUE"
                labels = paste(j, pam_r$clustering, sep = "_")
                labels = cbind(as.character(labels), as.character(exemplars))
                rownames(labels) = rownames(sm2)
                colnames(labels) = c("labels", "exemplars")
                labels = as.data.frame(labels)
                Label = rbind(Label, labels)
            }
        }
        else {
            labels <- data.frame(labels = rep(1:dim(sm2)[1],
                1), exemplars = rep("TRUE", dim(sm2)[1]), row.names = rownames(sm2))
            labels$labels = paste(j, labels$labels, sep = "_")
            Label <- rbind(Label, labels)
        }
        Label
    }
    stopCluster(cn1)
    label <- r2[rownames(s_seurat@meta.data), ]
    s_seurat@meta.data <- cbind(s_seurat@meta.data, label)
    s_seurat@meta.data$pca_k = mdim
    return(s_seurat@meta.data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
