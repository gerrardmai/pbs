#
#
#
#
library(Seurat)
library(apcluster)
library(foreach)
library(doParallel)
library(flexclust)
library(dummies)
library(RColorBrewer)
library(ggplot2)
library(fastcluster)
library(cluster)
library(igraph)


#####function

select_cell_multiday=function(count_data=NULL,count_logdata=NULL,pca_k=NULL,cell_anno=NULL, pca_method="varianceplot",cluster_method="kmeans",center_method="snn_search",k1=NULL,clustern=4)
{

  ####detect input
  if(is.null(count_data) & is.null(count_logdata))
  {
    stop("empty input")
  }
  if(is.null(cell_anno))
  {
    stop("empty cell_anno")
  }

  if(length(unique(cell_anno$Day))<2)
  {
    stop("use other function names select_cell_singleday")
  }

  cn1 <- makeCluster(clustern)
  registerDoParallel(cn1)
  r2=foreach(i=unique(names(table(cell_anno$Day))),
             .combine = rbind,
             .packages=c("Seurat","apcluster","flexclust","dummies","fastcluster","cluster","igraph")) %dopar% {

               cn <- rownames(cell_anno)[cell_anno$Day==i]

               ###CREAT OBJECT
               if( ! is.null(count_data))
               {
                 sm <- count_data[,cn]
                 ######seurat with logdata counts
                 s_seurat <- CreateSeuratObject(counts = sm, project = i)
                 s_seurat <- NormalizeData(s_seurat)
               }else if (is.null(count_data) & ! is.null(count_logdata))
               {
                 sm_log <- count_logdata[,cn]
                 s_seurat <- CreateSeuratObject(counts = sm_log, project = i)
               }

               ###hvg
               s_seurat <- FindVariableFeatures(s_seurat, selection.method = "vst")

               #### run pca
               s_seurat <- ScaleData(s_seurat)
               s_seurat <- RunPCA(s_seurat, features = VariableFeatures(object = s_seurat))

               ####select dim pc
               if(is.null(pca_k) & pca_method=="varianceplot")
               {
                 ep <- ElbowPlot(s_seurat,ndims = 50)
                 mdim <- sum(abs(diff(ep$data$stdev))>0.01)+1
                 if(mdim<10)
                 {
                   mdim <- 10
                 }

               }else if(is.null(pca_k) & pca_method=="jackstraw")
               {
                 s_seurat <- JackStraw(s_seurat, num.replicate = 100)
                 s_seurat <- ScoreJackStraw(s_seurat, dims = 1:20)

                 js_p <- s_seurat@reductions$pca@jackstraw$overall.p.values
                 mdim <- max(which(js_p[,2] < 1e-100))
                 if(mdim<10)
                 {
                   mdim <- 10
                 }
               }else{
                 mdim=pca_k
               }

               ####find neigh and cluster

               s_seurat <- FindNeighbors(s_seurat, dims = 1:mdim)
               s_seurat <- FindClusters(s_seurat, resolution = 0.1)

               ###check 1
               #DimPlot(s_seurat, reduction = "pca")

               ###save cluster
               s_seurat@meta.data$seurat_clusters=paste(i,s_seurat@meta.data$seurat_clusters,sep = "_")
               sub_cluster_lab=s_seurat@meta.data$seurat_clusters
               names(sub_cluster_lab) <- rownames(s_seurat@meta.data)



               Label=c()
               for(j in unique(sub_cluster_lab))
               {
                 cn2 <- names(sub_cluster_lab)[sub_cluster_lab== j]
                 sm2 <- s_seurat@reductions$pca@cell.embeddings[cn2,1:mdim]
                 snn2 <- s_seurat@graphs$RNA_snn[cn2,cn2]
                 l_snn2 <- snn2
                 diag(l_snn2)=0

                 ####caculate the k
                 if(is.null(k1))
                 {
                   k1=round(sqrt(dim(sm2)[1]))
                   if(k1<50){
                     k1=50
                   }
                 }

                 if(k1<dim(sm2)[1])
                 {
                   if(cluster_method=="lpc")
                   {
                     g_r <- graph_from_adjacency_matrix(l_snn2,mode = "undirected",weighted=TRUE)
                     lpc_r <- label.propagation.community(g_r)

                     cls <- communities(lpc_r)
                     cell_num <- sapply(cls,length)
                     lpc_l=unlist(lapply(communities(lpc_r), function(cl) length(cl)))
                     idx <- rep(1:length(cls), cell_num)
                     names(idx)=as.character(unlist(cls))
                     idx=idx[rownames(sm2)]

                     #####snn_search
                     gcsample=cls
                     exemplars=rep("FALSE",nrow(sm2))
                     exemplars[which(rownames(sm2) %in% unlist(lapply(gcsample, function(cl) snn_search(cl,s_seurat))))]="TRUE"


                     labels=paste(j,idx,sep = "_")
                     labels=cbind(as.character(labels),as.character(exemplars))
                     rownames(labels)=rownames(sm2)
                     colnames(labels)=c("labels","exemplars")
                     labels=as.data.frame(labels)
                     Label=rbind(Label,labels)

                   }

                   if(cluster_method=="apcluster")
                   {
                     ####apcluster
                     apk_r <- apclusterK(negDistMat(r=2),x=sm2,K=k1)


                     cls <- apk_r@clusters
                     cell_num <- sapply(cls,length)
                     idx <- rep(1:length(apk_r@exemplars), cell_num)
                     labels <- data.frame(labels = idx,
                                          exemplars = names(unlist(apk_r@clusters))%in%names(unlist(apk_r@exemplars)),
                                          row.names = names(unlist(apk_r@clusters)))
                     labels$labels=paste(j,labels$labels,sep = "_")
                     Label <- rbind(Label,labels)
                   }

                   if(cluster_method=="kmeans")
                   {
                     set.seed(1)
                     km_r=kcca(sm2,k=k1,family=kccaFamily("kmeans"),control =list(iter.max=1000))


                     if(center_method=="snn_search")
                     {
                       gcsample=lapply(unique(km_r@cluster), function(cl) names(km_r@cluster)[km_r@cluster ==cl])

                       exemplars=rep("FALSE",nrow(sm2))
                       exemplars[which(rownames(sm2) %in% unlist(lapply(gcsample, function(cl) snn_search(cl,s_seurat))))]="TRUE"
                     }


                     if(center_method=="nearpoint_search")
                     {
                       dsm=rbind(km_r@centers,sm2)
                       dsm_r=as.matrix(dist(dsm))[-c(1:k1),1:k1]

                       exemplars=rep("FALSE",nrow(sm2))
                       exemplars[apply(dsm_r, 2, function(cl) which(cl==min(cl)))]="TRUE"
                     }
                     labels=paste(j,km_r@cluster,sep = "_")
                     labels=cbind(as.character(labels),as.character(exemplars))
                     rownames(labels)=rownames(sm2)
                     colnames(labels)=c("labels","exemplars")
                     labels=as.data.frame(labels)
                     Label=rbind(Label,labels)


                   }

                   if(cluster_method=="hierarchy")
                   {
                     hr=fastcluster::hclust(as.dist(1-snn2)^2,method='ward.D')
                     hr_cutree=Cutree(hr,k1)
                     names(hr_cutree)=hr$labels

                     gcsample=lapply(unique(hr_cutree), function(cl) names(hr_cutree)[hr_cutree ==cl])

                     exemplars=rep("FALSE",nrow(sm2))
                     exemplars[which(rownames(sm2) %in% unlist(lapply(gcsample, function(cl) snn_search(cl,s_seurat))))]="TRUE"

                     labels=paste(j,hr_cutree,sep = "_")
                     labels=cbind(as.character(labels),as.character(exemplars))
                     rownames(labels)=rownames(sm2)
                     colnames(labels)=c("labels","exemplars")
                     labels=as.data.frame(labels)
                     Label=rbind(Label,labels)

                   }

                   if(cluster_method=="pam")
                   {
                     pam_r <- pam(sm2,k=k1)

                     exemplars=rep("FALSE",nrow(sm2))
                     exemplars[which(rownames(sm2) %in% rownames(pam_r$medoids))]="TRUE"

                     labels=paste(j,pam_r$clustering,sep = "_")
                     labels=cbind(as.character(labels),as.character(exemplars))
                     rownames(labels)=rownames(sm2)
                     colnames(labels)=c("labels","exemplars")
                     labels=as.data.frame(labels)
                     Label=rbind(Label,labels)
                   }
                 }else{
                   labels <- data.frame(labels = rep(1:dim(sm2)[1], 1),
                                        exemplars = rep("TRUE",dim(sm2)[1]),
                                        row.names = rownames(sm2))
                   labels$labels=paste(j,labels$labels,sep = "_")
                   Label <- rbind(Label,labels)
                 }

               }

               ###sort Label
               Label <- Label[rownames(s_seurat@meta.data),]
               s_seurat@meta.data <- cbind(s_seurat@meta.data,Label)
               s_seurat@meta.data$pca_k=mdim

               s_seurat@meta.data

             }
  stopCluster(cn1)
  return(r2)
}





