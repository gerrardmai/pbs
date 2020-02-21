#####function
snn_search=function(cl, s_seurat)
{
  snn_table=as.matrix(s_seurat@graphs$RNA_snn[cl,cl])
  colnames(snn_table)=cl
  snn_table[snn_table>0]=1
  center_name=names(colSums(snn_table))[which(colSums(snn_table)==max(colSums(snn_table)))]

  if(length(center_name)>1)
  {
    snn_colsum=colSums(as.matrix(s_seurat@graphs$RNA_snn[cl,center_name]))
    center_name=names(snn_colsum)[which(snn_colsum==max(snn_colsum))][1]
  }
  center_name
}
