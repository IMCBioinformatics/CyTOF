

library(Rtsne)
library(flowCore)
library(cytofkit)
library(dplyr)
source("../../code/helper_functions.R")




ff<-read.FCS("../../data/See_Science_2017/normalized_data/PBMC1.fcs")
#ff<-cytof_exprsExtract(fcsFile = "../../data/See_Science_2017/raw/raw files/Run1.fcs",transformMethod = 'none')
#concatFCS(ff)
#ffn<-normCytof(x=ff, y="dvs",remove_beads=FALSE)

ff2<-randomize_zeros(ff)
fft<-transform(fcsFile = ff2,transformMethod = "logicle",verbose = TRUE)

fixedNum=5000
ffn<-sample_n(tbl = as.data.frame(fft),size = fixedNum,replace = ifelse(nrow(fft) < fixedNum, TRUE, FALSE))

#select markers
ffn<-ffn[,c(1,9,19,20,33,23,24,48,41,44)]


##############################################
# fft<-transform(fcsFile = ff,transformMethod = "logicle",verbose = TRUE)
# 
# fixedNum=5000
# ffn<-sample_n(tbl = as.data.frame(fft),size = fixedNum,replace = ifelse(nrow(fft) < fixedNum, TRUE, FALSE))
# 
# #select markers
# ffn<-ffn[,c(1,9,19,20,33,23,24,48,41,44)]
# 
# plot_marker_exp(as.matrix(ffn),prefix = "transformed_norand")
###################################################

ffnt<-Rtsne(ffn)
ggsave(ggplot(,aes(x=ffnt$Y[,1],y=ffnt$Y[,2]))+
         ylab("tSNE-2")+xlab("tSNE-1") +
         geom_point(aes(colour =ffn[,1]),size=0.5)+
         scale_colour_gradientn(colours=rev(rainbow(4)),name="CD45")
       ,filename = "img/tsne1.png")



ffnt<-Rtsne(ffn)
ggsave(ggplot(,aes(x=ffnt$Y[,1],y=ffnt$Y[,2]))+
         ylab("tSNE-2")+xlab("tSNE-1") +
         geom_point(aes(colour =ffn[,1]),size=0.5)+
         scale_colour_gradientn(colours=rev(rainbow(4)),name="CD45")
       ,filename = "img/tsne2.png")

## cluster
FS<-cytof_cluster(ffnt$Y,ffn,method = "FlowSOM",FlowSOM_k = 20)
data.frame(dim1=ffnt$Y[,1],dim2=ffnt$Y[,2],FLOWSOM=FS)->data
p<-cytof_clusterPlot(data,xlab="dim1",ylab="dim2",cluster = "FLOWSOM")
ggsave(plot = p,filename = "img/FS_tsne.png",device ="png" )


CX<-cytof_cluster(ffnt$Y,ffn,method = "ClusterX")
data.frame(dim1=ffnt$Y[,1],dim2=ffnt$Y[,2],ClusterX=CX)->data
p<-cytof_clusterPlot(data,xlab="dim1",ylab="dim2",cluster = "ClusterX")
ggsave(plot = p,filename = "img/CX_tsne.png",device ="png" )



plot_marker_exp(as.matrix(ffn),prefix = "transformed")
library(pheatmap)
pheatmap(mat = ffn)

plot_tsne_marker(ffn,ffnt$Y,prefix = "transformed_tsne")

#PCA
pca<-prcomp(ffn)
plot(pca$x[,1],pca$x[,2])

ggsave(ggplot(,aes(pca$x[,1],pca$x[,2]))+ 
         geom_point(aes(colour =ffn[,1]),size=0.5)+
         scale_colour_gradientn(colours=rev(rainbow(4)),name="CD45")+
         ylab("PCA2")+xlab("PCA1"),filename = "img/pca_CD45.png")

