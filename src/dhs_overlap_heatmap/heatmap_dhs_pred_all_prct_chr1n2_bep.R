# A script for making heatmaps of DHS peaks for observed and predicted Heart DHSs (Figure 2)
#
# by Dongwon Lee

rm(list=ls())
library('gplots')
library('RColorBrewer')
source('my.vioplot.R')

set.seed(2)

plot_heatmap<-function(pred, br.heatmap, colors, col_colors=NULL, maintxt="", n=5000, dist_met="euclidean", clust_met="average") {
    custom_lmat<-rbind(c(0,4,4),
                       c(0,1,1),
                       c(3,2,2),
                       c(0,0,5))
    mat<-pred

    isdhs<-ifelse(rowSums(mat) > 0, T, F)

    print(length(isdhs))
    print(length(which(isdhs)))

    ind1<-sample(which(isdhs), n)

    heatmap.2(mat[ind1, ], 
              distfun=function(x) dist(x, method=dist_met), 
              hclustfun=function(x) hclust(x, method=clust_met),
              dendrogram="col",
              labRow="",
              labCol="", 
              main=maintxt,
              lmat=custom_lmat,
              lhei=c(0.2, 0.1, 2, 0.65),
              lwid=c(0.1, 2, 2),
              margin=c(1, 1),
              ColSideColors=col_colors,
              trace="none", 
              key=TRUE,
              keysize=1,
              key.title=NA,
              key.xlab="Percentile",
              density.info="none",
              breaks=br.heatmap, 
              col=colors) 

    mat[ind1, ]
}

pred1<-read.table("overlap.all_heart_dhs_summits.e300.all.prct.chr1n2.txt", header=T)
pred2<-read.table("overlap.allhrt_p0_avg50_preds.rec50.all.nodhs3.all.prct.chr1n2.txt", header=T)
files<-read.table("heart_related_dhs_files.txt", stringsAsFactors=F)
rdhs_ids<-gsub("-", ".", gsub("_summits.bed$", "", basename(files[[1]])))

n1<-nrow(pred1)
n2<-nrow(pred2)

intcol1<-which(colnames(pred1) %in% rdhs_ids)
intcol2<-which(colnames(pred2) %in% rdhs_ids)

pred1.int<-pred1[, intcol1]
pred2.int<-pred2[, intcol2]

preds<-as.matrix(rbind(pred1.int, pred2.int))

#read cell type classes
celltype_classes<-read.table('heart_related_dhs_files_celltypes.txt', sep="\t", stringsAsFactors=F)
colnames(celltype_classes)<-c("expid", "celltype", "misc")
celltype_classes$expid<-gsub("-", ".", gsub("_summits.bed$", "", basename(celltype_classes$expid)))

expid_celltypeind<-data.frame(expid=colnames(preds), ind=NA)
for(i in 1:nrow(expid_celltypeind)) {
    ind<-grep(expid_celltypeind$expid[i], celltype_classes$expid)
    if (length(ind)>0) expid_celltypeind$ind[i]<-ind[1]
}

#celltype_classes$celltype[expid_celltypeind$ind]

cols<-brewer.pal(9, "Set1")
cols2<-brewer.pal(8, "Set2")
celltype_cols<-rep(cols[9], ncol(preds)) #grey by default

#Levels: 
#heart - red(1)
#microvascular endothelial cells - blue(2)
#muscle - green(3)
#smooth muscle cells - light green (col2:1)
#myoblast - purple (4)
#myocyte - orange (5)
#fibroblast - yallow (6)
#myosatellite - brown (7)
celltype_cols[grep("heart", celltype_classes$celltype[expid_celltypeind$ind])]<-cols[1] 
celltype_cols[grep("endothelial cells", celltype_classes$celltype[expid_celltypeind$ind])]<-cols[2]
celltype_cols[grep("muscle", celltype_classes$celltype[expid_celltypeind$ind])]<-cols[3]
celltype_cols[grep("myoblast", celltype_classes$celltype[expid_celltypeind$ind])]<-cols[4]
celltype_cols[grep("myosatellite", celltype_classes$celltype[expid_celltypeind$ind])]<-cols[5]
celltype_cols[grep("myocyte", celltype_classes$celltype[expid_celltypeind$ind])]<-cols[6]
celltype_cols[grep("fibroblast", celltype_classes$celltype[expid_celltypeind$ind])]<-cols[7]
celltype_cols[grep("smooth muscle", celltype_classes$celltype[expid_celltypeind$ind])]<-cols2[1]

br.heatmap<-seq(-1, 101, 2)

pdf("heatmap_heart_dhs_prct_r5k_bep.pdf", width=3, height=4.5)
par(ps=7)
heatmap_cols1<-colorRampPalette(brewer.pal(8,"Reds")[-1])(length(br.heatmap)-1)
mat1<-plot_heatmap(preds[1:n1, ], br.heatmap, heatmap_cols1, col_colors=celltype_cols, maintxt=NULL, n=500, dist_met="binary", clust_met="average")
dev.off()

pdf("heatmap_heart_pred_prct_r5k_bep.pdf", width=3, height=4.5)
par(ps=7, mar=c(0,0,0,0)+0.2)
heatmap_cols2<-colorRampPalette(brewer.pal(8,"Blues")[-1])(length(br.heatmap)-1)
mat2<-plot_heatmap(preds[(n1+1):(n1+n2), ], br.heatmap, heatmap_cols2, col_colors=celltype_cols, maintxt=NULL, n=500, dist_met="binary", clust_met="average")
dev.off()

pdf("heatmap_heart_labels_prct_r5k_bep.pdf", width=3, height=5)
par(ps=7, cex=1, cex.lab=1)
plot(1, type="n", axes=F, xlab="", ylab="")
legend("topleft",
       legend=c("Fetal Heart tissues",
                "Vascular endothelial cells",
                "Muscle tissues",
                "Myoblasts",
                "Myosatellites",
                "Myocytes",
                "Fibroblasts",
                "Smooth muscle cells"),
       fill=c(cols[1], cols[2], cols[3], cols[4], cols[5], cols[6], cols[7], cols2[1]),
       bty="n")
dev.off()

prct1<-unlist(mat1[mat1>0])
prct2<-unlist(mat2[mat2>0])

# signal strengh 
pdf("vioplot_prct_dhs_pred_bep.pdf", width=4.5, height=5)
my.vioplot(prct1, prct2, names=c("Observed CREs", "Predicted CREs"), col=c(cols[1], cols[2]), wex=scale(sqrt(c(length(prct1), length(prct2))), center=F))
title("DHS signal strength", ylab="Percentile")
dev.off()

# number of peaks histogram
pdf("hist_npeaks_exp_bep.pdf", width=8, height=5)
npeaks1<-apply(mat1, 1, function(x) length(which(x>0))/136)
npeaks2<-apply(mat2, 1, function(x) length(which(x>0))/136)

npeaks.br<-seq(0, 1, 0.1)
h1<-hist(npeaks1, br=npeaks.br, plot=F)
h2<-hist(npeaks2, br=npeaks.br, plot=F)
barplot(rbind(h1$density/10, h2$density/10), 
        col=c(cols[1], cols[2]),
        names.arg=paste("-", seq(10,100, 10), "%", sep=""), 
        beside=T,
        xlab="Proportion of experiments with heart relevant DHSs",
        ylab="Fraction of regions")
dev.off()
