loading<-T

gtexf<-"./GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz"

#1. reading files
if (loading==T) {
    library("gplots")
    library("RColorBrewer")
    library("edgeR")
    
    gtex<-read.table(gtexf, header=T, sep="\t", skip=2, stringsAsFactors=F)
    tissue_names<-gsub(" $", "", gsub("\\.+", " ", colnames(gtex)[3:55]))

    # Tissue name for figure. 
    tissue_desc<-c("Adipose-subcutaneous", "Adipose-visceral", "Adrenal gland", "Artery-aorta",
                   "Artery-coronary", "Artery-tibial", "Bladder", "Brain-amygdala",
                   "Brain-ACC", "Brain-caudate", "Brain-cerebellar hemisphere", "Brain-cerebellum",
                   "Brain-cortex", "Brain-frontal cortex", "Brain-hippocampus", "Brain-hypothalamus",
                   "Brain-NAcc", "Brain-putamen", "Brain-spinal cord", "Brain-SN",
                   "Breast-mammary tissue", "EBV transformed lymphocytes",
                   "Transformed fibroblasts", "Cervix-ectocervix",
                   "Cervix-endocervix", "Colon-sigmoid", "Colon-transverse",
                   "Esophagus-GEJ", "Esophagus-mucosa", "Esophagus-muscularis",
                   "Fallopian tube", "Heart-atrial appendage", "Heart-left ventricle",
                   "Kidney-cortex", "Liver", "Lung", "Minor salivary gland", "Muscle-skeletal",
                   "Nerve-tibial","Ovary","Pancreas", "Pituitary", "Prostate",
                   "Skin-not sun exposed", "Skin-sun exposed",
                   "Small intestine-terminal ileum", "Spleen", "Stomach", "Testis", "Thyroid",
                   "Uterus", "Vagina", "Whole blood")

    colnames(gtex)<-substr(gsub("\\.+", ".", colnames(gtex)), 0, 16)

    gtex.lib.size<-colSums(gtex[,-c(1,2)])
    gtex.lib.size<-gtex.lib.size/median(gtex.lib.size)
    gtex.norm.factors<-calcNormFactors(as.matrix(gtex[, -c(1,2)]))
    gtex.norm.factors<-gtex.norm.factors*gtex.lib.size

    gtex<-data.frame(Name=gtex$Name, Description=gtex$Description, sweep(as.matrix(gtex[,-c(1,2)]), 2, gtex.norm.factors, `/`), stringsAsFactors=F)
    write.table(gtex, "normalized_expr.txt", col.names=T, row.names=F, sep="\t", quote=F)
    excl_ind<-grep("ENSGR", gtex$Name) #exclude PAR genes....
    gtex<-gtex[-excl_ind, ]

    cisbp<-read.table("cisbp_TF_Information_refined.txt", header=T, sep="\t", stringsAsFactors=F)

    c2h2<-read.table("c2h2_motif_tf.txt", header=T, sep="\t", stringsAsFactors=F)
    c2h2.motids<-unlist(lapply(strsplit(c2h2$motif, "_"), function(x) x[1]))
    c2h2.tfs<-unlist(lapply(strsplit(c2h2$motif, "_"), function(x) x[2]))
    c2h2<-data.frame(Motif_ID=c2h2.motids, TF_Name=c2h2.tfs, stringsAsFactors=F)
    
    all_tfs<-unique(c(cisbp$TF_Name, c2h2.tfs))
    length(unique(cisbp$TF_Name)) # 811
    length(unique(c2h2.tfs)) # 94
    length(all_tfs) # 905

    gtex.cisbp<-subset(gtex, Description %in% all_tfs)
    rownames(gtex.cisbp)<-gtex.cisbp$Description
    nrow(gtex.cisbp) # 904 <- due to the ZNF8 

    fimo<-read.table("fimo.cisbpv2_c2h2zf.11mers_nr.uniq.allhrt_p0_avg50_svmw.kmercnt_per_motif.5p.txt", header=T, stringsAsFactors=F)
    cisbp5p_motids<-unique(sort(fimo$motif[fimo$padj_5p<0.01]))
    print(length(cisbp5p_motids))

    # store motif-TF table for future analysis
    tf_mot<-rbind(cisbp[, c("Motif_ID", "TF_Name")], c2h2)
    write.table(tf_mot, "motid_tfname_cisbpv2_c2h2_zfs.txt", row.names=F, col.names=T, sep="\t", quote=F)

    cisbp5p_tfs<-unique(sort(tf_mot$TF_Name[tf_mot$Motif_ID %in% cisbp5p_motids]))
}

##################################################################################
##################################################################################

fpkm_cutoff<-1

#       | no Expr |  Expr
#-------+---------+--------
# no Mot|  m00    |   m10
#    Mot|  m01    |   m11

make_cmat<-function(ntotal, tfset, motset) {
    m11<-length(which(tfset %in% motset))
    m10<-length(tfset) - m11
    m01<-length(motset) - m11
    m00<-ntotal - (m11+m10+m01)

    matrix(c(m00, m01, m10, m11), ncol=2, nrow=2) 
}

find_exp_gene_not_enriched<-function(tfset, motset) {
    tfset[which(!(tfset %in% motset))]
}

res<-c()
for(i in 3:ncol(gtex.cisbp)) {
    gtex.cisbp.other<-subset(gtex.cisbp, gtex.cisbp[,i]>=fpkm_cutoff)

    cmat<-make_cmat(nrow(gtex.cisbp), gtex.cisbp.other$Description, cisbp5p_tfs)
    res<-rbind(res, c(cmat[1,1], cmat[2,1], cmat[1,2], cmat[2,2], fisher.test(cmat, alternative="g")$p.value))

    if (substr(colnames(gtex.cisbp)[i], 0, 5) == "Heart") {
        expgenes_not_enriched <- find_exp_gene_not_enriched(gtex.cisbp.other$Description, cisbp5p_tfs)
        write.table(sort(expgenes_not_enriched), 
                    paste(colnames(gtex.cisbp)[i], "exp_gene_not_enriched.all.txt", sep="."), 
                    sep="\t", col.names=F, row.names=F, quote=F)
    }
}

res<-data.frame(tissue=tissue_desc, res)
colnames(res)[2:6]<-c("m00", "m01", "m10", "m11", "fisher_pval")
#res$padj<-ifelse(res$fisher_pval*53>1, 1, res$fisher_pval*53)
write.table(res, "TF_enrichment_test_all_v3.txt", sep="\t", col.names=T, row.names=F, quote=F)

col.set1<-brewer.pal(9, "Set1")
tissue_colors<-rep(col.set1[9], nrow(res))
heart_ind<-grep("Heart", tissue_desc)
tissue_colors[heart_ind]<-col.set1[1]

pval_ordered<-order(log10(res$fisher_pval))

find_gene_expr_not_enriched
print 
############ plots
linewd <-0.5
wd <- 5
ht <- 2.25
fig.nrows <- 1 
fig.ncols <- 1
pt <- 6
cex.general <- 1 
cex.lab <- 1
cex.axis <- 1
cex.main <- 1
cex.legend <- 1

pdf("barplot_TF_enrichment_test_v3_alltfs.pdf", width=wd, height=ht)
par(mar=c(4.5,2.5,0.5,0.5)+0.1, mgp=c(1.3,0.3,0))
par(cex=cex.general, cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main, ps=pt, lwd=linewd, bty="n", tck=-0.03)
par(font.main=1)

nlogpval<-log10(res$fisher_pval)*(-1)

mp<-barplot(nlogpval[pval_ordered], col=tissue_colors[pval_ordered], ylab=expression(-log[10](italic(P)[Fisher])), main="", space=0.5, border=NA)
text(mp, par("usr")[3], labels=tissue_desc[pval_ordered], srt=45, adj=1.05, xpd=TRUE, cex=5/6)
abline(h=-1*log10(0.01/nrow(res)), lty=2)

mtext("c", 3, line=-1, outer=T, adj=0.02, cex=(4/3), font=2)

dev.off()

###########################
# removing ubiquitous TFs

ubiq_tf_cutoff<-53*0.9
ubiq_tf_inds<-which(apply(gtex.cisbp[,-c(1,2)], 1, function(x) length(which(x>=fpkm_cutoff)))>ubiq_tf_cutoff)
gtex.cisbp.spec<-gtex.cisbp[-c(ubiq_tf_inds), ]
print(nrow(gtex.cisbp))
print(nrow(gtex.cisbp.spec))
cisbp5p.spec<-cisbp5p_tfs[!(cisbp5p_tfs %in% names(ubiq_tf_inds))]

res2<-c()
for(i in 3:ncol(gtex.cisbp.spec)) {
    gtex.cisbp.other<-subset(gtex.cisbp.spec, gtex.cisbp.spec[,i]>=fpkm_cutoff)

    cmat<-make_cmat(nrow(gtex.cisbp.spec), gtex.cisbp.other$Description, cisbp5p.spec)
    res2<-rbind(res2, c(cmat[1,1], cmat[2,1], cmat[1,2], cmat[2,2], fisher.test(cmat, alternative="g")$p.value))

    if (substr(colnames(gtex.cisbp)[i], 0, 5) == "Heart") {
        expgenes_not_enriched <- find_exp_gene_not_enriched(gtex.cisbp.other$Description, cisbp5p_tfs)
        write.table(sort(expgenes_not_enriched), 
                    paste(colnames(gtex.cisbp)[i], "exp_gene_not_enriched.spec90.txt", sep="."), 
                    sep="\t", col.names=F, row.names=F, quote=F)
    }
}
res2<-data.frame(tissue=tissue_desc, res2)
colnames(res2)[2:6]<-c("m00", "m01", "m10", "m11", "fisher_pval")
#res2$padj<-ifelse(res2$fisher_pval*53>1, 1, res2$fisher_pval*53)
write.table(res2, "TF_enrichment_test_spec90_v3.txt", sep="\t", col.names=T, row.names=F, quote=F)

# use the same order for better presentation
#pval_ordered<-order(log10(res2$fisher_pval))

pdf("barplot_TF_enrichment_test_v3_spectfs.pdf", width=wd, height=ht)
par(mar=c(4.5,2.5,0.5,0.5)+0.1, mgp=c(1.3,0.3,0))
par(cex=cex.general, cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main, ps=pt, lwd=linewd, bty="n", tck=-0.03)
par(font.main=1)

nlogpval2<-log10(res2$fisher_pval)*(-1)

mp<-barplot(nlogpval2[pval_ordered], col=tissue_colors[pval_ordered], ylab=expression(-log[10](italic(P)[Fisher])), main="", space=0.5, border=NA)

text(mp, par("usr")[3], labels=tissue_desc[pval_ordered], srt=45, adj=1.05, xpd=TRUE, cex=5/6)
abline(h=-1*log10(0.01/nrow(res2)), lty=2)

mtext("d", 3, line=-1, outer=T, adj=0.02, cex=(4/3), font=2)

dev.off()

###########################
# extract TFs
extract_tfs<-function(fpkmdata, ind, fpkm_cutoff, tfset.pred) {
    fpkmdata.ind<-subset(fpkmdata, fpkmdata[,ind]>=fpkm_cutoff)
    tfset.pred.expr<-fpkmdata.ind$Description[fpkmdata.ind$Description %in% tfset.pred]

    tfset.pred.expr<-subset(fpkmdata.ind, Description %in% tfset.pred.expr)[, c(2, ind)]

    tfset.pred.expr<-tfset.pred.expr[order(tfset.pred.expr[[2]], decreasing=T),]
    colnames(tfset.pred.expr)[1]<-"TF_Name"
    tfset.pred.expr
}

hrt.ind<-grep("Heart.Left", colnames(gtex.cisbp))

gtex.cisbp.hrt.all<-extract_tfs(gtex.cisbp, hrt.ind, fpkm_cutoff, cisbp5p_tfs)
gtex.cisbp.hrt.spec<-extract_tfs(gtex.cisbp.spec, hrt.ind, fpkm_cutoff, cisbp5p.spec)

write.table(gtex.cisbp.hrt.all, "expr_pred_tfs_all_v3.txt", quote=F, row.names=F, col.names=F)
write.table(gtex.cisbp.hrt.spec, "expr_pred_tfs_spec90_v3.txt", quote=F, row.names=F, col.names=F)
