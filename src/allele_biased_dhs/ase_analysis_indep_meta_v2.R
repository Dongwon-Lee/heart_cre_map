##########################################
# Meta analysis of QuASAR
##########################################
args<-commandArgs(trailingOnly=T)

snpsetid<-args[1] 
min.cov<-as.numeric(args[2]) #5

expids<-c("chaklab_dnaseseq_ap.8.2.mapped.mapq30",
          "chaklab_dnaseseq_tp.8.2.mapped.mapq30",
          "ENCFF000SPN.cut.8.2.mapped.mapq30",
          "ENCFF000SPP.cut.8.2.mapped.mapq30",
          "roadmap_uw_fhrt_r2",
          "roadmap_uw_fhrt_r3",
          "roadmap_uw_fhrt_r4",
          "roadmap_uw_fhrt_r5",
          "roadmap_uw_fhrt_r6",
          "roadmap_uw_fhrt_r7",
          "roadmap_uw_fhrt_r8",
          "roadmap_uw_fhrt_r9",
          "roadmap_uw_fhrt_r10",
          "roadmap_uw_fhrt_r11",
          "roadmap_uw_fhrt_r12",
          "roadmap_uw_fhrt_r13",
          "roadmap_uw_hrt")

fileNames<-paste(expids, snpsetid, min.cov, "quasar.ase.txt", sep=".")
nexp<-length(fileNames)

res<-c()
for(ii in c(1:length(fileNames))) {
    fn<-fileNames[ii]
    cat(paste("reading", fn, "\n"))

    x<-read.table(fileNames[ii], header=T, stringsAsFactors=F)
    colnames(x)[4:9]<-paste(colnames(x)[4:9], ii, sep="")
    if (ii == 1) {
        res <- x
    } else {
        res <- merge(res, x, all=T)
    }
}

#reordering
res<-res[order(res$chr, res$pos0), ]

#calculate 
aux <- t(sapply(1:nrow(res), function(ii) {
            cat(paste(ii, "\n"))
            pvals <- res[ii, paste("pval", c(1:nexp), sep="")]
            As <- sum(res[ii, paste("A", c(1:nexp), sep="")], na.rm=T)
            Rs <- sum(res[ii, paste("R", c(1:nexp), sep="")], na.rm=T)
            betas <- res[ii, paste("beta", c(1:nexp), sep="")]
            betas <- betas[!is.na(betas)]
            inconsistent <- any(betas * betas[1] < 0, na.rm=T)
            avg_rho <- mean(t(res[ii, paste("rho", c(1:nexp), sep="")]), na.rm=T)
            chisqval <- (-2)*sum(log(pvals), na.rm=T)
            nsamples <- sum(!is.na(pvals))
            mpval <- pchisq(chisqval, 2*nsamples, lower.tail=F)
            c(As, Rs, avg_rho, chisqval, nsamples, inconsistent, mpval) 
        }))

bedinfo <- cbind(res[[2]], res[[3]], res[[3]]+1, res[[1]])
colnames(bedinfo)<-c("chr", "pos0", "pos", "rsID")
colnames(aux)<-c("Asum", "Rsum", "rhoAvg", "chisq", "nsamples", "inconsistent", "metaPval")
write.table(cbind(bedinfo, aux), paste("ase_analysis_indep_meta_v2", snpsetid, min.cov, "txt", sep="."), quote=F, row.names=F, col.names=T, sep="\t")
