##################################################################
## modified QuASAR Workflow
##################################################################
args<-commandArgs(trailingOnly=T)

snpsetid<-args[1]
min.cov<-as.numeric(args[2]) #5

library('QuASAR')
source('aseInferenceSingle.R')

fileNames<-list.files(".", "*.quasar.in.gz")

##################################################################    
## ASE analysis for each sample
##################################################################
for(ii in c(1:length(fileNames))) {
    fn <- fileNames[ii]
    expid<-gsub(".quasar.in.gz", "", fn)
    print(expid)

    ase.dat.single <- UnionExtractFields(fn, combine=TRUE)
    ase.dat.single.gt <- PrepForGenotyping(ase.dat.single, min.coverage=min.cov)
    ase.single <- fitAseNull(ase.dat.single.gt$ref, ase.dat.single.gt$alt, log.gmat=log(ase.dat.single.gt$gmat))
    str(ase.single)
    head(ase.single$gt)

    ## Saving the output genotype probabilities
    out_dat <- data.frame(ase.dat.single.gt$annotations[, -5], map=ase.single$gt)
    write.table(out_dat, file=paste(expid, snpsetid, min.cov, 'genotypes.txt', sep="."), row.names=F, col.names=F, quote=F, sep="\t")

    ase.single.inf <- aseInferenceSingle(gts=ase.single$gt, eps.vect=ase.single$eps, priors=ase.dat.single.gt$gmat, ref.mat=ase.dat.single.gt$ref, alt.mat=ase.dat.single.gt$alt, min.cov=min.cov, sample.name=expid, annos=ase.dat.single.gt$annotations)

    write.table(ase.single.inf, paste(expid, snpsetid, min.cov, "quasar.ase.txt", sep="."), quote=F, row.names=F, col.names=T, sep="\t")
}
