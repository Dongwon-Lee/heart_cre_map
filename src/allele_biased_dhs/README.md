## Scripts for allele biased DHS analysis

  * Main R scripts:
    * [ase_analysis_indep.R](./ase_analysis_indep.R) - This should run first
    * [ase_analysis_indep_meta_v2.R](./ase_analysis_indep_meta_v2.R)

  * Input files for the main analysis can be made from bam files using batch_run_quasar_preproc.sh.
  However, bam files remapped by [GSNAP](http://research-pub.gene.com/gmap/) is too big to upload. So, we instead provide the processed files:
  [link](https://drive.google.com/drive/folders/14pNK6Cyaw-tgkT9oJca1xcZHK1oMu4Fx?usp=sharing) 
    * These files should be used as inputs for the main script.

  * Please note that we adopted [QuASAR method](https://github.com/piquelab/QuASAR), 
  and slightly modified it so that we can analyze DNase-seq data sets.
