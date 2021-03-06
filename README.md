# Heart CRE map
Supplementary data for:

**"Human cardiac cis-regulatory elements, their cognate transcription factors and regulatory DNA sequence variants"**
Dongwon Lee<sup>\*</sup>, Ashish Kapoor, Alexias Safi, Lingyun Song, Marc K. Halushka, Gregory E. Crawford, Aravinda Chakravarti<sup>\*</sup>

<sup>\*</sup> Correspondence should be addressed to: Aravinda Chakravarti (aravinda.chakravarti AT nyumc DOT org) or Dongwon Lee (dongwon.lee AT nyumc DOT org)

*Notes:*
  * Most big files are separately stored in google drive folders, and only links are provided in this README file.

## DHS peak sets called by MACS2 for the five adult heart samples
  * Chaklab1:
  [narrowPeak](https://drive.google.com/uc?id=18CpfhX_B5Vfm3ff_YXRzLvrm-LJLYvFl&export=download),
  [bigWig](https://drive.google.com/uc?id=1N9nf0Vz9waE8mYBQp2r-6XBGmqcPpjiw&export=download),
  [extended](https://drive.google.com/uc?id=1afEsOu_xFAJ_9kaNRM4IdzUsY-Z2H2Ke&export=download)
  * Chaklab2:
  [narrowPeak](https://drive.google.com/uc?id=1GqB-iQ0fbiiC_sOCatM-vijwcy_6DFR6&export=download),
  [bigWig](https://drive.google.com/uc?id=189fB_VCfCnwTNJSPeW8XFlytKyLx7W_i&export=download),
  [extended](https://drive.google.com/uc?id=1BgurZ6rq3WnC381gLNIQL4VBmBoAfrTW&export=download)
  * ENCODE1:
  [narrowPeak](https://drive.google.com/uc?id=1TvwX0yB0Z3BO-VCwTM8eFGsWozRPS9Gv&export=download),
  [bigWig](https://drive.google.com/uc?id=17oZTZQxK0AzT7yVNqwPasM48u66oIAZ_&export=download),
  [extended](https://drive.google.com/uc?id=1GYWespB_00ZMcWBYxTDhNqPE6j2UnvMW&export=download)
  * ENCODE2:
  [narrowPeak](https://drive.google.com/uc?id=1P28RAOzZTN3blqy2OUNN0L87QEft6O7J&export=download),
  [bigWig](https://drive.google.com/uc?id=1G__924LPbTMhfdBACsCEsc2JYixkaCFe&export=download),
  [extended](https://drive.google.com/uc?id=1ZPAL3CxVteZzYybC3Q-bl-yjTiJ21heC&export=download)
  * Roadmap:
  [narrowPeak](https://drive.google.com/uc?id=1lzJP-2MvxCwBTcC1lveExgOVumS-hSs-&export=download),
  [bigWig](https://drive.google.com/uc?id=1OjPRLFHfrkoTZC_ET6ZlzZc9PMZCgccH&export=download),
  [extended](https://drive.google.com/uc?id=17CNwE1jmc_xSm6uwQpiVC3w5bKqbq9nW&export=download)
  * Combined Heart DHSs: [DHSs/all_heart_dhs_summits.e300.bed](./DHSs/all_heart_dhs_summits.e300.bed)

*Notes:*
  * Click [here](https://drive.google.com/drive/folders/10N8sbZ5TKVrAnGJou7PuCcbljkBmW8WY) to View all files
  * Raw data are available in the GEO database (accession no. GSE104989)
  * Extended peaks are defined as 600bp peaks centered at the MACS2 summits. Any overlapping peaks were then merged. Peaks in chrY and chrM were excluded.

## LS-GKM training
  * Generic models:
    * Sequence files in FASTA format for training and test:
    [Chaklab1](https://drive.google.com/uc?id=1G4m_l-HS-2cEpRE-owrGypsfHqojFoUZ&export=download),
    [Chaklab2](https://drive.google.com/uc?id=1-GVbjA5G97rhCSH1m7_ED01VBzpkA9jW&export=download),
    [ENCODE1](https://drive.google.com/uc?id=1ClTtlZ3trQQqWYBEblOvUniXMNVMzGfz&export=download),
    [ENCODE2](https://drive.google.com/uc?id=1piMkzC1OJxLvic9-vjH3qQE78-1fqYf0&export=download),
    [Roadmap](https://drive.google.com/uc?id=1kHLojiEKlofbEgYp5wpoK9bwXEbuXT9m&export=download)
    ([View All](https://drive.google.com/drive/folders/1tI9r-d-gEmIP1bYLnuCvsXa-pN1b-SvL))
    * [trained model files](https://drive.google.com/uc?id=1akku984BPNM8xoINxZx0kI0ebxXtwjeu&export=download)
  * Specific models:
    * Sequence files in FASTA format for training and test:
    [Chaklab1](https://drive.google.com/uc?id=1T6AzmOOns4pS1qRZLrhT4KOebmiucBlx&export=download),
    [Chaklab2](https://drive.google.com/uc?id=1yMPqsAQCCqG1KOXTYSVSvg-DHBx2uv-Z&export=download),
    [ENCODE1](https://drive.google.com/uc?id=1lQtBYyEA3hVoWyyP6Mrei0gHbVIYKF4d&export=download),
    [ENCODE2](https://drive.google.com/uc?id=1yKtMaOCe5l_TGARLTQPHTzelLyKC1o_B&export=download),
    [Roadmap](https://drive.google.com/uc?id=13EEpGTabssDQ57NMt4aG629m4a-FlvOo&export=download)
    ([View All](https://drive.google.com/drive/folders/1ckZ9V44wDg0kFg_BAErxWrarIaalmdGM))
    * [trained model files](https://drive.google.com/uc?id=1Q276l31PJi9VZqzuCBrTQIxUmT-Xy9SX&export=download)

## Genome-wide prediction of cardiac DHSs
  * Scores of all 600bp windows:
    [Chaklab1](https://drive.google.com/uc?id=1Pnf83Cba7-53j4Vwyt18iSa6qCoEX3G7&export=download),
    [Chaklab2](https://drive.google.com/uc?id=1NljkE6MdTYhjotwpPeAbDEvfGh5lpxzP&export=download),
    [ENCODE1](https://drive.google.com/uc?id=15WQl6q5NPo05spY_nEIkL4drVOujYEkY&export=download),
    [ENCODE2](https://drive.google.com/uc?id=1KIDCU2cE7DbqJJKmmlpRb9ThNktilEpQ&export=download),
    [Roadmap](https://drive.google.com/uc?id=1bQrMBU0fujllXgg-FO4g2VKSH05b1nZD&export=download)
    ([View All](https://drive.google.com/drive/folders/1Gj97jIgukILKFkL_Xcm4c8jMStpM3unR))
  * Predicted DHSs
    * all predicted DHSs:
    [DHSs/allhrt_p0_avg50_preds.bep.all.bed](./DHSs/allhrt_p0_avg50_preds.bep.all.bed)
    * Predicted DHSs after removing DHS overlaps:
    [DHSs/allhrt_p0_avg50_preds.bep.nodhs3.bed](./DHSs/allhrt_p0_avg50_preds.bep.nodhs3.bed)

## Human transcription factors (TFs) and their position weight matrices (PWMs)
  * CIS-BP database
    * [TF_pwms/TF_Information_refined_v2.txt](./TF_pwms/TF_Information_refined_v2.txt):
    Information about TFs and their selected PWMs from the [CIS-BP database](http://cisbp.ccbr.utoronto.ca/) (Weirauch et al., Cell, 2014)
    * [TF_pwms/CISBP_Human_refined_v2_nr.meme](./TF_pwms/CISBP_Human_refined_v2_nr.meme):
    775 non-redundant PWMs in MEME format
    * [TF_pwms/CISBP_Human_refined_v2_nr_8_11.meme](./TF_pwms/CISBP_Human_refined_v2_nr_8_11.meme):
    submatrix of PWMs with 8-11 bp lengths used for 11-mer scanning with FIMO

  * C2H2 Zinc Finger Transcription Factors
    * [TF_pwms/c2h2_zfs.meme](./TF_pwms/c2h2_zfs.meme):
    94 C2H2 Zinc Finger PWMs (Schmitges et al., Genome Res, 2016). Note that ZNF8 was later excluded from this study (no expression value in GTEx)
    * [TF_pwms/c2h2_zfs_8_11.meme](./TF_pwms/c2h2_zfs_8_11.meme):
    submatrix of PWMs with 8-11 bp lengths used for 11-mer scanning

## Predicted cardiac TFs
  * [TF_analysis/fimo.cisbpv2_c2h2zf.11mers_nr.uniq.allhrt_p0_avg50_svmw.kmercnt_per_motif.5p.txt](./TF_analysis/fimo.cisbpv2_c2h2zf.11mers_nr.uniq.allhrt_p0_avg50_svmw.kmercnt_per_motif.5p.txt):
  Enrichment tests of PWM matched 11-mers in the top 5th percentile of the SVM weight distribution
    * **motif**: Motif ID
    * **total**: total number of 11-mers matching the PWM
    * **kmer5p**: number of top 5th percentile scoring 11-mers that match the PWM
    * **kmer5p_pcrt**: (kmer5p)/(total)
    * **poisson_pval_5p**: Poisson test P-value (1-side) with the null rate 0.05
    * **padj_5p**: Bonferroni adjusted P-value
  * [TF_analysis/TF_enrichment_test_all_v3.txt](./TF_analysis/TF_enrichment_test_all_v3.txt):
  Association tests of predicted TFs and expressed TFs for multiple human tissues
    * **tissue**: GTEx human tissue
    * **m00**: number of genes that are *NOT expressed* in the given tissue and *NOT enriched* in heart DHSs
    * **m01**: number of genes that are *NOT expressed* in the given tissue but *enriched* in heart DHSs
    * **m10**: number of genes that are *expressed* in the given tissue but *NOT enriched* in heart DHSs
    * **m11**: number of genes that are *expressed* in the given tissue and *enriched* in heart DHSs
    * **fisher_pval**: Fisher's exact test P-value (1-side) using 2-by-2 contingency table
  * [GTEx V6p normalized gene expression](https://drive.google.com/uc?id=1xiki7TWi_p4A_G0dcgSwfh92SnKttKZS&export=download)
  * [TF_analysis/TF_enrichment_test_all_v3.txt](./TF_analysis/TF_enrichment_test_spec90_v3.txt):
  Association tests of predicted TFs and expressed TFs for multiple human tissues after excluding commonly expressed TFs (>90% of tissues). The format is the same as above.
  * [TF_analysis/cisbpv2_c2h2zf_motif_clus.0.98.txt](./TF_analysis/cisbpv2_c2h2zf_motif_clus.0.98.txt):
  Clustering PWMs based on the top 5th percentile 11-mer matching profile ([TF_analysis/cisbpv2_motif_kmer_heatmap_data.txt.gz](./TF_analysis/cisbpv2_motif_kmer_heatmap_data.txt.gz)). UPGMA Hierarchical cluster analysis followed by tree cuts (cutoff=0.98) was applied to find clusters.  "Asymmetric binary" is used to calculate distance between PWMs.
    * **clust**: cluster ID
    * **motid**: Motif ID
    * **center**: Exemplary motif ID, defined as the motif that has the minimum average distance to other motifs within its cluster
    * **avgdist**: Average distance of the examplary motif ID to others

## Predicted cardiac regulatory variants
  * deltaSVM scores of common variants (>1% MAF in EUR) using the cardiac generic model:
  [gkmscore.allhrt_avg50.1kgp3.hg19.eur.maf1.e10.dsvm.bed.gz](https://drive.google.com/uc?id=1OU5I3sqqMHSSoH3ZYOkBLy1SubpKtWMb&export=download)
    * cutoff: abs(deltaSVM_generic) > 0.894
  * deltaSVM scores of common variants (>1% MAF in EUR) using the cardiac specific model:
  [gkmscore.allhrt_spec30_avg50.1kgp3.hg19.eur.maf1.e10.dsvm.bed.gz](https://drive.google.com/uc?id=19sKVqCzxExxrPSD5ftszWG_tm4P8D56F&export=download)
    * cutoff: abs(deltaSVM_specific) > 0.586
    * Columns
      1. chr
      2. start
      3. end
      4. RSID
      5. ref. allele
      6. alt. allele
      7. alt. allele frequency in EUR
      8. deltaSVM
      9. ref. SVM score
      10. alt. SVM score

  * generic deltaSVM scores of SNPs overlapping heart CREs:
  [gkmscore.allhrt_avg50.1kgp3.hg19.eur.maf1.e10.dsvm.allhrtcre.bed.gz](https://drive.google.com/uc?id=1GQ2jj1O5GLBBr0zsQKG52uceSzLe2-Od&export=download)
  * specific deltaSVM scores of SNPs overlapping heart CREs:
  [gkmscore.allhrt_spec30_avg50.1kgp3.hg19.eur.maf1.e10.dsvm.allhrtcre.bed.gz](https://drive.google.com/uc?id=1sSzCOTrdfDx5-Q7azPTpmpcbwy0KryRn&export=download)

  * *significant* generic deltaSVM SNPs overlapping heart CREs:
  [gkmscore.allhrt_avg50.1kgp3.hg19.eur.maf1.e10.dsvm.allhrtcre.pos.bed.gz](https://drive.google.com/uc?id=1I1hS5e3x__7XbQvTfFs_u2JKjbnXpy5A&export=download)
  * *significant* specific deltaSVM SNPs overlapping heart CREs:
  [gkmscore.allhrt_spec30_avg50.1kgp3.hg19.eur.maf1.e10.dsvm.allhrtcre.pos.bed.gz](https://drive.google.com/uc?id=1_w-kgYJ-jILCG0Cbxeh6ODHS-1USx5hd&export=download)

## Locuszoom plots for QT-interval associated variants
  * [qtigc_locuszoom](./qtigc_locuszoom):
  Visualization of the 35 independent loci associated with QT-interval (Arking et al., Nature Genet, 2014) using LocusZoom (Pruim et al., Bioinformatics, 2010).

*last updated: 06/12/2018*
