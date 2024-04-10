# Reproducing this

```
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/10k_PBMC_Multiome_nextgem_Chromium_X/10k_PBMC_Multiome_nextgem_Chromium_X_atac_possorted_bam.bam
samtools index 10k_PBMC_Multiome_nextgem_Chromium_X_atac_possorted_bam.bam
samtools view 10k_PBMC_Multiome_nextgem_Chromium_X_atac_possorted_bam.bam chr22 -b -o 10k_PBMC_Multiome_nextgem_Chromium_X_atac.chr22.bam


```

