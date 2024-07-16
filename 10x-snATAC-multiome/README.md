# Reproducing this

### Analyze chromosome 22 -- only use one chromosome b/c its huge

```
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/10k_PBMC_Multiome_nextgem_Chromium_X/10k_PBMC_Multiome_nextgem_Chromium_X_atac_possorted_bam.bam
samtools index 10k_PBMC_Multiome_nextgem_Chromium_X_atac_possorted_bam.bam
samtools view 10k_PBMC_Multiome_nextgem_Chromium_X_atac_possorted_bam.bam chr22 -b -o 10k_PBMC_Multiome_nextgem_Chromium_X_atac.chr22.bam

samtools view 10k_PBMC_Multiome_nextgem_Chromium_X_atac_possorted_bam.bam | cut -f 2,11,13 | grep "MD" | grep -e A -e C -e T -e G > 10k_PBMC_Multiome_nextgem_Chromium_X_atac_chr22.tsv
```

### remap bulk atac from original mito paper
```
bwa mem -t 8 $idx SRR7245894_1.fastq SRR7245894_2.fastq | samtools sort -@ 8 -o SRR7245894.bam
samtools view SRR7245894.bam  | cut -f 2,11,13 | grep "MD" | awk 'length($3) > 7 {print $0}' > SRR7245894_bulkATAC.MD.tsv
```