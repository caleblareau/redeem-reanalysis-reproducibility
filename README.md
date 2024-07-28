# ReDeeM Reanalyses

Repository for reproducing followup analyses of the ReDeeM manuscript. 
This repository assumes the input of the processed data from the authors (as well as other
public datasets) and can be used with the code here to reproduce the summary statistics
and figure panels in the _bioRxiv_ response paper. 


# Setup

## Clone the repository and organize according the folder hierarchy below:

```
git clone https://github.com/caleblareau/redeem-reanalysis-reproducibility
```

Place the repository folder at the same level as another folder called `redeem-downloaded`
that will be populated with data from the authors as we show below:

```
redeem-downloaded
├── mito_data_redeem
│   ├── GSE261078_downloaded
│   ├── Old1.BMMC.Consensus.final
│   ├── Old1.HSPC.Consensus.final
│   ├── Old2.BMMC.Consensus.final
│   ├── Old2.HSPC.Consensus.final
│   ├── Youn2.BMMC.Consensus.final
│   ├── Youn2.HPC.Consensus.final
│   ├── Youn2.HSC.Consensus.final
│   ├── Young1.T1.BMMC.Consensus.final
│   ├── Young1.T1.HPC.Consensus.final
│   ├── Young1.T1.HSC.Consensus.final
│   ├── Young1.T2.BMMC.Consensus.final
│   ├── Young1.T2.HPC.Consensus.final
│   └── Young1.T2.HSC.Consensus.final
└── seurat_data_redeem
    ├── Old1.BMMC_HSPC.Seurat.RDS
    ├── Old2.BMMC_HSPC.Seurat.RDS
    ├── Young1.All.T1.Seurat.RDS
    ├── Young1.All.T2.Seurat.RDS
    ├── Young1.HSC.T1T2.Seurat.RDS
    ├── Young2.All.Seurat.RDS
    └── Young2.HSC.Seurat.RDS
redeem-reanalysis-reproducibility
├── 10x-snATAC-multiome
├── Additional-donors_GSE261078
├── README.md
├── code
├── data
├── final_plots
├── mouse-data_GSE259284
└── output
```


## Downloading ReDeeM data

The folder structure above specifies the data organization. Follow these instructions to populate the data:

- For the Seurat cell state data, download the individual `.rds` files from [this FigShare link](https://figshare.com/articles/dataset/Annotated_Seurat_objects/23290004/1) and 
place the files in `redeem-downloaded/seurat_data_redeem`.

- For the mitochondrial data, download the `.zip` files from [this FigShare link](https://figshare.com/articles/dataset/ReDeeM_raw_mutation_calling/24418966/1) and 
place the files in `redeem-downloaded/mito_data_redeem`.

- For the extended donors, download the processed data from [this GEO accession](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261078)
and place the folders (with the .txt files gzipped) in the `redeem-downloaded/mito_data_redeem/GSE261078_downloaded`

# Key scripts for reanalysis

### Connectivity quantifications

The proportions of "1s" in the ReDeeM matrix is quantified in `code/02_howMany1s.R`.
The impact of the connectivity under different thresholds is determined at `code/03_quantify_connectedness_1vs2.R`.

### LMHC calling

We use thresholds of 20 cells and a mean of 1.3 eUMIs with the reference mismatch to identify LMHCs. 

The annotation of the LMHCs occurs in `code/04_meta_profile.R`. Pileups of individual variants
can be generated using `04a_position_bias_examples.R`.

### Tree construction and MRCA analyses

Our attempts at reconstructing the ReDeeM phylogenetic trees with different inclusion/exclusion criteria
for Figure 2 are contained here `code/20_split_tree_make.R`. We quantify the degree of overlap 
by running a most recent common ancestor (MRCA) analysis here: `code/21_mrca_distance_1vs2p.R`. Note, 
this script takes a long time to run but can be modified with fewer iterations in the loop. The
visualization of the tree pre/post is available here: `code/22_viz_tree.R`.

### KS testing for position bias 

The script for computing this metric per variant is contained here: `code/06_biased_rate-split.R`. 

## Edge accumulations

The enrichment of transversions at the edge of molecules as a function of eUMI depth is found at
`code/08_visualize_edge_classification.R`  and the overall mismatches as a function of mutation type
are available at `code/09_what_nucleotide_changes.R`. 

## MQuad

We adapted the [MQuad](https://www.nature.com/articles/s41467-022-28845-0) workflow for 
compatibility with the ReDeeM data structure. Code to generate the input AD and DP matrices
is available `code/11_export_mquad.R` and for visualizing the MQuad results is `code/12_viz_mquad_pb.R`.


### Mouse and extended donor analyses

The author processed the murine and extended human donors with a distinct set of ReDeeM parameters
(total instead of sensitive), so these analyses proceed differently. The processed data was 
curated from GEO (as it was not available on FigShare). The relevant GEO accession numbers
are contained in the folder names: mouse (`mouse-data_GSE259284`) and additional hashed human donors
 (`Additional-donors_GSE261078`).

## Public data

To assess whether the edge mismatch accumulation occurs in stanadrd ATAC-seq, public
bulk ATAC and single-cell multiome analyses occur in `10x-snATAC-multiome`. Further, we
quantify the edge mismatches from the ATAC library upstream of the ReDeeM enrichment 
in `code/05a_variant_position_atac.R`.



## Questions?

[Email Caleb](mailto:lareauc@mskcc.org)

<br>
