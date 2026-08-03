# README

## Description

Evaluation experiment of Digital Cytometry (or transcriptomic deconvolution) methods using: 

-   Gene expression profiles (GEPs) of six pure cell lines as reference and ground truth.
-   Bulk data from our own mixtures generated *in vitro* mixtures using the six cell lines in different proportions (ground truth).

Publication: _currently under review; reference available once it is published_

**Evaluated methods:**

-   CIBERSORTx (regression-based): Newman, A.M., Steen, C.B., Liu, C.L. et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. Nat Biotechnol 37, 773–782 (2019). <https://doi.org/10.1038/s41587-019-0114-2>
-   DWLS (regression-based): Tsoucas, D., Dong, R., Chen, H. et al. Accurate estimation of cell-type composition from gene expression data. Nat Commun 10, 2975 (2019). <https://doi.org/10.1038/s41467-019-10802-z>
-   BayesPrism (Bayesian/Probabilistic): Chu, T., Wang, Z., Pe’er, D. et al. Cell type and gene expression deconvolution with BayesPrism enables Bayesian integrative analysis across bulk and single-cell RNA sequencing in oncology. Nat Cancer 3, 505–517 (2022). <https://doi.org/10.1038/s43018-022-00356-3>
-   InstaPrism (Fast implementation of BayesPrism): Mengying Hu, Maria Chikina, InstaPrism: an R package for fast implementation of BayesPrism, Bioinformatics, Volume 40, Issue 7 (2024). <https://doi.org/10.1093/bioinformatics/btae440>
-   (BLADE (Probabilistic)): Andrade Barbosa, B., van Asten, S.D., Oh, J.W. et al. Bayesian log-normal deconvolution for enhanced in silico microdissection of bulk gene expression data. Nat Commun 12, 6106 (2021). https://doi.org/10.1038/s41467-021-26328-2 

**Cell lines:**

-   UPCI-SCC-040 (HNSCC tumor cells)
-   UPCI-SCC-131 (HNSCC tumor cells)
-   OKF6 (immortalized keratinocytes)
-   Dermal primary fibroblasts
-   Jurkat T-cells (immune cells)
-   Natural Killer (NK-92) cells (immune cells)

**Mixtures composition:**

see data/InVitroMixtures.tsv

## Data

Available in GEO with accession number: GSE324026.

## Code

- Running deconvolution: see code/01_RunningDigitalCytometry 

- Modifying reference profiles for Use Case 3: see code/02_NoiseAdditionReferenc

- Output processing: see code/03_OutputProcessing (include ??)

- Distances between output and ground truth: see code/04_DistancesFromGroundTruth

- Recovery of the biological signal: see code/05_DownstreamAnalysis (might be modified for re-submission)

## Authors

-   TMO, IDC, Helmholtz Munich.
