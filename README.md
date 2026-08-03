# README

## Description

Evaluation experiment of Digital Cytometry (or transcriptomic deconvolution) methods using:

-   Gene expression profiles (GEPs) of six pure cell lines as reference and ground truth.
-   Bulk data from our own mixtures generated *in vitro* mixtures using the six cell lines in different proportions (ground truth).

![Schematic visualization of the experimental workflow for digital cytometry evaluation. This figure was created with Biorender.com](imgs/Fig1_experiment_workflow.png)

Publication: *currently under review; reference available once it is published#*

**Evaluated tools:**

| Name | Method | Estimation GEPs | Reference |
|------------------|------------------|------------------|------------------|
| **CIBERSORTx** | Regression | Yes | <https://doi.org/10.1038/s41587-019-0114-2> |
| **DWLS** | Regression | No | <https://doi.org/10.1038/s41467-019-10802-z> |
| **BayesPrism** | Bayesian | Yes | <https://doi.org/10.1038/s43018-022-00356-3> |
| **InstaPrism** | Fast Implementation of BayesPrism | Yes | <https://doi.org/10.1093/bioinformatics/btae440> |
| **BLADE** | Bayesian | Yes | <https://doi.org/10.1038/s41467-021-26328-2> |

## Data

Available in GEO with accession number: [GSE324026](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324026)

**Cell lines:**

-   UPCI-SCC-040 (HNSCC tumor cells)
-   UPCI-SCC-131 (HNSCC tumor cells)
-   OKF6 (immortalized keratinocytes)
-   Dermal primary fibroblasts
-   Jurkat T-cells (immune cells)
-   Natural Killer (NK-92) cells (immune cells)

**Mixtures composition:**

see data/InVitroMixtures.tsv

## Code

Main functions to perform evaluation experiment:

-   Running deconvolution: see code/01_RunningDigitalCytometry

-   Modifying reference profiles for Use Case 3: see code/02_NoiseAdditionReferenc

-   Output processing: see code/03_OutputProcessing (include ??)

-   Distances between output and ground truth: see code/04_DistancesFromGroundTruth

-   Recovery of the biological signal: see code/05_DownstreamAnalysis (might be modified for re-submission)

## Authors

-   TMO, IDC, Helmholtz Munich.
