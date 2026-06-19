# CONCISE

*Spatial co-expression and cell-cell communication inference from spatially resolved transcriptomics with CONCISE*

A statistically principled, robust, and fast method for spatial ligand-receptor interaction inference and gene-gene network construction.

![CONCISE\_overview](overview/overview.jpg)

CONCISE jointly models key features of ST count data, including spatial autocorrelation, heterogeneous total molecular counts and measurement errors, all of which can lead to spurious results if ignored. By explicitly accounting for these sources of confounding, CONCISE enables reliable inference with well-calibrated false-positive control and higher detection power. With efficient moment-based parameter estimation and analytically derived hypothesis testing, CONCISE avoids computationally intensive permutation procedures and restrictive distributional assumptions, achieving rigorous inference while maintaining high robustness and computational efficiency.

Extensive real-data permutation and biologically motivated negative-control studies have demonstrated its reliablility. Applications of CONCISE have uncovered distinct communication patterns between inflammation-associated fibroblasts and normal fibroblasts during intestinal inflammation, and have elucidated complex tumor-immune and tumor-stromal signaling networks within the tumor microenvironment.

