# CONCISE

A statistically principled, robust, and efficient method for spatial gene-gene co-expression and ligand-receptor interaction inference

![CONCISE\_overview](overview/overview.jpg)

CONCISE is a statistically principled approach for inferring spatially constrained gene co-expression and cell-cell communication (CCC) from spatial transcriptomics (ST) data. It jointly models key features of ST data, including spatial autocorrelation, the count-based nature of the data, heterogeneous total molecular counts and measurement errors, all of which can confound the spatial analyses and generate spurious discoveries if ignored. By explicitly accounting for these factors, CONCISE enables reliable inference with well-calibrated false-positive control and high detection power.

Built upon a flexible model with efficient moment-based parameter estimation and analytically derived hypothesis testing, CONCISE avoids computationally intensive permutation procedures and restrictive distributional assumptions. As a result, it delivers statistically rigorous inference while maintaining high robustness and computational efficiency.

The package supports both spatial ligand–receptor interaction (LRI) inference and gene-gene network construction. It facilitates the discovery of biologically meaningful signaling interactions and spatial gene modules across diverse biological systems and contexts. Applications of CONCISE have uncovered distinct communication patterns between inflammation-associated fibroblasts and normal fibroblasts during intestinal inflammation, and have elucidated complex tumor-immune and tumor-stromal signaling networks within the tumor microenvironment.

