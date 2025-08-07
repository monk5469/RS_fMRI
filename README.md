# DCM-PEB Scripts for fMRI Repetition Suppression Study

This repository contains analysis scripts used for Dynamic Causal Modelling (DCM) and Parametric Empirical Bayes (PEB) applied to a publicly available fMRI dataset. These scripts were adapted based on the methods described in:
   
    Lee, Sung-Mu, et al. "Effects of face repetition on ventral visual stream connectivity using dynamic causal modelling of fMRI data." NeuroImage 264 (2022): 119708.

## Key Modifications
Customized F-contrasts to identify condition-relevant activation

An optimized ROI selection approach based on group-level voxel-wise peak coverage, designed to improve inter-subject consistency

Design of B matrix input variables to model condition-specific modulatory effects

## Dataset
The dataset used in this study is publicly available from OpenNeuro: https://openneuro.org/datasets/ds000117/versions/1.1.0

## Usage
These scripts are intended for use with SPM12 in MATLAB. For detailed usage instructions, please refer to the original publication by Lee et al. (2022), along with the notes included in each script file.

