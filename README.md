# BIOME Saguenay Sediment Analysis (2022-2024)

This project utilizes the **G2Sd** R package (Fournier et al., 2014) for the statistical analysis of unconsolidated sediments. It provides a specialized pipeline to process, validate, and standardize granulometric data collected within the Saguenay Fjord.

## Overview

The analytical core of this repository is built upon a customized implementation of `G2Sd`. To ensure the highest level of scientific rigor and reproducibility, several internal functions have been refined to strictly adhere to the mathematical standards of Blott & Pye (2001) and Folk & Ward (1957).

## Core Features

* **Mathematical Corrections**: Custom R functions that override standard `G2Sd` internal methods to ensure strict mathematical compliance with original sedimentological equations.
* **Dual-Scale Statistics**: Automated calculation of textural parameters (Mean, Sorting, Skewness, Kurtosis) in both metric (µm) and logarithmic (Phi) units.

---

## Function Reference

To ensure reproducibility and precision, this project utilizes several modular functions located in the `R/` directory:

* **`calc_percentile.R`**: 
    * Contains `.percentile()`.
    * Uses linear interpolation between sieve/laser mesh sizes to calculate exact percentiles (D5 to D95). 
    * Includes a boundary correction for the 0 µm fraction to allow for Phi scale conversion.

* **`calc_stats.R`**:
    * **`.moment.arith()`**: Computes the arithmetic Method of Moments statistics based on true midpoints between size classes.
    * **`.moment.geom()`**: Computes the geometric Method of Moments. It utilizes natural logarithms (ln) to ensure higher precision than standard base-2 implementations.
    * **`.fowa.stat()`**: Implements the Folk and Ward (1957) graphical method, producing both numerical indices and qualitative textural descriptions (e.g., "Poorly Sorted", "Fine Skewed").

* **`sed_classes.R`**:
    * **`sed_classes()`**: Categorizes sediment samples into 17 detailed grain-size classes (from boulder to clay) following the GRADISTAT scale (Blott & Pye, 2001).

* **`process_sediment_data.R`**:
    * **`process_sediment_data()`**: The primary orchestrator. It runs the analysis pipeline, manages data merging, and renames all variables to match the project's official data dictionary.

---

## Repository Structure

* `R/`: Modular R scripts and custom functions.
* `data/raw/`: Original laboratory exports (Excel format).
* `data/processed/`: Final, cleaned, and validated datasets in CSV format.
* `gradistat.R`: The main execution script used to generate the published dataset.
