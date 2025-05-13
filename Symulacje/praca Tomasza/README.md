# Time Series Detrending and Spectral Analysis Pipeline

This repository provides an R-based workflow to:
- Simulate noisy time series with oscillatory and trend components
- Compare multiple detrending methods (linear, rolling, low/high-pass, FIR, polynomial, etc.)
- Apply Hanning window tapering and zero-padding
- Compute Fourier transforms and aggregate spectral power
- Perform permutation-based confidence intervals on spectral estimates
- Run the same analyses on real experimental data

## Repository Structure
- **functions.r**, **functions2.r**: Core utilities for data generation, preprocessing, detrending, padding, and FFT
- **simulations.r**, **symart2.r**: Example scripts to generate synthetic data, apply all detrending methods, and aggregate FFT results
- **analiza2.r**: Real‐data analysis script that reads CSVs, preprocesses, detrends, computes spectra, runs permutation tests, and plots results
- **orig.png**, **with_mov_Avg.png**: Illustrative figures showing raw vs. moving‐average detrending
- **podsumowanie.pdf**, **symulacje_i_analiza.odt**: Written reports and summaries of simulation and analysis results

## Requirements
- R (>= 3.5)
- CRAN packages: `imputeTS`, `plyr`, `reshape2`, `signal`, `zoo`, `ggplot2`, `dplyr`, `tidyverse`, `fs`, `progressr`, `future`, `furrr`
- GitHub package: `detrendeR` by Bartek Kroczek
  ```r
  devtools::install_github("bartekkroczek/detrendeR", build_vignettes = TRUE)
  ```

## Installation
1. Clone the repository:
   ```bash
   git clone <this-repo-url>
   cd <repo-directory>
   ```
2. Install required R packages (via CRAN and GitHub as above).

## Usage
1. Open R or RStudio in the project directory.
2. Source the core functions:
   ```r
   source("functions2.r")   # or source("functions.r") for the original version
   ```
3. Run synthetic simulations:
   ```r
   source("simulations.r")  # baseline pipeline
   # or
   source("symart2.r")      # alternate/demo version with additional methods
   ```
   Results will be aggregated FFT magnitudes per detrending method and frequency.
4. Analyse real experiment data:
   - Place your `.csv` files in `../../Experiment_1/Data/` (or adjust the path in analiza2.r).
   - Source and run:
     ```r
     source("analiza2.r")
     ```
   - Use the provided functions `fftwhole()`, `makedistr()`, `plotci()`, and `pvals()` to compute and visualize spectra with confidence intervals.

## Examples
- In `analiza2.r`, after computing the permutation matrix `ci`, you can plot confidence bands around the empirical spectrum:
  ```r
  est <- fftwhole(d1a)
  plotci(ci, est)
  ```

## License
Please see the LICENSE file or contact the author for licensing terms.