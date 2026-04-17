# Ensaya:An ensemble age model for prediction of chronological age in adolescents and young adults

This repository provides an R-based workflow to run **Ensaya**, an ensamble epigenetic age estimator.

**Ensaya** combines three established epigenetic clocks:

* cAge (Bernabeu E, McCartney DL, Gadd DA, Hillary RF, Lu AT, Murphy L, Wrobel N, Campbell A, Harris SE, Liewald D, Hayward C, Sudlow C, Cox SR, Evans KL, Horvath S, McIntosh AM, Robinson MR, Vallejos CA, Marioni RE. Refining epigenetic prediction of chronological and biological age. Genome Med. 2023 Feb 28;15(1):12. doi: 10.1186/s13073-023-01161-y. PMID: 36855161; PMCID: PMC9976489.)
* Garma young (Garma LD, Quintela-Fandino M. Applicability of epigenetic age models to next-generation methylation arrays. Genome Med. 2024 Oct 7;16(1):116. doi: 10.1186/s13073-024-01387-4. Erratum in: Genome Med. 2024 Dec 3;16(1):142. doi: 10.1186/s13073-024-01402-8. PMID: 39375688; PMCID: PMC11460231.)
* PAYA (Aanes, H., Bleka, Ø., Dahlberg, P.S. et al. A new blood based epigenetic age predictor for adolescents and young adults. Sci Rep 13, 2303 (2023). https://doi.org/10.1038/s41598-023-29381-7)

The script runs all three clocks internally and returns a single output table with individual predictions and the age Ensaya estimate.

---

## Input format

The input file should contain:

* **Rows:** CpG sites
* **Columns:** Samples
* **Values:** Beta values

Supported formats:

* `.rds`
* `.csv`
* `.tsv`

---

## Usage

Only one function is needed:

```r
run_clocks(input_file, out_file)
```

### Example

```r
source("ensaya.R")

run_clocks(
  input_file = "data/example_data.csv",
  out_file = "results/ensaya_results.csv"
)
```

---

## Output

The output is a `.csv` file:

| Sample | cAge | Garma | PAYA | Ensaya |
| ------ | ---- | ----- | ---- | ------ |

* **Ensaya**: Mean of the three clocks (primary output)
* **cAge, Garma, PAYA**: Individual clock predictions 

## Requirements

```r
tidyverse
```

---

