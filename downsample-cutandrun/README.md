# CUT&RUN Downsampling Pipeline

## 1  Overview

This workflow implements downsampling in CUT&RUN samples:

- **Input:** paired‑end BAMs + a metadata CSV (conditions, IP vs IgG, replicates)
- **Steps:** down‑sample to matched read depth, merge replicates, peak‑call with **MACS2** (narrow/broad), generate **bigWig** tracks (raw & RPGC), and collect QC (flagstat, benchmarks, logs)
- **Output:**
  - `results/peaks/` — `*_peaks.{narrowPeak|broadPeak}` per replicate and pooled
  - `results/bigwig/` — RAW & RPGC bigWigs for individual and merged BAMs
  - `results/merged/` — replicate‑merged BAM + BAI
  - `results/flagstats/`, `results/logs/`, `results/benchmarks/` 

---

## 2  Quick Start

```bash
# 1) clone your fork or the original repo
$ git clone git@github.com:rlbc/rlbc_compbio.git
$ cd rlbc_compbio/downsample-cutandrun

# 2) (optional) create a base environment with snakemake
$ mamba create -n snakemake -c conda-forge -c bioconda snakemake
$ conda activate snakemake

# 3) dry‑run the DAG
$ snakemake -np

# 4) execute with Conda integration on 4 cores
$ snakemake --use-conda -j 4
```

---

## 3  Project Structure

```
CUTRUN_pipeline/
├── Snakefile
├── config.yaml              # workflow parameters
├── samples_example.csv       # template metadata file
├── envs/                     # conda YAMLs (auto‑generated/pinned)
└── results/                  # created at runtime
```

---

## 4  Preparing the Metadata CSV

Create a `samples_<mark>.csv` with **exactly** these columns (order unimportant):

```csv
condition,type,replicate,sample,bam
WT,IP,1,WT_IP_R1,/path/to/WT_IP_R1.bam
WT,IP,2,WT_IP_R2,/path/to/WT_IP_R2.bam
WT,IgG,1,WT_IgG_R1,/path/to/WT_IgG_R1.bam
KO,IP,1,KO_IP_R1,/data/KO_IP_R1.bam
KO,IgG,1,KO_IgG_R1,/data/KO_IgG_R1.bam
...
```

Example:

```csv
K562,IP,1,H3K27ac_K562_R1,K562_H3K27ac_ENCFF423WAK.bam
K562,IP,2,H3K27ac_K562_R2,K562_H3K27ac_ENCFF867JTP.bam
K562,IgG,1,IgG_K562_R1,K562_INPUT_ENCFF139VMV.bam
K562,IgG,2,IgG_K562_R2,K562_INPUT_ENCFF423WGF.bam
```

**Field definitions**

| Column      | Allowed values / format                             |
| ----------- | --------------------------------------------------- |
| `condition` | Biological condition (e.g. `WT`, `KO`, `cell_line`) |
| `type`      | `IP` or `IgG` (case‑sensitive)                      |
| `replicate` | Integer **≥1**; used to match IP ↔ IgG              |
| `sample`    | Unique identifier, becomes file prefixes            |
| `bam`       | Absolute path to **coordinate‑sorted** BAM          |

---

## 5  Outputs

| Path                  | Contents                                                  |
| --------------------- | --------------------------------------------------------- |
| `results/peaks/`      | `sample_peaks.narrowPeak`, `cond_pooled_peaks.narrowPeak` |
| `results/bigwig/`     | `*.original.bw` + `*.RPGC.bw`                             |
| `results/merged/`     | `cond_IP.bam`, `cond_IgG.bam` + `.bai`                    |
| `results/flagstats/`  | Raw read counts per BAM                                   |
| `results/benchmarks/` | Runtime+RAM TSVs from Snakemake                           |
| `results/logs/`       | Stdout/err of each rule                                   |

