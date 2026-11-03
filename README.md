# 🧬 ProbioLens: Genome-Based Prediction of Probiotic Potential

ProbioLens is a command-line toolkit for predicting **probiotic potential** from bacterial genome sequences using a codon-usage-based machine-learning model. It supports assembled isolate genomes, individual MAGs, batches of MAG FASTA files, and read-based workflows that assemble sequencing reads with MEGAHIT.

The current implementation uses **one probiotic-potential classifier (ProbioLens)**. It does not provide the former embedding-based prediction workflow or a separate ProbioLens-S safety model.

---

## 🔍 Overview

ProbioLens provides one prediction task:

- **Probiotic Potential:** predicts the probability that an input genome has probiotic potential and assigns `Probiotic` or `Non-Probiotic` using a fixed probability threshold of **0.8**.

The model uses **61 codon-usage features** extracted from coding sequences predicted by Prodigal:

- 59 Relative Synonymous Codon Usage (RSCU) features
- Effective Number of Codons (ENC)
- Codon Adaptation Index (CAI)

The classifier is loaded from `model/probiolens_p_classifier.pkl`, located beside the script.

---

## ✨ Features

- **Assembled genome prediction:** predict from one genome FASTA.
- **Batch MAG prediction:** process every supported FASTA file in a directory, generating one result per file.
- **Read-to-isolate workflow:** assemble reads using MEGAHIT and predict directly from the resulting contigs.
- **Read-to-MAG workflow:** assemble reads, map reads to contigs, bin with MetaBAT2, assess bin quality with CheckM2, and predict each passing MAG separately.
- **Configurable MAG quality criteria:** default minimum completeness is 20% and maximum contamination is 5%.
- **CSV reporting:** prediction probability, label, input/assembly mode, quality metrics (where applicable), status, and errors are recorded.

---

## ⚙️ Requirements

### Python dependencies

The script uses:

- Python 3
- NumPy
- pandas
- scikit-learn

Install the Python dependencies in your environment, for example:

```bash
python -m pip install numpy pandas scikit-learn
```

### External bioinformatics tools

Install the tools required for the workflow you intend to run and make their executables available on `PATH`.

| Tool | Required for |
|---|---|
| **Prodigal** (`prodigal`) | All prediction workflows; CDS prediction |
| **MEGAHIT** (`megahit`) | Read-based assembly (`--assembly`) |
| **Bowtie2** (`bowtie2`, `bowtie2-build`) | Read mapping in `--meta` workflow |
| **samtools** (`samtools`) | BAM sorting and indexing in `--meta` workflow |
| **MetaBAT2** (`metabat2`) | MAG binning in `--meta` workflow |
| **CheckM2** (`checkm2`) | MAG quality assessment in `--meta` workflow |

CheckM2 must already be configured to locate its database (for example, through its `CHECKM2DB` setup). The pipeline uses **8 threads by default** for external tools.

---

## 📦 Repository setup

Clone the repository and enter it:

```bash
git clone https://github.com/shindebpratik/probiolens.git
cd probiolens
```

Ensure the repository contains the script and model file in this layout:

```text
probiolens/
├── probiolens.py
└── model/
    └── probiolens_p_classifier.pkl
```

The model file must be available at the exact relative location shown above (the script resolves it relative to its own location).

Run the script directly:

```bash
python probiolens.py --help
```

> Use the actual script filename present in your checkout if it differs from `probiolens.py`.

---

## 🧭 Input modes and command-line options

ProbioLens has two main input routes.

### A. Already assembled FASTA

Use `-i` / `--input` for:

- one assembled genome FASTA file; or
- a directory containing multiple MAG FASTA files.

Supported FASTA extensions: `.fa`, `.fasta`, `.fna`, and `.fas`.

### B. Sequencing reads requiring assembly

Use `--assembly` (alias `--assemble`) with `-i1` / `--input1` and optionally `-i2` / `--input2`.

With assembly, select exactly one workflow:

- `--isolate`: predict the MEGAHIT assembly directly (no binning or CheckM2).
- `--meta`: map reads, bin contigs, assess MAG quality, and predict each passing MAG independently.

### Main options

| Option | Description |
|---|---|
| `-i`, `--input` | Assembled genome FASTA or directory of MAG FASTA files |
| `-i1`, `--input1` | First read file for assembly; single-end if `-i2` is omitted |
| `-i2`, `--input2` | Second read file for paired-end sequencing; requires `-i1` and `--assembly` |
| `-o`, `--output` | **Required.** Output CSV path |
| `--assembly`, `--assemble` | Assemble sequencing reads with MEGAHIT |
| `--isolate` | With assembly, predict the assembly directly |
| `--meta` | With assembly, run binning/quality assessment and predict passing MAGs |
| `--complete` | Minimum CheckM2 completeness percentage; default `20` |
| `--contam` | Maximum CheckM2 contamination percentage; default `5` |
| `-h`, `--help` | Show help |

The `--isolate` and `--meta` options are mutually exclusive and are only used with `--assembly`. For read input, `-i1` alone specifies single-end reads; `-i1` together with `-i2` specifies paired-end reads.

---

## 🚀 Usage examples

### 1. Predict from one assembled genome

```bash
python probiolens.py \
    -i genome.fasta \
    -o genome_prediction.csv
```

### 2. Predict from one MAG

```bash
python probiolens.py \
    -i bin.fasta \
    -o mag_prediction.csv
```

### 3. Batch prediction for a directory of MAGs

```bash
python probiolens.py \
    -i /path/to/mags/ \
    -o mag_predictions.csv
```

Each FASTA file is processed independently and receives its own prediction row.

### 4. Assemble single-end isolate reads and predict

```bash
python probiolens.py \
    -i1 reads.fastq.gz \
    --assembly \
    --isolate \
    -o isolate_prediction.csv
```

### 5. Assemble paired-end isolate reads and predict

```bash
python probiolens.py \
    -i1 reads_R1.fastq.gz \
    -i2 reads_R2.fastq.gz \
    --assembly \
    --isolate \
    -o isolate_prediction.csv
```

### 6. Assemble single-end metagenomic reads, bin, and predict MAGs

```bash
python probiolens.py \
    -i1 metagenomic_reads.fastq.gz \
    --assembly \
    --meta \
    -o mag_predictions.csv
```

### 7. Assemble paired-end metagenomic reads with custom quality cutoffs

```bash
python probiolens.py \
    -i1 reads_R1.fastq.gz \
    -i2 reads_R2.fastq.gz \
    --assembly \
    --meta \
    --complete 50 \
    --contam 3 \
    -o mag_predictions.csv
```

### 8. Display command-line help

```bash
python probiolens.py --help
```

---

## 🧬 Prediction workflow

### Assembled genome or MAG

1. Read the input FASTA.
2. Run Prodigal to predict coding sequences and produce a nucleotide CDS FASTA (`.ffn`).
3. Retain usable CDSs (length at least 60 nt, divisible by 3, and without internal stop codons).
4. Require at least **50 usable CDSs** for codon-feature extraction.
5. Calculate 59 RSCU features, ENC, and CAI.
6. Apply the 61-feature probiotic classifier.
7. Report the probiotic probability and threshold-based label.

### Read-based isolate workflow (`--assembly --isolate`)

1. Assemble reads using MEGAHIT.
2. Run Prodigal and extract codon-usage features from the assembly.
3. Predict probiotic potential for the assembly.

### Read-based metagenome workflow (`--assembly --meta`)

1. Assemble reads using MEGAHIT.
2. Build a Bowtie2 index and map reads to the assembly.
3. Sort and index the BAM file with samtools.
4. Calculate contig depths using `jgi_summarize_bam_contig_depths`.
5. Bin contigs using MetaBAT2.
6. Assess bins using CheckM2.
7. Retain bins meeting the configured quality criteria (default: completeness ≥ 20%; contamination ≤ 5%).
8. Run an independent probiotic prediction for each passing MAG.
9. Include quality-filtered-out bins in the results with their filtering status.

---

## 📁 Output and intermediate files

The required `-o` argument specifies the prediction CSV. Pipeline intermediate files are automatically created beside it; no separate pipeline-directory option is needed.

For assembly workflows, the pipeline directory is:

```text
<output_stem>_pipeline/
├── assembly/
├── mapping/                 # --meta workflow
├── bins/                    # --meta workflow
├── checkm2/                 # --meta workflow
├── prodigal/
├── bin_quality_filter.csv   # --meta workflow
└── logs (*.stderr.log)
```

For non-assembly prediction, Prodigal output is stored in:

```text
<output_stem>_prodigal/
```

The script avoids overwriting non-empty pipeline output directories. Choose a new output filename or move previous intermediate results before rerunning.

### Prediction CSV fields

The CSV includes fields such as:

| Column | Description |
|---|---|
| `genome` | Genome, assembly, or bin name |
| `prediction_task` | `Probiotic potential` |
| `input_mode` | `complete` or `mag` |
| `assembly_mode` | Empty for direct FASTA input; `isolate` or `meta` for assembly workflows |
| `probability_probiotic` | Predicted probability for the probiotic class |
| `threshold` | Classification threshold (`0.8`) |
| `predicted_label` | `Probiotic`, `Non-Probiotic`, or a processing/filter status |
| `completeness_percent` | CheckM2 completeness for MAG workflow; otherwise blank/NA |
| `contamination_percent` | CheckM2 contamination for MAG workflow; otherwise blank/NA |
| `quality_filter_passed` | Whether the MAG passed CheckM2 criteria; `not_assessed` where not applicable |
| `contigs_total` | Number of FASTA sequences/contigs |
| `contigs_used` | Number of contigs recorded as used |
| `contigs_filtered` | Number of contigs filtered (zero for current direct codon-feature processing) |
| `status` | Processing status, e.g. `success`, `failed`, or `filtered_out_quality` |
| `error` | Error message, if applicable |

For the metagenome workflow, bins that fail CheckM2 criteria are included in the CSV as `filtered_out_quality` and are not passed to the classifier. A bin that passes quality filtering but fails prediction is recorded as `failed`.

---

## ⚠️ Important considerations

- **Prediction is not experimental confirmation:** a `Probiotic` result does not establish probiotic efficacy or clinical safety.
- **No safety classifier:** this version predicts probiotic potential only; it does not predict safety or replace safety assessment.
- **Fixed threshold:** the current script uses a probability threshold of **0.8**; the threshold is not exposed as a CLI option.
- **MAG quality is not a prediction guarantee:** completeness and contamination are used as bin filters in the `--meta` workflow, not as proof of biological quality.
- **Minimum CDS requirement:** genomes with fewer than 50 usable CDSs cannot be processed for codon features.
- **Feature/model compatibility:** the classifier must expect exactly 61 features in the same order (59 RSCU features, ENC, CAI).
- **Model loading:** only load the pickle model from a trusted source; serialized Python models can execute code during loading.
- **Probability interpretation:** the output is a model probability and should not automatically be interpreted as a calibrated clinical probability.
- **External tools:** read workflows require the relevant tools and databases to be installed/configured correctly.

ProbioLens is intended as a computational screening and prioritisation tool. Interpret results with biological context and experimental evidence.

---

## 📖 Citation

If you use ProbioLens in a publication, thesis, or other academic work, cite the associated ProbioLens publication when available.

**ProbioLens publication:** *To be added.*

Until a publication is available, cite the software release and repository URL where appropriate.

---

## 📜 License

ProbioLens is distributed under the license specified in the repository's `LICENSE` file, if present. Review the license and the redistribution terms for the bundled classifier before redistributing or modifying the software.

---

## 👨‍💻 Support and contributions

For bug reports, installation issues, or feature requests, please use the repository issue tracker.

- **GitHub:** https://github.com/shindebpratik/probiolens
- **Contact:** Pratik Shinde — pratiks@iiitd.ac.in

When reporting an issue, include the command used, operating system, Python and dependency versions, relevant external-tool versions, and the complete error message. Do not upload confidential genome data to public issues.
