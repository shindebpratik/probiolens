# 🧬 ProbioLens: Genome-Based Prediction of Probiotic Potential and Safety

ProbioLens is a genome-based prediction toolkit designed to evaluate **probiotic potential** and **safety potential** from bacterial genome sequences. It uses pretrained 8-mer Word2Vec embeddings and trained machine-learning classifiers to generate probability-based predictions from complete genomes and metagenome-assembled genomes (MAGs).

The toolkit provides a command-line interface for analysing individual genomes or processing multiple MAG FASTA files in a directory. Prediction models and their associated word-vector representations are bundled with the Python package for convenient installation and execution.

---

## 🔍 Overview

ProbioLens provides two complementary prediction tasks:

* **ProbioLens-P - Probiotic Potential:** Predicts whether a genome is classified as having probiotic potential.
* **ProbioLens-S - Safety Potential:** Predicts whether a genome is classified as having safety potential.

Both prediction tasks use genome-derived 8-mer embeddings and their corresponding trained classifiers. The selected classifier produces a positive-class probability, which is compared with a user-defined threshold to determine the final prediction label.

ProbioLens supports two input modes:

1. **Complete Genome Prediction:** Analyse a FASTA file representing a complete genome.
2. **MAG Prediction:** Analyse an individual MAG or multiple MAG FASTA files from a directory, with optional contig-length filtering.

The package is designed for command-line use and produces structured CSV output for downstream analysis.

---

## 🧠 Core Functionalities

* 🧬 **Genome-based prediction:** Evaluate probiotic potential and safety potential from nucleotide sequences.
* 🔤 **8-mer Word2Vec embeddings:** Represent genomes using occurrence-weighted mean embeddings of valid 8-mers.
* 🤖 **Machine-learning classification:** Apply the trained classifier corresponding to the selected prediction task.
* 🧫 **Complete-genome analysis:** Process a single genome FASTA file.
* 🧩 **MAG analysis:** Process individual MAGs or batches of MAG FASTA files.
* 📏 **Contig-length filtering:** Exclude short contigs in MAG mode using a configurable length cutoff.
* 🎚️ **Adjustable classification threshold:** Modify the probability threshold used to assign prediction labels.
* 📊 **CSV output:** Save prediction probabilities, labels, contig statistics, and processing status.
* 📦 **Portable Python package:** Install ProbioLens using pip, with the prediction models included in the package.

---


## ⚙️ Installation

### Option 1 - Install using pip

The recommended installation method for routine prediction is pip.

```bash
pip install probiolens
```

Check that the command-line interface is available:

```bash
probiolens --help
```

### Option 2 - Install in a dedicated Python environment

Using an isolated Python environment can help avoid dependency conflicts with other bioinformatics tools.

```bash
python3 -m venv probiolens-env
source probiolens-env/bin/activate
```

Install ProbioLens:

```bash
python -m pip install --upgrade pip
python -m pip install probiolens
```

Verify the installation:

```bash
probiolens --help
```

On Windows, activate the environment using:

```powershell
probiolens-env\Scripts\activate
```

Then install ProbioLens using the same pip command.

### Option 3 - Install from GitHub using Conda

This option allows you to run ProbioLens directly from the cloned GitHub repository using the standalone Python script and bundled model files.

**1. Clone the repository**

Replace `<YOUR_GITHUB_REPOSITORY_URL>` with the URL of your ProbioLens repository.

```bash
git clone <YOUR_GITHUB_REPOSITORY_URL>
cd ProbioLens
```

**2. Create the Conda environment**

Make sure Conda or Miniforge is installed.

```bash
conda env create -f probiolens.yml
conda activate probiolens
```

**3. Verify the installation**

Display the command-line help:

```bash
python probiolens.py --help
```

Ensure that the `model/` directory contains all four pretrained model files before running a prediction.

**4. Run a test prediction**

Use the included `test.fasta` file to check that the script and model files are accessible. For example, if the script accepts the arguments shown below:

```bash
python probiolens.py \
    --task P \
    --mode complete \
    --input test.fasta \
    --output test_predictions.csv
```

Check the script's help output for the exact argument names supported by your version.

**Note:** The saved classifiers use scikit-learn serialization. Compatibility with the dependency versions used during model training is important. If a model-loading or version-compatibility error occurs, report the error along with your Python, scikit-learn, and gensim versions.


---

## 🧭 Usage Overview

ProbioLens uses a unified command-line interface:

```bash
probiolens -i <input> -o <output.csv> -m <mode> -p <prediction> [options]
```

### Required arguments

| Argument | Description |
|---|---|
| `-i`, `--input` | Input FASTA file or directory containing MAG FASTA files |
| `-o`, `--output` | Output CSV file |
| `-m`, `--mode` | Input mode: `complete` or `mag` |
| `-p`, `--prediction` | Prediction task: `P` or `S` |

### Optional arguments

| Argument | Description | Default |
|---|---|---|
| `-t`, `--threshold` | Probability threshold for the positive prediction | `0.5` |
| `--min-contig` | Contig-length cutoff in base pairs for MAG mode | `5000` |
| `-h`, `--help` | Display the help message | - |

The probability threshold must be between 0 and 1. In MAG mode, only contigs **strictly longer than** the selected cutoff are retained. The default cutoff of 5,000 bp therefore excludes contigs of exactly 5,000 bp.

---

## 🧩 Detailed Usage Examples

### 1. ProbioLens-P: Probiotic Potential Prediction

Predict probiotic potential from a complete genome:

```bash
probiolens \
    -i genome.fasta \
    -o probiotic_predictions.csv \
    -m complete \
    -p P
```

The output includes the positive-class probability and the corresponding prediction label:

* `Probiotic`
* `Non-Probiotic`

The classification label is assigned by comparing the model probability with the specified threshold.

### 2. ProbioLens-S: Safety Potential Prediction

Predict safety potential from a complete genome:

```bash
probiolens \
    -i genome.fasta \
    -o safety_predictions.csv \
    -m complete \
    -p S
```

The output labels are:

* `Safe`
* `Not Safe`

These labels represent the classifier's prediction of safety potential. They should not be interpreted as definitive evidence of clinical safety or pathogenicity.

### 3. Predict Probiotic Potential from a MAG

Analyse an individual MAG FASTA file:

```bash
probiolens \
    -i bin.fasta \
    -o mag_probiotic_predictions.csv \
    -m mag \
    -p P
```

By default, only contigs longer than 5,000 bp contribute to the genome representation.

### 4. Predict Safety Potential from a MAG

```bash
probiolens \
    -i bin.fasta \
    -o mag_safety_predictions.csv \
    -m mag \
    -p S
```

This applies the safety-potential classifier to the retained contigs and saves the resulting prediction in CSV format.

### 5. Batch Prediction for Multiple MAGs

ProbioLens can process FASTA files contained in a directory:

```bash
probiolens \
    -i /path/to/MAGs/ \
    -o mag_predictions.csv \
    -m mag \
    -p P
```

The input directory should contain supported FASTA files, including `.fa`, `.fasta`, `.fna`, or `.fas` files. Files are processed individually, and their results are combined into a single output CSV.

For safety-potential prediction, use `-p S` instead of `-p P`.

### 6. Change the Probability Threshold

The default classification threshold is 0.5. To use a threshold of 0.7:

```bash
probiolens \
    -i genome.fasta \
    -o predictions_threshold_07.csv \
    -m complete \
    -p P \
    -t 0.7
```

Increasing the threshold makes a positive prediction more stringent. The appropriate threshold should be selected based on the intended application and relevant validation results.

### 7. Change the MAG Contig-Length Cutoff

To retain only contigs strictly longer than 10,000 bp:

```bash
probiolens \
    -i bin.fasta \
    -o mag_predictions.csv \
    -m mag \
    -p P \
    --min-contig 10000
```

This option changes the contigs used to calculate the genome embedding. Different cutoffs can change predictions, particularly for incomplete or fragmented MAGs.

### 8. Display Help

```bash
probiolens --help
```

This displays the command-line options and examples.

---

## 🧬 Prediction Workflow

ProbioLens follows a sequence-based prediction workflow:

1. **Input processing:** Read the supplied complete-genome or MAG FASTA file.
2. **Contig filtering:** In MAG mode, discard contigs that do not exceed the selected length cutoff.
3. **8-mer extraction:** Identify valid nucleotide 8-mers composed of A, C, G, and T.
4. **Embedding generation:** Retrieve the corresponding pretrained Word2Vec vectors and calculate an occurrence-weighted mean embedding.
5. **Embedding centering:** Centre the mean embedding using the vocabulary-mean vector representation.
6. **Model inference:** Apply the classifier associated with ProbioLens-P or ProbioLens-S.
7. **Threshold-based classification:** Compare the positive-class probability with the selected threshold.
8. **Result reporting:** Save predictions, probabilities, contig statistics, and processing status to a CSV file.

The package contains separate classifier and word-vector files for each prediction task. The selected model pair is loaded according to the `-p` argument.

---

## 📊 Output Description

ProbioLens writes a CSV file containing the following columns:

| Column | Description |
|---|---|
| `genome` | Input FASTA filename |
| `prediction_task` | Probiotic potential or safety potential |
| `mode` | Input mode (`complete` or `mag`) |
| `probability_positive` | Predicted probability for the positive class |
| `threshold` | Threshold used to assign the label |
| `predicted_label` | Final predicted class |
| `contigs_total` | Total number of contigs in the input |
| `contigs_used` | Number of contigs retained for analysis |
| `contigs_filtered` | Number of contigs excluded by the MAG length filter |
| `status` | Processing status (`success` or `failed`) |
| `error` | Error description if processing fails |

### Example output

The following values are illustrative examples of the output format, not actual model results.

```csv
genome,prediction_task,mode,probability_positive,threshold,predicted_label,contigs_total,contigs_used,contigs_filtered,status,error
example_genome.fasta,Probiotic potential,complete,0.91,0.5,Probiotic,120,120,0,success,
example_bin.fasta,Safety potential,mag,0.84,0.5,Safe,350,142,208,success,
```

For a failed prediction, the output includes a failure status and an error message to help identify the problem.

---

## ⚠️ Important Considerations

* **Prediction versus confirmation:** A positive prediction does not establish that a strain is an effective probiotic or clinically safe. Experimental validation and appropriate safety assessment remain necessary.
* **Safety interpretation:** The `Safe` and `Not Safe` labels are model outputs, not definitive regulatory or clinical safety determinations.
* **MAG completeness:** Fragmentation and incomplete genome recovery can affect the representation generated from a MAG. The contig-length filter does not itself measure genome completeness.
* **Probability interpretation:** The reported probability is a model output. It should not automatically be interpreted as a calibrated clinical probability.
* **Threshold selection:** A threshold of 0.5 is the default, not necessarily the optimal operating point for every application.
* **Model compatibility:** Saved classifiers may depend on the versions of Python libraries used during model development.
* **Input requirements:** Input sequences must be supplied in FASTA format. Only valid A/C/G/T 8-mers contribute to the embedding.

ProbioLens is intended as a computational screening tool to support research and prioritisation. Predictions should be interpreted alongside genomic context, relevant biological evidence, and experimental findings.

---

## 📖 Citation

If you use ProbioLens in a publication, thesis, or other academic work, please cite the associated ProbioLens publication when available.

**ProbioLens publication:**


Until a publication is available, you may cite the software release using its version number and a stable repository or archival record, if available.

---

## 📜 License

ProbioLens is distributed under the license specified in the accompanying `LICENSE` file and package metadata.

Please review the license terms before redistributing or modifying the software. The redistribution terms must also cover the bundled classifiers and word-vector files.

---

## 👨‍💻 Support and Contributions

For bug reports, installation issues, or feature requests, please use the project's issue tracker or contact the corresponding maintainer.

* **GitHub:** https://github.com/shindebpratik/probiolens
* **Contact:** Pratik Shinde - pratiks@iiitd.ac.in

When reporting an issue, include the ProbioLens version, Python version, relevant dependency versions, the command used, and the complete error message. Avoid sharing confidential genome data or other sensitive information in public issue reports.
