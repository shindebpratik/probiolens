#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ProbioLens: genome-based prediction of probiotic potential.

Input options:
  * -i/--input: an already assembled FASTA file or directory of MAG FASTAs.
  * --assembly --isolate: assemble reads with MEGAHIT and predict the assembly.
  * --assembly --meta: MEGAHIT -> read mapping -> MetaBAT2 -> CheckM2 ->
    quality filtering -> independent probiotic prediction for each passing MAG.

For metagenome-derived bins, default CheckM2 thresholds are completeness >= 20%
and contamination <= 5%; users can adjust them with --complete and --contam.
Tools including Prodigal are expected to be available on PATH. CheckM2 must already be configured
to find its database (for example through its CHECKM2DB setup).
Model files are loaded from ./model/ relative to this script.
"""

import argparse
import os
import pickle
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd


# ============================================================================
# CONFIGURATION
# ============================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
MODEL_FILES = {
    "classifier": SCRIPT_DIR / "model/probiolens_p_classifier.pkl",
}

CODON_FEATURE_COUNT = 61  # 59 RSCU + ENC + CAI
FASTA_EXTENSIONS = {".fa", ".fasta", ".fna", ".fas"}

# User-specified CheckM2 filtering rules. These are the ONLY bin-quality cutoffs.
MIN_BIN_COMPLETENESS = 20.0
MAX_BIN_CONTAMINATION = 5.0
DEFAULT_THREADS = 8
DEFAULT_PROBABILITY_THRESHOLD = 0.5



# ============================================================================
# GENERAL HELPERS
# ============================================================================

def resolve_executable(executable: str, tool_name: str) -> str:
    """Resolve an executable from PATH or an explicit path."""
    resolved = shutil.which(executable)
    if resolved:
        return resolved
    candidate = Path(executable).expanduser()
    if candidate.is_file() and os.access(candidate, os.X_OK):
        return str(candidate.resolve())
    raise FileNotFoundError(
        f"Could not find {tool_name} executable '{executable}'. Install it and "
        f"make sure it is available on PATH."
    )


def run_command(command: List[str], description: str, log_path: Optional[Path] = None) -> None:
    """Run a command, optionally sending stderr to a log file."""
    print(f"\n{description}:")
    print(" ".join(command))
    try:
        if log_path is None:
            subprocess.run(command, check=True)
        else:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            with log_path.open("w", encoding="utf-8") as log_handle:
                subprocess.run(command, check=True, stderr=log_handle)
    except subprocess.CalledProcessError as exc:
        extra = f" See log: {log_path}" if log_path else " See the tool output above."
        raise RuntimeError(
            f"{description} failed with exit code {exc.returncode}.{extra}"
        ) from exc


def ensure_new_or_empty_directory(path: Path, description: str) -> None:
    """Avoid accidentally overwriting an earlier workflow run."""
    if path.exists():
        if not path.is_dir():
            raise FileExistsError(f"{description} path exists but is not a directory: {path}")
        if any(path.iterdir()):
            raise FileExistsError(
                f"{description} directory already exists and is not empty: {path}\n"
                "Choose another path or move the previous output before rerunning."
            )
    else:
        path.mkdir(parents=True, exist_ok=True)


def parse_fasta(filepath: Path) -> Iterable[Tuple[str, str]]:
    """Yield (header, uppercase_sequence) pairs from a FASTA file."""
    header = None
    parts: List[str] = []
    with filepath.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                text = line[1:].strip()
                header = text.split()[0] if text else "unnamed"
                parts = []
            else:
                if header is None:
                    raise ValueError(f"Sequence encountered before FASTA header: {filepath}")
                parts.append(line.upper())
    if header is not None:
        yield header, "".join(parts)


def validate_nonempty_file(path: Path, description: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{description} was not created: {path}")
    if path.stat().st_size == 0:
        raise RuntimeError(f"{description} is empty: {path}")


# ============================================================================
# MODEL LOADING
# ============================================================================


def load_classifier(path: Path):
    if not path.is_file():
        raise FileNotFoundError(f"Classifier file not found: {path}")
    # Only load pickle files from trusted sources.
    with path.open("rb") as handle:
        model = pickle.load(handle)
    if not hasattr(model, "predict_proba"):
        raise TypeError("Classifier must support predict_proba().")
    if not hasattr(model, "classes_"):
        raise TypeError("Classifier has no classes_ attribute.")
    if not hasattr(model, "n_features_in_"):
        raise TypeError("Classifier has no n_features_in_ attribute.")
    if 1 not in list(model.classes_):
        raise ValueError(
            f"Class 1 is not present in model labels: {list(model.classes_)}. "
            "The model must encode probiotic as class 1."
        )
    return model


def positive_class_index(classifier) -> int:
    return list(classifier.classes_).index(1)


# ============================================================================
# PRODIGAL AND CODON-USAGE FEATURES
# ============================================================================

# Standard genetic code, grouped by amino acid. These groups reproduce the
# feature definitions used by the supplied codon-feature extraction script.
AA_TO_CODONS = {
    "F": ["TTT", "TTC"],
    "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
    "I": ["ATT", "ATC", "ATA"],
    "M": ["ATG"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "Y": ["TAT", "TAC"],
    "H": ["CAT", "CAC"],
    "Q": ["CAA", "CAG"],
    "N": ["AAT", "AAC"],
    "K": ["AAA", "AAG"],
    "D": ["GAT", "GAC"],
    "E": ["GAA", "GAG"],
    "C": ["TGT", "TGC"],
    "W": ["TGG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
}
STOP_CODONS = {"TAA", "TAG", "TGA"}
CODON_TO_AA = {
    codon: aa for aa, codons in AA_TO_CODONS.items() for codon in codons
}
SYNONYMOUS_AAS = {
    aa: codons for aa, codons in AA_TO_CODONS.items() if len(codons) > 1
}
RSCU_CODONS = [
    codon
    for aa in sorted(SYNONYMOUS_AAS.keys())
    for codon in sorted(SYNONYMOUS_AAS[aa])
]
CODON_FEATURE_NAMES = [f"RSCU_{codon}" for codon in RSCU_CODONS] + ["ENC", "CAI"]


def run_prodigal(
    fasta_path: Path,
    output_root: Path,
    display_name: str,
    prodigal_executable: str = "prodigal",
) -> Path:
    """Predict CDSs with Prodigal and return the nucleotide CDS FASTA (.ffn)."""
    executable = resolve_executable(prodigal_executable, "Prodigal")
    safe_stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", Path(display_name).stem).strip("._")
    if not safe_stem:
        safe_stem = "genome"
    sample_dir = output_root / safe_stem
    sample_dir.mkdir(parents=True, exist_ok=True)

    prodigal_out = sample_dir / f"{safe_stem}.out"
    ffn_path = sample_dir / f"{safe_stem}.ffn"
    command = [
        executable,
        "-i", str(fasta_path.resolve()),
        "-o", str(prodigal_out.resolve()),
        "-d", str(ffn_path.resolve()),
    ]
    run_command(command, f"Prodigal CDS prediction for {display_name}")
    validate_nonempty_file(prodigal_out, "Prodigal gene-calling output (.out)")
    validate_nonempty_file(ffn_path, "Prodigal CDS nucleotide FASTA (.ffn)")
    return ffn_path


def read_cds(ffn_path: Path) -> List[str]:
    """Read and filter CDS sequences using the supplied codon-feature logic."""
    sequences: List[str] = []
    for _, raw_sequence in parse_fasta(ffn_path):
        sequence = raw_sequence.upper().replace("U", "T")
        if len(sequence) % 3 != 0 or len(sequence) < 60:
            continue
        # Ignore the terminal codon when checking for internal stop codons.
        codons = [sequence[i:i + 3] for i in range(0, len(sequence) - 3, 3)]
        if any(codon in STOP_CODONS for codon in codons):
            continue
        sequences.append(sequence)
    return sequences


def count_codons(cds_list: List[str]) -> Counter:
    """Count codons over CDSs, excluding the terminal codon as in the supplied code."""
    counts = Counter()
    for sequence in cds_list:
        for i in range(0, len(sequence) - 3, 3):
            codon = sequence[i:i + 3]
            if "N" not in codon and codon not in STOP_CODONS:
                counts[codon] += 1
    return counts


def compute_rscu(codon_counts: Counter) -> Dict[str, float]:
    """Compute relative synonymous codon usage (RSCU) for each synonymous codon."""
    rscu: Dict[str, float] = {}
    for aa, codons in SYNONYMOUS_AAS.items():
        total = sum(codon_counts.get(codon, 0) for codon in codons)
        expected = total / len(codons) if total > 0 else 0
        for codon in codons:
            observed = codon_counts.get(codon, 0)
            rscu[codon] = observed / expected if expected > 0 else 0.0
    return rscu


def compute_enc(codon_counts: Counter) -> float:
    """Compute Effective Number of Codons using the supplied Wright-style formula."""
    by_fold = defaultdict(list)
    for aa, codons in AA_TO_CODONS.items():
        if aa in ("M", "W"):
            continue
        counts = [codon_counts.get(codon, 0) for codon in codons]
        total = sum(counts)
        if total < 5:
            continue
        homozygosity = sum((count / total) ** 2 for count in counts)
        by_fold[len(codons)].append(homozygosity)

    def mean_f(fold):
        values = by_fold.get(fold, [])
        return float(np.mean(values)) if values else None

    f2, f3, f4, f6 = (mean_f(fold) for fold in (2, 3, 4, 6))
    enc = 2.0
    enc += 9.0 / f2 if f2 and f2 > 0 else 9.0 / (1.0 / 2)
    enc += 1.0 / f3 if f3 and f3 > 0 else 3.0
    enc += 5.0 / f4 if f4 and f4 > 0 else 5.0 * (1.0 / 4) / (1.0 / 4)
    if f6 and f6 > 0:
        enc += 3.0 / f6
    return float(min(max(enc, 20.0), 61.0))


def build_reference_from_highly_expressed(
    cds_list: List[str], top_fraction: float = 0.05
) -> Dict[str, float]:
    """Build the self-reference from the genome's top 5% longest CDSs."""
    sorted_cds = sorted(cds_list, key=len, reverse=True)
    n_reference = max(1, int(len(sorted_cds) * top_fraction))
    reference_counts = count_codons(sorted_cds[:n_reference])
    weights: Dict[str, float] = {}
    for _, codons in SYNONYMOUS_AAS.items():
        max_count = max(reference_counts.get(codon, 0) for codon in codons)
        for codon in codons:
            weights[codon] = (
                reference_counts.get(codon, 0) / max_count if max_count > 0 else 0.0
            )
    return weights


def compute_cai(cds_list: List[str], reference_weights: Dict[str, float]) -> float:
    """Compute the self-reference Codon Adaptation Index (CAI)."""
    log_weight_sum = 0.0
    count = 0
    for sequence in cds_list:
        for i in range(0, len(sequence) - 3, 3):
            codon = sequence[i:i + 3]
            if codon in STOP_CODONS or "N" in codon:
                continue
            aa = CODON_TO_AA.get(codon)
            if aa and aa in SYNONYMOUS_AAS:
                weight = reference_weights.get(codon, 1e-10)
                log_weight_sum += np.log(max(weight, 1e-10))
                count += 1
    return float(np.exp(log_weight_sum / count)) if count > 0 else 0.0


def extract_codon_features(ffn_path: Path) -> Dict[str, float]:
    """Extract the same 59 RSCU, ENC, and CAI features from Prodigal's .ffn."""
    cds_list = read_cds(ffn_path)
    if len(cds_list) < 50:
        raise ValueError(
            f"Only {len(cds_list)} usable CDSs were found in {ffn_path}; "
            "at least 50 are required for reliable codon-usage features."
        )

    codon_counts = count_codons(cds_list)
    rscu = compute_rscu(codon_counts)
    enc = compute_enc(codon_counts)
    reference_weights = build_reference_from_highly_expressed(cds_list)
    cai = compute_cai(cds_list, reference_weights)

    features = {f"RSCU_{codon}": float(rscu.get(codon, 0.0)) for codon in RSCU_CODONS}
    features["ENC"] = enc
    features["CAI"] = cai
    if len(features) != CODON_FEATURE_COUNT:
        raise RuntimeError(
            f"Internal codon feature count error: expected {CODON_FEATURE_COUNT}, "
            f"generated {len(features)}."
        )
    return features


# ============================================================================
# MEGAHIT ASSEMBLY
# ============================================================================

def run_megahit_assembly(
    reads1: Path,
    reads2: Optional[Path],
    assembly_dir: Path,
    megahit_executable: str,
    threads: int,
    memory_fraction: Optional[float],
) -> Path:
    """Assemble reads with MEGAHIT and return final.contigs.fa."""
    validate_nonempty_file(reads1, "Input read file")
    if reads2 is not None:
        validate_nonempty_file(reads2, "Second input read file")

    executable = resolve_executable(megahit_executable, "MEGAHIT")
    # MEGAHIT expects a new output directory (rather than one pre-created with a log).
    if assembly_dir.exists():
        if not assembly_dir.is_dir():
            raise FileExistsError(f"MEGAHIT output path is not a directory: {assembly_dir}")
        if any(assembly_dir.iterdir()):
            raise FileExistsError(
                f"MEGAHIT output directory is not empty: {assembly_dir}. "
                "Choose a new --assembly-dir or move the previous run."
            )
        assembly_dir.rmdir()
    assembly_dir.parent.mkdir(parents=True, exist_ok=True)

    command = [executable, "-t", str(threads), "-o", str(assembly_dir)]
    if reads2 is not None:
        command.extend(["-1", str(reads1.resolve()), "-2", str(reads2.resolve())])
    else:
        command.extend(["-r", str(reads1.resolve())])
    if memory_fraction is not None:
        command.extend(["-m", str(memory_fraction)])

    run_command(command, "MEGAHIT assembly", assembly_dir.parent / "megahit.stderr.log")
    contigs_fasta = assembly_dir / "final.contigs.fa"
    validate_nonempty_file(contigs_fasta, "MEGAHIT final contigs FASTA")
    print(f"MEGAHIT assembly FASTA: {contigs_fasta}")
    return contigs_fasta


# ============================================================================
# READ MAPPING AND METABAT2 BINNING
# ============================================================================

def map_reads_and_make_depth_file(
    contigs_fasta: Path,
    reads1: Path,
    reads2: Optional[Path],
    mapping_dir: Path,
    bowtie2_build_executable: str,
    bowtie2_executable: str,
    samtools_executable: str,
    depth_executable: str,
    threads: int,
) -> Path:
    """Map reads to contigs and create MetaBAT2's per-contig depth file."""
    bowtie2_build = resolve_executable(bowtie2_build_executable, "bowtie2-build")
    bowtie2 = resolve_executable(bowtie2_executable, "Bowtie2")
    samtools = resolve_executable(samtools_executable, "samtools")
    depth_tool = resolve_executable(depth_executable, "jgi_summarize_bam_contig_depths")

    mapping_dir.mkdir(parents=True, exist_ok=True)
    index_prefix = mapping_dir / "assembly_index"
    run_command(
        [bowtie2_build, str(contigs_fasta.resolve()), str(index_prefix)],
        "Building Bowtie2 index",
        mapping_dir / "bowtie2_build.stderr.log",
    )

    bam_path = mapping_dir / "reads.sorted.bam"
    bowtie_command = [bowtie2, "-x", str(index_prefix), "-p", str(threads)]
    if reads2 is not None:
        bowtie_command.extend(["-1", str(reads1.resolve()), "-2", str(reads2.resolve())])
    else:
        bowtie_command.extend(["-U", str(reads1.resolve())])

    sort_command = [
        samtools, "sort", "-@", str(threads), "-o", str(bam_path), "-"
    ]
    print("\nMapping reads to assembly with Bowtie2, then sorting with samtools:")
    print(" ".join(bowtie_command) + " | " + " ".join(sort_command))
    bowtie_log_path = mapping_dir / "bowtie2.stderr.log"
    sort_log_path = mapping_dir / "samtools_sort.stderr.log"
    with bowtie_log_path.open("w", encoding="utf-8") as bowtie_log, \
            sort_log_path.open("w", encoding="utf-8") as sort_log:
        bowtie_proc = subprocess.Popen(
            bowtie_command, stdout=subprocess.PIPE, stderr=bowtie_log
        )
        assert bowtie_proc.stdout is not None
        sort_error = None
        try:
            subprocess.run(
                sort_command,
                stdin=bowtie_proc.stdout,
                check=True,
                stderr=sort_log,
            )
        except subprocess.CalledProcessError as exc:
            sort_error = exc
        finally:
            bowtie_proc.stdout.close()
        bowtie_return_code = bowtie_proc.wait()
        if sort_error is not None:
            raise RuntimeError(
                f"samtools sort failed with exit code {sort_error.returncode}. "
                f"See log: {sort_log_path}"
            ) from sort_error
        if bowtie_return_code != 0:
            raise RuntimeError(
                f"Bowtie2 failed with exit code {bowtie_return_code}. "
                f"See log: {bowtie_log_path}"
            )

    validate_nonempty_file(bam_path, "Sorted BAM")
    run_command(
        [samtools, "index", str(bam_path)],
        "Indexing sorted BAM",
        mapping_dir / "samtools_index.stderr.log",
    )

    depth_file = mapping_dir / "contig_depth.txt"
    run_command(
        [depth_tool, "--outputDepth", str(depth_file), str(bam_path)],
        "Calculating per-contig depth for MetaBAT2",
        mapping_dir / "depth.stderr.log",
    )
    validate_nonempty_file(depth_file, "MetaBAT2 depth file")
    return depth_file


def run_metabat2(
    contigs_fasta: Path,
    depth_file: Path,
    bins_dir: Path,
    metabat2_executable: str,
    threads: int,
) -> List[Path]:
    """Run MetaBAT2 and return its bin FASTA files."""
    metabat2 = resolve_executable(metabat2_executable, "MetaBAT2")
    bins_dir.mkdir(parents=True, exist_ok=True)
    output_prefix = bins_dir / "bin"
    run_command(
        [
            metabat2,
            "-i", str(contigs_fasta.resolve()),
            "-a", str(depth_file.resolve()),
            "-o", str(output_prefix),
            "-t", str(threads),
        ],
        "MetaBAT2 binning",
        bins_dir / "metabat2.stderr.log",
    )

    bin_files = sorted(
        path.resolve()
        for path in bins_dir.iterdir()
        if path.is_file()
        and path.name.startswith("bin.")
        and path.suffix.lower() in FASTA_EXTENSIONS
        and path.stat().st_size > 0
    )
    if not bin_files:
        raise RuntimeError(
            "MetaBAT2 did not produce any non-empty bin FASTA files. "
            "The assembly may not have yielded bins suitable for MetaBAT2."
        )
    print(f"MetaBAT2 bins generated: {len(bin_files)}")
    return bin_files


# ============================================================================
# CHECKM2 BIN QUALITY ASSESSMENT AND FILTERING
# ============================================================================

def normalise_bin_name(value: object) -> str:
    """Strip FASTA extensions without stripping numeric bin suffixes (e.g. bin.1)."""
    name = Path(str(value).strip()).name
    # MetaBAT2 commonly names bins bin.1.fa, bin.2.fa; CheckM2 may report bin.1.
    # Do not use Path.suffix generically because it treats the .1 in bin.1 as an extension.
    lower_name = name.lower()
    for extension in sorted(FASTA_EXTENSIONS, key=len, reverse=True):
        if lower_name.endswith(extension):
            return name[:-len(extension)]
    return name


def run_checkm2_and_filter_bins(
    bin_files: List[Path],
    checkm2_dir: Path,
    quality_filter_path: Path,
    checkm2_executable: str,
    threads: int,
    min_completeness: float,
    max_contamination: float,
) -> Tuple[Dict[str, Dict], Path]:
    """Run CheckM2, make a filter report, and return quality metrics by bin name."""
    checkm2 = resolve_executable(checkm2_executable, "CheckM2")
    # Let CheckM2 create its output directory itself; do not leave a log file in it.
    if checkm2_dir.exists():
        if not checkm2_dir.is_dir():
            raise FileExistsError(f"CheckM2 output path is not a directory: {checkm2_dir}")
        if any(checkm2_dir.iterdir()):
            raise FileExistsError(f"CheckM2 output directory is not empty: {checkm2_dir}")
        checkm2_dir.rmdir()
    checkm2_dir.parent.mkdir(parents=True, exist_ok=True)

    command = [
        checkm2, "predict",
        "--threads", str(threads),
        "--input", str(bin_files[0].parent.resolve()),
        "-x", "fa",
        "--output-directory", str(checkm2_dir.resolve()),
    ]
    run_command(command, "CheckM2 bin quality assessment", checkm2_dir.parent / "checkm2.stderr.log")
    quality_report = checkm2_dir / "quality_report.tsv"
    validate_nonempty_file(quality_report, "CheckM2 quality_report.tsv")

    report = pd.read_csv(quality_report, sep="\t")
    columns_by_lower = {str(column).strip().lower(): column for column in report.columns}
    name_col = columns_by_lower.get("name")
    completeness_col = columns_by_lower.get("completeness")
    contamination_col = columns_by_lower.get("contamination")
    if name_col is None or completeness_col is None or contamination_col is None:
        raise ValueError(
            "Could not find Name, Completeness, and Contamination columns in "
            f"CheckM2 report. Found columns: {list(report.columns)}"
        )

    quality_by_name: Dict[str, Dict] = {}
    filter_rows = []
    bin_path_by_name = {normalise_bin_name(path.name): path for path in bin_files}

    for _, row in report.iterrows():
        name = normalise_bin_name(row[name_col])
        completeness = pd.to_numeric(row[completeness_col], errors="coerce")
        contamination = pd.to_numeric(row[contamination_col], errors="coerce")
        passed = bool(
            pd.notna(completeness)
            and pd.notna(contamination)
            and float(completeness) >= min_completeness
            and float(contamination) <= max_contamination
        )
        matching_path = bin_path_by_name.get(name)
        if matching_path is None:
            # CheckM2 sometimes strips/changes extensions, so try basename matching.
            matching_path = next(
                (p for p in bin_files if normalise_bin_name(p.name) == name), None
            )
        quality_by_name[name] = {
            "bin": matching_path.name if matching_path else name,
            "bin_path": matching_path,
            "completeness": float(completeness) if pd.notna(completeness) else np.nan,
            "contamination": float(contamination) if pd.notna(contamination) else np.nan,
            "quality_filter_passed": passed,
        }
        filter_rows.append({
            "bin": matching_path.name if matching_path else name,
            "completeness": float(completeness) if pd.notna(completeness) else np.nan,
            "contamination": float(contamination) if pd.notna(contamination) else np.nan,
            "minimum_completeness_required": min_completeness,
            "maximum_contamination_allowed": max_contamination,
            "quality_filter_passed": passed,
        })

    # Include any bin that CheckM2 did not report as a failed quality assessment.
    seen = set(quality_by_name)
    for path in bin_files:
        name = normalise_bin_name(path.name)
        if name not in seen:
            quality_by_name[name] = {
                "bin": path.name,
                "bin_path": path,
                "completeness": np.nan,
                "contamination": np.nan,
                "quality_filter_passed": False,
            }
            filter_rows.append({
                "bin": path.name,
                "completeness": np.nan,
                "contamination": np.nan,
                "minimum_completeness_required": min_completeness,
                "maximum_contamination_allowed": max_contamination,
                "quality_filter_passed": False,
            })

    pd.DataFrame(filter_rows).to_csv(quality_filter_path, index=False)
    accepted = sum(bool(item["quality_filter_passed"]) for item in quality_by_name.values())
    print(
        f"Bins passing CheckM2 filter: {accepted}/{len(bin_files)} "
        f"(completeness >= {min_completeness:.1f}%, "
        f"contamination <= {max_contamination:.1f}%)"
    )
    return quality_by_name, quality_report


# ============================================================================
# CODON-ONLY PREDICTION
# ============================================================================



def predict_one(
    fasta_path: Path,
    display_name: str,
    input_mode: str,
    assembly_mode: str,
    threshold: float,
    classifier,
    prodigal_output_dir: Path,
) -> Dict:
    """Predict probiotic potential from the 61 codon-usage features only."""
    sequences = list(parse_fasta(fasta_path))
    if not sequences:
        raise ValueError(f"No FASTA sequences found in {fasta_path}.")
    stats = {
        "contigs_total": len(sequences),
        "contigs_used": len(sequences),
        "contigs_filtered": 0,
    }

    # Predict CDS regions and derive codon features from Prodigal's .ffn output.
    ffn_path = run_prodigal(
        fasta_path=fasta_path,
        output_root=prodigal_output_dir,
        display_name=display_name,
    )
    codon_features = extract_codon_features(ffn_path)
    codon_vector = np.asarray(
        [codon_features[name] for name in CODON_FEATURE_NAMES], dtype=np.float32
    )

    if len(codon_vector) != classifier.n_features_in_:
        raise ValueError(
            f"Codon feature vector has {len(codon_vector)} features, but the "
            f"classifier expects {classifier.n_features_in_}. The codon-only model "
            "must be trained with the same feature order: 59 RSCU features, ENC, CAI."
        )

    X = pd.DataFrame([codon_vector], columns=CODON_FEATURE_NAMES)
    probability = float(
        classifier.predict_proba(X)[0, positive_class_index(classifier)]
    )
    label = "Probiotic" if probability >= threshold else "Non-Probiotic"
    return {
        "genome": display_name,
        "prediction_task": "Probiotic potential",
        "input_mode": input_mode,
        "assembly_mode": assembly_mode,
        "probability_probiotic": probability,
        "threshold": threshold,
        "predicted_label": label,
        **stats,
        "status": "success",
        "error": "",
    }


def discover_fastas(input_path: Path) -> Tuple[List[Path], str]:
    """Discover FASTAs. A single file is treated as a genome; a directory as MAGs."""
    path = input_path.expanduser()
    if not path.exists():
        raise FileNotFoundError(f"Input path does not exist: {path}")
    if path.is_file():
        if path.suffix.lower() not in FASTA_EXTENSIONS:
            raise ValueError(
                f"'{path}' is not a recognised FASTA file. For sequencing reads, "
                "use --assembly with -i1 (and optionally -i2)."
            )
        return [path.resolve()], "complete"
    files = sorted(
        p.resolve() for p in path.iterdir()
        if p.is_file() and p.suffix.lower() in FASTA_EXTENSIONS
    )
    if not files:
        raise ValueError(f"No FASTA files found in directory: {path}")
    return files, "mag"


# ============================================================================
# COMMAND-LINE INTERFACE
# ============================================================================

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="probiolens.py",
        description=(
            "ProbioLens predicts probiotic potential. Supply an assembled FASTA "
            "with -i, or use --assembly with -i1 and optional -i2 to assemble reads."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples
--------
1. Predict from one assembled genome or MAG FASTA:
   python3 probiolens.py -i genome.fasta -o results.csv

2. Predict every FASTA in a directory (one prediction per file):
   python3 probiolens.py -i /path/to/mags/ -o results.csv

3. Assemble paired-end isolate reads and predict the assembly:
   python3 probiolens.py -i1 reads_R1.fastq.gz -i2 reads_R2.fastq.gz --assembly --isolate -o results.csv

4. Assemble single-end metagenomic reads, bin, assess quality, and predict each MAG:
   python3 probiolens.py -i1 reads.fastq.gz --assembly --meta --complete 20 --contam 5 -o results.csv

5. Paired-end metagenomic reads with stricter thresholds:
   python3 probiolens.py --input1 R1.fq.gz --input2 R2.fq.gz --assembly --meta --complete 50 --contam 3 --output results.csv

Main options
------------
-i, --input       Assembled FASTA file or directory of MAG FASTA files
-i1, --input1     First read file; single-end if -i2 is omitted
-i2, --input2     Second read file for paired-end reads
-o, --output      Required output CSV path
--assembly        Assemble reads with MEGAHIT (alias: --assemble)
--isolate         With assembly, predict the assembly directly
--meta            With assembly, bin contigs and predict each passing MAG
--complete        Minimum CheckM2 completeness percentage (default: 20)
--contam          Maximum CheckM2 contamination percentage (default: 5)

Tools must be available on PATH. CheckM2 must already be configured to locate
its database. Model files are loaded from the model/ directory beside this script.
""",
    )
    input_group = parser.add_mutually_exclusive_group(required=False)
    input_group.add_argument(
        "-i", "--input",
        help="Already assembled genome/MAG FASTA file or directory of MAG FASTAs.",
    )
    input_group.add_argument(
        "-i1", "--input1",
        help="First sequencing-read file when --assembly is used (single-end if -i2 omitted).",
    )
    parser.add_argument(
        "-i2", "--input2",
        help="Second sequencing-read file for paired-end reads; requires -i1 and --assembly.",
    )
    parser.add_argument("-o", "--output", required=True, help="Output prediction CSV path.")
    parser.add_argument(
        "--assembly", "--assemble", dest="assembly", action="store_true",
        help="Assemble sequencing reads with MEGAHIT before prediction.",
    )
    assembly_group = parser.add_mutually_exclusive_group(required=False)
    assembly_group.add_argument("--isolate", action="store_true", help="With --assembly, predict the assembly directly.")
    assembly_group.add_argument("--meta", action="store_true", help="With --assembly, bin contigs and predict each passing MAG.")
    parser.add_argument(
        "--complete", type=float, default=20.0,
        help="Minimum CheckM2 completeness percentage for MAGs (default: 20).",
    )
    parser.add_argument(
        "--contam", type=float, default=5.0,
        help="Maximum CheckM2 contamination percentage for MAGs (default: 5).",
    )

    return parser


def _base_result(
    genome: str,
    assembly_mode: str,
    threshold: float,
    status: str,
    error: str = "",
    completeness: float = np.nan,
    contamination: float = np.nan,
    quality_filter_passed: object = "not_assessed",
) -> Dict:
    return {
        "genome": genome,
        "prediction_task": "Probiotic potential",
        "input_mode": "mag" if assembly_mode == "meta" else "complete",
        "assembly_mode": assembly_mode,
        "completeness_percent": completeness,
        "contamination_percent": contamination,
        "quality_filter_passed": quality_filter_passed,
        "probability_probiotic": np.nan,
        "threshold": threshold,
        "predicted_label": status,
        "contigs_total": np.nan,
        "contigs_used": np.nan,
        "contigs_filtered": np.nan,
        "status": status,
        "error": error,
    }


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()


    if not 0.0 <= args.complete <= 100.0:
        parser.error("--complete must be between 0 and 100.")
    if not 0.0 <= args.contam <= 100.0:
        parser.error("--contam must be between 0 and 100.")

    if args.assembly:
        if not args.input1:
            parser.error("With --assembly, provide reads using -i1/--input1 (and optionally -i2/--input2).")
        if args.input2 and not args.input1:
            parser.error("-i2/--input2 requires -i1/--input1.")
        if args.input is not None:
            parser.error("Use -i for assembled FASTA input OR -i1/-i2 for reads with --assembly, not both.")
        if not (args.isolate or args.meta):
            parser.error("With --assembly, choose exactly one mode: --isolate or --meta.")
    else:
        if not args.input:
            parser.error("Without --assembly, provide an assembled FASTA file or directory using -i/--input.")
        if args.input1 or args.input2:
            parser.error("-i1/-i2 are read inputs and require --assembly.")
        if args.isolate or args.meta:
            parser.error("--isolate and --meta are only used with --assembly.")

    output_path = Path(args.output).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pipeline_dir: Optional[Path] = None

    try:
        print("\n" + "=" * 72)
        print("ProbioLens — probiotic-potential prediction")
        print("=" * 72)
        print("Prediction model: probiotic potential only")
        print(f"Probability threshold: {DEFAULT_PROBABILITY_THRESHOLD}")
        if args.assembly:
            mode = "isolate" if args.isolate else "meta"
            print(f"Assembly mode: {mode}")
            if mode == "meta":
                print(
                    f"CheckM2 criteria: completeness >= {args.complete:.1f}%; "
                    f"contamination <= {args.contam:.1f}%"
                )

        # Load the codon-only probiotic classifier; embeddings are not used.
        classifier = load_classifier(MODEL_FILES["classifier"])
        if classifier.n_features_in_ != CODON_FEATURE_COUNT:
            raise ValueError(
                "Codon-only model dimension mismatch: classifier expects "
                f"{classifier.n_features_in_} features, but this pipeline supplies "
                f"{CODON_FEATURE_COUNT} codon-usage features (59 RSCU + ENC + CAI). "
                "Use the updated 61-feature classifier at the existing model path."
            )

        results: List[Dict] = []
        if args.assembly:
            reads1 = Path(args.input1).expanduser().resolve()
            reads2 = Path(args.input2).expanduser().resolve() if args.input2 else None
            validate_nonempty_file(reads1, "Input read file (-i1)")
            if reads2 is not None:
                validate_nonempty_file(reads2, "Second input read file (-i2)")

            # Pipeline files are stored beside the requested results CSV; no path flag is needed.
            pipeline_dir = output_path.parent / f"{output_path.stem}_pipeline"
            ensure_new_or_empty_directory(pipeline_dir, "Pipeline output")
            assembly_dir = pipeline_dir / "assembly"
            contigs_fasta = run_megahit_assembly(
                reads1=reads1,
                reads2=reads2,
                assembly_dir=assembly_dir,
                megahit_executable="megahit",
                threads=DEFAULT_THREADS,
                memory_fraction=None,
            )

            if args.isolate:
                print("Isolate mode: predicting the assembled contigs directly; skipping binning and CheckM2.")
                result = predict_one(
                    fasta_path=contigs_fasta,
                    display_name=f"{reads1.name}_assembly",
                    input_mode="complete",
                    assembly_mode="isolate",
                    threshold=DEFAULT_PROBABILITY_THRESHOLD,
                    classifier=classifier,
                    prodigal_output_dir=pipeline_dir / "prodigal",
                )
                result.update({
                    "completeness_percent": np.nan,
                    "contamination_percent": np.nan,
                    "quality_filter_passed": "not_assessed",
                })
                results.append(result)
            else:
                mapping_dir = pipeline_dir / "mapping"
                depth_file = map_reads_and_make_depth_file(
                    contigs_fasta=contigs_fasta,
                    reads1=reads1,
                    reads2=reads2,
                    mapping_dir=mapping_dir,
                    bowtie2_build_executable="bowtie2-build",
                    bowtie2_executable="bowtie2",
                    samtools_executable="samtools",
                    depth_executable="jgi_summarize_bam_contig_depths",
                    threads=DEFAULT_THREADS,
                )
                bins_dir = pipeline_dir / "bins"
                bin_files = run_metabat2(
                    contigs_fasta=contigs_fasta,
                    depth_file=depth_file,
                    bins_dir=bins_dir,
                    metabat2_executable="metabat2",
                    threads=DEFAULT_THREADS,
                )
                quality_filter_path = pipeline_dir / "bin_quality_filter.csv"
                quality_by_name, quality_report = run_checkm2_and_filter_bins(
                    bin_files=bin_files,
                    checkm2_dir=pipeline_dir / "checkm2",
                    quality_filter_path=quality_filter_path,
                    checkm2_executable="checkm2",
                    threads=DEFAULT_THREADS,
                    min_completeness=args.complete,
                    max_contamination=args.contam,
                )
                print(f"Full CheckM2 report: {quality_report}")
                print(f"Quality filter summary: {quality_filter_path}")

                for bin_path in bin_files:
                    bin_key = normalise_bin_name(bin_path.name)
                    quality = quality_by_name.get(bin_key, {
                        "completeness": np.nan,
                        "contamination": np.nan,
                        "quality_filter_passed": False,
                    })
                    completeness = quality["completeness"]
                    contamination = quality["contamination"]
                    passed = bool(quality["quality_filter_passed"])
                    if not passed:
                        results.append(_base_result(
                            genome=bin_path.name,
                            assembly_mode="meta",
                            threshold=DEFAULT_PROBABILITY_THRESHOLD,
                            status="filtered_out_quality",
                            error=(
                                f"Did not meet CheckM2 criteria: completeness >= {args.complete:g}% "
                                f"and contamination <= {args.contam:g}%."
                            ),
                            completeness=completeness,
                            contamination=contamination,
                            quality_filter_passed=False,
                        ))
                        continue

                    try:
                        result = predict_one(
                            fasta_path=bin_path,
                            display_name=bin_path.name,
                            input_mode="mag",
                            assembly_mode="meta",
                            threshold=DEFAULT_PROBABILITY_THRESHOLD,
                            classifier=classifier,
                            prodigal_output_dir=pipeline_dir / "prodigal",
                        )
                        result.update({
                            "completeness_percent": completeness,
                            "contamination_percent": contamination,
                            "quality_filter_passed": True,
                        })
                        results.append(result)
                        print(
                            f"{bin_path.name}: completeness={completeness:.2f}%, "
                            f"contamination={contamination:.2f}%, "
                            f"probability={result['probability_probiotic']:.4f}, "
                            f"label={result['predicted_label']}"
                        )
                    except Exception as exc:
                        print(f"ERROR predicting {bin_path.name}: {exc}", file=sys.stderr)
                        results.append(_base_result(
                            genome=bin_path.name,
                            assembly_mode="meta",
                            threshold=DEFAULT_PROBABILITY_THRESHOLD,
                            status="failed",
                            error=str(exc),
                            completeness=completeness,
                            contamination=contamination,
                            quality_filter_passed=True,
                        ))
        else:
            fasta_files, inferred_input_mode = discover_fastas(Path(args.input))
            print(f"FASTA files to predict: {len(fasta_files)}")
            for index, fasta_path in enumerate(fasta_files, start=1):
                print(f"\n[{index}/{len(fasta_files)}] Processing: {fasta_path.name}")
                try:
                    result = predict_one(
                        fasta_path=fasta_path,
                        display_name=fasta_path.name,
                        input_mode=inferred_input_mode,
                        assembly_mode="",
                        threshold=DEFAULT_PROBABILITY_THRESHOLD,
                        classifier=classifier,
                        prodigal_output_dir=output_path.parent / f"{output_path.stem}_prodigal",
                    )
                    result.update({
                        "completeness_percent": np.nan,
                        "contamination_percent": np.nan,
                        "quality_filter_passed": "not_assessed",
                    })
                    results.append(result)
                    print(f"Probability probiotic: {result['probability_probiotic']:.4f}")
                    print(f"Label: {result['predicted_label']}")
                    print(f"Contigs used: {result['contigs_used']}/{result['contigs_total']}")
                except Exception as exc:
                    print(f"ERROR: {exc}", file=sys.stderr)
                    failed = _base_result(
                        genome=fasta_path.name,
                        assembly_mode="",
                        threshold=DEFAULT_PROBABILITY_THRESHOLD,
                        status="failed",
                        error=str(exc),
                    )
                    failed["input_mode"] = inferred_input_mode
                    results.append(failed)

        columns = [
            "genome",
            "prediction_task",
            "input_mode",
            "assembly_mode",
            "completeness_percent",
            "contamination_percent",
            "quality_filter_passed",
            "probability_probiotic",
            "threshold",
            "predicted_label",
        ]
        pd.DataFrame(results, columns=columns).to_csv(output_path, index=False)

        succeeded = sum(row.get("status") == "success" for row in results)
        filtered = sum(row.get("status") == "filtered_out_quality" for row in results)
        failures = sum(row.get("status") == "failed" for row in results)
        print("\n" + "=" * 72)
        print(f"Results saved: {output_path}")
        print(f"Successful predictions: {succeeded}")
        print(f"Bins excluded by CheckM2 filter: {filtered}")
        print(f"Failed predictions: {failures}")
        if pipeline_dir is not None:
            print(f"Pipeline outputs: {pipeline_dir}")
        print("=" * 72)

        if failures:
            sys.exit(1)
    except Exception as exc:
        print(f"\nFATAL ERROR: {exc}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
