
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ProbioLens standalone prediction tool.

Prediction:
    P = Probiotic potential
    S = Safety potential

Modes:
    complete = one complete-genome FASTA
    mag      = one MAG FASTA or a directory containing MAG FASTA files

MAG filtering:
    Only contigs strictly longer than 5,000 bp are used by default.

Model files are loaded from the directory containing this script.
"""

import argparse
import pickle
import sys
import warnings
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from gensim.models import KeyedVectors, Word2Vec
from sklearn.exceptions import InconsistentVersionWarning


# ============================================================
# CONFIGURATION — EDIT HERE IF MODEL FILENAMES CHANGE
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent

MODEL_FILES = {
    "P": {
        "classifier": SCRIPT_DIR / "model/probiolens_p_classifier.pkl",
        "vectors": SCRIPT_DIR / "model/probiolens_p_vectors.kv",
    },
    "S": {
        "classifier": SCRIPT_DIR / "model/probiolens_s_classifier.pkl",
        "vectors": SCRIPT_DIR / "model/probiolens_s_vectors.kv",
    },
}

K = 8
DEFAULT_MIN_CONTIG_LENGTH = 5000

FASTA_EXTENSIONS = {".fa", ".fasta", ".fna", ".fas"}

# Suppress the scikit-learn version warning in console output.
# This does not fix incompatibilities between model versions.
warnings.filterwarnings("ignore", category=InconsistentVersionWarning)


# ============================================================
# FASTA PARSER
# ============================================================

def parse_fasta(filepath):
    """Yield FASTA header and uppercase sequence."""

    header = None
    parts = []

    with open(filepath, "r", encoding="utf-8") as handle:
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
                    raise ValueError(
                        f"Sequence encountered before FASTA header: {filepath}"
                    )
                parts.append(line.upper())

    if header is not None:
        yield header, "".join(parts)


# ============================================================
# LOAD MODELS
# ============================================================

def load_word_vectors(path):
    if not path.is_file():
        raise FileNotFoundError(f"Word-vector file not found: {path}")

    if path.suffix.lower() == ".kv":
        return KeyedVectors.load(str(path))

    return Word2Vec.load(str(path)).wv


def load_classifier(path):
    if not path.is_file():
        raise FileNotFoundError(f"Classifier file not found: {path}")

    # Only load pickle files from trusted sources.
    with open(path, "rb") as handle:
        model = pickle.load(handle)

    if not hasattr(model, "predict_proba"):
        raise TypeError("Classifier must support predict_proba().")

    if not hasattr(model, "classes_"):
        raise TypeError("Classifier has no classes_ attribute.")

    if not hasattr(model, "n_features_in_"):
        raise TypeError("Classifier has no n_features_in_ attribute.")

    return model


def positive_class_index(classifier):
    """Return the probability column corresponding to class 1."""

    classes = list(classifier.classes_)

    if 1 not in classes:
        raise ValueError(
            f"Class 1 is not present in the model labels: {classes}. "
            "Check the training label encoding."
        )

    return classes.index(1)


# ============================================================
# 8-MER EMBEDDING
# ============================================================

def valid_kmer(kmer):
    return (
        len(kmer) == K
        and all(base in "ACGT" for base in kmer)
    )


def embed_fasta(fasta_path, word_vectors, vocabulary_mean, mode, min_contig_length,):
    """Generate an occurrence-weighted mean 8-mer embedding."""

    sequences = list(parse_fasta(fasta_path))

    if not sequences:
        raise ValueError("No FASTA sequences found.")

    total_contigs = len(sequences)
    used_contigs = 0

    embedding_sum = np.zeros(
        word_vectors.vector_size,
        dtype=np.float64,
    )

    known_occurrences = 0

    for _, sequence in sequences:

        # In MAG mode, retain only contigs strictly longer than the cutoff.
        if mode == "mag" and len(sequence) <= min_contig_length:
            continue

        used_contigs += 1

        if len(sequence) < K:
            continue

        counts = Counter(
            sequence[i:i + K]
            for i in range(len(sequence) - K + 1)
            if valid_kmer(sequence[i:i + K])
        )

        for kmer, count in counts.items():
            if kmer in word_vectors:
                embedding_sum += (
                    np.asarray(word_vectors[kmer], dtype=np.float64)
                    * count
                )
                known_occurrences += count

    if used_contigs == 0:
        raise ValueError(
            f"No contigs passed the length filter (> {min_contig_length} bp)."
        )

    if known_occurrences == 0:
        raise ValueError(
            "No valid 8-mers from the retained sequences were found "
            "in the word-vector vocabulary."
        )

    mean_embedding = embedding_sum / known_occurrences

    # Match the centering procedure in the original prediction pipeline.
    centered_embedding = mean_embedding - vocabulary_mean

    stats = {
        "contigs_total": total_contigs,
        "contigs_used": used_contigs,
        "contigs_filtered": total_contigs - used_contigs,
    }

    return centered_embedding.astype(np.float32), stats


# ============================================================
# INPUT FILE DISCOVERY
# ============================================================

def discover_fastas(input_path, mode):
    path = Path(input_path)

    if not path.exists():
        raise FileNotFoundError(f"Input path does not exist: {path}")

    if path.is_file():
        if path.suffix.lower() not in FASTA_EXTENSIONS:
            raise ValueError(f"Unsupported FASTA extension: {path.suffix}")
        return [path]

    if mode == "complete":
        raise ValueError(
            "Complete mode requires one FASTA file, not a directory."
        )

    files = sorted(
        p for p in path.iterdir()
        if p.is_file() and p.suffix.lower() in FASTA_EXTENSIONS
    )

    if not files:
        raise ValueError(f"No FASTA files found in directory: {path}")

    return files


# ============================================================
# PREDICTION
# ============================================================

def predict_one(fasta_path, mode, prediction, threshold, word_vectors, classifier, vocabulary_mean, min_contig_length,):
    embedding, stats = embed_fasta(
        fasta_path=fasta_path,
        word_vectors=word_vectors,
        vocabulary_mean=vocabulary_mean,
        mode=mode,
        min_contig_length=min_contig_length,
    )

    if len(embedding) != classifier.n_features_in_:
        raise ValueError(
            f"Embedding has {len(embedding)} features, but the classifier "
            f"expects {classifier.n_features_in_}."
        )

    feature_names = [f"dim_{i}" for i in range(len(embedding))]
    X = pd.DataFrame([embedding], columns=feature_names)

    class_index = positive_class_index(classifier)
    probability = float(classifier.predict_proba(X)[0, class_index])

    is_positive = probability >= threshold

    if prediction == "P":
        label = "Probiotic" if is_positive else "Non-Probiotic"
        task = "Probiotic potential"
    else:
        label = "Safe" if is_positive else "Not Safe"
        task = "Safety potential"

    return {
        "genome": fasta_path.name,
        "prediction_task": task,
        "mode": mode,
        "probability_positive": probability,
        "threshold": threshold,
        "predicted_label": label,
        **stats,
        "status": "success",
        "error": "",
    }


# ============================================================
# MAIN
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        prog="probiolens.py",
        description="""
ProbioLens: Genome-based prediction of probiotic potential and safety.

Choose a prediction task:
  P  Predict probiotic potential (probiotic / non-probiotic)
  S  Predict safety potential (safe / not safe)

Choose an input mode:
  complete  One complete-genome FASTA file
  mag       One MAG FASTA file or a directory of MAG FASTA files

In MAG mode, only contigs strictly longer than 5,000 bp
are used by default.
""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
EXAMPLES

1. Probiotic prediction for a complete genome:
   python3 probiolens.py -i genome.fasta -o results.csv -m complete -p P

2. Safety prediction for a complete genome:
   python3 probiolens.py -i genome.fasta -o results.csv -m complete -p S

3. Probiotic prediction for a MAG:
   python3 probiolens.py -i bin.fasta -o results.csv -m mag -p P

4. Safety prediction for a MAG:
   python3 probiolens.py -i bin.fasta -o results.csv -m mag -p S -t 0.6

5. Process all MAG FASTA files in a directory:
   python3 probiolens.py -i /path/to/MAGs/ -o results.csv -m mag -p S

6. Change the MAG contig-length cutoff:
   python3 probiolens.py -i bin.fasta -o results.csv -m mag -p P --min-contig 10000

7. Display this help:
   python3 probiolens.py --help

NOTES
  - Default probability threshold: 0.5
  - Default MAG contig cutoff: 5,000 bp (strictly greater than)
  - Model paths are configured inside the script.
  - P and S load their respective saved model and word-vector files.
""",
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input FASTA file, or a MAG directory in MAG mode.",
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output CSV file path.",
    )

    parser.add_argument(
        "-m", "--mode",
        required=True,
        choices=["complete", "mag"],
        help="Input mode: complete genome or MAG.",
    )

    parser.add_argument(
        "-p", "--prediction",
        required=True,
        choices=["P", "S"],
        help="Prediction task: P = probiotic potential; S = safety potential.",
    )

    parser.add_argument(
        "-t", "--threshold",
        type=float,
        default=0.5,
        help="Probability threshold between 0 and 1 (default: 0.5).",
    )

    parser.add_argument(
        "--min-contig",
        type=int,
        default=5000,
        help="MAG mode: retain contigs strictly longer than this length in bp (default: 5000).",
    )

    args = parser.parse_args()

    if not 0.0 <= args.threshold <= 1.0:
        parser.error("--threshold must be between 0 and 1.")

    if args.min_contig < 0:
        parser.error("--min-contig cannot be negative.")

    # Keep the rest of your existing main() function below this point:
    # model loading, prediction, and CSV output.

    try:
        print("\n" + "=" * 60)
        print("ProbioLens prediction")
        print("=" * 60)
        print(f"Mode: {args.mode}")
        print(f"Prediction: ProbioLens-{args.prediction}")
        print(f"Threshold: {args.threshold}")

        if args.mode == "mag":
            print(f"Contig filter: > {args.min_contig} bp")

        # Load only the selected model pair.
        selected = MODEL_FILES[args.prediction]

        word_vectors = load_word_vectors(selected["vectors"])
        classifier = load_classifier(selected["classifier"])

        if classifier.n_features_in_ != word_vectors.vector_size:
            raise ValueError(
                "Model dimension mismatch: classifier expects "
                f"{classifier.n_features_in_} features but the vectors have "
                f"{word_vectors.vector_size} dimensions."
            )

        positive_class_index(classifier)

        # Match the vocabulary-mean centering procedure used during training.
        vocabulary_mean = np.asarray(
            word_vectors.vectors.mean(axis=0),
            dtype=np.float64,
        )

        fasta_files = discover_fastas(args.input, args.mode)
        print(f"FASTA files found: {len(fasta_files)}")

        results = []

        for index, fasta_path in enumerate(fasta_files, start=1):
            print(
                f"\n[{index}/{len(fasta_files)}] "
                f"Processing: {fasta_path.name}"
            )

            try:
                result = predict_one(
                    fasta_path=fasta_path,
                    mode=args.mode,
                    prediction=args.prediction,
                    threshold=args.threshold,
                    word_vectors=word_vectors,
                    classifier=classifier,
                    vocabulary_mean=vocabulary_mean,
                    min_contig_length=args.min_contig,
                )

                results.append(result)

                print(
                    f"Probability: {result['probability_positive']:.4f}"
                )
                print(f"Label: {result['predicted_label']}")
                print(
                    f"Contigs used: {result['contigs_used']}/"
                    f"{result['contigs_total']}"
                )

            except Exception as exc:
                print(f"ERROR: {exc}", file=sys.stderr)

                results.append({
                    "genome": fasta_path.name,
                    "prediction_task": (
                        "Probiotic potential"
                        if args.prediction == "P"
                        else "Safety potential"
                    ),
                    "mode": args.mode,
                    "probability_positive": np.nan,
                    "threshold": args.threshold,
                    "predicted_label": "prediction_failed",
                    "contigs_total": np.nan,
                    "contigs_used": np.nan,
                    "contigs_filtered": np.nan,
                    "status": "failed",
                    "error": str(exc),
                })

        # Keep only the requested, useful output columns.
        columns = [
            "genome",
            "prediction_task",
            "mode",
            "probability_positive",
            "threshold",
            "predicted_label",
            "contigs_total",
            "contigs_used",
            "contigs_filtered",
            "status",
            "error",
        ]

        output_path = Path(args.output).resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)

        pd.DataFrame(results, columns=columns).to_csv(
            output_path,
            index=False,
        )

        successes = sum(r["status"] == "success" for r in results)
        failures = len(results) - successes

        print("\n" + "=" * 60)
        print(f"Results saved: {output_path}")
        print(f"Successful: {successes}")
        print(f"Failed: {failures}")
        print("=" * 60)

        if failures:
            sys.exit(1)

    except Exception as exc:
        print(f"\nFATAL ERROR: {exc}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()

