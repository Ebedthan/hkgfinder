"""library of sub-program for hkgfinder."""

from __future__ import annotations

import bz2
import gzip
import importlib.resources
import logging
import lzma
import multiprocessing.pool
import shutil
import sys
import warnings
from collections import namedtuple
from pathlib import Path
from typing import Optional, Union

import psutil
import pyhmmer
import pyrodigal
from Bio import SeqIO
from Bio.Seq import reverse_complement, translate
from pyfastxcli import pyfastx

# Constants
HMMDESC = {
    "MnmE": "tRNA uridine-5-carboxymethylaminomethyl(34) synthesis GTPase MnmE",
    "DnaK": "Molecular chaperonne DnaK",
    "GyrB": "DNA topoisomerase (ATP-hydrolyzing) subunit B",
    "RecA": "recombinase RecA",
    "rpoB": "DNA-directed RNA polymerase subunit beta",
    "infB": "translation initiation factor IF-2",
    "atpD": "F0F1 ATP synthase subunit beta",
    "GroEL": "chaperonin GroEL",
    "fusA": "Elongation factor G",
    "ileS": "isoleucine--tRNA ligase",
    "lepA": "translation elongation factor 4",
    "leuS_bact": "leucine--tRNA ligase bacteria",
    "leuS_arch": "leucine--tRNA ligase archaea",
    "PyrG": "CTP synthase (glutamine hydrolyzing)",
    "recG": "ATP-dependent DNA helicase RecG",
    "rplB_bact": "50S ribosomal protein L2",
    "nifH": "nitrogenase iron protein",
    "nodC": "chitooligosaccharide synthase NodC",
}
BUFFER_SIZE = 1024 * 1024  # 1 MB buffer
SEQ_WIDTH = 60


def write_sequences(
    base_filename: str,
    filtered_results: list,
    seqs: dict,
    *,
    split: bool,
    is_prot: bool,
) -> None:
    """Write sequences to output files based on the filtering results.

    Args:
    ----
    base_filename: Base filename.
    filtered_results: List of filtered result objects.
    seqs: Dictionary containing protein sequences.
    split: Boolean indicating whether to split output into separate files.
    is_prot: Boolean indicating the type of sequence to write.

    """
    ext = "faa" if is_prot else "fna"

    def format_sequence(seq_id: str, hmm_id: str, sequence: str) -> str:
        chunks = [
            sequence[i : i + SEQ_WIDTH]
            for i in range(0, len(sequence), SEQ_WIDTH)
        ]
        return f">{seq_id}_gene={hmm_id}\n" + "\n".join(chunks) + "\n"

    if split:
        # Group results by hmm_id to minimize file opens
        grouped = {}
        for result in filtered_results:
            grouped.setdefault(result.hmm_id, []).append(result)

        for hmm_id, results in grouped.items():
            filepath = Path(f"{base_filename}_{hmm_id}.{ext}")
            with filepath.open("a", encoding="utf-8") as out:
                buffer = []
                for result in results:
                    seq = seqs[result.hit_id].seq
                    buffer.append(format_sequence(result.hit_id, hmm_id, seq))
                    if len(buffer) >= 100:  # Write in chunks of 100 sequences
                        out.write("".join(buffer))
                        buffer = []
                if buffer:
                    out.write("".join(buffer))
    else:
        filepath = Path(f"{base_filename}.{ext}")
        with filepath.open("w", encoding="utf-8") as out:
            buffer = []
            for result in filtered_results:
                seq = seqs[result.hit_id].seq
                buffer.append(
                    format_sequence(result.hit_id, result.hmm_id, seq)
                )
                if len(buffer) >= 100:
                    out.write("".join(buffer))
                    buffer = []
            if buffer:
                out.write("".join(buffer))


def write_results(
    output_path: Path,
    filtered_results: list,
    *,
    to_stdout: bool = False,
) -> None:
    """Write filtered results to the specified output file or stdout.

    Args:
    ----
    output_path: Path to the output file.
    filtered_results: List of filtered result objects.
    to_stdout: Flag to indicate whether to write to stdout. Default is False.

    """
    header = "seq_name\tpred_gene\tgene_desc\te-value\tbitscore"
    lines = [header]

    for result in filtered_results:
        lines.append(
            f"{result.hit_id}\t{result.hmm_id}\t{result.hmm_desc}\t"
            f"{result.evalue:.3f}\t{result.bitscore:.3f}",
        )
    output = "\n".join(lines)
    if to_stdout:
        print(output)
    else:
        output_path.write_text(output, encoding="utf-8")


def search_hmm(seqdata: Path, cpus: int) -> dict:
    """Search HMM and return best results."""
    results = []
    best_results = {}

    # Memory check optimization
    available_memory = psutil.virtual_memory().available
    target_size = seqdata.stat().st_size
    should_prefetch = target_size < available_memory * 0.1

    hmm_path = importlib.resources.files("hkgfinder").joinpath("hkgfinder.hmm")
    with (
        hmm_path.open("rb") as hmm,
        pyhmmer.plan7.HMMFile(
            hmm,  # type: ignore[type-supported]
        ) as hmm_file,
        pyhmmer.easel.SequenceFile(seqdata, digital=True) as seq_file,
    ):
        # Pre-fetch targets into memory if size is appropriate
        if should_prefetch:
            logging.info("Pre-fetching targets into memory")
            targets = seq_file.read_block()
            mem_kb = sum(sys.getsizeof(target) for target in targets) / 1024
            logging.info("Database in-memory size: %.1f KiB", mem_kb)
        else:
            targets = seq_file

        HMMResult = namedtuple(
            "HMMResult",
            ["hmm_id", "hmm_desc", "hit_id", "evalue", "bitscore"],
        )

        # Process hits in bulk
        for hits in pyhmmer.hmmer.hmmsearch(
            hmm_file,
            targets,  # type: ignore[type-supported]
            cpus=cpus,
            bit_cutoffs="gathering",  # type: ignore[type-supported]
        ):
            hmm_id = hits.query.name.decode("utf-8")
            hmm_desc = HMMDESC[hmm_id]

            for hit in hits:
                if hit.included:
                    results.append(
                        HMMResult(
                            hmm_id,  # type: ignore[type-supported]
                            hmm_desc,
                            hit.name.decode("utf-8"),
                            hit.evalue,
                            hit.score,
                        ),
                    )

        # Find best results based on bitscore
        for result in results:
            existing = best_results.get(result.hit_id)
            if existing is None or result.bitscore > existing.bitscore:
                best_results[result.hit_id] = result

    return best_results


def find_protein_coding_genes(
    infile: str,
    outdir: str,
    ncpus: int,
    *,
    is_meta: bool = False,
    save_both: bool,
) -> None:
    """Find protein coding genes using pyrodigal."""
    records = SeqIO.index(infile, "fasta")
    gene_finder = pyrodigal.GeneFinder(meta=is_meta)

    if not is_meta:
        gene_finder.train(*(bytes(record.seq) for record in records.values()))

    protein_path = Path(outdir, "my.proteins.faa")
    dna_path = Path(outdir, "my.dna.fna") if save_both else None

    with multiprocessing.pool.ThreadPool(processes=ncpus) as pool:
        sequences = (bytes(record.seq) for record in records.values())
        predictions = pool.map(gene_finder.find_genes, sequences)

        with protein_path.open("w", encoding="utf-8") as prot_out:
            if dna_path and save_both:
                with dna_path.open("w", encoding="utf-8") as dna_out:
                    for gene in predictions:
                        gene.write_translations(prot_out, sequence_id="gene")
                        gene.write_genes(dna_out, sequence_id="gene")
            else:
                for gene in predictions:
                    gene.write_translations(prot_out, sequence_id="gene")


def handle_cpus(asked_cpus: int, available_cpus: int | None) -> int:
    """Allocate the good number of CPUs based on asked cpus vs available cpus."""
    available = available_cpus or 1
    return min(max(1, asked_cpus), available) if asked_cpus > 0 else available


def detect_compression(file_path: Union[str, Path]) -> Optional[str]:
    """Detect compression type by reading magic numbers from file header.

    Returns:
        'gzip', 'bzip2', 'lzma', or None if uncompressed

    """
    magic_numbers = {
        b"\x1f\x8b\x08": "gzip",
        b"\x42\x5a\x68": "bzip2",
        b"\xfd\x37\x7a\x58\x5a\x00": "lzma",
    }
    with open(file_path, "rb") as f:
        file_start = f.read(6)  # Read enough bytes to detect all formats

    for magic, format_type in magic_numbers.items():
        if file_start.startswith(magic):
            return format_type
    return None


def decompress_file(
    input_file: Union[str, Path],
    output_dir: Union[str, Path],
    output_filename: str = "decompressed_output.fa",
    force_extension: bool = False,
) -> Path:
    """Enhanced decompression function with magic number detection.

    Args:
        input_file: Path to the compressed input file.
        output_dir: Directory where the decompressed file will be saved.
        output_filename: Name of the decompressed file.
        force_extension: If True, only use file extension for detection.

    Returns:
        Path to the decompressed file.

    """
    input_path = Path(input_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / output_filename

    # Determine compression type
    compression = None

    if force_extension:
        # Simple extension-based detection
        suffixes = input_path.suffixes
        if ".gz" in suffixes:
            compression = "gzip"
        elif ".bz2" in suffixes:
            compression = "bzip2"
        elif ".xz" in suffixes or ".lzma" in suffixes:
            compression = "lzma"
    else:
        # First try magic number detection
        compression = detect_compression(input_path)

        # Fall back to extension detection if magic number fails
        if compression is None:
            suffixes = input_path.suffixes
            if ".gz" in suffixes:
                compression = "gzip"
            elif ".bz2" in suffixes:
                compression = "bzip2"
            elif ".xz" in suffixes or ".lzma" in suffixes:
                compression = "lzma"

    if compression is None:
        raise ValueError(
            f"Cannot determine compression format for {input_path}"
        )

    # Decompress the file
    decompressors = {
        "gzip": gzip.open,
        "bzip2": bz2.open,
        "lzma": lzma.open,
    }

    try:
        with decompressors[compression](input_path, "rb") as f_in:
            with open(output_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
    except Exception as e:
        if output_path.exists():
            output_path.unlink()
        raise OSError(f"Failed to decompress {input_path}: {e!s}")

    return output_path


def do_translation(infile: str, outfile: Path, sw: int = 60) -> None:
    """Translate a DNA fasta file into a proteins fasta file.

    Args:
    ----
    infile: input file path as string.
    outfile: Output file as path.
    sw: Sequence width. Default: 60.

    Returns:
    -------
    Translated sequences

    """
    fa = pyfastx.Fasta(infile)
    buffer = []

    with outfile.open("w") as protfile:
        for sequence in fa:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                protseq = _translate_seq(sequence.seq)
                for idx, frame in enumerate(protseq, 1):
                    chunks = [
                        frame[i : i + sw] for i in range(0, len(frame), sw)
                    ]
                    buffer.append(
                        f">{sequence.name}_frame={idx}\n"
                        + "\n".join(chunks)
                        + "\n"
                    )
                    if len(buffer) >= 100:
                        protfile.write("".join(buffer))
                        buffer = []
        if buffer:
            protfile.write("".join(buffer))


def _translate_seq(seq: str) -> list[str]:
    """Translate DNA sequence to proteins in the six frames.

    Args:
    ----
    seq: DNA sequence to translate.

    Returns:
    -------
    A list of amino-acids in six frames

    """
    rc_seq = reverse_complement(seq)
    return [
        str(translate(seq)),
        str(translate(seq[1:])),
        str(translate(seq[2:])),
        str(translate(rc_seq)),
        str(translate(rc_seq[1:])),
        str(translate(rc_seq[2:])),
    ]
