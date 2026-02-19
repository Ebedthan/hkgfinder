"""library of sub-program for hkgfinder."""

from __future__ import annotations

import bz2
import datetime
import gzip
import logging
import lzma
import multiprocessing as mp
import os
import platform
import shutil
import sys
import warnings
from collections import deque
from functools import lru_cache
from pathlib import Path
from typing import BinaryIO, Iterator, Optional, Union

import pyrodigal
from Bio import SeqIO
from Bio.Seq import reverse_complement, translate
from pyfastxcli import pyfastx

from hkgfinder.models import RunMode

from .config import HKGFinderConfig
from .models import ProcessingContext, SequenceType
from .writers import GenBankWriter, SequenceWriter

# Constants for optimization
CHUNK_SIZE = 1000  # number of sequences to process in parallel batches
IO_BUFFER_SIZE = 65536  # 64KB for file I/O


def setup_logging(quiet: bool, debug: bool = False) -> None:
    """
    Configure logging for the application.

    Args:
        quiet: Suppress informational messages
        debug: Enable debug level logging
    """
    if debug:
        level = logging.DEBUG
    elif quiet:
        level = logging.CRITICAL
    else:
        level = logging.INFO

    logging.basicConfig(
        format="[%(asctime)s][%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
        level=level,
    )


@lru_cache(maxsize=1)
def get_username() -> str:
    """Get current username safely with caching."""
    return os.environ.get("USER", os.environ.get("USERNAME", "unknown user"))


def log_startup_info(config: HKGFinderConfig, run_mode: RunMode) -> None:
    """Log startup information."""
    logging.info("This is hkgfinder %s", config.VERSION)
    logging.info("Written by %s", config.AUTHOR)
    logging.info("Available at %s", config.URL)
    logging.info(
        "Localtime is %s",
        datetime.datetime.now(tz=datetime.timezone.utc).strftime("%H:%M:%S"),
    )
    logging.info("You are %s", get_username())
    logging.info("Operating system is %s", platform.system())
    logging.info("You are running hkgfinder in %s mode", run_mode.value)


def determine_run_mode(args) -> RunMode:
    """Determine execution mode from arguments."""
    if args.g:
        return RunMode.GENOME
    elif args.m:
        return RunMode.METAGENOME
    else:
        return RunMode.NORMAL


def determine_cpu_count(requested: int) -> int:
    """
    Determine number of CPUs to use with optimized defaults.

    Args:
        requested: Number of CPUs requested by user

    Returns:
        Optimal number of CPUs to use
    """
    available_cpus = os.cpu_count() or 1
    logging.info("System has %s cores", available_cpus)

    cpus = _handle_cpus(requested, available_cpus)
    logging.info("We will use maximum of %s cores", cpus)

    return cpus


def process_sequences(
    ctx: ProcessingContext,
    config: HKGFinderConfig,
    tmpdir: str,
) -> Path:
    """
    Process sequences based on run mode with optimizations.

    Args:
        ctx: Processing context
        config: Configuration object
        tmpdir: Temporary directory path

    Returns:
        Path to processed sequence file
    """
    # Genome or metagenome mode - predict genes
    if ctx.run_mode in (RunMode.GENOME, RunMode.METAGENOME):
        logging.info("Predicting protein-coding genes")
        find_protein_coding_genes(
            ctx.input_file,
            tmpdir,
            ctx.num_cpus,
            is_meta=(ctx.run_mode == RunMode.METAGENOME),
            save_both=ctx.save_dna,
        )
        return Path(tmpdir, "my.proteins.faa")

    # Normal mode with DNA - translate
    elif ctx.sequence_type == SequenceType.DNA:
        logging.info("Translating sequences into 6 frames")
        output_path = Path(tmpdir, "input_translate.fa")
        do_translation(ctx.input_file, output_path)

        # Validate translated sequences aren't too long
        fa = pyfastx.Fasta(str(output_path))
        if len(fa.longest) >= config.max_seq_length:
            logging.error(
                "Translated sequence length exceeds %s bp. Consider genome mode (-g)?",
                config.max_seq_length,
            )
            sys.exit(1)

        return output_path

    # Normal mode with protein - use directly
    else:
        return Path(ctx.input_file)


def write_output_sequences(
    args,
    filtered_results: list,
    tmpdir: str,
    ctx: ProcessingContext,
    original_fasta,
    config: HKGFinderConfig,
) -> None:
    """
    Write output sequences if requested.

    Args:
        args: Command line arguments
        filtered_results: List of filtered HMM results
        tmpdir: Temporary directory path
        ctx: Processing context
        original_fasta: Original FASTA file object
        config: Configuration object
    """
    if not (args.faa or args.fna or args.genbank):
        return

    # Load protein sequences - cache to avoid re-reading
    proteins = _load_sequences(
        ctx, tmpdir, is_protein=True, original_fasta=original_fasta
    )

    # Sort if requested (in-place for efficiency)
    if args.s:
        filtered_results.sort(key=lambda r: r.hit_id)

    # Write protein sequences
    if args.faa:
        logging.info("Writing predicted protein sequences")
        writer = SequenceWriter(
            os.path.splitext(args.faa)[0],
            "faa",
            args.s,
            config,
        )
        writer.write_sequences(filtered_results, proteins)

    # Write DNA sequences
    if args.fna:
        logging.info("Writing predicted DNA sequences")
        dna_sequences = _load_sequences(
            ctx, tmpdir, is_protein=False, original_fasta=original_fasta
        )

        writer = SequenceWriter(
            os.path.splitext(args.fna)[0],
            "fna",
            args.s,
            config,
        )
        writer.write_sequences(filtered_results, dna_sequences)

    # Write GenBank format
    if args.genbank:
        logging.info("Writing sequences in GenBank format")

        # Determine if we're writing proteins or DNA
        is_protein = args.faa or ctx.sequence_type == SequenceType.PROTEIN
        sequences = (
            proteins
            if is_protein
            else _load_sequences(
                ctx, tmpdir, is_protein=False, original_fasta=original_fasta
            )
        )

        gb_writer = GenBankWriter(
            os.path.splitext(args.genbank)[0],
            args.s,
            config,
        )
        gb_writer.write_genbank(
            filtered_results,
            sequences,
            is_protein=is_protein,
            organism=getattr(args, "organism", None),
        )


def _load_sequences(
    ctx: ProcessingContext,
    tmpdir: str,
    is_protein: bool,
    original_fasta,
) -> pyfastx.Fasta:
    """
    Load sequences based on context.

    Args:
        ctx: Processing context
        tmpdir: Temporary directory
        is_protein: Whether to load protein or DNA sequences
        original_fasta: Original FASTA file object

    Returns:
        pyfastx.Fasta object
    """
    if ctx.run_mode in (RunMode.GENOME, RunMode.METAGENOME):
        filename = "my.proteins.faa" if is_protein else "my.dna.fna"
        return pyfastx.Fasta(str(Path(tmpdir, filename)))
    elif is_protein and ctx.sequence_type == SequenceType.PROTEIN:
        return original_fasta
    elif is_protein:
        return pyfastx.Fasta(str(Path(tmpdir, "input_translate.fa")))
    else:
        return original_fasta


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


def find_protein_coding_genes(
    infile: str,
    outdir: str,
    ncpus: int,
    *,
    is_meta: bool = False,
    save_both: bool = False,
) -> None:
    """
    Find protein coding genes using pyrodigal with parallel processing.

    Args:
        infile: Input FASTA file path
        outdir: Output directory
        ncpus: Number of CPUs to use
        is_meta: Use metagenome mode
        save_both: Save both protein and DNA sequences
    """
    records = SeqIO.index(infile, "fasta")
    gene_finder = pyrodigal.GeneFinder(meta=is_meta)

    # Training phase (not needed for metagenome mode)
    if not is_meta:
        logging.info("Training gene finder on input sequences")
        # Convert to bytes efficiently
        training_sequences = (bytes(record.seq) for record in records.values())
        gene_finder.train(*training_sequences)

    protein_path = Path(outdir, "my.proteins.faa")
    dna_path = Path(outdir, "my.dna.fna") if save_both else None

    # Use process pool for CPU-intensive gene finding
    num_sequences = len(records)
    logging.info(
        f"Finding genes in {num_sequences} sequences using {ncpus} cores"
    )

    with mp.Pool(processes=ncpus) as pool:
        # Create sequence iterator
        sequences = (bytes(record.seq) for record in records.values())

        # Process in chunks for better memory management
        chunk_size = max(1, num_sequences // (ncpus * 4))
        predictions = pool.imap(
            gene_finder.find_genes, sequences, chunksize=chunk_size
        )

        # Write results with buffering
        with protein_path.open(
            "w", encoding="utf-8", buffering=IO_BUFFER_SIZE
        ) as prot_out:
            if dna_path and save_both:
                with dna_path.open(
                    "w", encoding="utf-8", buffering=IO_BUFFER_SIZE
                ) as dna_out:
                    _write_gene_predictions(predictions, prot_out, dna_out)  # type: ignore[type-supported]
            else:
                _write_gene_predictions(predictions, prot_out, None)  # type: ignore[type-supported]

    logging.info("Gene prediction complete")


def _write_gene_predictions(
    predictions: Iterator,
    protein_file: BinaryIO,
    dna_file: Optional[BinaryIO],
) -> None:
    """
    Write gene predictions to output files.

    Args:
        predictions: Iterator of gene predictions
        protein_file: Output file for proteins
        dna_file: Optional output file for DNA sequences
    """
    genes_found = 0

    for genes in predictions:
        if genes:
            genes_found += len(genes)
            genes.write_translations(
                protein_file,
                sequence_id="gene",
                include_stop=False,
                full_id=True,
            )
            if dna_file:
                genes.write_genes(
                    dna_file,
                    sequence_id="gene",
                    full_id=True,
                )

    logging.debug(f"Found {genes_found} genes total")


def _handle_cpus(asked_cpus: int, available_cpus: int | None) -> int:
    """
    Allocate optimal number of CPUs.

    Args:
        asked_cpus: Number of CPUs requested
        available_cpus: Number of CPUs available

    Returns:
        Number of CPUs to use
    """
    available = available_cpus or 1

    if asked_cpus <= 0:
        # Use all available, but leave one for system
        return max(1, available - 1)

    return min(max(1, asked_cpus), available)


@lru_cache(maxsize=128)
def _detect_compression(file_path: Union[str, Path]) -> Optional[str]:
    """
    Detect compression type by reading magic numbers from file header.

    Cached for performance when checking the same file multiple times.

    Args:
        file_path: Path to file

    Returns:
        'gzip', 'bzip2', 'lzma', or None if uncompressed
    """
    magic_numbers = {
        b"\x1f\x8b\x08": "gzip",
        b"\x42\x5a\x68": "bzip2",
        b"\xfd\x37\x7a\x58\x5a\x00": "lzma",
    }

    try:
        with open(file_path, "rb") as f:
            file_start = f.read(6)
    except (IOError, OSError) as e:
        logging.warning(f"Could not read file {file_path}: {e}")
        return None

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
    """
    Enhanced decompression function with magic number detection.

    Args:
        input_file: Path to the compressed input file
        output_dir: Directory where the decompressed file will be saved
        output_filename: Name of the decompressed file
        force_extension: If True, only use file extension for detection

    Returns:
        Path to the decompressed file

    Raises:
        ValueError: If compression format cannot be determined
        OSError: If decompression fails
    """
    input_path = Path(input_file)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / output_filename

    # Determine compression type
    compression = None

    if force_extension:
        compression = _detect_compression_by_extension(input_path)
    else:
        # First try magic number detection
        compression = _detect_compression(input_path)

        # Fall back to extension detection if magic number fails
        if compression is None:
            compression = _detect_compression_by_extension(input_path)

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

    logging.info(f"Decompressing {input_path} using {compression}")

    try:
        with decompressors[compression](input_path, "rb") as f_in:
            with open(output_path, "wb", buffering=IO_BUFFER_SIZE) as f_out:
                # Use efficient copying with large buffer
                shutil.copyfileobj(f_in, f_out, length=IO_BUFFER_SIZE)
    except Exception as e:
        if output_path.exists():
            output_path.unlink()
        raise OSError(f"Failed to decompress {input_path}: {e!s}") from e

    logging.info(f"Decompressed to {output_path}")
    return output_path


def _detect_compression_by_extension(file_path: Path) -> Optional[str]:
    """
    Detect compression from file extension.

    Args:
        file_path: Path to file

    Returns:
        Compression type or None
    """
    suffixes = file_path.suffixes

    if ".gz" in suffixes:
        return "gzip"
    elif ".bz2" in suffixes:
        return "bzip2"
    elif ".xz" in suffixes or ".lzma" in suffixes:
        return "lzma"

    return None


def do_translation(
    infile: str,
    outfile: Path,
    seq_width: int = 60,
) -> None:
    """
    Translate a DNA FASTA file into a proteins FASTA file with optimized I/O.

    Args:
        infile: Input file path
        outfile: Output file path
        seq_width: Sequence width for output (default: 60)
    """
    fa = pyfastx.Fasta(infile)
    num_sequences = len(fa)
    logging.info(f"Translating {num_sequences} DNA sequences")

    # Use deque for efficient buffering
    buffer = deque()
    buffer_size = 0
    max_buffer_size = 100  # Number of sequences to buffer

    sequences_processed = 0

    with outfile.open(
        "w", encoding="utf-8", buffering=IO_BUFFER_SIZE
    ) as protfile:
        for sequence in fa:
            # Suppress Biopython warnings for partial codons
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")

                # Translate to all 6 frames
                protseqs = _translate_seq(sequence.seq)

                for idx, frame in enumerate(protseqs, 1):
                    # Format sequence with line wrapping
                    formatted = _format_translated_sequence(
                        sequence.name,
                        idx,
                        frame,
                        seq_width,
                    )
                    buffer.append(formatted)
                    buffer_size += 1

                    # Write buffer when full
                    if buffer_size >= max_buffer_size:
                        protfile.write("".join(buffer))
                        buffer.clear()
                        buffer_size = 0

            sequences_processed += 1
            if sequences_processed % 1000 == 0:
                logging.debug(
                    f"Translated {sequences_processed}/{num_sequences} sequences"
                )

        # Write remaining buffer
        if buffer:
            protfile.write("".join(buffer))

    logging.info(f"Translation complete: {sequences_processed} sequences")


def _format_translated_sequence(
    seq_name: str,
    frame_num: int,
    sequence: str,
    seq_width: int,
) -> str:
    """
    Format a translated sequence with header and line wrapping.

    Args:
        seq_name: Sequence name
        frame_num: Frame number (1-6)
        sequence: Protein sequence
        seq_width: Width for line wrapping

    Returns:
        Formatted FASTA entry
    """
    # Use generator for memory efficiency on long sequences
    chunks = [
        sequence[i : i + seq_width] for i in range(0, len(sequence), seq_width)
    ]

    return f">{seq_name}_frame={frame_num}\n" + "\n".join(chunks) + "\n"


def _translate_seq(seq: str) -> list[str]:
    """
    Translate DNA sequence to proteins in six frames.

    Uses to_stop=False to handle partial codons gracefully.

    Args:
        seq: DNA sequence to translate

    Returns:
        List of amino acid sequences in six frames
    """
    # Calculate reverse complement once
    rc_seq = reverse_complement(seq)

    # Translate all frames efficiently
    # Using to_stop=False to avoid exceptions on partial codons
    return [
        str(translate(seq, to_stop=False, table="1")),
        str(translate(seq[1:], to_stop=False, table="1")),
        str(translate(seq[2:], to_stop=False, table="1")),
        str(translate(rc_seq, to_stop=False, table="1")),
        str(translate(rc_seq[1:], to_stop=False, table="1")),
        str(translate(rc_seq[2:], to_stop=False, table="1")),
    ]
