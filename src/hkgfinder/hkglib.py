"""library of sub-program for hkgfinder."""

from __future__ import annotations

import bz2
import datetime
import gzip
import logging
import lzma
import multiprocessing.pool
import os
import platform
import shutil
import sys
import warnings
from pathlib import Path
from typing import Optional, Union

import pyrodigal
from Bio import SeqIO
from Bio.Seq import reverse_complement, translate
from pyfastxcli import pyfastx

from hkgfinder.models import RunMode

from .config import HKGFinderConfig
from .models import ProcessingContext, SequenceType
from .writers import GenBankWriter, SequenceWriter


def setup_logging(quiet: bool) -> None:
    """Configure logging for the application."""
    level = logging.CRITICAL if quiet else logging.INFO
    logging.basicConfig(
        format="[%(asctime)s][%(levelname)s] %(message)s",
        datefmt="%H:%M:%S",
        level=level,
    )


def get_username() -> str:
    """Get current username safely."""
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
    """Determine number of CPUs to use."""
    available_cpus = os.cpu_count()
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
    Process sequences based on run mode.

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
    """Write output sequences if requested."""
    if not (args.faa or args.fna or args.genbank):
        return

    # Load protein sequences
    if ctx.run_mode in (RunMode.GENOME, RunMode.METAGENOME):
        proteins = pyfastx.Fasta(str(Path(tmpdir, "my.proteins.faa")))
    elif ctx.sequence_type == SequenceType.PROTEIN:
        proteins = original_fasta
    else:
        proteins = pyfastx.Fasta(str(Path(tmpdir, "input_translate.fa")))

    # Sort if requested
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
        if ctx.run_mode in (RunMode.GENOME, RunMode.METAGENOME):
            dna_sequences = pyfastx.Fasta(str(Path(tmpdir, "my.dna.fna")))
        else:
            dna_sequences = original_fasta

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
        if args.faa or ctx.sequence_type == SequenceType.PROTEIN:
            # Write protein GenBank
            gb_writer = GenBankWriter(
                os.path.splitext(args.genbank)[0],
                args.s,
                config,
            )
            gb_writer.write_genbank(
                filtered_results,
                proteins,
                is_protein=True,
                organism=args.organism,
            )
        else:
            # Write DNA GenBank
            if ctx.run_mode in (RunMode.GENOME, RunMode.METAGENOME):
                dna_sequences = pyfastx.Fasta(str(Path(tmpdir, "my.dna.fna")))
            else:
                dna_sequences = original_fasta

            gb_writer = GenBankWriter(
                os.path.splitext(args.genbank)[0],
                args.s,
                config,
            )
            gb_writer.write_genbank(
                filtered_results,
                dna_sequences,
                is_protein=False,
                organism=args.organism,
            )


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


def _handle_cpus(asked_cpus: int, available_cpus: int | None) -> int:
    """Allocate the good number of CPUs based on asked cpus vs available cpus."""
    available = available_cpus or 1
    return min(max(1, asked_cpus), available) if asked_cpus > 0 else available


def _detect_compression(file_path: Union[str, Path]) -> Optional[str]:
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
        compression = _detect_compression(input_path)

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
