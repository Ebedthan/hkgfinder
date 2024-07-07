"""library of sub-program for hkgfinder."""

from __future__ import annotations

import importlib.resources
import logging
import sys
import warnings
from collections import namedtuple
from pathlib import Path

import psutil
import pyfastx
import pyhmmer
import pyrodigal
import xphyle
from Bio import SeqIO
from Bio.Seq import reverse_complement, translate
from xphyle import paths

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
    newline = "\n"
    ext = "faa" if is_prot else "fna"

    if split:
        for result in filtered_results:
            filepath = Path(f"{base_filename}_{result.hmm_id}.{ext}")
            with filepath.open("a", encoding="utf-8") as out:
                seq = seqs[result.hit_id].seq
                formatted_seq = [seq[i : i + 60] for i in range(0, len(seq), 60)]
                out.write(
                    f">{result.hit_id}_gene={result.hmm_id}\n{newline.join(formatted_seq)}\n",
                )
    else:
        filepath = Path(f"{base_filename}.{ext}")
        with filepath.open("w", encoding="utf-8") as out:
            for result in filtered_results:
                seq = seqs[result.hit_id].seq
                formatted_seq = [seq[i : i + 60] for i in range(0, len(seq), 60)]
                out.write(
                    f">{result.hit_id}_gene={result.hmm_id}\n{newline.join(formatted_seq)}\n",
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

    def format_result(result) -> str:
        return f"{result.hit_id}\t{result.hmm_id}\t{result.hmm_desc}\t{result.evalue:.3f}\t{result.bitscore:.3f}"

    output_lines = [header] + [format_result(result) for result in filtered_results]

    if to_stdout:
        print("\n".join(output_lines))
    else:
        with output_path.open("w", encoding="utf-8") as ofile:
            ofile.write("\n".join(output_lines) + "\n")


def search_hmm(seqdata: Path, cpus: int) -> dict:
    """Search HMM and return best results."""
    results = []
    best_results = {}
    available_memory = psutil.virtual_memory().available
    target_size = seqdata.stat().st_size

    # Open HMM file
    hmm_path = importlib.resources.files("hkgfinder").joinpath("hkgfinder.hmm")
    with hmm_path.open("rb") as hmm, pyhmmer.plan7.HMMFile(
        hmm,  # type: ignore[type-supported]
    ) as hmm_file, pyhmmer.easel.SequenceFile(seqdata, digital=True) as seq_file:
        # Pre-fetch targets into memory if size is appropriate
        if target_size < available_memory * 0.1:
            logging.info("Pre-fetching targets into memory")
            targets = seq_file.read_block()
            mem = sys.getsizeof(targets) + sum(
                sys.getsizeof(target) for target in targets
            )
            mem_kb = mem / 1024
            logging.info("Database in-memory size: %s KiB", f"{mem_kb:.1f}")
        else:
            targets = seq_file

        HMMResult = namedtuple(
            "HMMResult",
            ["hmm_id", "hmm_desc", "hit_id", "evalue", "bitscore"],
        )

        # Perform hmmsearch
        for hits in pyhmmer.hmmer.hmmsearch(
            hmm_file,
            targets,  # type: ignore[type-supported]
            cpus=cpus,
            bit_cutoffs="gathering",  # type: ignore[type-supported]
        ):
            hmm_id = hits.query_name
            hmm_desc = HMMDESC[str(hmm_id, encoding="utf-8")]  # type: ignore[type-supported]

            for hit in hits:
                if hit.included:
                    results.append(
                        HMMResult(
                            str(hmm_id, encoding="utf-8"),  # type: ignore[type-supported]
                            hmm_desc,
                            str(hit.name, encoding="utf-8"),
                            hit.evalue,
                            hit.score,
                        ),
                    )

        # Find best results based on bitscore
        for result in results:
            if result.hit_id in best_results:
                previous_bitscore = best_results[result.hit_id].bitscore
                if result.bitscore > previous_bitscore:
                    best_results[result.hit_id] = result
            else:
                best_results[result.hit_id] = result

    return best_results


def find_protein_coding_genes(
    infile: str,
    outdir: str,
    *,
    is_meta: bool = False,
    save_both: bool,
) -> None:
    """Find protein coding genes using pyrodigal."""
    records = SeqIO.index(infile, "fasta")
    orf_finder = pyrodigal.OrfFinder(meta=is_meta)
    orf_finder.train(*(bytes(record.seq) for record in records.values()))
    with Path.open(
        Path(outdir, "my.proteins.faa"),
        "w+",
        encoding="utf-8",
    ) as out:
        for seqid, seq in records.items():
            genes = orf_finder.find_genes(bytes(seq.seq))
            genes.write_translations(out, sequence_id=seqid, width=60)
    if save_both:
        with Path.open(
            Path(outdir, "my.dna.fna"),
            "w+",
            encoding="utf-8",
        ) as out:
            for seqid, seq in records.items():
                genes = orf_finder.find_genes(bytes(seq.seq))
                genes.write_genes(out, sequence_id=seqid, width=60)


def handle_cpus(asked_cpus: int, available_cpus: int | None) -> int:
    """Allocate the good number of CPUs based on asked cpus vs available cpus."""
    cpus = 1
    if asked_cpus == 0 or asked_cpus > int(available_cpus or 1):
        cpus = available_cpus
    else:
        cpus = asked_cpus
    return int(cpus or 1)


def decompress_file(input_file: str, outdir: str) -> str:
    """Decompress a file and return the path to the uncompressed file."""
    output_path = Path(outdir, "hkgfinder_input.fa")
    # Handling input file supply
    if input_file == "<stdin>":
        with output_path.open("w+", encoding="utf-8") as binput, xphyle.xopen(
            paths.STDIN,
            context_wrapper=True,
        ) as infile:
            for line in infile:
                binput.write(line)  # type: ignore[type-supported]
    else:
        with output_path.open("w+", encoding="utf-8") as binput, xphyle.xopen(
            input_file,
            context_wrapper=True,
        ) as infile:
            for line in infile:
                binput.write(line)  # type: ignore[type-supported]
    return str(output_path)


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
    output_path = Path(outfile)
    # Use buffer to reduce I/O operations
    buffer_size = 1024 * 1024  # 1 MB buffer
    buffer = []
    with output_path.open("w") as protfile:
        for sequence in fa:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                protseq = _translate_seq(sequence.seq)
                for idx, frame in enumerate(protseq):
                    seq_letters = [frame[i : i + sw] for i in range(0, len(frame), sw)]
                    nl = "\n"
                    buffer.append(
                        f">{sequence.name}_frame={idx + 1}\n{nl.join(seq_letters)}\n",
                    )
                    if len(buffer) > buffer_size:
                        protfile.write("".join(buffer))
                        buffer = []
        # Write remaining buffer to file
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
