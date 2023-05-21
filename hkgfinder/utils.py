"""Module providing utils for hkgfinder."""

# Copyright 2022-2023 Anicet Ebou.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed except according
# to those terms.

import warnings
from datetime import datetime
from sys import stderr
import os
import pyfastx
from pathlib import Path

from Bio.Seq import reverse_complement, translate


def get_file_type(filename):
    """
    get_file_type get file compression type matching the magic bytes

    :filename: File name of the file to check
    """

    magic_dict = {
        b"\x1f\x8b": "gz",
        b"\x42\x5a": "bz2",
        b"\xfd\x37\x7a\x58\x5a": "lzma",
    }

    # The longer byte to read in at file start
    # is max(len(x) for x in magic_dict) which gives 7 as result
    with open(filename, "rb") as f:
        # Read at most 7 bytes at file start
        file_start = f.read(7)

        # Match bytes read with compression type
        for magic, filetype in magic_dict.items():
            if file_start.startswith(magic):
                return filetype

        # If no match, the compression is not known or not compressed
        return "unknown"


def get_best_match(mdict):
    result = {}
    for key in mdict.keys():
        sorted_x = dict(
            sorted(
                mdict[key].items(), key=lambda item: item[1][1], reverse=True
            )
        )
        fk, vk = next(iter(sorted_x.items()))
        result[key] = f"{fk}#{vk[0]}#{vk[1]}"

    return result


def print_install_tool(tool):
    """
    print_install_tool print useful installation
    instruction for required tools.
    """

    if tool == "hmmsearch":
        err("hmmsearch not found. Please visit https://hmmer3.org.")
    elif tool == "prodigal":
        err(
            "prodigal not found. Please visit "
            + "https://github.com/hyattpd/Prodigal"
        )


def elapsed_since(start):
    walltime = datetime.now() - start
    return walltime


def err(text):
    print(f"Error: {text}", file=stderr)


def msg(text, is_quiet):
    """
    msg produces nice message and info output on terminal.

    :text: Message to print to STDOUT.
    """
    if not is_quiet:
        print(f"[{datetime.now().strftime('%H:%M:%S')}][INFO] {text}")


def warn(text, is_quiet):
    if not is_quiet:
        print(f"[{datetime.now().strftime('%H:%M:%S')}][WARN] {text}")


def is_dna_or_aa(s):
    """
    is_dna_or_aa test if input sequence is DNA or proteins.

    :s: input sequence
    """

    dna = "ATCG"
    prot = "ABCDEFGHIKLMNPQRSTVWYZ*X"
    stype = ""

    if all(i in dna for i in s):
        stype = "DNA"
    elif all(i in prot for i in s):
        stype = "protein"
    else:
        stype = "unknown"

    return stype


def _translate_seq(seq):
    """
    _translate_seq translate DNA sequence to proteins in the six frames.

    :seq: DNA sequence to translate.
    """

    seqlist = []
    # frame 1
    seqlist.append(translate(seq))
    # frame 2
    seqlist.append(translate(seq[1:]))
    # frame 3
    seqlist.append(translate(seq[2:]))
    # frame 4
    seqlist.append(translate(reverse_complement(seq)))
    # frame 5
    seqlist.append(translate(reverse_complement(seq)[1:]))
    # frame 6
    seqlist.append(translate(reverse_complement(seq)[2:]))

    return seqlist


def do_translation(infile: str, outfile: Path, width=60):
    """
    do_translation translate a DNA fasta file into proteins
    fasta file.

    :infile: Pyfasta object.
    :outfile: Output file.
    :sw: Sequence width. Default: 60.
    """

    file = pyfastx.Fasta(infile)

    with open(outfile, "w", encoding="utf-8") as protfile:
        for sequence in file:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                protseq = _translate_seq(sequence.seq)
                for idx, frame in enumerate(protseq):
                    seq_letters = [
                        frame[i : i + width]
                        for i in range(0, len(frame), width)
                    ]
                    new_line = "\n"
                    protfile.write(
                        f">{sequence.name}_frame={idx + 1}\n"
                        + f"{new_line.join(map(str, seq_letters))}\n"
                    )


hmmdesc = {
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


def get_hmm_desc(hmm_id):
    return hmmdesc[hmm_id]


def parse_cpu(cpu, is_quiet):
    cpus = 0
    available_cpus = os.cpu_count()
    if cpu == 0:
        cpus = available_cpus
    elif cpu > available_cpus:
        warn(
            f"Option -t asked for {cpu} threads,"
            + f" but system has only {available_cpus}.",
            is_quiet,
        )
        cpus = available_cpus

    return cpus


def write_seq(file, seq_id, seq):
    seq = [seq[i : i + 60] for i in range(0, len(seq), 60)]
    newline = "\n"
    file.write(f">{seq_id}" + f"\n{newline.join(map(str, seq))}\n")
