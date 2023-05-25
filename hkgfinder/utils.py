"""Module providing utils for hkgfinder."""

# Copyright 2022-2023 Anicet Ebou.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed except according
# to those terms.

import typing
import warnings
from pathlib import Path

import pyfastx
from Bio.Seq import reverse_complement, translate


def _translate_seq(seq: str) -> typing.List[str]:
    """
    _translate_seq translate DNA sequence to proteins in the six frames.

    :seq: DNA sequence to translate.
    """

    seqlist = []

    rev_comp = reverse_complement(seq)

    # frame 1
    seqlist.append(translate(seq))
    # frame 2
    seqlist.append(translate(seq[1:]))
    # frame 3
    seqlist.append(translate(seq[2:]))
    # frame 4
    seqlist.append(translate(rev_comp))
    # frame 5
    seqlist.append(translate(rev_comp[1:]))
    # frame 6
    seqlist.append(translate(rev_comp[2:]))

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


def write_seq(file: typing.TextIO, seq_id: str, seq: str):
    """
    write_seq provide a convenient interface to write a sequence
    to an already opened file.
    """
    chunks = [seq[i : i + 60] for i in range(0, len(seq), 60)]
    newline = "\n"
    file.write(f">{seq_id}\n{newline.join(map(str, chunks))}\n")
