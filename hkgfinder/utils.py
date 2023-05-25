"""Module providing utils for hkgfinder."""

# Copyright 2022-2023 Anicet Ebou.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed except according
# to those terms.

import warnings
from datetime import datetime
import os
import pyfastx
from pathlib import Path
import typing

from Bio.Seq import reverse_complement, translate


def elapsed_since(start):
    walltime = datetime.now() - start
    return walltime


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


def write_seq(file, seq_id, seq):
    seq = [seq[i : i + 60] for i in range(0, len(seq), 60)]
    newline = "\n"
    file.write(f">{seq_id}" + f"\n{newline.join(map(str, seq))}\n")
