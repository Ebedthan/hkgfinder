# Copyright 2022-2023 Anicet Ebou.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed except according
# to those terms.

from datetime import datetime
import sys
from Bio.Seq import reverse_complement, translate
import warnings


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
    print(f"Error: {text}", file=sys.stderr)


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


def do_translation(infile, outfile, sw=60):
    """
    do_translation translate a DNA fasta file into proteins
    fasta file.

    :infile: Pyfasta object.
    :outfile: Output file.
    :sw: Sequence width. Default: 60.
    """

    with open(outfile, "w") as protfile:
        for sequence in infile:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                protseq = _translate_seq(sequence.seq)
                for idx, frame in enumerate(protseq):
                    # Rule E203 from flacke8 check for extraneous whitespace
                    # before a colon. But black follow PEP8 rules.
                    # A PR is open to resolve this issue:
                    # https://github.com/PyCQA/pycodestyle/pull/914
                    seq_letters = [
                        frame[i : i + sw]  # noqa: E203
                        for i in range(0, len(frame), sw)
                    ]
                    nl = "\n"
                    protfile.write(
                        f">{sequence.name}_frame={idx + 1}\n"
                        + f"{nl.join(map(str, seq_letters))}\n"
                    )
