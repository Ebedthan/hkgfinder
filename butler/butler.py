#!/usr/bin/env python3

# Copyright 2022 Anicet Ebou.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed except according
# to those terms.


import argparse
from Bio import SearchIO
from collections import defaultdict
from datetime import datetime
from distutils.spawn import find_executable
import fileinput
import os
from pathlib import Path
import platform
import pyfastx
from random import randrange
from subprocess import run
import sys
from tempfile import TemporaryDirectory


AUTHOR = "Anicet Ebou <anicet.ebou@gmail.com>"
URL = "https://github.com/Ebedthan/butler"
VERSION = "0.1"

# Define command-line arguments----------------------------------------------
parser = argparse.ArgumentParser(
    prog="butler",
    usage="butler [options] [<FILE>]",
    add_help=False,
)

parser.add_argument(
    "file",
    nargs="?",
    type=argparse.FileType("r"),
    default=sys.stdin,
    help=argparse.SUPPRESS,
)
parser.add_argument(
    "-o",
    nargs="?",
    type=argparse.FileType("w"),
    metavar="FILE",
    default=sys.stdout,
    help="output result to FILE [stdout]",
)
parser.add_argument(
    "-a", action="store_true", help="activate anonymous mode [false]"
)
parser.add_argument(
    "--faa",
    type=str,
    metavar="FILE",
    help="output matched proteins sequences to FILE",
)
parser.add_argument(
    "--fna",
    type=str,
    metavar="FILE",
    help="output matched DNA sequences to FILE",
)
parser.add_argument(
    "-t",
    type=int,
    metavar="INT",
    default=1,
    help="number of threads [1]",
)
parser.add_argument(
    "-q", action="store_true", help="decrease program verbosity"
)
parser.add_argument(
    "-v", "--version", action="version", version="%(prog)s " + f"{VERSION}"
)
parser.add_argument(
    "-h", "--help", action="help", help="show this help message and exit"
)
parser.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)

args = parser.parse_args()


def main():
    # Preparing the environment -----------------------------------------------
    # Record the program start time
    startime = datetime.now()

    # Create a temporary directory
    butler_temp = TemporaryDirectory()

    # Handling db directory path
    # Are we in a docker file ? If yes the ENV variable IS_DOCKER is True
    dbdir = ""
    try:
        # get env var telling if we are in docker
        is_docker = os.environ["IS_DOCKER"]
    except KeyError:
        is_docker = False

    if is_docker:
        # path to hmm db
        dbdir = "/usr/local/lib/python3.8/dist-packages/db"

    else:
        # find path to butler.py
        bindir = Path(__file__).resolve().parent
        bindir = bindir.parent

        try:
            dbdir = (
                # case when app is installed by cloned repo
                Path(bindir, "db")
                # case when path is set by env variable
                or os.environ["BUTLERDB"]
                # case when app is installed through pip
                or "/usr/local/lib/python3.8/dist-packages/db"
            )
        except KeyError:
            print(
                "Error: HMMs models not found in $PATH. "
                + "Please set BUTLERDB environment variable to the path "
                + "where models are stored. "
                + "Visit https://github.com/Ebedthan/butler "
                + "for more informations.",
                file=sys.stderr,
            )
            sys.exit(1)

    # Get current user name
    try:
        user = os.environ["USER"]
    except KeyError:
        user = "not telling me who you are"

    # Handling CLI arguments --------------------------------------------------
    # Handling input file supply
    from_stdin = False
    if args.file.name == "<stdin>":
        from_stdin = True
        binput = open(Path(butler_temp.name, "butler_input.fa"), "w")
        for line in fileinput.input(files="-"):
            binput.write(line)
        binput.close()
    else:
        infile_type = get_file_type(args.file.name)
        if infile_type not in ["gz", "unknown"]:
            err(
                f"Input file compression ({infile_type})"
                + " is not currently supported"
            )
            sys.exit(1)

    # Check Alphabet of supplied sequences
    fa = ""
    is_prot = False
    if args.file.name == "<stdin>":
        fa = pyfastx.Fasta(str(Path(butler_temp.name, "butler_input.fa")))
    else:
        fa = pyfastx.Fasta(args.file.name)

    if fa.type not in ["DNA", "protein"]:
        err("Supplied file should be a DNA or protein file")
        sys.exit(1)

    if fa.type == "protein":
        is_prot = True

    if is_prot and args.fna:
        err(
            "Cannot retrieve matching DNA sequences from protein file. "
            + "Please remove --fna option"
        )
        sys.exit(1)

    # Start program ---------------------------------------------------------
    msg(f"This is butler {VERSION}")
    msg(f"Written by {AUTHOR}")
    msg(f"Available at {URL}")
    msg(f"Localtime is {datetime.now().strftime('%H:%M:%S')}")
    msg(f"You are {user}")
    msg(f"Operating system is {platform.system()}")
    if args.a:
        msg("You are running butler in anonymous mode")
    else:
        msg("You are running butler in normal mode")

    # Handling number of threads -----------------------------------------------
    cpus = args.t
    available_cpus = os.cpu_count()
    msg(f"System has {available_cpus} cores")

    if args.t == 0:
        cpus = available_cpus
    elif args.t > available_cpus:
        warn(
            f"Option -t asked for {args.cpus} threads,"
            + f" but system has only {available_cpus}."
        )
        cpus = available_cpus
    msg(f"We will use maximum of {cpus} cores")

    # Verify presence of needed tools ---------------------------------------
    needed_tools = ("hmmsearch", "prodigal")
    for tool in needed_tools:
        if find_executable(tool) is not None:
            msg(f"Found {tool}")
        else:
            print_install_tool(tool)
            sys.exit(1)

    # Running tools -----------------------------------------------------------
    # Predict protein-coding genes
    if not is_prot:
        msg("Predicting protein-coding genes")
        if args.a:
            if from_stdin:
                run(
                    [
                        "prodigal",
                        "-q",
                        "-p",
                        "anon",
                        "-o",
                        Path(butler_temp.name, "mygenes"),
                        "-i",
                        Path(butler_temp.name, "butler_input.fa"),
                        "-a",
                        Path(butler_temp.name, "my.proteins.faa"),
                    ]
                )
            else:
                run(
                    [
                        "prodigal",
                        "-q",
                        "-p",
                        "anon",
                        "-o",
                        Path(butler_temp.name, "mygenes"),
                        "-i",
                        str(args.file.name),
                        "-a",
                        Path(butler_temp.name, "my.proteins.faa"),
                    ]
                )

        else:
            if from_stdin:
                run(
                    [
                        "prodigal",
                        "-q",
                        "-o",
                        Path(butler_temp.name, "mygenes"),
                        "-i",
                        Path(butler_temp.name, "butler_input.fa"),
                        "-a",
                        Path(butler_temp.name, "my.proteins.faa"),
                    ]
                )
            else:
                run(
                    [
                        "prodigal",
                        "-q",
                        "-o",
                        Path(butler_temp.name, "mygenes"),
                        "-i",
                        str(args.file.name),
                        "-a",
                        Path(butler_temp.name, "my.proteins.faa"),
                    ]
                )

    # Classifying sequences into housekkeping genes
    msg("classifying sequences into housekeeping genes")
    if not is_prot:
        run(
            [
                "hmmsearch",
                "--cut_ga",
                "--cpu",
                str(cpus),
                "--noali",
                "-o",
                Path(butler_temp.name, "butler.hmmsearch"),
                Path(dbdir, "butler.hmm"),
                Path(butler_temp.name, "my.proteins.faa"),
            ]
        )
    else:
        if from_stdin:
            run(
                [
                    "hmmsearch",
                    "--cut_ga",
                    "--cpu",
                    str(cpus),
                    "--noali",
                    "-o",
                    Path(butler_temp.name, "butler.hmmsearch"),
                    Path(dbdir, "butler.hmm"),
                    Path(butler_temp.name, "butler_input.fa"),
                ]
            )
        else:
            run(
                [
                    "hmmsearch",
                    "--cut_ga",
                    "--cpu",
                    str(cpus),
                    "--noali",
                    "-o",
                    Path(butler_temp.name, "butler.hmmsearch"),
                    Path(dbdir, "butler.hmm"),
                    str(args.file.name),
                ]
            )

    hmmdict = defaultdict(lambda: defaultdict(list))

    # Second iteration over output file to get evalues and hsps
    with open(Path(butler_temp.name, "butler.hmmsearch")) as hmmfile:
        for record in SearchIO.parse(hmmfile, "hmmer3-text"):
            hits = record.hits
            for hit in hits:
                for _ in hit.hsps:
                    hmmdict[hit.id][record.id].extend(
                        [hit.evalue, hit.bitscore]
                    )
    hmmfile.close()

    classif = get_best_match(hmmdict)

    # Write output ------------------------------------------------------------
    msg("Writing output")
    ofile = sys.stdout
    if args.o.name != "<stdout>":
        ofile = open(args.o.name, "w")

    print("seq\tgene\te-value\tbitscore", file=ofile)
    for k, v in classif.items():
        res = v.split("#")
        print(f"{k}\t{res[0]}\t{res[1]}\t{res[2]}", file=ofile)

    # Get matched proteins
    if args.faa or args.fna:
        prots = ""
        nl = "\n"

        if not is_prot:
            prots = pyfastx.Fasta(
                str(Path(butler_temp.name, "my.proteins.faa"))
            )
        else:
            prots = pyfastx.Fasta(args.file.name)

        if args.faa:
            msg("Writing out predicted proteins sequences")
            oprots = open(args.faa, "w")
            for k, v in classif.items():
                seq = [
                    prots[k].seq[i : i + 60]  # noqa: E203
                    for i in range(0, len(prots[k]), 60)
                ]

                oprots.write(
                    f">{k}_gene={v.split('#')[0]}\n{nl.join(map(str, seq))}\n"
                )

        if args.fna:
            msg("Writing out predicted DNA sequences")
            ofna = open(args.fna, "w")
            mg = ""

            if from_stdin:
                mg = pyfastx.Fasta(
                    str(Path(butler_temp.name, "butler_input.fa"))
                )
            else:
                mg = pyfastx.Fasta(args.file.name)

            for k, v in classif.items():
                data = prots[k].description.split("#")[:4]
                seq = []
                newk = "_".join(k.split("_")[:-1])
                ss = int(data[3].strip())
                if ss == 1:
                    subseq = mg.fetch(
                        newk, (int(data[1].strip()), int(data[2].strip()))
                    )
                    seq = [
                        subseq[i : i + 60]  # noqa: E203
                        for i in range(0, len(subseq), 60)
                    ]
                elif ss == -1:
                    subseq = mg.fetch(
                        newk,
                        (int(data[1].strip()), int(data[2].strip())),
                        strand="-",
                    )
                    seq = [
                        subseq[i : i + 60]  # noqa: E203
                        for i in range(0, len(subseq), 60)
                    ]

                ofna.write(
                    f">{k}_gene={v.split('#')[0]}\n{nl.join(map(str, seq))}\n"
                )

    # Cleaning around ---------------------------------------------------------
    if not from_stdin:
        try:
            os.remove(f"{args.file.name}.fxi")
        except OSError:
            pass

    msg("Task finished successfully")
    msg(f"Walltime used (hh:mm:ss.ms): {elapsed_since(startime)}")
    if randrange(0, 100000) % 2:
        msg("Nice to have you. Share, enjoy and come back!")
    else:
        msg("Thanks you, come again.")


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
            "prodigal not found. Please visit"
            + "https://github.com/hyattpd/Prodigal"
        )


def elapsed_since(start):
    walltime = datetime.now() - start
    return walltime


def err(text):
    print(f"Error: {text}", file=sys.stderr)


def msg(text):
    """
    msg produces nice message and info output on terminal.

    :text: Message to print to STDOUT.
    """
    if not args.q:
        print(f"[{datetime.now().strftime('%H:%M:%S')}][INFO] {text}")


def warn(text):
    if not args.q:
        print(f"[{datetime.now().strftime('%H:%M:%S')}][WARN] {text}")


def exception_handler(
    exception_type,
    exception,
    traceback,
    debug_hook=sys.excepthook,
):
    """
    exception_handler remove default debug info and traceback
    from python output on command line. Use program --debug
    option to re-enable default behaviour.
    """

    if args.debug:
        debug_hook(exception_type, exception, traceback)
    else:
        print(f"{exception_type.__name__}, {exception}")


sys.excepthook = exception_handler

if __name__ == "__main__":
    main()
