#!/usr/bin/env python3

# Copyright 2022-2023 Anicet Ebou.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed except according
# to those terms.

import argparse
from Bio import SearchIO
from collections import defaultdict
from datetime import datetime
from distutils.spawn import find_executable
import os
from pathlib import Path
import platform
import pyfastx
from random import randrange
from subprocess import run
import sys
from tempfile import TemporaryDirectory
from . import utils
import xphyle
import xphyle.paths


AUTHOR = "Anicet Ebou <anicet.ebou@gmail.com>"
URL = "https://github.com/Ebedthan/hkgfinder"
VERSION = "0.1"

# Define command-line arguments----------------------------------------------
parser = argparse.ArgumentParser(
    prog="hkgfinder",
    usage="hkgfinder [options] [<FILE>]",
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
group = parser.add_mutually_exclusive_group()
group.add_argument(
    "-g",
    action="store_true",
    help="activate genome mode [false]",
)
group.add_argument(
    "-m",
    action="store_true",
    help="activate metagenome mode [false]",
)
parser.add_argument(
    "--faa",
    type=str,
    metavar="FILE",
    help="output matched protein sequences to FILE [false]",
)
parser.add_argument(
    "--fna",
    type=str,
    metavar="FILE",
    help="output matched DNA sequences to FILE [false]",
)
parser.add_argument(
    "-s",
    action="store_true",
    help="output sequences in file by gene [false]",
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
    hkgfinder_temp = TemporaryDirectory()

    # Is hkgfinder ran in quiet mode?
    is_quiet = False
    if args.q:
        is_quiet = True

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
        # find path to hkgfinder.py
        bindir = Path(__file__).resolve().parent
        bindir = bindir.parent

        try:
            dbdir = (
                # case when app is installed by cloned repo
                Path(bindir, "db")
                # case when app is installed through pip
                or "/usr/local/lib/python3.8/dist-packages/db"
                # case when path is set by env variable
                or os.environ["HKGDB"]
            )
        except KeyError:
            print(
                "Error: HMMs models not found in $PATH. "
                + "Please set HKGDB environment variable to the path "
                + "where models are stored. "
                + "Visit https://github.com/Ebedthan/hkgfinder "
                + "for more informations.",
                file=sys.stderr,
            )
            sys.exit(1)

        # Get current user name
    try:
        user = os.environ["USER"]
    except KeyError:
        user = "not telling me who you are"

    # Handling CLI arguments ---------------------------------------------
    if args.s:
        if not args.faa and not args.fna:
            utils.warn(
                "Unecessary -s option in absence of --faa or --fna", is_quiet
            )

    # Handling input file supply
    if args.file.name == "<stdin>":
        binput = open(Path(hkgfinder_temp.name, "hkgfinder_input.fa"), "w")
        with xphyle.xopen(xphyle.paths.STDIN, context_wrapper=True) as infile:
            for line in infile:
                binput.write(line)  # type: ignore
    else:
        binput = open(Path(hkgfinder_temp.name, "hkgfinder_input.fa"), "w")
        with xphyle.xopen(str(args.file.name), context_wrapper=True) as infile:
            for line in infile:
                binput.write(line)  # type: ignore

    # Check Alphabet of supplied sequences
    fa = ""
    if args.file.name == "<stdin>":
        fa = pyfastx.Fasta(str(Path(hkgfinder_temp.name, "hkgfinder_input.fa")))
    else:
        fa = pyfastx.Fasta(args.file.name)

    fatype = utils.is_dna_or_aa(fa[0].seq)

    if args.g and fatype == "protein":
        utils.err("Cannot run genome mode on proteins")
        sys.exit(1)

    if args.m and fatype == "protein":
        utils.err("Cannot run metagenome mode on proteins")
        sys.exit(1)

    if fatype not in ["DNA", "protein"]:
        utils.err("Supplied file should be a DNA or protein file")
        sys.exit(1)

    if fatype == "protein" and args.fna:
        utils.err(
            "Cannot retrieve matching DNA sequences from protein file. "
            + "Please remove --fna option"
        )
        sys.exit(1)

    # Check if fasta file does not contain duplicate sequences
    # which would break hmmsearch
    ids = fa.keys()
    if len(ids) != len(set(ids)):
        utils.err(
            "Supplied FASTA file contains duplicate sequences. "
            + "Please remove them before launching hkgfinder."
        )
        try:
            os.remove(f"{args.file.name}.fxi")
        except OSError:
            pass
        sys.exit(1)

    # Start program ---------------------------------------------------------
    utils.msg(f"This is hkgfinder {VERSION}", is_quiet)
    utils.msg(f"Written by {AUTHOR}", is_quiet)
    utils.msg(f"Available at {URL}", is_quiet)
    utils.msg(f"Localtime is {datetime.now().strftime('%H:%M:%S')}", is_quiet)
    utils.msg(f"You are {user}", is_quiet)
    utils.msg(f"Operating system is {platform.system()}", is_quiet)
    if args.g:
        utils.msg("You are running hkgfinder in genome mode", is_quiet)
    elif args.m:
        utils.msg("You are running hkgfinder in metagenome mode", is_quiet)
    else:
        utils.msg("You are running kgfinder in normal mode", is_quiet)

    # Handling number of threads -----------------------------------------------
    cpus = args.t
    available_cpus = os.cpu_count()
    utils.msg(f"System has {available_cpus} cores", is_quiet)
    if args.t == 0:
        cpus = available_cpus
    elif args.t > available_cpus:
        utils.warn(
            f"Option -t asked for {args.cpus} threads,"
            + f" but system has only {available_cpus}.",
            is_quiet,
        )
        cpus = available_cpus
    utils.msg(f"We will use maximum of {cpus} cores", is_quiet)

    # Verify presence of needed tools ---------------------------------------
    needed_tools = ("hmmsearch", "prodigal")
    for tool in needed_tools:
        if find_executable(tool) is not None:
            utils.msg(f"Found {tool}", is_quiet)
        else:
            utils.print_install_tool(tool)
            sys.exit(1)

    # Running tools -----------------------------------------------------------
    # In genome mode ----------------------------------------------------------
    if args.g:
        utils.msg("Predicting protein-coding genes", is_quiet)
        run(
            [
                "prodigal",
                "-q",
                "-o",
                Path(hkgfinder_temp.name, "mygenes"),
                "-i",
                Path(hkgfinder_temp.name, "hkgfinder_input.fa"),
                "-a",
                Path(hkgfinder_temp.name, "my.proteins.faa"),
            ]
        )

    # In metagenome mode ------------------------------------------------------
    elif args.m:
        utils.msg("Predicting protein-coding genes", is_quiet)
        run(
            [
                "prodigal",
                "-q",
                "-p",
                "anon",
                "-o",
                Path(hkgfinder_temp.name, "mygenes"),
                "-i",
                Path(hkgfinder_temp.name, "hkgfinder_input.fa"),
                "-a",
                Path(hkgfinder_temp.name, "my.proteins.faa"),
            ]
        )

    # In normal mode ----------------------------------------------------------
    else:
        if fatype == "DNA":
            utils.msg("Translating sequences into 6 frames", is_quiet)
            fa = pyfastx.Fasta(
                str(Path(hkgfinder_temp.name, "hkgfinder_input.fa"))
            )

            # translate sequences
            utils.do_translation(
                fa, Path(hkgfinder_temp.name, "input_translate.fa")
            )

    # Classifying sequences into housekkeping genes
    utils.msg("Classifying sequences into housekeeping genes", is_quiet)
    if args.g or args.m:
        run(
            [
                "hmmsearch",
                "--cut_ga",
                "--cpu",
                str(cpus),
                "--noali",
                "-o",
                Path(hkgfinder_temp.name, "hkgfinder.hmmsearch"),
                Path(dbdir, "hkgfinder.hmm"),
                Path(hkgfinder_temp.name, "my.proteins.faa"),
            ]
        )
    else:
        if fatype == "protein":
            run(
                [
                    "hmmsearch",
                    "--cut_ga",
                    "--cpu",
                    str(cpus),
                    "--noali",
                    "-o",
                    Path(hkgfinder_temp.name, "hkgfinder.hmmsearch"),
                    Path(dbdir, "hkgfinder.hmm"),
                    Path(hkgfinder_temp.name, "hkgfinder_input.fa"),
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
                    Path(hkgfinder_temp.name, "hkgfinder.hmmsearch"),
                    Path(dbdir, "hkgfinder.hmm"),
                    Path(hkgfinder_temp.name, "input_translate.fa"),
                ]
            )

    hmmdict = defaultdict(lambda: defaultdict(list))

    # Second iteration over output file to get evalues and hsps
    utils.msg("Parsing HMM result", is_quiet)

    hmmresult = SearchIO.parse(
        Path(hkgfinder_temp.name, "hkgfinder.hmmsearch"), "hmmer3-text"
    )
    for record in hmmresult:
        hits = record.hits
        for hit in hits:
            for _ in hit.hsps:
                hmmdict[f"{hit.id}#{hit.description}"][
                    f"{record.id}#{record.description}"
                ].extend([hit.evalue, hit.bitscore])

    classif = utils.get_best_match(hmmdict)

    # Write output ------------------------------------------------------------
    utils.msg("Writing output", is_quiet)
    ofile = sys.stdout
    if args.o.name != "<stdout>":
        ofile = open(args.o.name, "w")
    print(
        "seq_name\tseq_desc\tpred_gene\tgene_desc\te-value\tbitscore",
        file=ofile,
    )
    for k, v in classif.items():
        res = v.split("#")
        seqinfo = k.split("#", 1)
        print(
            f"{seqinfo[0]}\t{seqinfo[1]}"
            + f"\t{res[0]}\t{res[1]}\t{res[2]}\t{res[3]}",
            file=ofile,
        )
    # Get matched proteins
    if args.faa or args.fna:
        genes = {}
        prots = ""
        nl = "\n"
        if args.p:
            prots = pyfastx.Fasta(
                str(Path(hkgfinder_temp.name, "my.proteins.faa"))
            )
        else:
            if fatype == "protein":
                prots = pyfastx.Fasta(args.file.name)
            else:
                prots = fa

        if args.s:
            for k, v in classif.items():
                mygene = v.split("#")[0]
                if mygene in genes:
                    genes[mygene].append(k)
                else:
                    genes[v.split("#")[0]] = [k]

        if args.faa:
            utils.msg("Writing out predicted proteins sequences", is_quiet)

            if args.s:
                for gene in genes.keys():
                    oprots = open(
                        f"{os.path.splitext(args.faa)}_{gene}.faa", "w"
                    )
                    for s in genes[gene]:
                        seq = [
                            prots[s].seq[i : i + 60]  # noqa: E203
                            for i in range(0, len(prots[s]), 60)
                        ]
                        oprots.write(
                            f">{s}_gene={gene}\n{nl.join(map(str, seq))}\n"
                        )
            else:
                oprots = open(f"{os.path.splitext(args.faa)}.faa", "w")
                for k, v in classif.items():
                    seq = [
                        prots[k].seq[i : i + 60]  # noqa: E203
                        for i in range(0, len(prots[k]), 60)
                    ]
                    oprots.write(
                        f">{k}_gene={v.split('#')[0]}\n"
                        + f"{nl.join(map(str, seq))}\n"
                    )

        if args.fna:
            utils.msg("Writing out predicted DNA sequences", is_quiet)

            ofna = open(args.fna, "w")
            mg = pyfastx.Fasta(
                str(Path(hkgfinder_temp.name, "hkgfinder_input.fa"))
            )

            for gene in genes.keys():
                ofna = open(f"{os.path.splitext(args.faa)}_{gene}.faa", "w")
                for s in genes[gene]:
                    data = prots[s].description.split("#")[:4]
                    seq = []
                    newk = "_".join(s.split("_")[:-1])
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

                    ofna.write(f">{s}_gene={gene}\n{nl.join(map(str, seq))}\n")

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
    try:
        os.remove(f"{args.file.name}.fxi")
    except OSError:
        pass
    utils.msg("Task finished successfully", is_quiet)
    utils.msg(
        f"Walltime used (hh:mm:ss.ms): {utils.elapsed_since(startime)}",
        is_quiet,
    )
    if randrange(0, 100000) % 2:
        utils.msg("Nice to have you. Share, enjoy and come back!", is_quiet)
    else:
        utils.msg("Thanks you, come again.", is_quiet)


if __name__ == "__main__":
    main()
