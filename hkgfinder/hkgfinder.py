#!/usr/bin/env python3

# Copyright 2022-2023 Anicet Ebou.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed except according
# to those terms.

import argparse
import os
import platform
import sys
from collections import namedtuple
from datetime import datetime
from operator import attrgetter
from pathlib import Path
from random import randrange
from tempfile import TemporaryDirectory

import psutil
import pyfastx
import pyhmmer
import pyrodigal
import xphyle
import xphyle.paths
from Bio import SeqIO

from . import utils

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
        binput = open(
            Path(hkgfinder_temp.name, "hkgfinder_input.fa"),
            "w",
            encoding="utf-8",
        )
        with xphyle.xopen(xphyle.paths.STDIN, context_wrapper=True) as infile:
            for line in infile:
                binput.write(line)  # type: ignore
    else:
        binput = open(
            Path(hkgfinder_temp.name, "hkgfinder_input.fa"),
            "w",
            encoding="utf-8",
        )
        with xphyle.xopen(str(args.file.name), context_wrapper=True) as infile:
            for line in infile:
                binput.write(line)  # type: ignore

    # Check Alphabet of supplied sequences
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

    if len(fa.longest) >= 10000 and not args.g:
        utils.err(
            "Sequence length greater than 10000 bp. Do you want genome mode?"
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
    available_cpus = os.cpu_count()
    utils.msg(f"System has {available_cpus} cores", is_quiet)
    cpus = utils.parse_cpu(args.t, is_quiet)

    utils.msg(f"We will use maximum of {cpus} cores", is_quiet)

    # Running tools -----------------------------------------------------------
    # In genome mode ----------------------------------------------------------
    if args.g:
        utils.msg("Predicting protein-coding genes", is_quiet)
        records = SeqIO.index(
            str(Path(hkgfinder_temp.name, "hkgfinder_input.fa")), "fasta"
        )
        orf_finder = pyrodigal.OrfFinder()
        orf_finder.train(*(bytes(record.seq) for record in records.values()))
        with open(
            Path(hkgfinder_temp.name, "my.proteins.faa"), "w+", encoding="utf-8"
        ) as prot_file:
            for seq in records.values():
                for i, prediction in enumerate(
                    orf_finder.find_genes(bytes(seq.seq))
                ):
                    prot = prediction.translate()
                    prot_file.write(f">{seq.id}_{i+1}\n{prot}\n")

        if args.fna:
            with open(
                Path(hkgfinder_temp.name, "my.dna.fna"), "w+", encoding="utf-8"
            ) as dna_file:
                for seq in records.values():
                    for i, prediction in enumerate(
                        orf_finder.find_genes(bytes(seq.seq))
                    ):
                        dna = prediction
                        dna_file.write(f">{seq.id}_{i+1}\n{dna}\n")

    # In metagenome mode ------------------------------------------------------
    elif args.m:
        utils.msg("Predicting protein-coding genes", is_quiet)
        records = SeqIO.index(
            str(Path(hkgfinder_temp.name, "hkgfinder_input.fa")), "fasta"
        )
        orf_finder = pyrodigal.OrfFinder(meta=True)
        with open(
            Path(hkgfinder_temp.name, "my.proteins.faa"), "w+", encoding="utf-8"
        ) as prot_file:
            for seq in records.values():
                for i, prediction in enumerate(
                    orf_finder.find_genes(bytes(seq.seq))
                ):
                    prot = prediction.translate()
                    prot_file.write(f">{seq.id}_{i+1}\n{prot}\n")

        if args.fna:
            with open(
                Path(hkgfinder_temp.name, "my.dna.fna"), "w+", encoding="utf-8"
            ) as dna_file:
                for seq in records.values():
                    for i, prediction in enumerate(
                        orf_finder.find_genes(bytes(seq.seq))
                    ):
                        dna = prediction
                        dna_file.write(f">{seq.id}_{i+1}\n{dna}\n")

    # In normal mode ----------------------------------------------------------
    else:
        if fatype == "DNA":
            utils.msg("Translating sequences into 6 frames", is_quiet)

            # translate sequences
            utils.do_translation(
                str(Path(hkgfinder_temp.name, "hkgfinder_input.fa")),
                Path(hkgfinder_temp.name, "input_translate.fa"),
            )
            fa = pyfastx.Fasta(
                str(Path(hkgfinder_temp.name, "input_translate.fa"))
            )
            if len(fa.longest) >= 10000 and not args.g:
                utils.err(
                    "Sequence length greater than 10000 bp. Do you want genome mode?"
                )
                sys.exit(1)

    # Classifying sequences into housekeeping genes
    utils.msg("Classifying sequences into housekeeping genes", is_quiet)
    if args.g or args.m:
        seqdata = Path(hkgfinder_temp.name, "my.proteins.faa")
    else:
        if fatype == "protein":
            seqdata = Path(hkgfinder_temp.name, "hkgfinder_input.fa")
        else:
            seqdata = Path(hkgfinder_temp.name, "input_translate.fa")

    results = []
    best_results = {}
    available_memory = psutil.virtual_memory().available
    target_size = os.stat(seqdata).st_size
    # importlib.resources.open_binary(__name__, "hkgfinder.hmm")
    bindir = Path(__file__).resolve().parent
    bindir = bindir.parent
    with pyhmmer.plan7.HMMFile(
        str(Path(bindir, "db", "hkgfinder.hmm"))
    ) as hmm_file:
        with pyhmmer.easel.SequenceFile(seqdata, digital=True) as seq_file:
            if target_size < available_memory * 0.1:
                utils.msg("Pre-fetching targets into memory", is_quiet)
                targets = seq_file.read_block()
                mem = sys.getsizeof(targets) + sum(
                    sys.getsizeof(target) for target in targets
                )
                mem = mem / 1024
                utils.msg(
                    "Database in-memory size: " + f"{mem:.1f} KiB",
                    is_quiet,
                )
            else:
                targets = seq_file

            HMMResult = namedtuple(
                "HMMResult",
                [
                    "hmm_id",
                    "hmm_desc",
                    "hit_id",
                    "evalue",
                    "bitscore",
                ],
            )

            for hits in pyhmmer.hmmer.hmmsearch(
                hmm_file,
                targets,  # type: ignore
                cpus=cpus,  # type: ignore
                bit_cutoffs="gathering",  # type: ignore
            ):
                hmm_id = hits.query_name
                hmm_desc = utils.get_hmm_desc(str(hmm_id, encoding="utf-8"))  # type: ignore
                for hit in hits:
                    if hit.included:
                        results.append(
                            HMMResult(
                                str(hmm_id, encoding="utf-8"),  # type: ignore
                                hmm_desc,
                                str(hit.name, encoding="utf-8"),
                                hit.evalue,  # type: ignore
                                hit.score,  # type: ignore
                            )
                        )

                # Second iteration over output file to get evalues and hsps
                # utils.msg("Parsing HMM result", is_quiet)
                for result in results:
                    if result.hit_id in best_results:
                        previous_bitscore = best_results[result.hit_id].bitscore
                        if result.bitscore > previous_bitscore:
                            best_results[result.hit_id] = result
                    else:
                        best_results[result.hit_id] = result

    # Write output ------------------------------------------------------------
    filtered_results = [best_results[k] for k in sorted(best_results)]

    utils.msg("Writing output", is_quiet)
    if args.o.name != "<stdout>":
        with open(args.o.name, "w", encoding="utf-8") as ofile:
            print(
                "seq_name\tpred_gene\tgene_desc\te-value\tbitscore",
                file=ofile,
            )
            for result in filtered_results:
                print(
                    result.hit_id,
                    result.hmm_id,
                    result.hmm_desc,
                    f"{result.evalue:.3f}",
                    f"{result.bitscore:.3f}",
                    sep="\t",
                    file=ofile,
                )
    else:
        print(
            "seq_name\tseq_desc\tpred_gene\tgene_desc\te-value\tbitscore",
        )
        for result in filtered_results:
            print(
                result.hit_id,
                result.hmm_id,
                result.hmm_desc,
                f"{result.evalue:.3f}",
                f"{result.bitscore:.3f}",
                sep="\t",
            )
    # Get matched proteins
    if args.faa or args.fna:
        newline = "\n"

        if args.g or args.m:
            prots = pyfastx.Fasta(
                str(Path(hkgfinder_temp.name, "my.proteins.faa"))
            )
        else:
            if fatype == "protein":
                prots = fa
            else:
                prots = pyfastx.Fasta(
                    str(Path(hkgfinder_temp.name, "input_translate.fa"))
                )

        if args.s:
            sorted(filtered_results, key=attrgetter("hit_id"))

        if args.faa:
            utils.msg("Writing out predicted proteins sequences", is_quiet)

            if args.s:
                for result in filtered_results:
                    with open(
                        f"{os.path.splitext(args.faa)[0]}_{result.hmm_id}.faa",
                        "a",
                        encoding="utf-8",
                    ) as out:
                        seq = [
                            prots[result.hit_id].seq[i : i + 60]
                            for i in range(0, len(prots[result.hit_id]), 60)
                        ]
                        out.write(
                            f">{result.hit_id}_gene={result.hmm_id}"
                            + f"\n{newline.join(map(str, seq))}\n"
                        )
            else:
                with open(
                    f"{os.path.splitext(args.faa)[0]}.faa",
                    "w",
                    encoding="utf-8",
                ) as out:
                    for result in filtered_results:
                        seq = [
                            prots[result.hit_id].seq[i : i + 60]
                            for i in range(0, len(prots[result.hit_id]), 60)
                        ]
                        out.write(
                            f">{result.hit_id}_gene={result.hmm_id}\n"
                            + f"{newline.join(map(str, seq))}\n"
                        )
        if args.fna:
            utils.msg("Writing out predicted DNA sequences", is_quiet)
            if args.g or args.m:
                dna_file = pyfastx.Fasta(
                    str(Path(hkgfinder_temp.name, "my.dna.fna"))
                )
            else:
                dna_file = pyfastx.Fasta(args.file)

            if args.s:
                for result in filtered_results:
                    with open(
                        f"{os.path.splitext(args.faa)[0]}_{result.hmm_id}.fna",
                        "a",
                        encoding="utf-8",
                    ) as out:
                        seq = [
                            dna_file[result.hit_id].seq[i : i + 60]
                            for i in range(0, len(dna_file[result.hit_id]), 60)
                        ]
                        out.write(
                            f">{result.hit_id}_gene={result.hmm_id}\n"
                            + f"{newline.join(map(str, seq))}\n"
                        )
            else:
                with open(
                    f"{os.path.splitext(args.faa)[0]}.fna",
                    "a",
                    encoding="utf-8",
                ) as out:
                    for result in filtered_results:
                        seq = [
                            dna_file[result.hit_id].seq[i : i + 60]
                            for i in range(0, len(dna_file[result.hit_id]), 60)
                        ]
                        out.write(
                            f">{result.hit_id}_gene={result.hmm_id}\n"
                            + f"{newline.join(map(str, seq))}\n"
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
