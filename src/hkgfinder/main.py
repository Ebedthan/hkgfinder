#!/usr/bin/env python3

"""Main hkgfinder program."""

# Copyright 2022-2024 Anicet Ebou.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed except according
# to those terms.

import contextlib
import datetime
import logging
import os
import platform
import sys
from operator import attrgetter
from pathlib import Path
from random import randrange
from tempfile import TemporaryDirectory

import pyfastx

from . import cli, hkglib

AUTHOR = "Anicet Ebou <anicet.ebou@gmail.com>"
URL = "https://github.com/Ebedthan/hkgfinder"
VERSION = "0.3"

args = cli.args

MAX_SEQ_LENGTH = 10000

QUIETNESS_LEVEL = logging.CRITICAL if args.q else logging.INFO

logging.basicConfig(
    format="[%(asctime)s][%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
    level=QUIETNESS_LEVEL,
)


def main() -> None:
    """Contains the main hkgfinder program."""
    # Preparing the environment -----------------------------------------------
    # Record the program start time
    startime = datetime.datetime.now(tz=datetime.timezone.utc)

    # Get current user name
    try:
        user = os.environ["USER"]
    except KeyError:
        user = "not telling me who you are"
    # Handling CLI arguments ---------------------------------------------
    if args.s and not args.faa and not args.fna:
        logging.warning("Unecessary -s option in absence of --faa or --fna")

    # Create a temporary directory
    with TemporaryDirectory() as tmpdir:
        uncompressed_input = hkglib.decompress_file(str(args.file.name), tmpdir)

        # Check Alphabet of supplied sequences
        fasta = pyfastx.Fasta(str(uncompressed_input))
        fatype = fasta.type

        if args.g and fatype == "protein":
            logging.error("Cannot run genome mode on proteins")
            sys.exit(1)

        if args.m and fatype == "protein":
            logging.error("Cannot run metagenome mode on proteins")
            sys.exit(1)

        if fatype not in ["DNA", "protein"]:
            logging.error("Supplied file should be a DNA or protein file")
            sys.exit(1)

        if fatype == "protein" and args.fna:
            logging.error(
                "Cannot retrieve matching DNA sequences from protein file. Remove --fna",
            )
            sys.exit(1)

        if len(fasta.longest) >= MAX_SEQ_LENGTH and not args.g:
            logging.error(
                "Sequence length greater than %s bp. Do you want genome mode (-g)?",
                MAX_SEQ_LENGTH,
            )
            sys.exit(1)

        # Check if fasta file does not contain duplicate sequences
        # which would break hmmsearch
        ids = fasta.keys()
        if len(ids) != len(set(ids)):
            logging.error("Supplied FASTA file contains duplicate sequences.")
            with contextlib.suppress(OSError):
                Path.unlink(Path(f"{args.file.name}.fxi"))
            sys.exit(1)

        # Start program ---------------------------------------------------------
        logging.info("This is hkgfinder %s", VERSION)
        logging.info("Written by %s", AUTHOR)
        logging.info("Available at %s", URL)
        logging.info(
            "Localtime is %s",
            datetime.datetime.now(tz=datetime.timezone.utc).strftime("%H:%M:%S"),
        )
        logging.info("You are %s", user)
        logging.info("Operating system is %s", platform.system())
        if args.g:
            logging.info("You are running hkgfinder in genome mode")
        elif args.m:
            logging.info("You are running hkgfinder in metagenome mode")
        else:
            logging.info("You are running kgfinder in normal mode")

        # Handling number of threads -----------------------------------------------
        available_cpus = os.cpu_count()
        logging.info("System has %s cores", available_cpus)

        cpus = hkglib.handle_cpus(args.t, available_cpus)
        logging.info("We will use maximum of %s cores", cpus)

        # Running tools -----------------------------------------------------------
        logging.info("Predicting protein-coding genes")
        # In genome mode ----------------------------------------------------------
        if args.g:
            hkglib.find_protein_coding_genes(
                str(uncompressed_input),
                tmpdir,
                save_both=bool(args.fna),
            )

        # In metagenome mode ------------------------------------------------------
        elif args.m:
            hkglib.find_protein_coding_genes(
                str(uncompressed_input),
                tmpdir,
                is_meta=True,
                save_both=bool(args.fna),
            )

        # In normal mode ----------------------------------------------------------
        elif fatype == "DNA":
            logging.info("Translating sequences into 6 frames")

            # translate sequences
            hkglib.do_translation(
                str(uncompressed_input),
                Path(tmpdir, "input_translate.fa"),
            )
            fa = pyfastx.Fasta(str(Path(tmpdir, "input_translate.fa")))
            if len(fa.longest) >= MAX_SEQ_LENGTH and not args.g:
                logging.error(
                    "Seq length greater than 10000 bp. Do you want genome mode?",
                )
                sys.exit(1)

        # Classifying sequences into housekeeping genes
        logging.info("Classifying sequences into housekeeping genes")
        if args.g or args.m:
            seqdata = Path(tmpdir, "my.proteins.faa")
        elif fatype == "protein":
            seqdata = uncompressed_input
        else:
            seqdata = Path(tmpdir, "input_translate.fa")

        # Run HMM
        best_results = hkglib.search_hmm(Path(seqdata), cpus)

        # Write output ------------------------------------------------------------
        filtered_results = [best_results[k] for k in sorted(best_results)]

        logging.info("Writing output")
        hkglib.write_results(
            Path(args.o.name),
            filtered_results,
            to_stdout=(args.o.name == "<stdout>"),
        )

        # Get matched proteins
        if args.faa or args.fna:
            if args.g or args.m:
                prots = pyfastx.Fasta(str(Path(tmpdir, "my.proteins.faa")))
            elif fatype == "protein":
                prots = fasta
            else:
                prots = pyfastx.Fasta(str(Path(tmpdir, "input_translate.fa")))

            if args.s:
                sorted(filtered_results, key=attrgetter("hit_id"))

            if args.faa:
                logging.info("Writing out predicted proteins sequences")
                hkglib.write_sequences(
                    os.path.splitext(args.faa)[0],
                    filtered_results,
                    prots,
                    split=args.s,
                    is_prot=True,
                )

            if args.fna:
                logging.info("Writing out predicted DNA sequences")
                if args.g or args.m:
                    dnas = pyfastx.Fasta(str(Path(tmpdir, "my.dna.fna")))
                else:
                    dnas = pyfastx.Fasta(args.file)

                hkglib.write_sequences(
                    os.path.splitext(args.fna)[0],
                    filtered_results,
                    dnas,
                    split=args.s,
                    is_prot=False,
                )

        # Cleaning around ---------------------------------------------------------
        with contextlib.suppress(OSError):
            Path.unlink(Path(f"{args.file.name}.fxi"))
        logging.info("Task finished successfully")
        logging.info(
            "Walltime used (hh:mm:ss.ms): %s",
            datetime.datetime.now(tz=datetime.timezone.utc) - startime,
        )
        if randrange(0, 100000) % 2:  # noqa: S311
            logging.info("Nice to have you. Share, enjoy and come back!")
        else:
            logging.info("Thanks you, come again.")


if __name__ == "__main__":
    main()
