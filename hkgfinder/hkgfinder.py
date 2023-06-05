#!/usr/bin/env python3

# Copyright 2022-2023 Anicet Ebou.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed except according
# to those terms.

import argparse
import logging
import os
import platform
import sys
import typing
import warnings
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
from Bio.Seq import reverse_complement, translate

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

if args.q:
    QUIETNESS_LEVEL = logging.CRITICAL
else:
    QUIETNESS_LEVEL = logging.INFO

logging.basicConfig(
    format="[%(asctime)s][%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
    level=QUIETNESS_LEVEL,
)

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


def main():
    # Preparing the environment -----------------------------------------------
    # Record the program start time
    startime = datetime.now()

    # Get current user name
    try:
        user = os.environ["USER"]
    except KeyError:
        user = "not telling me who you are"
    # Handling CLI arguments ---------------------------------------------
    if args.s:
        if not args.faa and not args.fna:
            logging.warning("Unecessary -s option in absence of --faa or --fna")

    # Create a temporary directory
    with TemporaryDirectory() as tmpdir:
        # Handling input file supply
        if args.file.name == "<stdin>":
            with open(
                Path(tmpdir, "hkgfinder_input.fa"),
                "w",
                encoding="utf-8",
            ) as binput:
                with xphyle.xopen(
                    xphyle.paths.STDIN, context_wrapper=True
                ) as infile:
                    for line in infile:
                        binput.write(line)  # type: ignore
        else:
            with open(
                Path(tmpdir, "hkgfinder_input.fa"),
                "w",
                encoding="utf-8",
            ) as binput:
                with xphyle.xopen(
                    str(args.file.name), context_wrapper=True
                ) as infile:
                    for line in infile:
                        binput.write(line)  # type: ignore

        # Check Alphabet of supplied sequences
        if args.file.name == "<stdin>":
            fasta = pyfastx.Fasta(str(Path(tmpdir, "hkgfinder_input.fa")))
        else:
            fasta = pyfastx.Fasta(args.file.name)

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
                "Cannot retrieve matching DNA sequences from protein file. Remove --fna"
            )
            sys.exit(1)

        if len(fasta.longest) >= 10000 and not args.g:
            logging.error(
                "Sequence length greater than 10000 bp. Do you want genome mode?"
            )
            sys.exit(1)

        # Check if fasta file does not contain duplicate sequences
        # which would break hmmsearch
        ids = fasta.keys()
        if len(ids) != len(set(ids)):
            logging.error("Supplied FASTA file contains duplicate sequences.")
            try:
                os.remove(f"{args.file.name}.fxi")
            except OSError:
                pass
            sys.exit(1)

        # Start program ---------------------------------------------------------
        logging.info("This is hkgfinder %s", VERSION)
        logging.info("Written by %s", AUTHOR)
        logging.info("Available at %s", URL)
        logging.info("Localtime is %s", datetime.now().strftime("%H:%M:%S"))
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
        cpus = 1
        available_cpus = os.cpu_count()
        if args.t == 0:
            cpus = available_cpus
        elif args.t > available_cpus:
            logging.warning(
                "Option -t asked for %s threads but system has only %s",
                args.t,
                available_cpus,
            )
            cpus = available_cpus
        else:
            cpus = args.t

        logging.info("We will use maximum of %s cores", cpus)

        # Running tools -----------------------------------------------------------
        # In genome mode ----------------------------------------------------------
        if args.g:
            logging.info("Predicting protein-coding genes")
            records = SeqIO.index(
                str(Path(tmpdir, "hkgfinder_input.fa")), "fasta"
            )
            orf_finder = pyrodigal.OrfFinder()
            orf_finder.train(
                *(bytes(record.seq) for record in records.values())
            )
            with open(
                Path(tmpdir, "my.proteins.faa"), "w+", encoding="utf-8"
            ) as prot_file:
                for seq in records.values():
                    prediction = orf_finder.find_genes(bytes(seq.seq))
                    prediction.write_translations(prot_file, seq.id)

            if args.fna:
                with open(
                    Path(tmpdir, "my.dna.fna"), "w+", encoding="utf-8"
                ) as dna_file:
                    for seq in records.values():
                        prediction = orf_finder.find_genes(bytes(seq.seq))
                        prediction.write_genes(dna_file, seq.id)

        # In metagenome mode ------------------------------------------------------
        elif args.m:
            logging.info("Predicting protein-coding genes")
            records = SeqIO.index(
                str(Path(tmpdir, "hkgfinder_input.fa")), "fasta"
            )
            orf_finder = pyrodigal.OrfFinder(meta=True)
            with open(
                Path(tmpdir, "my.proteins.faa"), "w+", encoding="utf-8"
            ) as prot_file:
                for seq in records.values():
                    for i, prediction in enumerate(
                        orf_finder.find_genes(bytes(seq.seq))
                    ):
                        prot = prediction.translate()
                        prot_file.write(f">{seq.id}_{i+1}\n{prot}\n")

            if args.fna:
                with open(
                    Path(tmpdir, "my.dna.fna"), "w+", encoding="utf-8"
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
                logging.info("Translating sequences into 6 frames")

                # translate sequences
                do_translation(
                    str(Path(tmpdir, "hkgfinder_input.fa")),
                    Path(tmpdir, "input_translate.fa"),
                )
                fa = pyfastx.Fasta(str(Path(tmpdir, "input_translate.fa")))
                if len(fa.longest) >= 10000 and not args.g:
                    logging.error(
                        "Seq length greater than 10000 bp. Do you want genome mode?"
                    )
                    sys.exit(1)

        # Classifying sequences into housekeeping genes
        logging.info("Classifying sequences into housekeeping genes")
        if args.g or args.m:
            seqdata = Path(tmpdir, "my.proteins.faa")
        else:
            if fatype == "protein":
                seqdata = Path(tmpdir, "hkgfinder_input.fa")
            else:
                seqdata = Path(tmpdir, "input_translate.fa")

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
                    logging.info("Pre-fetching targets into memory")
                    targets = seq_file.read_block()
                    mem = sys.getsizeof(targets) + sum(
                        sys.getsizeof(target) for target in targets
                    )
                    mem = mem / 1024
                    logging.info(
                        "Database in-memory size: %s KiB", f"{mem:.1f}"
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
                    hmm_desc = HMMDESC[str(hmm_id, encoding="utf-8")]  # type: ignore
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
                    # logging.info("Parsing HMM result", is_quiet)
                    for result in results:
                        if result.hit_id in best_results:
                            previous_bitscore = best_results[
                                result.hit_id
                            ].bitscore
                            if result.bitscore > previous_bitscore:
                                best_results[result.hit_id] = result
                        else:
                            best_results[result.hit_id] = result

        # Write output ------------------------------------------------------------
        filtered_results = [best_results[k] for k in sorted(best_results)]

        logging.info("Writing output")
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
                prots = pyfastx.Fasta(str(Path(tmpdir, "my.proteins.faa")))
            else:
                if fatype == "protein":
                    prots = fasta
                else:
                    prots = pyfastx.Fasta(
                        str(Path(tmpdir, "input_translate.fa"))
                    )

            if args.s:
                sorted(filtered_results, key=attrgetter("hit_id"))

            if args.faa:
                logging.info("Writing out predicted proteins sequences")

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
                logging.info("Writing out predicted DNA sequences")
                if args.g or args.m:
                    dna_file = pyfastx.Fasta(str(Path(tmpdir, "my.dna.fna")))
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
                                for i in range(
                                    0, len(dna_file[result.hit_id]), 60
                                )
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
                                for i in range(
                                    0, len(dna_file[result.hit_id]), 60
                                )
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
        logging.info("Task finished successfully")
        logging.info(
            "Walltime used (hh:mm:ss.ms): %s", datetime.now() - startime
        )
        if randrange(0, 100000) % 2:
            logging.info("Nice to have you. Share, enjoy and come back!")
        else:
            logging.info("Thanks you, come again.")


if __name__ == "__main__":
    main()
