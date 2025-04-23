"""CLI program."""

import argparse
import sys

VERSION = "0.3"

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
parser.add_argument("-q", action="store_true", help="decrease program verbosity")
parser.add_argument("-d", action="store_true", help="enable debug mode")
parser.add_argument(
    "-v",
    "--version",
    action="version",
    version="%(prog)s " + f"{VERSION}",
)
parser.add_argument(
    "-h",
    "--help",
    action="help",
    help="show this help message and exit",
)
parser.add_argument("--debug", action="store_true", help=argparse.SUPPRESS)

args = parser.parse_args()
