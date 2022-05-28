# Copyright 2022 Anicet Ebou.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed except according
# to those terms.

from datetime import datetime
import sys


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
