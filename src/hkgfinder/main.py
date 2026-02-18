#!/usr/bin/env python3

"""Main hkgfinder program."""

# Copyright 2022-2026 Anicet Ebou.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed except according
# to those terms.
import contextlib
import datetime
import logging
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

from pyfastxcli import pyfastx
from xopen import xopen

from hkgfinder.hmm_search import HMMSearcher
from hkgfinder.validation import InputValidator
from hkgfinder.writers import ResultWriter

from . import hkglib
from .cli import get_args
from .config import create_config


def main() -> None:
    """Main hkgfinder program."""
    args = get_args()

    # Create config from CLI args
    config = create_config(args)

    # Setup
    hkglib.setup_logging(args.q)
    startime = datetime.datetime.now(tz=datetime.timezone.utc)

    # determine run mode
    run_mode = hkglib.determine_run_mode(args)

    if args.s and not args.faa and not args.fna and not args.genbank:
        logging.warning(
            "Unecessary -s option in absence of --faa, --fna or --genbank"
        )

    try:
        with TemporaryDirectory() as tmpdir:
            with xopen(str(args.file.name)) as uncompressed_input:
                # validate input file
                validator = InputValidator()
                validation = validator.validate_fasta_file(
                    uncompressed_input.name, run_mode, config
                )

                if not validation.is_valid:
                    logging.error(validation.error_message)
                    sys.exit(1)

                if validation.warnings is not None:
                    for warning in validation.warnings:
                        logging.warning(warning)

                # Validate output compatibility
                output_validation = validator.validate_output_compatibility(
                    validation.sequence_type, bool(args.fna), run_mode
                )

                if not output_validation.is_valid:
                    logging.error(output_validation.error_message)
                    sys.exit(1)

                # llog startup info
                hkglib.log_startup_info(config, run_mode)

                # determine CPU count
                num_cpus = hkglib.determine_cpu_count(args.t)

                # create processing context
                ctx = hkglib.ProcessingContext(
                    input_file=uncompressed_input.name,
                    temp_dir=tmpdir,
                    num_cpus=num_cpus,
                    run_mode=run_mode,
                    sequence_type=validation.sequence_type,
                    save_dna=bool(args.fna),
                    save_protein=bool(args.faa),
                    split_output=bool(args.s),
                )

                # Process sequences
                processed_file = hkglib.process_sequences(ctx, config, tmpdir)

                # Run HMM search
                logging.info("Classifying sequences into housekeeping genes")
                searcher = HMMSearcher(num_cpus)
                best_results = searcher.search(processed_file, config)

                # sort and write main results
                filtered_results = [
                    best_results[k] for k in sorted(best_results)
                ]

                logging.info("Writing output")
                ResultWriter.write_results(
                    Path(args.o.name),
                    filtered_results,
                    to_stdout=(args.o.name == "<stdout>"),
                )

                # write sequence outputs if requested
                original_fasta = pyfastx.Fasta(uncompressed_input.name)
                hkglib.write_output_sequences(
                    args, filtered_results, tmpdir, ctx, original_fasta, config
                )
    except Exception as e:
        logging.error(f"Fatal error: {e}", exc_info=args.debug)
        sys.exit(1)

    finally:
        with contextlib.suppress(OSError):
            Path.unlink(Path(f"{args.file.name}.fxi"))

    # success message
    elapsed = datetime.datetime.now(tz=datetime.timezone.utc) - startime
    logging.info("Task finished successfully")
    logging.info(f"Walltime used (hh:mm:ss:ms): {elapsed}")


if __name__ == "__main__":
    main()
