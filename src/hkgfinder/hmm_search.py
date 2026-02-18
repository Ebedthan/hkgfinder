"""HMM search functionality with improved error handling."""

import importlib.resources
import logging
import sys
from pathlib import Path
from typing import Dict

import psutil
import pyhmmer

from .config import HKGFinderConfig
from .models import HMMResult


class HMMSearcher:
    """Performs HMM searches on sequences."""

    def __init__(self, num_cpus: int):
        """
        Initialize HMM searcher.

        Args:
            num_cpus: Number of CPUs to use
        """
        self.num_cpus = num_cpus

    def search(
        self, sequence_file: Path, config: HKGFinderConfig
    ) -> Dict[str, HMMResult]:
        """
        Search sequences against HMM database.

        Args:
            sequence_file: Path to sequence file

        Returns:
            Dictionary mapping hit_id to best HMMResult

        Raises:
            FileNotFoundError: If HMM or sequence file not found
            RuntimeError: If HMM search fails
        """
        if not sequence_file.exists():
            raise FileNotFoundError(f"Sequence file not found: {sequence_file}")

        results = []
        best_results = {}

        # Determine if we should prefetch into memory
        should_prefetch = self._should_prefetch_targets(sequence_file, config)

        try:
            hmm_path = importlib.resources.files("hkgfinder").joinpath(
                "hkgfinder.hmm"
            )

            with (
                hmm_path.open("rb") as hmm,
                pyhmmer.plan7.HMMFile(
                    hmm  # type: ignore[type-supported]
                ) as hmm_file,
                pyhmmer.easel.SequenceFile(
                    sequence_file, digital=True
                ) as seq_file,
            ):
                targets = self._prepare_targets(seq_file, should_prefetch)
                results = self._run_hmmsearch(hmm_file, targets, config)

        except FileNotFoundError as e:
            logging.error(f"HMM database file not found: {e}")
            raise
        except Exception as e:
            logging.error(f"HMM search failed: {e}")
            raise RuntimeError(f"HMM search failed: {e}") from e

        # Find best hit for each sequence
        for result in results:
            existing = best_results.get(result.hit_id)
            if existing is None or result.bitscore > existing.bitscore:
                best_results[result.hit_id] = result

        logging.info(f"Found {len(best_results)} sequences with HMM matches")
        return best_results

    def _should_prefetch_targets(
        self, sequence_file: Path, config: HKGFinderConfig
    ) -> bool:
        """Determine if targets should be prefetched into memory."""
        try:
            available_memory = psutil.virtual_memory().available
            file_size = sequence_file.stat().st_size
            return (
                file_size < available_memory * config.memory_prefetch_threshold
            )
        except Exception as e:
            logging.warning(f"Could not determine memory availability: {e}")
            return False

    def _prepare_targets(self, seq_file, should_prefetch: bool):
        """Prepare target sequences for HMM search."""
        if should_prefetch:
            logging.info("Pre-fetching targets into memory")
            targets = seq_file.read_block()
            mem_kb = sum(sys.getsizeof(target) for target in targets) / 1024
            logging.info(f"Database in-memory size: {mem_kb:.1f} KiB")
            return targets
        else:
            logging.info("Streaming targets from disk")
            return seq_file

    def _run_hmmsearch(
        self, hmm_file, targets, config: HKGFinderConfig
    ) -> list[HMMResult]:
        """Run HMM search and collect results."""
        results = []

        try:
            for hits in pyhmmer.hmmer.hmmsearch(
                hmm_file,
                targets,
                cpus=self.num_cpus,
                bit_cutoffs="gathering",
            ):
                hmm_id = hits.query.name.decode("utf-8")
                hmm_desc = config.HMM_DESCRIPTIONS.get(  # type: ignore[type-supported]
                    hmm_id,
                    f"Unknown HMM: {hmm_id}",
                )

                for hit in hits:
                    if hit.included:
                        results.append(
                            HMMResult(
                                hmm_id=hmm_id,
                                hmm_desc=hmm_desc,
                                hit_id=hit.name.decode("utf-8"),
                                evalue=hit.evalue,
                                bitscore=hit.score,
                            )
                        )
        except Exception as e:
            logging.error(f"Error during HMM search: {e}")
            raise

        logging.info(f"HMM search found {len(results)} total hits")
        return results
