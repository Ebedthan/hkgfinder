"""Input validation for hkgfinder."""

import logging
from typing import Set

from pyfastxcli import pyfastx

from hkgfinder.config import HKGFinderConfig

from .models import RunMode, SequenceType, ValidationResult


class InputValidator:
    """Validates input files and parameters."""

    @staticmethod
    def validate_fasta_file(
        file_path: str, run_mode: RunMode, config: HKGFinderConfig
    ) -> ValidationResult:
        """
        Validate FASTA file comprehensively.

        Args:
            file_path: Path to FASTA file
            run_mode: Execution mode

        Returns:
            ValidationResult with validation outcome
        """
        warnings = []

        try:
            fasta = pyfastx.Fasta(file_path)
        except Exception as e:
            return ValidationResult(
                is_valid=False,
                sequence_type=SequenceType.UNKNOW,
                error_message=f"Failed to parse FASTA file: {e}",
            )

        # Determine sequence type
        seq_type_str = fasta.type
        if seq_type_str == "DNA":
            seq_type = SequenceType.DNA
        elif seq_type_str == "protein":
            seq_type = SequenceType.PROTEIN
        else:
            return ValidationResult(
                is_valid=False,
                sequence_type=SequenceType.UNKNOW,
                error_message=f"Unknown sequence type: {seq_type_str}. Expected DNA or protein.",
            )

        # Validate mode compatibility
        if (
            run_mode in (RunMode.GENOME, RunMode.METAGENOME)
            and seq_type == SequenceType.PROTEIN
        ):
            return ValidationResult(
                is_valid=False,
                sequence_type=seq_type,
                error_message=f"Cannot run {run_mode.value} mode on protein sequences",
            )

        # Check for duplicate IDs
        ids = fasta.keys()
        if len(ids) != len(set(ids)):
            duplicates = InputValidator._find_duplicates(ids)
            return ValidationResult(
                is_valid=False,
                sequence_type=seq_type,
                error_message=f"FASTA file contains duplicate sequence IDs: {', '.join(list(duplicates)[:5])}",
            )

        # Check sequence length
        max_length = len(fasta.longest)
        if max_length >= config.max_seq_length and run_mode == RunMode.NORMAL:
            warnings.append(
                f"Sequences exceed {config.max_seq_length} bp (max: {max_length} bp). "
                "Consider using genome mode (-g) or tweak advanced parameters."
            )

        # Check if file is empty
        if len(ids) == 0:
            return ValidationResult(
                is_valid=False,
                sequence_type=seq_type,
                error_message="FASTA file contains no sequences",
            )

        logging.info(f"Validated {len(ids)} sequences (type: {seq_type.value})")

        return ValidationResult(
            is_valid=True,
            sequence_type=seq_type,
            warnings=warnings,
        )

    @staticmethod
    def _find_duplicates(items: list) -> Set[str]:
        """Find duplicate items in a list."""
        seen = set()
        duplicates = set()
        for item in items:
            if item in seen:
                duplicates.add(item)
            seen.add(item)
        return duplicates

    @staticmethod
    def validate_output_compatibility(
        seq_type: SequenceType,
        save_dna: bool,
        run_mode: RunMode,
    ) -> ValidationResult:
        """Validate output format compatibility."""
        if (
            seq_type == SequenceType.PROTEIN
            and save_dna
            and run_mode == RunMode.NORMAL
        ):
            return ValidationResult(
                is_valid=False,
                sequence_type=seq_type,
                error_message="Cannot retrieve DNA sequences from protein input in normal mode. Remove --fna.",
            )

        return ValidationResult(
            is_valid=True,
            sequence_type=seq_type,
        )
