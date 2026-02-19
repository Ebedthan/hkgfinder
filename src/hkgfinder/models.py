# Copyright 2022-2026 Anicet Ebou.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed except according
# to those terms.

"""Data models for hkgfinder."""

from dataclasses import dataclass
from enum import Enum
from typing import Optional


class SequenceType(Enum):
    """Supported sequence types."""

    DNA = "DNA"
    PROTEIN = "protein"
    UNKNOW = "unknown"


class RunMode(Enum):
    """Execution modes for hkgfinder."""

    NORMAL = "normal"
    GENOME = "genome"
    METAGENOME = "metagenome"


@dataclass(frozen=True)
class HMMResult:
    """Result from HMM search."""

    hmm_id: str
    hmm_desc: str
    hit_id: str
    evalue: float
    bitscore: float

    def __lt__(self, other: "HMMResult") -> bool:
        """Compare by bitscore for sorting."""
        return self.bitscore < other.bitscore

    def to_tsv_line(self) -> str:
        """Convert to TSV format."""
        return (
            f"{self.hit_id}\t{self.hmm_id}\t{self.hmm_desc}\t"
            f"{self.evalue:.3f}\t{self.bitscore:.3f}"
        )


@dataclass
class ValidationResult:
    """Result of input validation."""

    is_valid: bool
    sequence_type: SequenceType
    error_message: Optional[str] = None
    warnings: list[str] | None = None

    def __post_init__(self):
        if self.warnings is None:
            self.warnings = []


@dataclass
class ProcessingContext:
    """Context for processing pipeline."""

    input_file: str
    temp_dir: str
    num_cpus: int
    run_mode: RunMode
    sequence_type: SequenceType
    save_dna: bool = False
    save_protein: bool = False
    split_output: bool = False
