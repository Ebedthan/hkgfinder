"""Configuration management for hkgfinder."""

from dataclasses import dataclass
from typing import Dict


@dataclass
class HKGFinderConfig:
    """Central configuration for hkgfinder."""

    # Version and metadata
    VERSION: str = "0.4"
    AUTHOR: str = "Anicet Ebou <anicet dot ebou at gmail.com>"
    URL: str = "https://github.com/Ebedthan/hkgfinder"

    # Processing limits (can be overriden via CLI)
    max_seq_length: int = 10_000
    buffer_size: int = 1024 * 1024  # 1 MB
    seq_width: int = 60
    write_chunk_size: int = 100
    memory_prefetch_threshold: float = 0.1

    # HMM gene descriptions
    HMM_DESCRIPTIONS: Dict[str, str] | None = None

    def __post_init__(self):
        if self.HMM_DESCRIPTIONS is None:
            self.HMM_DESCRIPTIONS = {
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


def create_config(args) -> HKGFinderConfig:
    """Create config from CLI arguments."""
    # preset profiles
    profiles = {
        "default": {"max_seq_length": 10_000, "buffer_size": 1},
        "large-genome": {"max_seq_length": 100_000, "buffer_size": 4},
        "metagenome": {"max_seq_length": 50_000, "buffer_size": 2},
    }
    preset = profiles[args.profile]
    config = HKGFinderConfig(
        max_seq_length=preset["max_seq_length"],
        buffer_size=preset["buffer_size"] * 1024 * 1024,
    )

    # Overide with CLI arguments if provided
    if args.max_seq_length != 10000:
        config.max_seq_length = args.max_seq_length

    if args.buffer_size != 1024 * 1024:
        config.buffer_size = args.buffer_size * 1024 * 1024

    return config
