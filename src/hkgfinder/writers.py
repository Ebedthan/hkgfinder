# Copyright 2022-2026 Anicet Ebou.
# Licensed under the MIT license (http://opensource.org/licenses/MIT)
# This file may not be copied, modified, or distributed except according
# to those terms.

"""Efficient sequence and result writers."""

import logging
from datetime import datetime
from pathlib import Path
from typing import Optional, TextIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .config import HKGFinderConfig
from .models import HMMResult


class SequenceWriter:
    """Handles efficient writing of sequences to files."""

    def __init__(
        self,
        base_filename: str,
        extension: str,
        split: bool,
        config: HKGFinderConfig,
    ):
        """
        Initialize sequence writer.

        Args:
            base_filename: Base name for output files
            extension: File extension (faa or fna)
            split: Whether to split by gene
            config: Configuration object
        """
        self.base_filename = base_filename
        self.extension = extension
        self.split = split
        self.config = config
        self._file_handles: dict[str, TextIO] = {}

    def write_sequences(
        self,
        results: list[HMMResult],
        sequences: dict,
    ) -> None:
        """
        Write sequences to output files.

        Args:
            results: List of HMM results
            sequences: Dictionary of sequence objects
        """
        if self.split:
            self._write_split(results, sequences)
        else:
            self._write_combined(results, sequences)

    def _write_split(
        self,
        results: list[HMMResult],
        sequences: dict,
    ) -> None:
        """Write sequences to separate files by gene."""
        # Group by hmm_id
        grouped: dict[str, list[HMMResult]] = {}
        for result in results:
            grouped.setdefault(result.hmm_id, []).append(result)

        for hmm_id, gene_results in grouped.items():
            filepath = Path(f"{self.base_filename}_{hmm_id}.{self.extension}")
            try:
                with filepath.open("a", encoding="utf-8") as out:
                    self._write_batch(out, gene_results, sequences)
                logging.info(
                    f"Wrote {len(gene_results)} sequences to {filepath}"
                )
            except IOError as e:
                logging.error(f"Failed to write to {filepath}: {e}")
                raise

    def _write_combined(
        self,
        results: list[HMMResult],
        sequences: dict,
    ) -> None:
        """Write all sequences to a single file."""
        filepath = Path(f"{self.base_filename}.{self.extension}")
        try:
            with filepath.open("w", encoding="utf-8") as out:
                self._write_batch(out, results, sequences)
            logging.info(f"Wrote {len(results)} sequences to {filepath}")
        except IOError as e:
            logging.error(f"Failed to write to {filepath}: {e}")
            raise

    def _write_batch(
        self,
        file_handle: TextIO,
        results: list[HMMResult],
        sequences: dict,
    ) -> None:
        """Write a batch of sequences with buffering."""
        buffer = []

        for result in results:
            try:
                seq = sequences[result.hit_id].seq
                formatted = self._format_sequence(
                    result.hit_id, result.hmm_id, seq
                )
                buffer.append(formatted)

                if len(buffer) >= self.config.write_chunk_size:
                    file_handle.write("".join(buffer))
                    buffer = []
            except KeyError:
                logging.warning(
                    f"Sequence not found for hit_id: {result.hit_id}"
                )
                continue

        # Write remaining
        if buffer:
            file_handle.write("".join(buffer))

    def _format_sequence(self, seq_id: str, hmm_id: str, sequence: str) -> str:
        """Format sequence in FASTA format with proper line wrapping."""
        chunks = [
            sequence[i : i + self.config.seq_width]
            for i in range(0, len(sequence), self.config.seq_width)
        ]
        return f">{seq_id}_gene={hmm_id}\n" + "\n".join(chunks) + "\n"


class GenBankWriter:
    """Handles writing sequences to GenBank format with annotations."""

    def __init__(
        self, base_filename: str, split: bool, config: HKGFinderConfig
    ):
        """
        Initialize GenBank writer.

        Args:
            base_filename: Base name for output files
            split: Whether to split by gene
            config: Configuration object
        """
        self.base_filename = base_filename
        self.split = split
        self.config = config

    def write_genbank(
        self,
        results: list[HMMResult],
        sequences: dict,
        is_protein: bool = False,
        organism: Optional[str] = None,
    ) -> None:
        """
        Write sequences to GenBank format with feature annotations.

        Args:
            results: List of HMM results
            sequences: Dictionary of sequence objects
            is_protein: Whether sequences are proteins (True) or DNA (False)
            organism: Optional organism name
            taxonomy: Optional taxonomy lineage
        """
        if self.split:
            self._write_genbank_split(results, sequences, is_protein, organism)
        else:
            self._write_genbank_combined(
                results, sequences, is_protein, organism
            )

    def _write_genbank_split(
        self,
        results: list[HMMResult],
        sequences: dict,
        is_protein: bool,
        organism: Optional[str],
    ) -> None:
        """Write GenBank records to separate files by gene."""
        # Group by hmm_id
        grouped: dict[str, list[HMMResult]] = {}
        for result in results:
            grouped.setdefault(result.hmm_id, []).append(result)

        for hmm_id, gene_results in grouped.items():
            filepath = Path(f"{self.base_filename}_{hmm_id}.gb")
            try:
                records = self._create_genbank_records(
                    gene_results,
                    sequences,
                    is_protein,
                    organism,
                )
                SeqIO.write(records, filepath, "genbank")
                logging.info(
                    f"Wrote {len(gene_results)} GenBank records to {filepath}"
                )
            except Exception as e:
                logging.error(f"Failed to write GenBank file {filepath}: {e}")
                raise

    def _write_genbank_combined(
        self,
        results: list[HMMResult],
        sequences: dict,
        is_protein: bool,
        organism: Optional[str],
    ) -> None:
        """Write all GenBank records to a single file."""
        filepath = Path(f"{self.base_filename}.gb")
        try:
            records = self._create_genbank_records(
                results,
                sequences,
                is_protein,
                organism,
            )
            SeqIO.write(records, filepath, "genbank")
            logging.info(f"Wrote {len(results)} GenBank records to {filepath}")
        except Exception as e:
            logging.error(f"Failed to write GenBank file {filepath}: {e}")
            raise

    def _create_genbank_records(
        self,
        results: list[HMMResult],
        sequences: dict,
        is_protein: bool,
        organism: Optional[str],
    ) -> list[SeqRecord]:
        """
        Create GenBank SeqRecord objects with annotations.

        Args:
            results: List of HMM results
            sequences: Dictionary of sequence objects
            is_protein: Whether sequences are proteins
            organism: Organism name

        Returns:
            List of SeqRecord objects
        """
        records = []

        for result in results:
            try:
                seq_obj = sequences[result.hit_id]
                seq_str = str(seq_obj.seq)

                # Create SeqRecord
                record = SeqRecord(
                    Seq(seq_str),
                    id=result.hit_id,
                    name=result.hmm_id,
                    description=result.hmm_desc,
                )

                # Set molecule type
                record.annotations["molecule_type"] = (
                    "protein" if is_protein else "DNA"
                )

                # Add metadata
                record.annotations["topology"] = "linear"
                record.annotations["date"] = (
                    datetime.now().strftime("%d-%b-%Y").upper()
                )

                if organism:
                    record.annotations["organism"] = organism

                # Add source for hkgfinder
                record.annotations["source"] = "hkgfinder"
                record.annotations["comment"] = (
                    f"Identified by hkgfinder v{self.config.VERSION} using HMM search. "
                    f"E-value: {result.evalue:.3e}, Bit score: {result.bitscore:.1f}"
                )

                # Add features
                features = self._create_features(
                    result, len(seq_str), is_protein
                )
                record.features = features

                records.append(record)

            except KeyError:
                logging.warning(
                    f"Sequence not found for hit_id: {result.hit_id}"
                )
                continue
            except Exception as e:
                logging.warning(
                    f"Failed to create GenBank record for {result.hit_id}: {e}"
                )
                continue

        return records

    def _create_features(
        self,
        result: HMMResult,
        seq_length: int,
        is_protein: bool,
    ) -> list[SeqFeature]:
        """
        Create GenBank features for the sequence.

        Args:
            result: HMM result
            seq_length: Length of sequence
            is_protein: Whether sequence is protein

        Returns:
            List of SeqFeature objects
        """
        features = []

        # Create source feature
        source_feature = SeqFeature(
            FeatureLocation(0, seq_length),
            type="source",
            qualifiers={
                "mol_type": "protein" if is_protein else "genomic DNA",
            },
        )
        features.append(source_feature)

        # Create gene feature for DNA or CDS for protein
        if is_protein:
            # For proteins, add a CDS feature
            cds_feature = SeqFeature(
                FeatureLocation(0, seq_length),
                type="CDS",
                qualifiers={
                    "gene": result.hmm_id,
                    "product": result.hmm_desc,
                    "note": f"E-value: {result.evalue:.3e}",
                    "inference": "profile:hkgfinder:HMM",
                    "translation": "",  # Will be filled by BioPython
                },
            )
            features.append(cds_feature)
        else:
            # For DNA, add gene and CDS features
            gene_feature = SeqFeature(
                FeatureLocation(0, seq_length),
                type="gene",
                qualifiers={
                    "gene": result.hmm_id,
                    "note": f"E-value: {result.evalue:.3e}",
                },
            )
            features.append(gene_feature)

            cds_feature = SeqFeature(
                FeatureLocation(0, seq_length),
                type="CDS",
                qualifiers={
                    "gene": result.hmm_id,
                    "product": result.hmm_desc,
                    "inference": "profile:hkgfinder:HMM",
                },
            )
            features.append(cds_feature)

        return features


class ResultWriter:
    """Handles writing of tabular results."""

    @staticmethod
    def write_results(
        output_path: Path,
        results: list[HMMResult],
        to_stdout: bool = False,
    ) -> None:
        """
        Write results to TSV format.

        Args:
            output_path: Output file path
            results: List of HMM results
            to_stdout: Write to stdout instead of file
        """
        header = "seq_name\tpred_gene\tgene_desc\te-value\tbitscore"
        lines = [header] + [result.to_tsv_line() for result in results]
        output = "\n".join(lines)

        if to_stdout:
            print(output)
        else:
            try:
                output_path.write_text(output, encoding="utf-8")
                logging.info(
                    f"Wrote results for {len(results)} sequences to {output_path}"
                )
            except IOError as e:
                logging.error(f"Failed to write results to {output_path}: {e}")
                raise
