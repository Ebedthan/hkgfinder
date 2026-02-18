"""Unit tests for writers.py module."""

import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from Bio import SeqIO

from hkgfinder.config import HKGFinderConfig
from hkgfinder.models import HMMResult

# Assuming your package structure
from hkgfinder.writers import GenBankWriter, ResultWriter, SequenceWriter


class MockSequence:
    """Mock sequence object to simulate pyfastx behavior."""

    def __init__(self, seq_string):
        self.seq = seq_string


class TestSequenceWriter(unittest.TestCase):
    """Test cases for SequenceWriter class."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = TemporaryDirectory()
        self.config = HKGFinderConfig()

        # Create mock HMM results
        self.results = [
            HMMResult(
                hmm_id="DnaK",
                hmm_desc="Molecular chaperonne DnaK",
                hit_id="seq_001",
                evalue=1.5e-50,
                bitscore=168.5,
            ),
            HMMResult(
                hmm_id="RecA",
                hmm_desc="recombinase RecA",
                hit_id="seq_002",
                evalue=2.3e-45,
                bitscore=155.2,
            ),
            HMMResult(
                hmm_id="DnaK",
                hmm_desc="Molecular chaperonne DnaK",
                hit_id="seq_003",
                evalue=3.1e-48,
                bitscore=162.1,
            ),
        ]

        # Create mock sequences
        self.sequences = {
            "seq_001": MockSequence("ATGCGATCGATCGATCGATCG"),
            "seq_002": MockSequence("GCTAGCTAGCTAGCTAGCTAG"),
            "seq_003": MockSequence("TTAACCGGTTAACCGGTTAAC"),
        }

    def tearDown(self):
        """Clean up test fixtures."""
        self.temp_dir.cleanup()

    def test_init(self):
        """Test SequenceWriter initialization."""
        writer = SequenceWriter(
            "test_output",
            "faa",
            split=False,
            config=self.config,
        )

        self.assertEqual(writer.base_filename, "test_output")
        self.assertEqual(writer.extension, "faa")
        self.assertFalse(writer.split)
        self.assertEqual(writer.config, self.config)

    def test_write_combined_sequences(self):
        """Test writing all sequences to a single file."""
        output_path = Path(self.temp_dir.name) / "output"

        writer = SequenceWriter(
            str(output_path),
            "faa",
            split=False,
            config=self.config,
        )

        writer.write_sequences(self.results, self.sequences)

        # Check file was created
        expected_file = Path(f"{output_path}.faa")
        self.assertTrue(expected_file.exists())

        # Check content
        content = expected_file.read_text()
        self.assertIn(">seq_001_gene=DnaK", content)
        self.assertIn(">seq_002_gene=RecA", content)
        self.assertIn(">seq_003_gene=DnaK", content)
        self.assertIn("ATGCGATCGATCGATCGATCG", content)

    def test_write_split_sequences(self):
        """Test writing sequences to separate files by gene."""
        output_path = Path(self.temp_dir.name) / "output"

        writer = SequenceWriter(
            str(output_path),
            "fna",
            split=True,
            config=self.config,
        )

        writer.write_sequences(self.results, self.sequences)

        # Check files were created
        dnak_file = Path(f"{output_path}_DnaK.fna")
        reca_file = Path(f"{output_path}_RecA.fna")

        self.assertTrue(dnak_file.exists())
        self.assertTrue(reca_file.exists())

        # Check DnaK file has 2 sequences
        dnak_content = dnak_file.read_text()
        self.assertIn(">seq_001_gene=DnaK", dnak_content)
        self.assertIn(">seq_003_gene=DnaK", dnak_content)
        self.assertNotIn(">seq_002_gene=RecA", dnak_content)

        # Check RecA file has 1 sequence
        reca_content = reca_file.read_text()
        self.assertIn(">seq_002_gene=RecA", reca_content)
        self.assertNotIn("DnaK", reca_content)

    def test_format_sequence_line_wrapping(self):
        """Test sequence formatting with line wrapping."""
        writer = SequenceWriter(
            "test",
            "faa",
            split=False,
            config=self.config,
        )

        # Create a long sequence that should wrap
        long_seq = "A" * 200
        formatted = writer._format_sequence("test_id", "test_gene", long_seq)

        lines = formatted.strip().split("\n")

        # First line should be header
        self.assertEqual(lines[0], ">test_id_gene=test_gene")

        # Remaining lines should be max seq_width long (except possibly last)
        for line in lines[1:-1]:
            self.assertEqual(len(line), self.config.seq_width)

    def test_write_with_missing_sequence(self):
        """Test handling of missing sequences."""
        # Add a result for a non-existent sequence
        results_with_missing = self.results + [
            HMMResult(
                hmm_id="GyrB",
                hmm_desc="DNA topoisomerase",
                hit_id="missing_seq",
                evalue=1e-40,
                bitscore=140.0,
            )
        ]

        output_path = Path(self.temp_dir.name) / "output"

        writer = SequenceWriter(
            str(output_path),
            "faa",
            split=False,
            config=self.config,
        )

        # Should not raise exception, just log warning
        with patch("hkgfinder.writers.logging") as mock_logging:
            writer.write_sequences(results_with_missing, self.sequences)

            # Check warning was logged
            mock_logging.warning.assert_called()

        # Check file was still created with valid sequences
        expected_file = Path(f"{output_path}.faa")
        self.assertTrue(expected_file.exists())

        content = expected_file.read_text()
        self.assertNotIn("missing_seq", content)

    def test_write_io_error(self):
        """Test handling of I/O errors."""
        writer = SequenceWriter(
            "/invalid/path/output",
            "faa",
            split=False,
            config=self.config,
        )

        with self.assertRaises(IOError):
            writer.write_sequences(self.results, self.sequences)


class TestGenBankWriter(unittest.TestCase):
    """Test cases for GenBankWriter class."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = TemporaryDirectory()
        self.config = HKGFinderConfig()

        self.results = [
            HMMResult(
                hmm_id="DnaK",
                hmm_desc="Molecular chaperonne DnaK",
                hit_id="seq_001",
                evalue=1.5e-50,
                bitscore=168.5,
            ),
            HMMResult(
                hmm_id="RecA",
                hmm_desc="recombinase RecA",
                hit_id="seq_002",
                evalue=2.3e-45,
                bitscore=155.2,
            ),
        ]

        # Protein sequences
        self.protein_sequences = {
            "seq_001": MockSequence("MKKLAVAGALALAGSANA"),
            "seq_002": MockSequence("MTELLNVLQAVRKTG"),
        }

        # DNA sequences
        self.dna_sequences = {
            "seq_001": MockSequence("ATGAAAAAACTGGCAGTTGCAGGCGCACTG"),
            "seq_002": MockSequence("ATGACAGAACTGCTGAATGTGCTGCAGGCA"),
        }

    def tearDown(self):
        """Clean up test fixtures."""
        self.temp_dir.cleanup()

    def test_init(self):
        """Test GenBankWriter initialization."""
        writer = GenBankWriter(
            "test_output",
            split=False,
            config=self.config,
        )

        self.assertEqual(writer.base_filename, "test_output")
        self.assertFalse(writer.split)
        self.assertEqual(writer.config, self.config)

    def test_write_protein_genbank_combined(self):
        """Test writing protein sequences to GenBank format."""
        output_path = Path(self.temp_dir.name) / "output"

        writer = GenBankWriter(
            str(output_path),
            split=False,
            config=self.config,
        )

        writer.write_genbank(
            self.results,
            self.protein_sequences,
            is_protein=True,
            organism="Escherichia coli",
        )

        # Check file was created
        gb_file = Path(f"{output_path}.gb")
        self.assertTrue(gb_file.exists())

        # Parse GenBank file
        records = list(SeqIO.parse(gb_file, "genbank"))

        self.assertEqual(len(records), 2)

        # Check first record
        record1 = records[0]
        self.assertEqual(record1.id, "seq_001")
        self.assertEqual(record1.name, "DnaK")
        self.assertIn("Molecular chaperonne DnaK", record1.description)
        self.assertEqual(record1.annotations["molecule_type"], "protein")
        self.assertEqual(record1.annotations["organism"], "Escherichia coli")

        # Check features
        self.assertTrue(len(record1.features) >= 2)  # At least source and CDS

        # Check for CDS feature
        cds_features = [f for f in record1.features if f.type == "CDS"]
        self.assertEqual(len(cds_features), 1)
        self.assertEqual(cds_features[0].qualifiers["gene"][0], "DnaK")
        self.assertEqual(
            cds_features[0].qualifiers["product"][0],
            "Molecular chaperonne DnaK",
        )

    def test_write_dna_genbank_combined(self):
        """Test writing DNA sequences to GenBank format."""
        output_path = Path(self.temp_dir.name) / "output"

        writer = GenBankWriter(
            str(output_path),
            split=False,
            config=self.config,
        )

        writer.write_genbank(
            self.results,
            self.dna_sequences,
            is_protein=False,
            organism="Bacillus subtilis",
        )

        gb_file = Path(f"{output_path}.gb")
        self.assertTrue(gb_file.exists())

        records = list(SeqIO.parse(gb_file, "genbank"))
        self.assertEqual(len(records), 2)

        # Check DNA-specific features
        record1 = records[0]
        self.assertEqual(record1.annotations["molecule_type"], "DNA")

        # DNA should have both gene and CDS features
        gene_features = [f for f in record1.features if f.type == "gene"]
        cds_features = [f for f in record1.features if f.type == "CDS"]

        self.assertEqual(len(gene_features), 1)
        self.assertEqual(len(cds_features), 1)

    def test_write_genbank_split(self):
        """Test writing GenBank files split by gene."""
        # Add another DnaK result
        results_with_duplicate = self.results + [
            HMMResult(
                hmm_id="DnaK",
                hmm_desc="Molecular chaperonne DnaK",
                hit_id="seq_003",
                evalue=1.0e-45,
                bitscore=160.0,
            )
        ]

        sequences_with_extra = {
            **self.protein_sequences,
            "seq_003": MockSequence("MKKLAVAGALALAGSANAAA"),
        }

        output_path = Path(self.temp_dir.name) / "output"

        writer = GenBankWriter(
            str(output_path),
            split=True,
            config=self.config,
        )

        writer.write_genbank(
            results_with_duplicate,
            sequences_with_extra,
            is_protein=True,
        )

        # Check separate files were created
        dnak_file = Path(f"{output_path}_DnaK.gb")
        reca_file = Path(f"{output_path}_RecA.gb")

        self.assertTrue(dnak_file.exists())
        self.assertTrue(reca_file.exists())

        # Check DnaK file has 2 records
        dnak_records = list(SeqIO.parse(dnak_file, "genbank"))
        self.assertEqual(len(dnak_records), 2)

        # Check RecA file has 1 record
        reca_records = list(SeqIO.parse(reca_file, "genbank"))
        self.assertEqual(len(reca_records), 1)

    def test_genbank_annotations(self):
        """Test GenBank annotations are properly set."""
        output_path = Path(self.temp_dir.name) / "output"

        writer = GenBankWriter(
            str(output_path),
            split=False,
            config=self.config,
        )

        writer.write_genbank(
            [self.results[0]],
            self.protein_sequences,
            is_protein=True,
            organism="Test organism",
        )

        gb_file = Path(f"{output_path}.gb")
        records = list(SeqIO.parse(gb_file, "genbank"))
        record = records[0]

        # Check all required annotations
        self.assertEqual(record.annotations["source"], "hkgfinder")
        self.assertIn("hkgfinder", record.annotations["comment"])
        self.assertIn("E-value", record.annotations["comment"])
        self.assertIn("Bit score", record.annotations["comment"])
        self.assertEqual(record.annotations["topology"], "linear")

    def test_genbank_with_missing_sequence(self):
        """Test GenBank writing with missing sequences."""
        results_with_missing = self.results + [
            HMMResult(
                hmm_id="GyrB",
                hmm_desc="DNA topoisomerase",
                hit_id="missing_seq",
                evalue=1e-40,
                bitscore=140.0,
            )
        ]

        output_path = Path(self.temp_dir.name) / "output"

        writer = GenBankWriter(
            str(output_path),
            split=False,
            config=self.config,
        )

        with patch("hkgfinder.writers.logging") as mock_logging:
            writer.write_genbank(
                results_with_missing,
                self.protein_sequences,
                is_protein=True,
            )

            # Check warning was logged
            mock_logging.warning.assert_called()

        # Check file was created with valid sequences only
        gb_file = Path(f"{output_path}.gb")
        records = list(SeqIO.parse(gb_file, "genbank"))

        # Should only have 2 records, not 3
        self.assertEqual(len(records), 2)

    def test_create_features_protein(self):
        """Test feature creation for protein sequences."""
        writer = GenBankWriter(
            "test",
            split=False,
            config=self.config,
        )

        features = writer._create_features(
            self.results[0],
            seq_length=100,
            is_protein=True,
        )

        # Should have source and CDS
        self.assertEqual(len(features), 2)

        # Check source feature
        source = features[0]
        self.assertEqual(source.type, "source")
        self.assertEqual(source.qualifiers["mol_type"], "protein")

        # Check CDS feature
        cds = features[1]
        self.assertEqual(cds.type, "CDS")
        self.assertEqual(cds.qualifiers["gene"], "DnaK")
        self.assertIn("E-value", cds.qualifiers["note"])

    def test_create_features_dna(self):
        """Test feature creation for DNA sequences."""
        writer = GenBankWriter(
            "test",
            split=False,
            config=self.config,
        )

        features = writer._create_features(
            self.results[0],
            seq_length=300,
            is_protein=False,
        )

        # Should have source, gene, and CDS
        self.assertEqual(len(features), 3)

        # Check types
        types = [f.type for f in features]
        self.assertIn("source", types)
        self.assertIn("gene", types)
        self.assertIn("CDS", types)


class TestResultWriter(unittest.TestCase):
    """Test cases for ResultWriter class."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = TemporaryDirectory()

        self.results = [
            HMMResult(
                hmm_id="DnaK",
                hmm_desc="Molecular chaperonne DnaK",
                hit_id="seq_001",
                evalue=1.5e-50,
                bitscore=168.5,
            ),
            HMMResult(
                hmm_id="RecA",
                hmm_desc="recombinase RecA",
                hit_id="seq_002",
                evalue=2.3e-45,
                bitscore=155.2,
            ),
        ]

    def tearDown(self):
        """Clean up test fixtures."""
        self.temp_dir.cleanup()

    def test_write_results_to_file(self):
        """Test writing results to a file."""
        output_path = Path(self.temp_dir.name) / "results.tsv"

        ResultWriter.write_results(output_path, self.results, to_stdout=False)

        self.assertTrue(output_path.exists())

        content = output_path.read_text()
        lines = content.strip().split("\n")

        # Check header
        self.assertEqual(
            lines[0], "seq_name\tpred_gene\tgene_desc\te-value\tbitscore"
        )

        # Check first result
        self.assertIn("seq_001", lines[1])
        self.assertIn("DnaK", lines[1])
        self.assertIn("Molecular chaperonne DnaK", lines[1])

        # Check second result
        self.assertIn("seq_002", lines[2])
        self.assertIn("RecA", lines[2])

    def test_write_results_to_stdout(self):
        """Test writing results to stdout."""
        output_path = Path(self.temp_dir.name) / "dummy.tsv"

        with patch("builtins.print") as mock_print:
            ResultWriter.write_results(
                output_path, self.results, to_stdout=True
            )

            # Check print was called
            mock_print.assert_called_once()

            # Get the printed content
            printed_content = mock_print.call_args[0][0]
            self.assertIn("seq_name\tpred_gene", printed_content)
            self.assertIn("seq_001", printed_content)

    def test_write_results_formatting(self):
        """Test TSV formatting of results."""
        output_path = Path(self.temp_dir.name) / "results.tsv"

        ResultWriter.write_results(output_path, self.results, to_stdout=False)

        content = output_path.read_text()
        lines = content.strip().split("\n")

        # Check tab separation
        for line in lines[1:]:  # Skip header
            fields = line.split("\t")
            self.assertEqual(len(fields), 5)  # 5 columns

            # Check e-value and bitscore formatting
            evalue = fields[3]
            bitscore = fields[4]

            # E-value should have 3 decimal places
            self.assertRegex(evalue, r"\d+\.\d{3}")

            # Bitscore should have 3 decimal places
            self.assertRegex(bitscore, r"\d+\.\d{3}")

    def test_write_empty_results(self):
        """Test writing empty results."""
        output_path = Path(self.temp_dir.name) / "results.tsv"

        ResultWriter.write_results(output_path, [], to_stdout=False)

        self.assertTrue(output_path.exists())

        content = output_path.read_text()
        lines = content.strip().split("\n")

        # Should only have header
        self.assertEqual(len(lines), 1)
        self.assertIn("seq_name", lines[0])

    def test_write_results_io_error(self):
        """Test handling of I/O errors."""
        invalid_path = Path("/invalid/path/results.tsv")

        with self.assertRaises(IOError):
            ResultWriter.write_results(
                invalid_path, self.results, to_stdout=False
            )


class TestIntegration(unittest.TestCase):
    """Integration tests for writers."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = TemporaryDirectory()
        self.config = HKGFinderConfig()

        self.results = [
            HMMResult(
                hmm_id="DnaK",
                hmm_desc="Molecular chaperonne DnaK",
                hit_id="seq_001",
                evalue=1.5e-50,
                bitscore=168.5,
            ),
            HMMResult(
                hmm_id="RecA",
                hmm_desc="recombinase RecA",
                hit_id="seq_002",
                evalue=2.3e-45,
                bitscore=155.2,
            ),
            HMMResult(
                hmm_id="DnaK",
                hmm_desc="Molecular chaperonne DnaK",
                hit_id="seq_003",
                evalue=3.1e-48,
                bitscore=162.1,
            ),
        ]

        self.sequences = {
            "seq_001": MockSequence("MKKLAVAGALALAGSANA"),
            "seq_002": MockSequence("MTELLNVLQAVRKTG"),
            "seq_003": MockSequence("MKKLAVAGALALAGS"),
        }

    def tearDown(self):
        """Clean up test fixtures."""
        self.temp_dir.cleanup()

    def test_write_all_formats(self):
        """Test writing FASTA, GenBank, and TSV outputs together."""
        output_base = Path(self.temp_dir.name) / "output"

        # Write FASTA
        seq_writer = SequenceWriter(
            str(output_base),
            "faa",
            split=False,
            config=self.config,
        )
        seq_writer.write_sequences(self.results, self.sequences)

        # Write GenBank
        gb_writer = GenBankWriter(
            str(output_base),
            split=False,
            config=self.config,
        )
        gb_writer.write_genbank(
            self.results,
            self.sequences,
            is_protein=True,
            organism="Test organism",
        )

        # Write TSV
        ResultWriter.write_results(
            Path(f"{output_base}.tsv"),
            self.results,
        )

        # Check all files exist
        self.assertTrue(Path(f"{output_base}.faa").exists())
        self.assertTrue(Path(f"{output_base}.gb").exists())
        self.assertTrue(Path(f"{output_base}.tsv").exists())

        # Verify content consistency
        fasta_content = Path(f"{output_base}.faa").read_text()
        self.assertEqual(fasta_content.count(">seq_"), 3)

        gb_records = list(SeqIO.parse(f"{output_base}.gb", "genbank"))
        self.assertEqual(len(gb_records), 3)

        tsv_lines = Path(f"{output_base}.tsv").read_text().strip().split("\n")
        self.assertEqual(len(tsv_lines), 4)  # Header + 3 results


if __name__ == "__main__":
    unittest.main()
