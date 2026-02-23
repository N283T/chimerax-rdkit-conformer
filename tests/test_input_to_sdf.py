import pytest


class TestParseInput:
    """Tests for parse_input() dispatcher function."""

    def test_parse_smiles(self, script_module):
        mol = script_module.parse_input("CCO", "smiles")
        assert mol is not None
        assert mol.GetNumAtoms() == 3  # C, C, O (no explicit H)

    def test_parse_inchi(self, script_module):
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        mol = script_module.parse_input(inchi, "inchi")
        assert mol is not None
        assert mol.GetNumAtoms() == 3  # C, C, O

    def test_parse_sequence_protein(self, script_module):
        mol = script_module.parse_input("GGG", "sequence")
        assert mol is not None
        assert mol.GetNumAtoms() > 0

    def test_parse_dna(self, script_module):
        mol = script_module.parse_input("ACGT", "dna")
        assert mol is not None
        assert mol.GetNumAtoms() > 0

    def test_parse_rna(self, script_module):
        mol = script_module.parse_input("ACGU", "rna")
        assert mol is not None
        assert mol.GetNumAtoms() > 0

    def test_parse_fasta(self, script_module):
        fasta = ">test\nGGG"
        mol = script_module.parse_input(fasta, "fasta")
        assert mol is not None
        assert mol.GetNumAtoms() > 0

    def test_parse_helm(self, script_module):
        helm = "PEPTIDE1{G.G.G}$$$$"
        mol = script_module.parse_input(helm, "helm")
        assert mol is not None
        assert mol.GetNumAtoms() > 0

    def test_invalid_smiles_raises(self, script_module):
        with pytest.raises(ValueError, match="Invalid smiles input"):
            script_module.parse_input("not_a_smiles_$$$", "smiles")

    def test_unsupported_format_raises(self, script_module):
        with pytest.raises(ValueError, match="Unsupported format"):
            script_module.parse_input("CCO", "xyz")
