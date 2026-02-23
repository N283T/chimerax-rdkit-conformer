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


class TestInputTo3d:
    """Integration tests for input_to_3d() — full parse+embed+optimize pipeline."""

    def test_smiles_generates_sdf(self, script_module, tmp_path):
        sdf_path = script_module.input_to_3d("CCO", "smiles", tmp_path / "out.sdf")
        assert sdf_path.exists()
        assert sdf_path.suffix == ".sdf"
        content = sdf_path.read_text()
        assert "V2000" in content or "V3000" in content

    def test_inchi_generates_sdf(self, script_module, tmp_path):
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        sdf_path = script_module.input_to_3d(inchi, "inchi", tmp_path / "out.sdf")
        assert sdf_path.exists()
        content = sdf_path.read_text()
        assert "V2000" in content or "V3000" in content

    def test_sequence_generates_sdf(self, script_module, tmp_path):
        sdf_path = script_module.input_to_3d("GGG", "sequence", tmp_path / "out.sdf")
        assert sdf_path.exists()

    def test_temp_file_created_when_no_output(self, script_module):
        sdf_path = script_module.input_to_3d("CCO", "smiles")
        assert sdf_path.exists()
        assert sdf_path.suffix == ".sdf"
        sdf_path.unlink()  # cleanup

    def test_invalid_input_raises(self, script_module, tmp_path):
        with pytest.raises(ValueError, match="Invalid smiles input"):
            script_module.input_to_3d("not_valid_$$$", "smiles", tmp_path / "out.sdf")
