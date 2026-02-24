import warnings
from unittest.mock import patch

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


class TestInputToJson:
    """Integration tests for input_to_json() — full parse+embed+optimize pipeline."""

    def test_smiles_returns_json(self, script_module):
        result = script_module.input_to_json("CCO", "smiles")
        assert "atoms" in result
        assert "bonds" in result
        assert len(result["atoms"]) > 0
        assert len(result["bonds"]) > 0

    def test_inchi_returns_json(self, script_module):
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        result = script_module.input_to_json(inchi, "inchi")
        assert len(result["atoms"]) > 0

    def test_sequence_returns_json(self, script_module):
        result = script_module.input_to_json("GGG", "sequence")
        assert len(result["atoms"]) > 0

    def test_invalid_input_raises(self, script_module):
        with pytest.raises(ValueError, match="Invalid smiles input"):
            script_module.input_to_json("not_valid_$$$", "smiles")

    def test_json_round_trip(self, script_module):
        """Verify JSON output is actually JSON-serializable."""
        import json

        result = script_module.input_to_json("CCO", "smiles")
        json_str = json.dumps(result)
        parsed = json.loads(json_str)
        assert parsed["atoms"] == result["atoms"]
        assert parsed["bonds"] == result["bonds"]

    def test_embed_failure_raises_runtime_error(self, script_module):
        """Verify RuntimeError when 3D embedding fails."""
        with patch("rdkit.Chem.AllChem.EmbedMolecule", return_value=-1):
            with pytest.raises(RuntimeError, match="Failed to generate 3D coordinates"):
                script_module.input_to_json("CCO", "smiles")

    def test_mmff_failure_warns_but_returns_result(self, script_module):
        """When MMFF optimization fails, a warning is issued but valid JSON is still produced."""
        with patch("rdkit.Chem.AllChem.MMFFOptimizeMolecule", return_value=-1):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = script_module.input_to_json("CCO", "smiles")
                assert len(w) == 1
                assert "MMFF force field setup failed" in str(w[0].message)
        assert "atoms" in result
        assert len(result["atoms"]) > 0

    def test_mmff_non_convergence_warns(self, script_module):
        """When MMFF does not converge, a warning is issued."""
        with patch("rdkit.Chem.AllChem.MMFFOptimizeMolecule", return_value=1):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = script_module.input_to_json("CCO", "smiles")
                assert len(w) == 1
                assert "did not converge" in str(w[0].message)
        assert "atoms" in result


class TestMolToJson:
    """Tests for mol_to_json() — RDKit Mol to JSON-serializable dict."""

    def test_returns_atoms_and_bonds_keys(self, script_module):
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.AddHs(Chem.MolFromSmiles("C"))
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        result = script_module.mol_to_json(mol)
        assert "atoms" in result
        assert "bonds" in result

    def test_atom_has_element_and_coords(self, script_module):
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.AddHs(Chem.MolFromSmiles("C"))
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        result = script_module.mol_to_json(mol)
        atom = result["atoms"][0]
        assert "element" in atom
        assert "x" in atom
        assert "y" in atom
        assert "z" in atom
        assert isinstance(atom["x"], float)

    def test_ethanol_atom_count(self, script_module):
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        result = script_module.mol_to_json(mol)
        # CCO = 3 heavy + 8 H = 9 atoms total (with explicit H)
        assert len(result["atoms"]) == 9

    def test_bond_has_begin_end_order(self, script_module):
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.AddHs(Chem.MolFromSmiles("C"))
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        result = script_module.mol_to_json(mol)
        bond = result["bonds"][0]
        assert "begin" in bond
        assert "end" in bond
        assert "order" in bond

    def test_ethanol_bond_orders(self, script_module):
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        result = script_module.mol_to_json(mol)
        orders = [b["order"] for b in result["bonds"]]
        # All bonds in ethanol are single bonds
        assert all(o == 1.0 for o in orders)

    def test_benzene_has_aromatic_bonds(self, script_module):
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        result = script_module.mol_to_json(mol)
        orders = [b["order"] for b in result["bonds"]]
        # Benzene has 6 aromatic C-C bonds (1.5) + 6 C-H single bonds (1.0)
        aromatic_count = sum(1 for o in orders if o == 1.5)
        assert aromatic_count == 6

    def test_double_bond_order(self, script_module):
        """Ethylene (C=C) should have a double bond with order 2.0."""
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.AddHs(Chem.MolFromSmiles("C=C"))
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        result = script_module.mol_to_json(mol)
        orders = [b["order"] for b in result["bonds"]]
        assert 2.0 in orders

    def test_triple_bond_order(self, script_module):
        """Acetylene (C#C) should have a triple bond with order 3.0."""
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.AddHs(Chem.MolFromSmiles("C#C"))
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        result = script_module.mol_to_json(mol)
        orders = [b["order"] for b in result["bonds"]]
        assert 3.0 in orders

    def test_bond_indices_within_atom_range(self, script_module):
        """Bond begin/end indices must reference valid atom positions."""
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        result = script_module.mol_to_json(mol)
        num_atoms = len(result["atoms"])
        for bond in result["bonds"]:
            assert 0 <= bond["begin"] < num_atoms
            assert 0 <= bond["end"] < num_atoms

    def test_no_conformer_raises_value_error(self, script_module):
        """mol_to_json raises ValueError when mol has no conformer."""
        from rdkit import Chem

        mol = Chem.MolFromSmiles("C")
        with pytest.raises(ValueError, match="mol_to_json requires a Mol with 3D"):
            script_module.mol_to_json(mol)
