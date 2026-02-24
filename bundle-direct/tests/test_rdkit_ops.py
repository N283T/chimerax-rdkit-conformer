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
    """Integration tests for input_to_json() -- full parse+embed+optimize pipeline."""

    def test_smiles_returns_json(self, script_module):
        result = script_module.input_to_json("CCO", "smiles")
        assert isinstance(result, list)
        assert len(result) == 1
        assert "atoms" in result[0]
        assert "bonds" in result[0]
        assert len(result[0]["atoms"]) > 0
        assert len(result[0]["bonds"]) > 0

    def test_inchi_returns_json(self, script_module):
        inchi = "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
        result = script_module.input_to_json(inchi, "inchi")
        assert len(result[0]["atoms"]) > 0

    def test_sequence_returns_json(self, script_module):
        result = script_module.input_to_json("GGG", "sequence")
        assert len(result[0]["atoms"]) > 0

    def test_invalid_input_raises(self, script_module):
        with pytest.raises(ValueError, match="Invalid smiles input"):
            script_module.input_to_json("not_valid_$$$", "smiles")

    def test_json_round_trip(self, script_module):
        """Verify JSON output is actually JSON-serializable."""
        import json

        result = script_module.input_to_json("CCO", "smiles")
        json_str = json.dumps(result)
        parsed = json.loads(json_str)
        assert parsed[0]["atoms"] == result[0]["atoms"]
        assert parsed[0]["bonds"] == result[0]["bonds"]

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
                result = script_module.input_to_json("CCO", "smiles", optimize=True)
                assert len(w) == 1
                assert "MMFF force field setup failed" in str(w[0].message)
        assert "atoms" in result[0]
        assert len(result[0]["atoms"]) > 0

    def test_mmff_non_convergence_warns(self, script_module):
        """When MMFF does not converge, a warning is issued."""
        with patch("rdkit.Chem.AllChem.MMFFOptimizeMolecule", return_value=1):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = script_module.input_to_json("CCO", "smiles", optimize=True)
                assert len(w) == 1
                assert "did not converge" in str(w[0].message)
        assert "atoms" in result[0]


class TestMolToJson:
    """Tests for mol_to_json() -- RDKit Mol to JSON-serializable dict."""

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

    def test_unknown_bond_type_warns_and_defaults_to_single(self, script_module):
        """Unknown bond types emit a warning and default to order 1.0."""
        from unittest.mock import MagicMock
        from rdkit import Chem
        from rdkit.Chem import AllChem

        mol = Chem.AddHs(Chem.MolFromSmiles("CC"))
        AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

        # Patch GetBondType on the first bond to return an unrecognized type
        original_get_bonds = mol.GetBonds

        def patched_get_bonds():
            bonds = list(original_get_bonds())
            fake_bond = MagicMock(wraps=bonds[0])
            fake_type = MagicMock()
            fake_type.name = "FAKE"
            fake_bond.GetBondType = MagicMock(return_value=fake_type)
            return [fake_bond] + bonds[1:]

        mol.GetBonds = patched_get_bonds

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = script_module.mol_to_json(mol)
            unknown_warnings = [x for x in w if "Unknown bond type" in str(x.message)]
            assert len(unknown_warnings) == 1
            assert "defaulting to 1.0" in str(unknown_warnings[0].message)
        assert result["bonds"][0]["order"] == 1.0


class TestMultiConformer:
    """Tests for multi-conformer generation via input_to_json()."""

    def test_single_conformer_returns_list_of_one(self, script_module):
        """Default num_confs=1 returns a list with one conformer dict."""
        result = script_module.input_to_json("CCO", "smiles")
        assert isinstance(result, list)
        assert len(result) == 1
        assert "atoms" in result[0]
        assert "bonds" in result[0]

    def test_multiple_conformers_returns_list(self, script_module):
        """num_confs=3 returns a list of conformer dicts."""
        result = script_module.input_to_json("CCO", "smiles", num_confs=3)
        assert isinstance(result, list)
        assert len(result) >= 1
        assert len(result) <= 3
        for conf in result:
            assert "atoms" in conf
            assert "bonds" in conf

    def test_all_conformers_have_same_atom_count(self, script_module):
        """All conformers for same molecule must have identical atom count."""
        result = script_module.input_to_json("c1ccccc1", "smiles", num_confs=5)
        atom_counts = [len(conf["atoms"]) for conf in result]
        assert len(set(atom_counts)) == 1  # all same

    def test_all_conformers_have_same_bond_count(self, script_module):
        """All conformers must have identical bond topology."""
        result = script_module.input_to_json("c1ccccc1", "smiles", num_confs=5)
        bond_counts = [len(conf["bonds"]) for conf in result]
        assert len(set(bond_counts)) == 1  # all same

    def test_conformers_have_different_coordinates(self, script_module):
        """Different conformers should have different 3D coordinates."""
        # n-butane is flexible enough to produce distinct conformers
        result = script_module.input_to_json("CCCC", "smiles", num_confs=5)
        if len(result) < 2:
            pytest.skip("RMS pruning left only 1 conformer")
        coords_0 = [(a["x"], a["y"], a["z"]) for a in result[0]["atoms"]]
        coords_1 = [(a["x"], a["y"], a["z"]) for a in result[1]["atoms"]]
        assert coords_0 != coords_1

    def test_pruning_may_return_fewer_conformers(self, script_module):
        """Requesting many conformers of a rigid molecule may yield fewer due to RMS pruning."""
        # Methane is extremely rigid -- most conformers will be pruned
        result = script_module.input_to_json("C", "smiles", num_confs=10)
        assert isinstance(result, list)
        assert len(result) >= 1
        assert len(result) <= 10

    def test_zero_conformers_raises(self, script_module):
        """num_confs=0 should raise ValueError."""
        with pytest.raises(ValueError, match="conformers must be between 1 and 50"):
            script_module.input_to_json("CCO", "smiles", num_confs=0)

    def test_over_max_conformers_raises(self, script_module):
        """num_confs=51 should raise ValueError."""
        with pytest.raises(ValueError, match="conformers must be between 1 and 50"):
            script_module.input_to_json("CCO", "smiles", num_confs=51)

    def test_embed_multiple_failure_raises_runtime_error(self, script_module):
        """RuntimeError when EmbedMultipleConfs produces zero conformers."""
        with patch("rdkit.Chem.AllChem.EmbedMultipleConfs", return_value=[]):
            with pytest.raises(
                RuntimeError, match="Failed to generate any 3D conformers"
            ):
                script_module.input_to_json("CCO", "smiles", num_confs=3)


class TestOptimizeFlag:
    """Tests for the optimize parameter in input_to_json()."""

    def test_optimize_false_skips_mmff(self, script_module):
        """When optimize=False, MMFF should not be called."""
        with patch("rdkit.Chem.AllChem.MMFFOptimizeMolecule") as mock_mmff:
            result = script_module.input_to_json("CCO", "smiles", optimize=False)
            mock_mmff.assert_not_called()
        assert len(result[0]["atoms"]) > 0

    def test_optimize_true_calls_mmff(self, script_module):
        """When optimize=True, MMFF should be called."""
        with patch(
            "rdkit.Chem.AllChem.MMFFOptimizeMolecule", return_value=0
        ) as mock_mmff:
            script_module.input_to_json("CCO", "smiles", optimize=True)
            mock_mmff.assert_called_once()

    def test_optimize_default_is_false(self, script_module):
        """Default optimize=False means no MMFF."""
        with patch("rdkit.Chem.AllChem.MMFFOptimizeMolecule") as mock_mmff:
            script_module.input_to_json("CCO", "smiles")
            mock_mmff.assert_not_called()

    def test_optimize_true_multi_conformer(self, script_module):
        """optimize=True with multi-conformer calls MMFFOptimizeMoleculeConfs."""
        with patch(
            "rdkit.Chem.AllChem.MMFFOptimizeMoleculeConfs",
            return_value=[(0, 1.0), (0, 2.0)],
        ) as mock_mmff:
            script_module.input_to_json("CCO", "smiles", num_confs=3, optimize=True)
            mock_mmff.assert_called_once()

    def test_optimize_false_multi_conformer_skips_mmff(self, script_module):
        """optimize=False with multi-conformer skips MMFFOptimizeMoleculeConfs."""
        with patch("rdkit.Chem.AllChem.MMFFOptimizeMoleculeConfs") as mock_mmff:
            script_module.input_to_json("CCO", "smiles", num_confs=3, optimize=False)
            mock_mmff.assert_not_called()

    def test_multi_conformer_mmff_failure_warns(self, script_module):
        """MMFF failure on a conformer emits a warning but still returns results."""
        with patch(
            "rdkit.Chem.AllChem.MMFFOptimizeMoleculeConfs",
            return_value=[(-1, 0.0), (0, 1.0)],
        ):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = script_module.input_to_json(
                    "CCO", "smiles", num_confs=3, optimize=True
                )
                mmff_warnings = [x for x in w if "MMFF" in str(x.message)]
                assert len(mmff_warnings) == 1
                assert "force field setup failed" in str(mmff_warnings[0].message)
        assert len(result) >= 1

    def test_multi_conformer_mmff_non_convergence_warns(self, script_module):
        """MMFF non-convergence on a conformer emits a warning."""
        with patch(
            "rdkit.Chem.AllChem.MMFFOptimizeMoleculeConfs",
            return_value=[(0, 1.0), (1, 2.0)],
        ):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                result = script_module.input_to_json(
                    "CCO", "smiles", num_confs=3, optimize=True
                )
                mmff_warnings = [x for x in w if "MMFF" in str(x.message)]
                assert len(mmff_warnings) == 1
                assert "did not converge" in str(mmff_warnings[0].message)
        assert len(result) >= 1
