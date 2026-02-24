# Direct Model Creation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Eliminate the SDF file intermediary by building ChimeraX `AtomicStructure` directly from RDKit coordinate data passed as JSON over stdout.

**Architecture:** The subprocess script (`input_to_json.py`) serializes RDKit Mol objects to JSON (atoms with coordinates, bonds with orders). `cmd.py` parses the JSON and constructs an `AtomicStructure` using ChimeraX's `struct_edit` API (`add_atom`, `add_bond`). The `output` keyword is removed entirely.

**Tech Stack:** RDKit (subprocess), ChimeraX AtomicStructure API (`chimerax.atomic`), JSON

---

### Background

**Current pipeline:**
```
subprocess (RDKit) → temp SDF file → ChimeraX `open` command → delete temp file
```

**New pipeline:**
```
subprocess (RDKit) → JSON stdout → cmd.py parses → AtomicStructure in memory
```

**Key ChimeraX API references** (from ChimeraX source):
```python
from chimerax.atomic import AtomicStructure, Element, Bond
from chimerax.atomic.struct_edit import add_atom, add_bond
from numpy import array, float64

# Bond.register_attr(session, "order", "rdkit_conformer", attr_type=float)
# s = AtomicStructure(session, name="EtOH")
# r = s.new_residue("UNL", " ", 1)
# a = add_atom("C1", "C", r, array([x, y, z], dtype=float64))
# b = add_bond(a1, a2)
# b.order = 1.0
# session.models.add([s])
```

**Test command:**
```bash
uv run --no-project --with rdkit --with pytest pytest tests/ -v
```

---

### Task 1: Add `mol_to_json()` Tests (RED)

**Files:**
- Modify: `tests/test_input_to_sdf.py` (will be renamed in Task 3)

**Step 1: Write failing tests for `mol_to_json()`**

Add the following test class at the end of `tests/test_input_to_sdf.py`:

```python
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
```

**Step 2: Run tests to verify they fail**

Run: `uv run --no-project --with rdkit --with pytest pytest tests/test_input_to_sdf.py::TestMolToJson -v`

Expected: FAIL with `AttributeError: module 'input_to_sdf' has no attribute 'mol_to_json'`

**Step 3: Commit**

```bash
git add tests/test_input_to_sdf.py
git commit -m "test: add mol_to_json tests (RED)"
```

---

### Task 2: Implement `mol_to_json()` and Rename Script (GREEN)

**Files:**
- Delete: `src/input_to_sdf.py`
- Create: `src/input_to_json.py`
- Modify: `tests/conftest.py:10` (update script path)

**Step 1: Create `src/input_to_json.py`**

Copy `src/input_to_sdf.py` to `src/input_to_json.py` and make these changes:

1. Update module docstring to `"""Convert molecular notation to 3D JSON using RDKit ETKDGv3."""`
2. Add `import json` to imports
3. Remove `import tempfile`
4. Add `mol_to_json()` function after `parse_input()`:

```python
def mol_to_json(mol) -> dict:
    """Serialize an RDKit Mol with 3D coordinates to a JSON-compatible dict.

    Args:
        mol: RDKit Mol object with embedded 3D coordinates and explicit Hs.

    Returns:
        Dict with 'atoms' and 'bonds' keys.
    """
    conf = mol.GetConformer()
    atoms = []
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        atoms.append({
            "element": mol.GetAtomWithIdx(i).GetSymbol(),
            "x": pos.x,
            "y": pos.y,
            "z": pos.z,
        })

    bond_type_map = {
        Chem.BondType.SINGLE: 1.0,
        Chem.BondType.DOUBLE: 2.0,
        Chem.BondType.TRIPLE: 3.0,
        Chem.BondType.AROMATIC: 1.5,
    }
    bonds = []
    for bond in mol.GetBonds():
        bonds.append({
            "begin": bond.GetBeginAtomIdx(),
            "end": bond.GetEndAtomIdx(),
            "order": bond_type_map.get(bond.GetBondType(), 1.0),
        })

    return {"atoms": atoms, "bonds": bonds}
```

5. Replace `input_to_3d()` with `input_to_json()`:

```python
def input_to_json(input_str: str, fmt: str = "smiles") -> dict:
    """Convert molecular notation to a 3D JSON dict using ETKDGv3.

    Args:
        input_str: Molecular notation string.
        fmt: Input format (default: smiles).

    Returns:
        Dict with 'atoms' and 'bonds' keys.
    """
    mol = parse_input(input_str, fmt)
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        raise RuntimeError(f"Failed to generate 3D coordinates for: {input_str}")

    opt_result = AllChem.MMFFOptimizeMolecule(mol)
    if opt_result == -1:
        import warnings

        warnings.warn(
            f"MMFF force field setup failed for: {input_str}. "
            "3D coordinates are unoptimized.",
            stacklevel=2,
        )

    return mol_to_json(mol)
```

6. Replace `main()`:

```python
def main() -> None:
    parser = argparse.ArgumentParser(description="Convert molecular notation to 3D JSON")
    parser.add_argument("input", help="Molecular notation string")
    parser.add_argument(
        "-f",
        "--format",
        default="smiles",
        choices=VALID_FORMATS,
        help="Input format (default: smiles)",
    )
    args = parser.parse_args()

    try:
        data = input_to_json(args.input, args.format)
    except (ValueError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {type(e).__name__}: {e}", file=sys.stderr)
        sys.exit(2)

    print(json.dumps(data))
```

7. Remove unused imports: `tempfile`, `Path` (no longer needed)

The full file should have these imports:
```python
import argparse
import json
import sys

from rdkit import Chem
from rdkit.Chem import AllChem
```

**Step 2: Update `tests/conftest.py` to reference new script**

Change line 10 from:
```python
    script_path = Path(__file__).parent.parent.joinpath("src", "input_to_sdf.py")
```
to:
```python
    script_path = Path(__file__).parent.parent.joinpath("src", "input_to_json.py")
```

Update the fixture docstring:
```python
    """Import input_to_json.py as a module (bypassing PEP 723 metadata)."""
```

Update `spec_from_file_location` module name:
```python
    spec = importlib.util.spec_from_file_location("input_to_json", str(script_path))
```

**Step 3: Delete `src/input_to_sdf.py`**

```bash
git rm src/input_to_sdf.py
```

**Step 4: Run `mol_to_json` tests to verify they pass**

Run: `uv run --no-project --with rdkit --with pytest pytest tests/test_input_to_sdf.py::TestMolToJson -v`

Expected: All 6 tests PASS

**Step 5: Commit**

```bash
git add src/input_to_json.py tests/conftest.py
git commit -m "feat: add mol_to_json and rename script to input_to_json.py"
```

---

### Task 3: Update Integration Tests for JSON Output

**Files:**
- Rename: `tests/test_input_to_sdf.py` → `tests/test_input_to_json.py`
- Modify: `tests/test_input_to_json.py`

**Step 1: Rename test file**

```bash
git mv tests/test_input_to_sdf.py tests/test_input_to_json.py
```

**Step 2: Update `TestInputTo3d` → `TestInputToJson`**

Replace the entire `TestInputTo3d` class with:

```python
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
```

**Step 3: Run all tests**

Run: `uv run --no-project --with rdkit --with pytest pytest tests/ -v`

Expected: All tests PASS (9 parse_input + 5 input_to_json + 6 mol_to_json = 20 tests)

**Step 4: Commit**

```bash
git add tests/test_input_to_json.py
git commit -m "test: update integration tests for JSON output"
```

---

### Task 4: Update `cmd.py` with Direct Model Building

**Files:**
- Modify: `src/cmd.py`

This is the core change. Replace the SDF-based pipeline with JSON parsing and `AtomicStructure` construction.

**Step 1: Update imports**

Replace the imports section (lines 1-10) with:

```python
"""Command implementation for rdkconf."""

import json
import re
import shutil
import subprocess
from pathlib import Path

from chimerax.atomic import AtomicStructure, Bond
from chimerax.atomic.struct_edit import add_atom, add_bond
from chimerax.core.commands import CmdDesc, StringArg, BoolArg, run
from chimerax.core.errors import UserError
from numpy import array, float64
```

**Step 2: Update `_find_script()`**

Change the script filename from `input_to_sdf.py` to `input_to_json.py`:

```python
def _find_script() -> Path:
    """Locate the bundled input_to_json.py script."""
    return Path(__file__).parent.joinpath("input_to_json.py")
```

**Step 3: Remove `_validate_sdf_path()`**

Delete the entire `_validate_sdf_path()` function (lines 35-42 in the current file).

**Step 4: Add `_build_model()` function**

Add this after `_detect_format()`:

```python
def _build_model(session, mol_data, name="UNL"):
    """Build an AtomicStructure from parsed JSON molecule data.

    Parameters
    ----------
    session : chimerax.core.session.Session
    mol_data : dict
        JSON-parsed dict with 'atoms' and 'bonds' keys.
    name : str
        Model and residue name.

    Returns
    -------
    AtomicStructure
    """
    Bond.register_attr(session, "order", "rdkit_conformer", attr_type=float)

    s = AtomicStructure(session, name=name)
    r = s.new_residue(name, " ", 1)

    element_count = {}
    atoms = []
    for atom_info in mol_data["atoms"]:
        elem = atom_info["element"]
        n = element_count.get(elem, 0) + 1
        element_count[elem] = n
        atom_name = f"{elem}{n}"
        xyz = array([atom_info["x"], atom_info["y"], atom_info["z"]], dtype=float64)
        a = add_atom(atom_name, elem, r, xyz)
        atoms.append(a)

    for bond_info in mol_data["bonds"]:
        b = add_bond(atoms[bond_info["begin"]], atoms[bond_info["end"]])
        b.order = bond_info.get("order", 1.0)

    return s
```

**Step 5: Replace `rdkconf()` function**

Replace the entire `rdkconf()` function with:

```python
def rdkconf(session, input_str, format=None, name="UNL", hydrogen=True):
    """Generate 3D conformer from molecular notation using RDKit ETKDGv3.

    Parameters
    ----------
    session : chimerax.core.session.Session
    input_str : str
        Molecular notation string (SMILES by default).
    format : str, optional
        Input format. Auto-detects InChI. Default: smiles.
    name : str
        Residue name (default: UNL).
    hydrogen : bool
        Show hydrogens (default: True).
    """
    fmt = _detect_format(input_str, format)
    name = _validate_name(name)
    script_path = _find_script()
    uv_path = _find_uv()

    cmd_args = [
        uv_path,
        "run",
        "--script",
        str(script_path),
        input_str,
        "--format",
        fmt,
    ]

    try:
        result = subprocess.run(
            cmd_args,
            capture_output=True,
            text=True,
            timeout=60,
        )
    except subprocess.TimeoutExpired:
        raise UserError("RDKit 3D generation timed out (60s)")

    if result.returncode != 0:
        raise UserError(f"RDKit error: {result.stderr.strip()}")

    try:
        mol_data = json.loads(result.stdout)
    except json.JSONDecodeError as e:
        raise UserError(f"Failed to parse RDKit output: {e}")

    model = _build_model(session, mol_data, name)
    session.models.add([model])

    if not hydrogen:
        model_spec = f"#{model.id_string}"
        run(session, f"hide {model_spec} & H")

    session.logger.info(f"Generated 3D conformer ({fmt}) from: {input_str}")
```

**Step 6: Update `rdkconf_desc`**

Remove `output` from the keyword list and remove `SaveFileNameArg` import:

```python
rdkconf_desc = CmdDesc(
    required=[("input_str", StringArg)],
    keyword=[
        ("format", StringArg),
        ("name", StringArg),
        ("hydrogen", BoolArg),
    ],
    synopsis="Generate 3D conformer from molecular notation using RDKit ETKDGv3",
)
```

**Step 7: Remove unused imports**

Remove `os` from imports (was used for `os.unlink`). The `import os` line should be deleted.

**Step 8: Run unit tests**

Run: `uv run --no-project --with rdkit --with pytest pytest tests/ -v`

Expected: All 20 tests PASS (script-level tests don't depend on cmd.py)

**Step 9: Commit**

```bash
git add src/cmd.py
git commit -m "feat: build AtomicStructure directly from JSON, remove SDF intermediary"
```

---

### Task 5: Update pyproject.toml, Smoke Test, and README

**Files:**
- Modify: `pyproject.toml`
- Modify: `scripts/smoke.cxc`
- Modify: `README.md`

**Step 1: Update `pyproject.toml`**

Change version from `0.3.0` to `0.4.0`:
```toml
version = "0.4.0"
```

Change package-data from `input_to_sdf.py` to `input_to_json.py`:
```toml
[chimerax.package-data]
"src/" = ["input_to_json.py"]
```

**Step 2: Update `scripts/smoke.cxc`**

Replace contents with:

```
# Smoke test for ChimeraX-RDKitConformer
# SMILES (default)
rdkconf c1ccccc1
# InChI (auto-detect)
rdkconf "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
# Sequence with explicit format
rdkconf GGG format sequence
# Custom name
rdkconf CCO name EtOH
# Hide hydrogens
rdkconf CCO hydrogen false
log show
```

**Step 3: Update `README.md`**

In the Usage section, remove the `output` example line:
```
rdkconf CCO output ~/ethanol.sdf                   # Save SDF file
```

In the "How It Works" section, update the description to reflect direct model construction:

Replace:
```
ChimeraX's built-in SMILES support (`open smiles:CCO`) depends on NCI's web service
(requires internet, no control over 3D generation quality). This bundle uses RDKit
ETKDGv3 locally via `uv run --script` subprocess for high-quality 3D coordinate
generation without network dependency.
```

With:
```
ChimeraX's built-in SMILES support (`open smiles:CCO`) depends on NCI's web service
(requires internet, no control over 3D generation quality). This bundle uses RDKit
ETKDGv3 locally via `uv run --script` subprocess for high-quality 3D coordinate
generation without network dependency. The 3D structure is built directly in memory
using ChimeraX's AtomicStructure API — no intermediate files are written.
```

**Step 4: Run tests**

Run: `uv run --no-project --with rdkit --with pytest pytest tests/ -v`

Expected: All 20 tests PASS

**Step 5: Commit**

```bash
git add pyproject.toml scripts/smoke.cxc README.md
git commit -m "chore: update package config, smoke test, and README for v0.4.0"
```

---

### Task 6: End-to-End Verification in ChimeraX

This task cannot be automated in unit tests (requires ChimeraX runtime). Perform manually.

**Step 1: Install bundle in ChimeraX**

```
ChimeraX --nogui --exit --cmd 'devel install /path/to/chimerax-rdkit-conformer'
```

Or via MCP:
```
chimerax_run: devel install /path/to/chimerax-rdkit-conformer
```

**Step 2: Verify version**

Check log output shows `v0.4.0` and `input_to_json.py` is bundled.

**Step 3: Test commands**

```
rdkconf c1ccccc1
rdkconf "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
rdkconf GGG format sequence
rdkconf ACGT format dna
rdkconf CCO name EtOH
rdkconf CCO hydrogen false
```

Each should generate a model directly (no temp file flash in output).

**Step 4: Run smoke test**

```
chimerax_run: open /path/to/chimerax-rdkit-conformer/scripts/smoke.cxc
```

**Step 5: Verify no regressions**

- Models appear correctly in 3D view
- Atom names are correct (C1, C2, O1, H1, etc.)
- Bond connectivity looks correct
- Hydrogen hiding works
- Custom name appears in model panel

**Step 6: Commit (if any fixes needed)**

```bash
git commit -m "fix: address issues found during E2E verification"
```
