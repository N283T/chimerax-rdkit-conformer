# Multi-Format Input Support Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extend `rdkconf` to accept InChI, FASTA, Sequence, HELM, DNA, and RNA inputs in addition to SMILES, with InChI auto-detection.

**Architecture:** Add `parse_input()` dispatcher to the standalone RDKit script (renamed `input_to_sdf.py`), accept `--format` CLI flag. In `cmd.py`, add `format` keyword with auto-detect logic for InChI. All formats follow the same pipeline: parse -> AddHs -> ETKDGv3 embed -> MMFF optimize -> SDF.

**Tech Stack:** Python 3.12, RDKit (MolFromSmiles, MolFromInchi, MolFromFASTA, MolFromSequence, MolFromHELM), ChimeraX BundleBuilder, pytest

**Design doc:** `docs/plans/2026-02-23-multi-format-input-design.md`

---

### Task 1: Set up test infrastructure and write parse_input tests

**Files:**
- Create: `tests/conftest.py`
- Create: `tests/test_input_to_sdf.py`

**Step 1: Create test fixture for importing the script module**

The script uses PEP 723 metadata and chimerax-incompatible imports, so use `importlib.util` to load it directly.

```python
# tests/conftest.py
import importlib.util

import pytest
from pathlib import Path


@pytest.fixture
def script_module():
    """Import input_to_sdf.py as a module (bypassing PEP 723 metadata)."""
    script_path = Path(__file__).parent.parent.joinpath("src", "input_to_sdf.py")
    spec = importlib.util.spec_from_file_location("input_to_sdf", str(script_path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module
```

**Step 2: Write parse_input tests**

```python
# tests/test_input_to_sdf.py
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
```

**Step 3: Run tests to verify they fail**

Run: `uv run --with rdkit --with pytest pytest tests/ -v`
Expected: FAIL — `input_to_sdf.py` does not exist yet, `parse_input` not defined

**Step 4: Commit test scaffolding**

```bash
git add tests/conftest.py tests/test_input_to_sdf.py
git commit -m "test: add parse_input tests for multi-format support (RED)"
```

---

### Task 2: Create input_to_sdf.py with parse_input

**Files:**
- Delete: `src/smiles_to_sdf.py`
- Create: `src/input_to_sdf.py`

**Step 1: Create input_to_sdf.py with parse_input**

Rename `smiles_to_sdf.py` to `input_to_sdf.py`. Add `parse_input()` dispatcher. Keep `smiles_to_3d` as `input_to_3d` (generalized). All formats need `AddHs()` before `EmbedMolecule()`.

```python
#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "rdkit",
# ]
# ///
"""Convert molecular notation to 3D SDF using RDKit ETKDGv3."""

import argparse
import sys
import tempfile
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

VALID_FORMATS = ("smiles", "inchi", "fasta", "sequence", "helm", "dna", "rna")


def parse_input(input_str: str, fmt: str = "smiles") -> Chem.Mol:
    """Parse molecular notation string into an RDKit Mol object.

    Args:
        input_str: Molecular notation string.
        fmt: Input format. One of: smiles, inchi, fasta, sequence, helm, dna, rna.

    Returns:
        RDKit Mol object.

    Raises:
        ValueError: If format is unsupported or input is invalid.
    """
    parsers = {
        "smiles": lambda s: Chem.MolFromSmiles(s),
        "inchi": lambda s: Chem.MolFromInchi(s),
        "fasta": lambda s: Chem.MolFromFASTA(s),
        "sequence": lambda s: Chem.MolFromSequence(s),
        "helm": lambda s: Chem.MolFromHELM(s),
        "dna": lambda s: Chem.MolFromSequence(s, flavor=6),
        "rna": lambda s: Chem.MolFromSequence(s, flavor=2),
    }
    parser = parsers.get(fmt)
    if parser is None:
        raise ValueError(f"Unsupported format: {fmt}")
    mol = parser(input_str)
    if mol is None:
        raise ValueError(f"Invalid {fmt} input: {input_str}")
    return mol


def input_to_3d(input_str: str, fmt: str = "smiles", output: Path | None = None) -> Path:
    """Convert molecular notation to a 3D SDF file using ETKDGv3.

    Args:
        input_str: Molecular notation string.
        fmt: Input format (default: smiles).
        output: Output file path. If None, creates a temp file.

    Returns:
        Path to the generated SDF file.
    """
    mol = parse_input(input_str, fmt)
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        raise RuntimeError(f"Failed to generate 3D coordinates for: {input_str}")

    AllChem.MMFFOptimizeMolecule(mol)

    if output is None:
        tmp = tempfile.NamedTemporaryFile(suffix=".sdf", delete=False)
        output = Path(tmp.name)
        tmp.close()

    with Chem.SDWriter(str(output)) as writer:
        writer.write(mol)

    return output


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert molecular notation to 3D SDF")
    parser.add_argument("input", help="Molecular notation string")
    parser.add_argument(
        "-f", "--format", default="smiles", choices=VALID_FORMATS,
        help="Input format (default: smiles)",
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=None, help="Output SDF file path"
    )
    args = parser.parse_args()

    try:
        path = input_to_3d(args.input, args.format, args.output)
    except (ValueError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    print(str(path))


if __name__ == "__main__":
    main()
```

**Step 2: Run parse_input tests to verify they pass**

Run: `uv run --with rdkit --with pytest pytest tests/test_input_to_sdf.py::TestParseInput -v`
Expected: All 9 tests PASS

**Step 3: Commit**

```bash
git rm src/smiles_to_sdf.py
git add src/input_to_sdf.py
git commit -m "feat: add parse_input dispatcher with multi-format support"
```

---

### Task 3: Write input_to_3d integration tests

**Files:**
- Modify: `tests/test_input_to_sdf.py`

**Step 1: Add integration tests for 3D generation**

Append to `tests/test_input_to_sdf.py`:

```python
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
```

**Step 2: Run tests to verify they pass**

Run: `uv run --with rdkit --with pytest pytest tests/test_input_to_sdf.py -v`
Expected: All tests PASS (input_to_3d already implemented in Task 2)

**Step 3: Commit**

```bash
git add tests/test_input_to_sdf.py
git commit -m "test: add input_to_3d integration tests"
```

---

### Task 4: Update cmd.py with format support and auto-detect

**Files:**
- Modify: `src/cmd.py:13-15` (script reference)
- Modify: `src/cmd.py:45-97` (rdkconf function)
- Modify: `src/cmd.py:100-108` (CmdDesc)

**Step 1: Update _find_script to reference new filename**

Change `src/cmd.py:13-15`:

```python
def _find_script() -> Path:
    """Locate the bundled input_to_sdf.py script."""
    return Path(__file__).parent.joinpath("input_to_sdf.py")
```

**Step 2: Add _detect_format function**

Add after `_validate_sdf_path` (after line 42):

```python
_VALID_FORMATS = {"smiles", "inchi", "fasta", "sequence", "helm", "dna", "rna"}


def _detect_format(input_str: str, user_format: str | None) -> str:
    """Detect or validate input format.

    Auto-detects InChI (starts with 'InChI='). Defaults to SMILES.
    """
    if user_format is not None:
        fmt = user_format.lower()
        if fmt not in _VALID_FORMATS:
            raise UserError(
                f"Invalid format: {user_format!r}. "
                f"Valid formats: {', '.join(sorted(_VALID_FORMATS))}"
            )
        return fmt
    if input_str.startswith("InChI="):
        return "inchi"
    return "smiles"
```

**Step 3: Update rdkconf function to use format**

Replace `rdkconf` function:

```python
def rdkconf(session, input_str, format=None, output=None, name="UNL", hydrogen=True):
    """Generate 3D conformer from molecular notation using RDKit ETKDGv3.

    Parameters
    ----------
    session : chimerax.core.session.Session
    input_str : str
        Molecular notation string (SMILES by default).
    format : str, optional
        Input format. Auto-detects InChI. Default: smiles.
    output : str, optional
        Save SDF to this path.
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
        uv_path, "run", "--script", str(script_path),
        input_str, "--format", fmt,
    ]
    if output:
        cmd_args.extend(["-o", output])

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

    sdf_path = _validate_sdf_path(result.stdout.strip())

    try:
        open_cmd = f'open "{sdf_path}" name "{name}"'
        models = run(session, open_cmd)

        if not hydrogen and models:
            model_spec = f"#{models[0].id_string}"
            run(session, f"hide {model_spec} & H")
    finally:
        if output is None:
            try:
                os.unlink(sdf_path)
            except OSError:
                pass

    session.logger.info(f"Generated 3D conformer ({fmt}) from: {input_str}")
```

**Step 4: Update CmdDesc to include format keyword**

```python
rdkconf_desc = CmdDesc(
    required=[("input_str", StringArg)],
    keyword=[
        ("format", StringArg),
        ("output", SaveFileNameArg),
        ("name", StringArg),
        ("hydrogen", BoolArg),
    ],
    synopsis="Generate 3D conformer from molecular notation using RDKit ETKDGv3",
)
```

**Step 5: Commit**

```bash
git add src/cmd.py
git commit -m "feat: add format keyword with InChI auto-detect to rdkconf command"
```

---

### Task 5: Update pyproject.toml and smoke test

**Files:**
- Modify: `pyproject.toml:7` (version)
- Modify: `pyproject.toml:25` (package-data)
- Modify: `scripts/smoke.cxc`

**Step 1: Update pyproject.toml**

Change version from `0.2.0` to `0.3.0`.

Change package-data from `"src/" = ["smiles_to_sdf.py"]` to `"src/" = ["input_to_sdf.py"]`.

```toml
[project]
name = "ChimeraX-RDKitConformer"
version = "0.3.0"
```

```toml
[chimerax.package-data]
"src/" = ["input_to_sdf.py"]
```

**Step 2: Update smoke test**

```cxc
# Smoke test for ChimeraX-RDKitConformer
# SMILES (default)
rdkconf c1ccccc1
# InChI (auto-detect)
rdkconf "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"
# Sequence with explicit format
rdkconf GGG format sequence
log show
```

**Step 3: Commit**

```bash
git add pyproject.toml scripts/smoke.cxc
git commit -m "chore: update package-data for renamed script, bump version to 0.3.0"
```

---

### Task 6: Update README

**Files:**
- Modify: `README.md`

**Step 1: Update Usage section**

Replace the Usage section with:

````markdown
## Usage

```
rdkconf CCO                                        # SMILES (default)
rdkconf c1ccccc1                                   # Benzene
rdkconf "CC(=O)Oc1ccccc1C(=O)O"                   # Aspirin
rdkconf "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"      # InChI (auto-detected)
rdkconf GGG format sequence                        # Peptide sequence
rdkconf ACGT format dna                            # DNA sequence
rdkconf ACGU format rna                            # RNA sequence
rdkconf ">seq\nGGG" format fasta                   # FASTA format
rdkconf "PEPTIDE1{G.G.G}$$$$" format helm          # HELM notation
rdkconf CCO name EtOH                              # Custom residue name
rdkconf CCO hydrogen false                         # Hide hydrogens
rdkconf CCO output ~/ethanol.sdf                   # Save SDF file
```

### Supported Formats

| Format | `format` keyword | Auto-detected |
|--------|-----------------|---------------|
| SMILES | `smiles` (default) | Yes (fallback) |
| InChI | `inchi` | Yes (`InChI=` prefix) |
| Protein sequence | `sequence` | No |
| DNA sequence | `dna` | No |
| RNA sequence | `rna` | No |
| FASTA | `fasta` | No |
| HELM | `helm` | No |
````

**Step 2: Commit**

```bash
git add README.md
git commit -m "docs: update README with multi-format usage examples"
```

---

### Post-Implementation Checklist

- [ ] All unit tests pass: `uv run --with rdkit --with pytest pytest tests/ -v`
- [ ] Smoke test passes in ChimeraX: `echidna run --script scripts/smoke.cxc`
- [ ] README updated with all formats
- [ ] Version bumped to 0.3.0
- [ ] Create PR to `feature/chimerax-bundle`
