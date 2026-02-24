# Multi-Conformer Generation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add `conformers N` keyword to `rdkconf` command to generate multiple 3D conformers as separate ChimeraX models.

**Architecture:** `input_to_json.py` gains `num_confs` param using `EmbedMultipleConfs` + `MMFFOptimizeMoleculeConfs`. Output changes from a single dict to a JSON array. `cmd.py` loops over the array calling `_build_model()` per conformer, naming models `{name}_1`, `{name}_2`, etc.

**Tech Stack:** RDKit (EmbedMultipleConfs, MMFFOptimizeMoleculeConfs), ChimeraX AtomicStructure API, pytest

---

### Task 1: Add multi-conformer tests for `input_to_json()`

**Files:**
- Modify: `tests/test_input_to_json.py`

**Context:** `input_to_json()` currently returns a single dict `{"atoms":[], "bonds":[]}`. After this feature, it will return a list of such dicts. We need to update existing tests and add new ones.

**Step 1: Write failing tests for multi-conformer `input_to_json()`**

Add a new test class `TestMultiConformer` at the end of `tests/test_input_to_json.py`:

```python
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
        result = script_module.input_to_json("c1ccccc1", "smiles", num_confs=5)
        if len(result) < 2:
            pytest.skip("RMS pruning left only 1 conformer")
        coords_0 = [(a["x"], a["y"], a["z"]) for a in result[0]["atoms"]]
        coords_1 = [(a["x"], a["y"], a["z"]) for a in result[1]["atoms"]]
        assert coords_0 != coords_1

    def test_pruning_may_return_fewer_conformers(self, script_module):
        """Requesting many conformers of a rigid molecule may yield fewer due to RMS pruning."""
        # Methane is extremely rigid — most conformers will be pruned
        result = script_module.input_to_json("C", "smiles", num_confs=10)
        assert isinstance(result, list)
        assert len(result) >= 1
        assert len(result) <= 10
```

**Step 2: Run tests to verify they fail**

Run: `uv run --no-project --with rdkit --with pytest pytest tests/test_input_to_json.py::TestMultiConformer -v`
Expected: FAIL — `input_to_json()` returns a dict, not a list. `TypeError` on `isinstance(result, list)` or `num_confs` unexpected keyword.

**Step 3: Commit failing tests**

```bash
git add tests/test_input_to_json.py
git commit -m "test: add multi-conformer tests for input_to_json (RED)"
```

---

### Task 2: Update existing tests for JSON array format

**Files:**
- Modify: `tests/test_input_to_json.py`

**Context:** The return type of `input_to_json()` is changing from `dict` to `list[dict]`. All existing `TestInputToJson` tests need to index into `result[0]` instead of using `result` directly.

**Step 1: Update all TestInputToJson tests**

In `TestInputToJson`, update each test that calls `input_to_json()`:

- `test_smiles_returns_json`: Change `result` assertions to `result[0]`
- `test_inchi_returns_json`: Change `result` assertions to `result[0]`
- `test_sequence_returns_json`: Change `result` assertions to `result[0]`
- `test_invalid_input_raises`: No change (exception before return)
- `test_json_round_trip`: Update to use `result[0]`
- `test_embed_failure_raises_runtime_error`: No change (exception before return)
- `test_mmff_failure_warns_but_returns_result`: Update to use `result[0]`
- `test_mmff_non_convergence_warns`: Update to use `result[0]`

Specific changes:

```python
# test_smiles_returns_json
result = script_module.input_to_json("CCO", "smiles")
assert isinstance(result, list)
assert len(result) == 1
assert "atoms" in result[0]
assert "bonds" in result[0]
assert len(result[0]["atoms"]) > 0
assert len(result[0]["bonds"]) > 0

# test_inchi_returns_json
result = script_module.input_to_json(inchi, "inchi")
assert len(result[0]["atoms"]) > 0

# test_sequence_returns_json
result = script_module.input_to_json("GGG", "sequence")
assert len(result[0]["atoms"]) > 0

# test_json_round_trip
result = script_module.input_to_json("CCO", "smiles")
json_str = json.dumps(result)
parsed = json.loads(json_str)
assert parsed[0]["atoms"] == result[0]["atoms"]
assert parsed[0]["bonds"] == result[0]["bonds"]

# test_mmff_failure_warns_but_returns_result
result = script_module.input_to_json("CCO", "smiles")
# ...
assert "atoms" in result[0]
assert len(result[0]["atoms"]) > 0

# test_mmff_non_convergence_warns
result = script_module.input_to_json("CCO", "smiles")
# ...
assert "atoms" in result[0]
```

**Step 2: Run tests to verify they still fail**

Run: `uv run --no-project --with rdkit --with pytest pytest tests/test_input_to_json.py::TestInputToJson -v`
Expected: FAIL — `input_to_json()` still returns dict, not list.

**Step 3: Commit**

```bash
git add tests/test_input_to_json.py
git commit -m "test: update existing tests for JSON array return type (RED)"
```

---

### Task 3: Implement multi-conformer in `input_to_json.py`

**Files:**
- Modify: `src/input_to_json.py:114-157` (input_to_json function)
- Modify: `src/input_to_json.py:160-183` (main function)

**Context:** `input_to_json()` must gain a `num_confs` param. For `num_confs=1`, keep using `EmbedMolecule` (deterministic with seed 42). For `num_confs>1`, use `EmbedMultipleConfs` with `pruneRmsThresh=0.5`. Always return a list of dicts. `main()` needs `--conformers` CLI arg.

**Step 1: Update `input_to_json()` signature and implementation**

Replace the `input_to_json` function (lines 114-157) with:

```python
_MAX_CONFORMERS = 50


def input_to_json(input_str: str, fmt: str = "smiles", num_confs: int = 1) -> list[dict]:
    """Convert molecular notation to 3D JSON dicts using ETKDGv3.

    Generates 3D coordinates via ETKDGv3 embedding, then optimizes geometry
    using MMFF force field. If MMFF optimization fails or does not converge,
    a warning is emitted and coordinates are returned as-is.

    Args:
        input_str: Molecular notation string.
        fmt: Input format (default: smiles).
        num_confs: Number of conformers to generate (default: 1, max: 50).

    Returns:
        List of dicts, each with 'atoms' and 'bonds' keys.

    Raises:
        ValueError: If format is unsupported, input is invalid, or num_confs is out of range.
        RuntimeError: If 3D coordinate generation fails.
    """
    if num_confs < 1 or num_confs > _MAX_CONFORMERS:
        raise ValueError(
            f"conformers must be between 1 and {_MAX_CONFORMERS}, got {num_confs}"
        )

    mol = parse_input(input_str, fmt)
    mol = Chem.AddHs(mol)

    import warnings

    params = AllChem.ETKDGv3()
    params.randomSeed = 42

    if num_confs == 1:
        result = AllChem.EmbedMolecule(mol, params)
        if result == -1:
            raise RuntimeError(
                f"Failed to generate 3D coordinates for: {input_str}"
            )

        opt_result = AllChem.MMFFOptimizeMolecule(mol)
        if opt_result == -1:
            warnings.warn(
                f"MMFF force field setup failed for: {input_str}. "
                "3D coordinates are unoptimized.",
                stacklevel=2,
            )
        elif opt_result == 1:
            warnings.warn(
                f"MMFF optimization did not converge for: {input_str}. "
                "3D coordinates may have suboptimal geometry.",
                stacklevel=2,
            )

        return [mol_to_json(mol)]
    else:
        params.pruneRmsThresh = 0.5
        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
        if len(conf_ids) == 0:
            raise RuntimeError(
                f"Failed to generate any 3D conformers for: {input_str}"
            )

        opt_results = AllChem.MMFFOptimizeMoleculeConfs(mol)
        for i, (converged, _energy) in enumerate(opt_results):
            if converged == -1:
                warnings.warn(
                    f"MMFF force field setup failed for conformer {i} of: "
                    f"{input_str}. 3D coordinates are unoptimized.",
                    stacklevel=2,
                )
            elif converged == 1:
                warnings.warn(
                    f"MMFF optimization did not converge for conformer {i} of: "
                    f"{input_str}. 3D coordinates may have suboptimal geometry.",
                    stacklevel=2,
                )

        conformers = []
        for conf_id in conf_ids:
            conformers.append(mol_to_json(mol, conf_id=conf_id))
        return conformers
```

**Step 2: Update `mol_to_json()` to accept optional `conf_id`**

Change the `mol_to_json` signature and the `GetConformer` call (line 71):

```python
def mol_to_json(mol: Chem.Mol, conf_id: int = -1) -> dict:
    """Serialize an RDKit Mol with 3D coordinates to a JSON-compatible dict.

    Args:
        mol: RDKit Mol object. Must have at least one conformer
            (embedded 3D coordinates). Typically with explicit Hs
            added via Chem.AddHs().
        conf_id: Conformer ID to serialize (default: -1, meaning first/only).

    Returns:
        Dict with 'atoms' (list of element/x/y/z) and 'bonds' (list of begin/end/order).

    Raises:
        ValueError: If mol has no conformer.
    """
    if mol.GetNumConformers() == 0:
        raise ValueError(
            "mol_to_json requires a Mol with 3D coordinates. "
            "Call AllChem.EmbedMolecule first."
        )

    conf = mol.GetConformer(conf_id)
    # ... rest unchanged
```

**Step 3: Update `main()` to add `--conformers` arg**

Add the argument after `--format` and pass it through:

```python
def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert molecular notation to 3D JSON"
    )
    parser.add_argument("input", help="Molecular notation string")
    parser.add_argument(
        "-f",
        "--format",
        default="smiles",
        choices=VALID_FORMATS,
        help="Input format (default: smiles)",
    )
    parser.add_argument(
        "-c",
        "--conformers",
        type=int,
        default=1,
        help="Number of conformers to generate (default: 1, max: 50)",
    )
    args = parser.parse_args()

    try:
        data = input_to_json(args.input, args.format, num_confs=args.conformers)
    except (ValueError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {type(e).__name__}: {e}", file=sys.stderr)
        sys.exit(2)

    print(json.dumps(data))
```

**Step 4: Run all tests**

Run: `uv run --no-project --with rdkit --with pytest pytest tests/test_input_to_json.py -v`
Expected: ALL PASS (27 existing updated + 6 new = 33 total)

**Step 5: Commit**

```bash
git add src/input_to_json.py
git commit -m "feat: add multi-conformer support to input_to_json()"
```

---

### Task 4: Update `cmd.py` for multi-conformer support

**Files:**
- Modify: `src/cmd.py:16` (add IntArg import)
- Modify: `src/cmd.py:153-218` (rdkconf function)
- Modify: `src/cmd.py:221-229` (rdkconf_desc)

**Context:** `rdkconf()` needs a `conformers` keyword (default=1, max=50). It must pass `--conformers N` to the subprocess, parse the JSON array, loop `_build_model()` per entry, and name models `{name}_1`, `{name}_2`... for multiple.

**Step 1: Add `IntArg` import**

Change line 16:
```python
from chimerax.core.commands import CmdDesc, StringArg, BoolArg, IntArg, run
```

**Step 2: Add `_MAX_CONFORMERS` constant after `_VALID_FORMATS`**

```python
_MAX_CONFORMERS = 50
```

**Step 3: Rewrite `rdkconf()` function**

Replace the `rdkconf` function (lines 153-218) with:

```python
def rdkconf(session, input_str, format=None, name="UNL", hydrogen=True, conformers=1):
    """Generate 3D conformer(s) from molecular notation using RDKit ETKDGv3.

    Parameters
    ----------
    session : chimerax.core.session.Session
    input_str : str
        Molecular notation string (SMILES by default).
    format : str, optional
        Input format. Auto-detects InChI. Default: smiles.
    name : str
        Model and residue name (default: UNL).
    hydrogen : bool
        Show hydrogens (default: True).
    conformers : int
        Number of conformers to generate (default: 1, max: 50).

    Raises
    ------
    UserError
        If uv is not found, input format is invalid, name is invalid,
        conformers is out of range, the RDKit subprocess fails or times out,
        or the output is unparseable.
    """
    fmt = _detect_format(input_str, format)
    name = _validate_name(name)

    if conformers < 1 or conformers > _MAX_CONFORMERS:
        raise UserError(
            f"conformers must be between 1 and {_MAX_CONFORMERS}, got {conformers}"
        )

    script_path = _find_script()
    uv_path = _find_uv()

    cmd_args = [
        uv_path,
        "run",
        "--no-project",
        "--script",
        str(script_path),
        input_str,
        "--format",
        fmt,
        "--conformers",
        str(conformers),
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

    if result.stderr.strip():
        session.logger.warning(f"RDKit warning: {result.stderr.strip()}")

    try:
        conformer_list = json.loads(result.stdout)
    except json.JSONDecodeError as e:
        raise UserError(f"Failed to parse RDKit output: {e}")

    if not isinstance(conformer_list, list) or len(conformer_list) == 0:
        raise UserError("RDKit output contains no conformers")

    models = []
    for i, mol_data in enumerate(conformer_list):
        if len(conformer_list) == 1:
            model_name = name
        else:
            model_name = f"{name}_{i + 1}"
        model = _build_model(session, mol_data, model_name)
        models.append(model)

    session.models.add(models)

    if not hydrogen:
        for model in models:
            model_spec = f"#{model.id_string}"
            run(session, f"hide {model_spec} & H")

    actual_count = len(conformer_list)
    if conformers > 1 and actual_count < conformers:
        session.logger.info(
            f"Generated {actual_count} of {conformers} requested conformers "
            f"({fmt}) from: {input_str} (duplicates pruned by RMS threshold)"
        )
    else:
        session.logger.info(
            f"Generated {actual_count} conformer(s) ({fmt}) from: {input_str}"
        )
```

**Step 4: Update `rdkconf_desc`**

Replace (lines 221-229):

```python
rdkconf_desc = CmdDesc(
    required=[("input_str", StringArg)],
    keyword=[
        ("format", StringArg),
        ("name", StringArg),
        ("hydrogen", BoolArg),
        ("conformers", IntArg),
    ],
    synopsis="Generate 3D conformer(s) from molecular notation using RDKit ETKDGv3",
)
```

**Step 5: Run tests**

Run: `uv run --no-project --with rdkit --with pytest pytest tests/test_input_to_json.py -v`
Expected: ALL PASS

**Step 6: Commit**

```bash
git add src/cmd.py
git commit -m "feat: add conformers keyword to rdkconf command"
```

---

### Task 5: Update docs, smoke test, and version

**Files:**
- Modify: `pyproject.toml:7` (version bump to 0.5.0)
- Modify: `scripts/smoke.cxc` (add multi-conformer test)
- Modify: `README.md` (add conformers usage example)

**Step 1: Bump version in `pyproject.toml`**

Change line 7:
```toml
version = "0.5.0"
```

**Step 2: Add multi-conformer smoke test**

Add before `log show` in `scripts/smoke.cxc`:
```
# Multiple conformers
rdkconf c1ccccc1 conformers 3
```

**Step 3: Update README.md usage section**

Add to the usage code block, after the `hydrogen false` example:
```
rdkconf c1ccccc1 conformers 5                     # Multiple conformers
rdkconf CCO conformers 10 name EtOH               # Named multi-conformers
```

**Step 4: Commit**

```bash
git add pyproject.toml scripts/smoke.cxc README.md
git commit -m "docs: update README and smoke test for multi-conformer, bump to v0.5.0"
```

---

### Task 6: E2E verification in ChimeraX

**Context:** Build and install the updated bundle in ChimeraX, then run the smoke test to verify multi-conformer generation works end-to-end.

**Step 1: Build and install**

```bash
echidna install
```

If ChimeraX is running, stop and restart it to clear module cache.

**Step 2: Run smoke test commands**

In ChimeraX (or via MCP):
```
rdkconf c1ccccc1 conformers 3
```

Expected: Multiple models appear (e.g., #1 `UNL_1`, #2 `UNL_2`, #3 `UNL_3` — or fewer due to RMS pruning).

**Step 3: Verify single conformer still works**

```
rdkconf CCO
```

Expected: Single model named `UNL` (no `_1` suffix).

**Step 4: Verify named multi-conformers**

```
rdkconf CCO conformers 5 name EtOH
```

Expected: Models named `EtOH_1`, `EtOH_2`, etc.

**Step 5: Verify error handling**

```
rdkconf CCO conformers 0
rdkconf CCO conformers 51
```

Expected: Both should show UserError with range guidance.

**Step 6: Commit any fixes found during E2E**

```bash
git add -A
git commit -m "fix: address issues found during E2E verification"
```

Only commit this step if fixes were needed.
