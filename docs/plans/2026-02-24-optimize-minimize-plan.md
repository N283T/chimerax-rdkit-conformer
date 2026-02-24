# Optimize & Minimize Options Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add `optimize` (RDKit MMFF) and `minimize` (ChimeraX AMBER) as independent boolean flags to the `rdkconf` command, both defaulting to false.

**Architecture:** Two independent optimization stages. `optimize` controls MMFF in the RDKit subprocess (`input_to_json.py`), toggled via a `--optimize` CLI flag. `minimize` runs ChimeraX's `minimize` command (1.11+) after model creation in `cmd.py`, with best-effort fallback on older versions. Breaking change: MMFF is no longer always-on.

**Tech Stack:** RDKit MMFF94, ChimeraX AMBER14/OpenMM (via `minimize` command), Python argparse, pytest

---

### Task 1: Add `optimize` parameter tests to input_to_json (RED)

**Files:**
- Modify: `tests/test_input_to_json.py`

**Step 1: Write the failing tests**

Add a new test class `TestOptimizeFlag` at the end of the file. Also update existing `TestInputToJson` tests that assume MMFF always runs.

```python
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
        with patch("rdkit.Chem.AllChem.MMFFOptimizeMolecule", return_value=0) as mock_mmff:
            script_module.input_to_json("CCO", "smiles", optimize=True)
            mock_mmff.assert_called_once()

    def test_optimize_default_is_false(self, script_module):
        """Default optimize=False means no MMFF."""
        with patch("rdkit.Chem.AllChem.MMFFOptimizeMolecule") as mock_mmff:
            script_module.input_to_json("CCO", "smiles")
            mock_mmff.assert_not_called()

    def test_optimize_true_multi_conformer(self, script_module):
        """optimize=True with multi-conformer calls MMFFOptimizeMoleculeConfs."""
        with patch("rdkit.Chem.AllChem.MMFFOptimizeMoleculeConfs", return_value=[(0, 1.0), (0, 2.0)]) as mock_mmff:
            script_module.input_to_json("CCO", "smiles", num_confs=3, optimize=True)
            mock_mmff.assert_called_once()

    def test_optimize_false_multi_conformer_skips_mmff(self, script_module):
        """optimize=False with multi-conformer skips MMFFOptimizeMoleculeConfs."""
        with patch("rdkit.Chem.AllChem.MMFFOptimizeMoleculeConfs") as mock_mmff:
            script_module.input_to_json("CCO", "smiles", num_confs=3, optimize=False)
            mock_mmff.assert_not_called()
```

Also update these existing tests in `TestInputToJson` that rely on MMFF being called by default. They must now pass `optimize=True`:

- `test_mmff_failure_warns_but_returns_result`: change `script_module.input_to_json("CCO", "smiles")` to `script_module.input_to_json("CCO", "smiles", optimize=True)`
- `test_mmff_non_convergence_warns`: same change — add `optimize=True`
- `test_embed_failure_raises_runtime_error`: no change needed (embedding, not MMFF)

**Step 2: Run tests to verify they fail**

Run: `uv run --with pytest --with rdkit --no-project -- pytest tests/test_input_to_json.py::TestOptimizeFlag -v`
Expected: FAIL — `input_to_json()` doesn't accept `optimize` parameter yet.

**Step 3: Commit**

```bash
git add tests/test_input_to_json.py
git commit -m "test: add optimize flag tests for input_to_json (RED)"
```

---

### Task 2: Implement `optimize` parameter in input_to_json.py (GREEN)

**Files:**
- Modify: `src/input_to_json.py:118-191` (function `input_to_json`) and `src/input_to_json.py:194-228` (function `main`)

**Step 1: Update `input_to_json()` function signature and body**

Change the function signature at line 118-120:

```python
def input_to_json(
    input_str: str, fmt: str = "smiles", num_confs: int = 1, optimize: bool = False
) -> list[dict]:
```

Update the docstring to document the new `optimize` parameter:

```python
    """Convert molecular notation to 3D JSON dicts using ETKDGv3.

    Generates 3D coordinates via ETKDGv3 embedding. Optionally optimizes
    geometry using MMFF force field when optimize=True.

    Args:
        input_str: Molecular notation string.
        fmt: Input format (default: smiles).
        num_confs: Number of conformers to generate (default: 1, max: 50).
        optimize: Run MMFF force field optimization (default: False).

    Returns:
        List of dicts, each with 'atoms' and 'bonds' keys.

    Raises:
        ValueError: If format is unsupported, input is invalid, or num_confs is out of range.
        RuntimeError: If 3D coordinate generation fails.
    """
```

In the single-conformer path (around line 150-169), wrap the MMFF block with `if optimize:`:

```python
    if num_confs == 1:
        result = AllChem.EmbedMolecule(mol, params)
        if result == -1:
            raise RuntimeError(f"Failed to generate 3D coordinates for: {input_str}")

        if optimize:
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
```

In the multi-conformer path (around line 170-191), wrap the MMFF block with `if optimize:`:

```python
    else:
        params.pruneRmsThresh = 0.5
        conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
        if len(conf_ids) == 0:
            raise RuntimeError(f"Failed to generate any 3D conformers for: {input_str}")

        if optimize:
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

        return [mol_to_json(mol, conf_id=conf_id) for conf_id in conf_ids]
```

**Step 2: Update `main()` to add `--optimize` CLI flag**

Add after the `--conformers` argument (around line 212):

```python
    parser.add_argument(
        "--optimize",
        action="store_true",
        default=False,
        help="Run MMFF force field optimization (default: off)",
    )
```

Update the `input_to_json()` call in `main()` (around line 216):

```python
        data = input_to_json(args.input, args.format, num_confs=args.conformers, optimize=args.optimize)
```

**Step 3: Run all tests**

Run: `uv run --with pytest --with rdkit --no-project -- pytest tests/test_input_to_json.py -v`
Expected: ALL PASS (35+ tests including new TestOptimizeFlag tests)

**Step 4: Commit**

```bash
git add src/input_to_json.py tests/test_input_to_json.py
git commit -m "feat: add optimize parameter to input_to_json()"
```

---

### Task 3: Add `optimize` and `minimize` keywords to cmd.py (with subprocess wiring)

**Files:**
- Modify: `src/cmd.py:156-265`

**Step 1: Add `optimize` parameter to `rdkconf()` and wire to subprocess**

Update the function signature at line 156:

```python
def rdkconf(session, input_str, format=None, name="UNL", hydrogen=True, conformers=1, optimize=False, minimize=False):
```

Update the docstring to document both new parameters.

In the `cmd_args` list (around line 196-207), conditionally add `--optimize`:

```python
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
    if optimize:
        cmd_args.append("--optimize")
```

**Step 2: Add `minimize` post-processing after model creation**

After `session.models.add(models)` (around line 247), add minimize logic:

```python
    if minimize:
        for model in models:
            model_spec = f"#{model.id_string}"
            try:
                run(session, f"minimize {model_spec} maxSteps 1000")
            except Exception:
                session.logger.warning(
                    f"minimize command failed for {model_spec} "
                    "(requires ChimeraX 1.11+); skipping AMBER optimization"
                )
```

**Step 3: Update `rdkconf_desc` to include new keywords**

Update the keyword list (around line 266-273):

```python
rdkconf_desc = CmdDesc(
    required=[("input_str", StringArg)],
    keyword=[
        ("format", StringArg),
        ("name", StringArg),
        ("hydrogen", BoolArg),
        ("conformers", IntArg),
        ("optimize", BoolArg),
        ("minimize", BoolArg),
    ],
    synopsis="Generate 3D conformer(s) from molecular notation using RDKit ETKDGv3",
)
```

**Step 4: Run tests to verify no regressions**

Run: `uv run --with pytest --with rdkit --no-project -- pytest tests/ -v`
Expected: ALL PASS

**Step 5: Commit**

```bash
git add src/cmd.py
git commit -m "feat: add optimize and minimize keywords to rdkconf command"
```

---

### Task 4: Update docs, smoke test, and version bump

**Files:**
- Modify: `README.md`
- Modify: `scripts/smoke.cxc`
- Modify: `pyproject.toml` (version bump)

**Step 1: Update README.md usage section**

Add optimize/minimize examples to the usage block:

```
rdkconf CCO optimize true                  # RDKit MMFF optimization
rdkconf CCO minimize true                  # ChimeraX AMBER minimization (1.11+)
rdkconf CCO optimize true minimize true    # Both: MMFF + AMBER
```

Update the "How It Works" section to mention the optional optimization stages.

Update the Requirements section: ChimeraX 1.6+ (1.11+ for `minimize` keyword).

**Step 2: Update smoke test**

Add to `scripts/smoke.cxc`:

```
# Optimize with MMFF
rdkconf CCO optimize true
# Minimize with AMBER (ChimeraX 1.11+ only)
rdkconf CCO minimize true
```

**Step 3: Bump version**

In `pyproject.toml`, bump version from `0.5.0` to `0.6.0`.

**Step 4: Commit**

```bash
git add README.md scripts/smoke.cxc pyproject.toml
git commit -m "docs: update README and smoke test for optimize/minimize, bump to v0.6.0"
```

---

### Task 5: E2E verification in ChimeraX

**Files:** None (manual verification)

**Step 1: Build and install**

```bash
echidna build && echidna install --chimerax /Applications/ChimeraX-1.11.app/Contents/bin/ChimeraX
```

Or if using ChimeraX 1.10:

```bash
echidna build && echidna install --chimerax /Applications/ChimeraX-1.10.app/Contents/bin/ChimeraX
```

**Step 2: Start ChimeraX and test**

Start ChimeraX fresh (full restart if previously running).

Test these commands in sequence:

```
# Default: ETKDGv3 only (no optimization)
rdkconf CCO

# MMFF optimization
rdkconf CCO optimize true name EtOH_mmff

# AMBER minimization (1.11+ only; expect warning on 1.10)
rdkconf CCO minimize true name EtOH_amber

# Both optimizations
rdkconf CCO optimize true minimize true name EtOH_both

# Multi-conformer with minimize
rdkconf c1ccccc1 conformers 3 minimize true name BEN

# Verify error cases still work
rdkconf CCO conformers 0
```

**Step 3: Verify expected behavior**

- Default: model created, no optimization messages
- `optimize true`: model created, MMFF runs in subprocess
- `minimize true` on 1.11+: model created, AMBER minimization runs
- `minimize true` on 1.10: model created, warning logged about ChimeraX version
- Both: model created, MMFF then AMBER both run
- Multi-conformer + minimize: each conformer model minimized independently
- Error cases: UserError with range guidance
