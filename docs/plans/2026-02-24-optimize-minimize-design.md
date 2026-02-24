# Optimize & Minimize Options Design

## Goal

Add two independent optimization flags to `rdkconf`: `optimize` (RDKit MMFF in subprocess) and `minimize` (ChimeraX AMBER post-processing). Both default to false, making ETKDGv3-only the new default.

## Command Interface

```
rdkconf CCO                                # ETKDGv3 only (new default)
rdkconf CCO optimize true                  # ETKDGv3 + RDKit MMFF
rdkconf CCO minimize true                  # ETKDGv3 + ChimeraX AMBER
rdkconf CCO optimize true minimize true    # ETKDGv3 + MMFF + AMBER
rdkconf CCO minimize true conformers 5    # Multi-conformer + AMBER
```

New keywords: `optimize` (BoolArg, default False), `minimize` (BoolArg, default False).

Breaking change: MMFF was previously always-on. Now requires explicit `optimize true`.

## Architecture

### Data Flow

```
rdkconf CCO optimize true minimize true conformers 3
  |
  v  subprocess (input_to_json.py)
  input_to_json.py "CCO" --format smiles --conformers 3 --optimize
  |  (--optimize flag controls MMFF; omitted = skip MMFF)
  |
  v  stdout: JSON array
  [{atoms, bonds}, {atoms, bonds}, {atoms, bonds}]
  |
  v  cmd.py: _build_model() per conformer → session.models.add()
  |
  v  if minimize: for each model → run(session, "minimize #N maxSteps 1000")
  |  (ChimeraX 1.11+ required; older versions get a warning and skip)
  |
  Done
```

Two optimization steps are fully independent:
- **optimize** (MMFF): runs in the RDKit subprocess, before JSON output
- **minimize** (AMBER): runs in ChimeraX after model creation

### Changes

#### input_to_json.py

- `input_to_json()` gains `optimize: bool = False` parameter.
- `optimize=False`: skip MMFF entirely, return raw ETKDGv3 coordinates.
- `optimize=True`: current MMFF behavior (single: `MMFFOptimizeMolecule`, multi: `MMFFOptimizeMoleculeConfs`).
- CLI gains `--optimize` flag (`store_true`, default `False`).

#### cmd.py

- `optimize` keyword added to `rdkconf_desc` (BoolArg, default False).
- `minimize` keyword added to `rdkconf_desc` (BoolArg, default False).
- When `optimize=True`: pass `--optimize` to subprocess.
- After model creation, if `minimize=True`:
  - Run `minimize #{model_id} maxSteps 1000` per model.
  - On older ChimeraX (no `minimize` command): log warning, keep model as-is.

## Error Handling

| Scenario | Behavior |
|----------|----------|
| `optimize true`, MMFF fails | Warning (current behavior, unchanged) |
| `minimize true`, ChimeraX < 1.11 | Warning: "minimize requires ChimeraX 1.11+; skipping" |
| `minimize true`, minimize command fails | Warning + keep model as-is (best-effort) |
| `minimize true conformers 5` | Minimize each model independently |

## Version Detection for minimize

Use try/except around `run(session, "minimize ...")`:

```python
try:
    run(session, f"minimize #{model.id_string} maxSteps 1000")
except Exception:
    session.logger.warning(
        "minimize command not available (requires ChimeraX 1.11+); skipping"
    )
```

This is simpler and more robust than version-checking APIs.

## Design Decisions

- **Both default false**: ETKDGv3 coordinates are adequate for visualization. Optimization is opt-in.
- **Independent flags, not enum**: The two optimizations are different processes (subprocess vs ChimeraX) and can be combined.
- **Best-effort minimize**: If ChimeraX version doesn't support `minimize`, warn and continue rather than error.
- **maxSteps 1000**: Reasonable default. Not exposed as a user parameter (YAGNI).
- **Breaking change accepted**: MMFF defaulting to off is intentional. ETKDGv3 alone produces reasonable 3D coordinates.

## ChimeraX minimize Command (Reference)

Available since ChimeraX 1.11 (December 2025). Uses OpenMM with AMBER14 force field. Handles small molecules via GAFF/Antechamber automatically when `dockPrep true` (default).

```
minimize [atomic-model-spec] [dockPrep true|false] [maxSteps N]
```
