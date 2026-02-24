# Multi-Conformer Generation Design

## Goal

Allow generating multiple 3D conformers from a single molecular input, each as a separate ChimeraX model.

## Command Interface

```
rdkconf c1ccccc1                          # 1 conformer (default)
rdkconf c1ccccc1 conformers 10            # 10 conformers as separate models
rdkconf c1ccccc1 conformers 10 name LIG   # each model named LIG_1, LIG_2, ...
```

New keyword: `conformers` (IntArg, default=1, max=50).

## Architecture

### Data Flow

```
rdkconf CCO conformers 5
  |
  v  subprocess
input_to_json.py "CCO" --format smiles --conformers 5
  |
  v  RDKit EmbedMultipleConfs(numConfs=5, pruneRmsThresh=0.5)
  v  MMFFOptimizeMoleculeConfs
  |
  v  stdout: JSON array
[{"atoms":[...], "bonds":[...]}, {"atoms":[...], "bonds":[...]}, ...]
  |
  v  cmd.py: loop over array, _build_model() per conformer
5 separate AtomicStructure models added to session
```

### JSON Output Format

Always a JSON array, even for single conformers:

```json
[{"atoms": [...], "bonds": [...]}]
```

This simplifies parsing in cmd.py (no conditional format handling). No external consumers of the JSON format exist.

## Changes

### input_to_json.py

- `input_to_json()` gains `num_confs` parameter (default=1).
- `num_confs=1`: Uses `EmbedMolecule` (current single-conformer behavior).
- `num_confs>1`: Uses `EmbedMultipleConfs` + `MMFFOptimizeMoleculeConfs`.
- `mol_to_json()` unchanged — called once per conformer.
- `pruneRmsThresh=0.5` hardcoded (not exposed to user). YAGNI.
- Output is always a JSON array of conformer dicts.
- CLI gains `--conformers N` argument (default=1).
- Actual conformer count may be less than requested due to RMS pruning. This is normal, not an error.

### cmd.py

- `conformers` keyword added to `rdkconf_desc` (IntArg).
- Validated: 1 <= conformers <= 50. Out-of-range raises UserError.
- Passed as `--conformers N` to subprocess.
- JSON array parsed, `_build_model()` called per entry.
- Model naming: single conformer uses `name` as-is. Multiple conformers use `{name}_1`, `{name}_2`, etc.
- Log message: "Generated N conformers (fmt) from: input".

## Error Handling

- `conformers 0` or `conformers 51` -> UserError with range guidance.
- RMS pruning returns fewer than requested -> info log ("Generated 3 of 5 requested conformers"), not an error.
- Empty result array (all conformers pruned) -> UserError.

## Design Decisions

- **Separate models, not coordsets**: Each conformer is an independent AtomicStructure. Easier to individually manipulate, hide/show, color differently.
- **Max 50 conformers**: Prevents accidental resource exhaustion. Most docking workflows use 10-50.
- **pruneRmsThresh not exposed**: Internal default of 0.5 A. Keep command simple.
- **Always JSON array**: Simplifies cmd.py parsing. Breaking change to JSON format but no external consumers.
