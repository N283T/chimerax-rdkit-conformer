# Direct Model Creation Design

Date: 2026-02-24

## Goal

Eliminate the SDF file intermediary. Build ChimeraX `AtomicStructure` directly from RDKit coordinate data passed as JSON over stdout.

## Current Pipeline

```
subprocess (RDKit) → temp SDF file → ChimeraX `open` command → delete temp file
```

## New Pipeline

```
subprocess (RDKit) → JSON stdout → cmd.py parses JSON → AtomicStructure in memory
```

## Benefits

- No temp file I/O (faster, no cleanup needed)
- No file path handling or injection surface
- More control over atom names, residue names, bond orders
- Simpler error handling (no file existence checks)

## JSON Format

```json
{
  "atoms": [
    {"element": "C", "x": 0.0, "y": 0.0, "z": 0.0},
    {"element": "O", "x": 1.54, "y": 1.2, "z": 0.0}
  ],
  "bonds": [
    {"begin": 0, "end": 1, "order": 1}
  ]
}
```

- Atoms indexed by array position
- Bond orders included (nice to have for display)
- No formal charges (YAGNI)

## UX

```
rdkconf CCO                                   # Direct model (no file)
rdkconf CCO name EtOH                         # Custom residue name
rdkconf CCO hydrogen false                    # Hide hydrogens
rdkconf "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3" # InChI auto-detect
rdkconf GGG format sequence                   # Peptide sequence
```

The `output` keyword is removed. Users can export from ChimeraX if needed.

## Changes

### 1. `input_to_sdf.py` -> `input_to_json.py` (rename)

- Default output: JSON to stdout
- Remove SDF writing, `--output` flag, tempfile logic
- New function: `mol_to_json(mol) -> dict` serializes RDKit Mol to JSON dict
- `main()` prints `json.dumps(mol_to_json(mol))`

### 2. `cmd.py`

- Remove: `output` keyword, `_validate_sdf_path()`, SDF open/cleanup
- Add: `_build_model(session, mol_data, name)` using AtomicStructure API
- Add: JSON parsing of subprocess stdout
- Bond orders: set if ChimeraX API supports it

### 3. `pyproject.toml`

- Update `package-data` script name
- Bump version

### 4. `__init__.py`

- No changes (command registration unchanged)

### 5. `README.md`

- Remove `output` examples
- Update "How It Works" section

## What's Removed

- `output` keyword from `rdkconf` command
- Temp file creation/cleanup in both script and cmd.py
- `_validate_sdf_path()` helper
- SDF writing in subprocess script
- `SaveFileNameArg` import

## Error Handling

- JSON parse failure -> `UserError` with raw stderr
- Missing atoms/bonds keys -> `UserError`
- Subprocess timeout/error -> same as current

## Decisions

- JSON over stdout (not pickle, not binary) for safety and debuggability
- Bond orders nice-to-have, not essential
- No formal charges in JSON (YAGNI)
- `output` keyword dropped entirely (export from ChimeraX instead)
- Script renamed to `input_to_json.py` to reflect new output format
