# Phase 1: ChimeraX `rdksmiles` Bundle Implementation

## Goal

Build a ChimeraX bundle (`ChimeraX-RDKitSMILES`) that registers an `rdksmiles` command. The command accepts a SMILES string, shells out to `uv run --script smiles_to_sdf.py` (using RDKit in an isolated environment), and opens the resulting SDF file in ChimeraX.

## Architecture

**Key decision**: Subprocess bridge. RDKit is unavailable inside ChimeraX's Python, so we delegate 3D coordinate generation to an external PEP 723 script run by `uv`.

**End-to-end flow**:
1. User types `rdksmiles "c1ccccc1"` in ChimeraX
2. ChimeraX calls `cmd.rdksmiles(session, "c1ccccc1")`
3. `_find_uv()` locates the `uv` binary (with PATH fallback for macOS Finder launch)
4. `_get_script_path()` locates the bundled `smiles_to_sdf.py` via `importlib.resources`
5. `subprocess.run(["uv", "run", "--script", script_path, smiles])` executes
6. Script prints SDF file path to stdout
7. ChimeraX opens the SDF file, then deletes the temp file

## Project Structure

```
chimerax-rdkit-smiles/
├── pyproject.toml              # Bundle metadata and build config
├── src/
│   ├── __init__.py             # BundleAPI registration
│   ├── cmd.py                  # Command implementation (rdksmiles)
│   └── smiles_to_sdf.py       # PEP 723 script (bundled as package data)
├── scripts/
│   └── smoke.cxc               # Smoke test for echidna run
├── tests/
│   └── test_cmd.py             # Unit tests
└── README.md
```

## Tasks

### 1. Scaffold with echidna
- [ ] Run `echidna init --name rdksmiles --bundle-name ChimeraX-RDKitSMILES .`
- [ ] Verify generated skeleton

### 2. Create `pyproject.toml`
- [ ] Set `package = "chimerax.rdkitsmiles"`, `pure = true`
- [ ] Declare `[chimerax.command.rdksmiles]` with category "Structure Generation"
- [ ] Include `smiles_to_sdf.py` via `[chimerax.package-data]`
- [ ] Add ruff config (target Python 3.11 for ChimeraX compatibility)

### 3. Create `src/__init__.py` (BundleAPI)
- [ ] Define `_RdksmilesAPI(BundleAPI)` with `api_version = 1`
- [ ] Implement `register_command()` that registers the `rdksmiles` command
- [ ] Set `bundle_api = _RdksmilesAPI()`

### 4. Create `src/cmd.py` (Core implementation)
- [ ] `_find_uv()` -- locate uv binary with PATH + common fallback locations
- [ ] `_get_script_path()` -- use `importlib.resources.files()` to find bundled script
- [ ] `rdksmiles()` function -- subprocess bridge with 60s timeout, error handling
- [ ] `rdksmiles_desc` -- `CmdDesc` with required `smiles` (StringArg) and keyword args: `output_file` (SaveFileNameArg), `name` (StringArg), `optimize` (BoolArg)
- [ ] All errors raised as `chimerax.core.errors.UserError`
- [ ] Clean up temp SDF file after opening

### 5. Copy `smiles_to_sdf.py` into `src/`
- [ ] Copy existing script as package data (no modifications needed for v1)

### 6. Create `scripts/smoke.cxc`
- [ ] Test `rdksmiles "CCO"` (ethanol)
- [ ] Verify model opens

### 7. Create `tests/test_cmd.py`
- [ ] Test `_find_uv()` with mocked `shutil.which`
- [ ] Test `_find_uv()` raises when uv not found
- [ ] Integration test: valid SMILES produces SDF via subprocess
- [ ] Integration test: invalid SMILES returns error

### 8. Build and test
- [ ] `echidna run --script scripts/smoke.cxc` -- build, install, launch with smoke test
- [ ] `echidna test` -- run unit tests
- [ ] Verify `rdksmiles "c1ccccc1"` opens benzene in ChimeraX

## Error Handling

| Failure Mode | Detection | User Message |
|---|---|---|
| `uv` not installed | `_find_uv()` fallback exhausted | "Cannot find 'uv' command. Install from..." |
| Invalid SMILES | Script exits code 1, stderr | "Failed to generate 3D structure: Invalid SMILES: ..." |
| 3D embedding fails | Script exits code 1, stderr | "Failed to generate 3D structure: Failed to generate 3D coordinates..." |
| Script timeout (60s) | `TimeoutExpired` caught | "Timed out generating 3D structure for: ..." |
| SDF file missing | Path check post-subprocess | "Script did not produce a valid SDF file." |

## Gotchas

1. **PATH isolation**: ChimeraX launched from macOS Finder has minimal PATH. `_find_uv()` must check `~/.local/bin/uv`, `~/.cargo/bin/uv`, `~/.nix-profile/bin/uv` as fallbacks.
2. **First-run latency**: First `uv run --script` downloads and caches RDKit (10-30s). Subsequent calls are fast (<2s).
3. **`importlib.resources` vs `__file__`**: Using `importlib.resources.files()` is more robust for locating package data than `__file__`.
4. **Package-data syntax**: `[chimerax.package-data]` uses ChimeraX bundle builder's own mechanism. If it fails, fall back to `[tool.setuptools.package-data]`.
5. **`optimize` keyword**: Declared in CmdDesc but not wired to script in v1 (MMFF always runs). Low-risk v2 enhancement.
6. **Temp file cleanup**: ChimeraX reads SDF into memory during `open`, so deleting afterward is safe.

## Prerequisites

- UCSF ChimeraX installed
- `echidna` CLI installed
- `uv` installed (already at `~/.nix-profile/bin/uv`)

---
- [ ] **DONE** - Phase complete
