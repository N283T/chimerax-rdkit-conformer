# Phase 1: ChimeraX RDKit SMILES Bundle

## Context

ChimeraX has built-in SMILES support (`open smiles:CCO`) but it depends on NCI's web service (online only, no control over 3D generation quality). We already have a working PEP 723 script (`smiles_to_sdf.py`) that uses RDKit ETKDGv3 for high-quality 3D coordinate generation locally.

**Goal**: Package this into a ChimeraX bundle so users can type `rdksmiles CCO` in ChimeraX's command line.

## Architecture

```
rdksmiles CCO
    │
    ▼ (ChimeraX command handler)
subprocess: uv run --script smiles_to_sdf.py "CCO"
    │
    ▼ (stdout: /tmp/xxx.sdf)
run(session, "open /tmp/xxx.sdf")
    │
    ▼ (cleanup temp file)
molecule displayed in ChimeraX
```

**Why subprocess?** RDKit is NOT available in ChimeraX's isolated Python. Using `uv run --script` cleanly manages the RDKit dependency outside ChimeraX.

## Project Structure

```
chimerax-rdkit-smiles/
├── pyproject.toml
├── src/
│   ├── __init__.py          # BundleAPI registration
│   ├── cmd.py               # rdksmiles command implementation
│   └── smiles_to_sdf.py     # PEP 723 script (bundled as package data)
├── scripts/
│   └── smoke.cxc            # Smoke test
└── README.md
```

## Step-by-Step Implementation

### Step 1: Install echidna

```bash
ghq get N283T/echidna
cd "$(ghq root)/github.com/N283T/echidna"
cargo install --path .
```

### Step 2: Scaffold with echidna

```bash
cd /Users/nagaet/chimerax-smiles
echidna init --name rdkit-smiles .
```

This generates the command template. We then customize the generated files.

### Step 3: pyproject.toml

```toml
[build-system]
requires = ["ChimeraX-BundleBuilder"]
build-backend = "chimerax.bundle_builder.cx_pep517"

[project]
name = "ChimeraX-RDKitSMILES"
version = "0.1.0"
description = "Generate 3D molecules from SMILES using RDKit ETKDGv3"
license = { text = "MIT" }
dynamic = ["classifiers", "requires-python"]

[chimerax]
package = "chimerax.rdkit_smiles"
pure = true
min-session-version = 1
max-session-version = 1
categories = ["Structure Generation"]
classifiers = []

[chimerax.command.rdksmiles]
category = "Structure Generation"
description = "Generate 3D structure from SMILES using RDKit ETKDGv3"

[chimerax.package-data]
"src/" = ["smiles_to_sdf.py"]
```

### Step 4: src/__init__.py

```python
"""ChimeraX-RDKitSMILES - Generate 3D molecules from SMILES using RDKit ETKDGv3"""

from chimerax.core.toolshed import BundleAPI


class _RdksmilesAPI(BundleAPI):
    api_version = 1

    @staticmethod
    def register_command(bundle_info, command_info, logger):
        from . import cmd
        if command_info.name == "rdksmiles":
            from chimerax.core.commands import register
            register(command_info.name, cmd.rdksmiles_desc, cmd.rdksmiles)

bundle_api = _RdksmilesAPI()
```

### Step 5: src/cmd.py

```python
"""Command implementation for rdksmiles."""

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

from chimerax.core.commands import CmdDesc, StringArg, BoolArg, SaveFileNameArg, run
from chimerax.core.errors import UserError


def _find_script() -> Path:
    """Locate the bundled smiles_to_sdf.py script."""
    return Path(__file__).parent / "smiles_to_sdf.py"


def _find_uv() -> str:
    """Find the uv executable."""
    uv_path = shutil.which("uv")
    if uv_path is None:
        raise UserError(
            "uv not found in PATH. Install uv: https://docs.astral.sh/uv/"
        )
    return uv_path


def rdksmiles(session, smiles, output=None, name="UNL", hydrogen=True):
    """Generate 3D structure from SMILES using RDKit ETKDGv3.

    Parameters
    ----------
    session : chimerax.core.session.Session
    smiles : str
        SMILES string.
    output : str, optional
        Save SDF to this path.
    name : str
        Residue name (default: UNL).
    hydrogen : bool
        Show hydrogens (default: True).
    """
    script_path = _find_script()
    uv_path = _find_uv()

    # Build command
    cmd_args = [uv_path, "run", "--script", str(script_path), smiles]
    if output:
        cmd_args.extend(["-o", output])

    # Run external RDKit process
    try:
        result = subprocess.run(
            cmd_args,
            capture_output=True,
            text=True,
            timeout=60,
            env={**os.environ, "UV_CACHE_DIR": os.environ.get("UV_CACHE_DIR", "")},
        )
    except subprocess.TimeoutExpired:
        raise UserError("RDKit 3D generation timed out (60s)")

    if result.returncode != 0:
        raise UserError(f"RDKit error: {result.stderr.strip()}")

    sdf_path = result.stdout.strip()
    if not Path(sdf_path).exists():
        raise UserError(f"SDF file not generated: {sdf_path}")

    # Open in ChimeraX
    open_cmd = f'open "{sdf_path}" name "{name}"'
    run(session, open_cmd)

    if not hydrogen:
        run(session, "hide H")

    # Clean up temp file (only if no explicit output)
    if output is None:
        try:
            os.unlink(sdf_path)
        except OSError:
            pass

    session.logger.info(f"Generated 3D structure from SMILES: {smiles}")


rdksmiles_desc = CmdDesc(
    required=[("smiles", StringArg)],
    keyword=[
        ("output", SaveFileNameArg),
        ("name", StringArg),
        ("hydrogen", BoolArg),
    ],
    synopsis="Generate 3D structure from SMILES using RDKit ETKDGv3",
)
```

### Step 6: src/smiles_to_sdf.py

Move existing `/Users/nagaet/chimerax-smiles/smiles_to_sdf.py` into `src/`.

### Step 7: scripts/smoke.cxc

```
# Smoke test for ChimeraX-RDKitSMILES
rdksmiles c1ccccc1
log show
```

### Step 8: Build and test

```bash
echidna run   # Build, install, launch ChimeraX
# In ChimeraX: rdksmiles CCO
# In ChimeraX: rdksmiles "CC(=O)Oc1ccccc1C(=O)O"
```

## Verification

1. `echidna run` succeeds (build + install)
2. In ChimeraX: `rdksmiles c1ccccc1` opens benzene
3. In ChimeraX: `rdksmiles "CC(=O)Oc1ccccc1C(=O)O"` opens aspirin
4. In ChimeraX: `rdksmiles CCO name EtOH` opens ethanol with custom name
5. Error handling: `rdksmiles INVALID_SMILES` shows a user-friendly error

---
- [ ] **DONE** - Phase complete
