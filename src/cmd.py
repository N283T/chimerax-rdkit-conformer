"""Command implementation for rdksmiles."""

import os
import shutil
import subprocess
from pathlib import Path

from chimerax.core.commands import CmdDesc, StringArg, BoolArg, SaveFileNameArg, run
from chimerax.core.errors import UserError


def _find_script() -> Path:
    """Locate the bundled smiles_to_sdf.py script."""
    return Path(__file__).parent.joinpath("smiles_to_sdf.py")


def _find_uv() -> str:
    """Find the uv executable."""
    uv_path = shutil.which("uv")
    if uv_path is None:
        raise UserError("uv not found in PATH. Install uv: https://docs.astral.sh/uv/")
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

    cmd_args = [uv_path, "run", "--script", str(script_path), smiles]
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

    sdf_path = result.stdout.strip()
    if not Path(sdf_path).exists():
        raise UserError(f"SDF file not generated: {sdf_path}")

    open_cmd = f'open "{sdf_path}" name "{name}"'
    run(session, open_cmd)

    if not hydrogen:
        run(session, "hide H")

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
