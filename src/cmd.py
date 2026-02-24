"""Command implementation for rdkconf."""

import os
import re
import shutil
import subprocess
from pathlib import Path

from chimerax.core.commands import CmdDesc, StringArg, BoolArg, SaveFileNameArg, run
from chimerax.core.errors import UserError


def _find_script() -> Path:
    """Locate the bundled input_to_sdf.py script."""
    return Path(__file__).parent.joinpath("input_to_sdf.py")


def _find_uv() -> str:
    """Find the uv executable."""
    uv_path = shutil.which("uv")
    if uv_path is None:
        raise UserError("uv not found in PATH. Install uv: https://docs.astral.sh/uv/")
    return uv_path


def _validate_name(name: str) -> str:
    """Validate residue name is safe for ChimeraX command strings."""
    if not re.match(r"^[A-Za-z0-9_]+$", name):
        raise UserError(
            f"Invalid residue name: {name!r} (alphanumeric and underscore only)"
        )
    return name


def _validate_sdf_path(raw_path: str) -> Path:
    """Validate subprocess output looks like a valid SDF path."""
    path = Path(raw_path)
    if not path.is_absolute() or path.suffix != ".sdf":
        raise UserError(f"Unexpected output from RDKit script: {raw_path!r}")
    if not path.exists():
        raise UserError(f"SDF file not generated: {raw_path}")
    return path


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
        uv_path,
        "run",
        "--script",
        str(script_path),
        input_str,
        "--format",
        fmt,
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
