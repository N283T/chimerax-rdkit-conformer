"""Command implementation for rdkconf."""

import json
import re
import shutil
import subprocess
from pathlib import Path

from chimerax.atomic import AtomicStructure, Bond
from chimerax.atomic.struct_edit import add_atom, add_bond
from chimerax.core.commands import CmdDesc, StringArg, BoolArg, run
from chimerax.core.errors import UserError
from numpy import array, float64


def _find_script() -> Path:
    """Locate the bundled input_to_json.py script."""
    return Path(__file__).parent.joinpath("input_to_json.py")


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


def _build_model(session, mol_data, name="UNL"):
    """Build an AtomicStructure from parsed JSON molecule data.

    Parameters
    ----------
    session : chimerax.core.session.Session
    mol_data : dict
        JSON-parsed dict with 'atoms' and 'bonds' keys.
    name : str
        Model and residue name.

    Returns
    -------
    AtomicStructure
    """
    if not mol_data.get("atoms"):
        raise UserError("RDKit output contains no atoms")
    if "bonds" not in mol_data:
        raise UserError("RDKit output is missing bonds data")

    try:
        Bond.register_attr(session, "order", "rdkit_conformer", attr_type=float)
    except ValueError:
        pass

    s = AtomicStructure(session, name=name)
    r = s.new_residue(name, " ", 1)

    element_count = {}
    atoms = []
    for atom_info in mol_data["atoms"]:
        elem = atom_info["element"]
        n = element_count.get(elem, 0) + 1
        element_count[elem] = n
        atom_name = f"{elem}{n}"
        xyz = array([atom_info["x"], atom_info["y"], atom_info["z"]], dtype=float64)
        a = add_atom(atom_name, elem, r, xyz)
        atoms.append(a)

    for bond_info in mol_data["bonds"]:
        b = add_bond(atoms[bond_info["begin"]], atoms[bond_info["end"]])
        b.order = bond_info.get("order", 1.0)

    return s


def rdkconf(session, input_str, format=None, name="UNL", hydrogen=True):
    """Generate 3D conformer from molecular notation using RDKit ETKDGv3.

    Parameters
    ----------
    session : chimerax.core.session.Session
    input_str : str
        Molecular notation string (SMILES by default).
    format : str, optional
        Input format. Auto-detects InChI. Default: smiles.
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
        "--no-project",
        "--script",
        str(script_path),
        input_str,
        "--format",
        fmt,
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

    try:
        mol_data = json.loads(result.stdout)
    except json.JSONDecodeError as e:
        raise UserError(f"Failed to parse RDKit output: {e}")

    model = _build_model(session, mol_data, name)
    session.models.add([model])

    if not hydrogen:
        model_spec = f"#{model.id_string}"
        run(session, f"hide {model_spec} & H")

    session.logger.info(f"Generated 3D conformer ({fmt}) from: {input_str}")


rdkconf_desc = CmdDesc(
    required=[("input_str", StringArg)],
    keyword=[
        ("format", StringArg),
        ("name", StringArg),
        ("hydrogen", BoolArg),
    ],
    synopsis="Generate 3D conformer from molecular notation using RDKit ETKDGv3",
)
