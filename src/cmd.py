"""ChimeraX command implementation for rdkconf.

Generates 3D conformers from molecular notations by running RDKit
in a uv-managed subprocess and building AtomicStructure models
directly from the JSON output.
"""

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

    When user_format is provided, validates it against _VALID_FORMATS.
    Otherwise, auto-detects InChI (starts with 'InChI=') and defaults to SMILES.

    Raises
    ------
    UserError
        If user_format is not a recognized format.
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
        Model and residue name. The residue is placed in chain ' '
        with sequence number 1.

    Returns
    -------
    AtomicStructure

    Raises
    ------
    UserError
        If mol_data contains no atoms, is missing 'bonds' data,
        or contains malformed atom/bond entries.

    Notes
    -----
    Registers a custom 'order' attribute on Bond if not already registered.
    """
    if not mol_data.get("atoms"):
        raise UserError("RDKit output contains no atoms")
    if "bonds" not in mol_data:
        raise UserError("RDKit output is missing bonds data")

    try:
        Bond.register_attr(session, "order", "rdkit_conformer", attr_type=float)
    except ValueError as e:
        if "already" not in str(e).lower():
            raise UserError(f"Failed to register bond order attribute: {e}") from e

    s = AtomicStructure(session, name=name)
    try:
        r = s.new_residue(name, " ", 1)

        element_count = {}
        atoms = []
        for i, atom_info in enumerate(mol_data["atoms"]):
            try:
                elem = atom_info["element"]
                xyz = array(
                    [atom_info["x"], atom_info["y"], atom_info["z"]], dtype=float64
                )
            except (KeyError, TypeError) as e:
                raise UserError(f"Malformed atom data at index {i}: {e}") from e
            n = element_count.get(elem, 0) + 1
            element_count[elem] = n
            atom_name = f"{elem}{n}"
            a = add_atom(atom_name, elem, r, xyz)
            atoms.append(a)

        num_atoms = len(atoms)
        for i, bond_info in enumerate(mol_data["bonds"]):
            try:
                begin_idx = bond_info["begin"]
                end_idx = bond_info["end"]
            except KeyError as e:
                raise UserError(f"Malformed bond data at index {i}: {e}") from e
            if (
                begin_idx < 0
                or begin_idx >= num_atoms
                or end_idx < 0
                or end_idx >= num_atoms
            ):
                raise UserError(
                    f"Invalid bond indices ({begin_idx}, {end_idx}) "
                    f"for {num_atoms} atoms"
                )
            b = add_bond(atoms[begin_idx], atoms[end_idx])
            b.order = bond_info.get("order", 1.0)
    except Exception:
        s.delete()
        raise

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
        Model and residue name (default: UNL).
    hydrogen : bool
        Show hydrogens (default: True).

    Raises
    ------
    UserError
        If uv is not found, input format is invalid, name is invalid,
        the RDKit subprocess fails or times out, or the output is unparseable.
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

    if result.stderr.strip():
        session.logger.warning(f"RDKit warning: {result.stderr.strip()}")

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
