#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "rdkit",
# ]
# ///
"""Convert SMILES to 3D SDF using RDKit ETKDGv3."""

import argparse
import sys
import tempfile
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem


def smiles_to_3d(smiles: str, output: Path | None = None) -> Path:
    """Convert a SMILES string to a 3D SDF file using ETKDGv3.

    Args:
        smiles: SMILES string to convert.
        output: Output file path. If None, creates a temp file.

    Returns:
        Path to the generated SDF file.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        raise RuntimeError(f"Failed to generate 3D coordinates for: {smiles}")

    AllChem.MMFFOptimizeMolecule(mol)

    if output is None:
        tmp = tempfile.NamedTemporaryFile(suffix=".sdf", delete=False)
        output = Path(tmp.name)
        tmp.close()

    with Chem.SDWriter(str(output)) as writer:
        writer.write(mol)

    return output


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert SMILES to 3D SDF")
    parser.add_argument("smiles", help="SMILES string")
    parser.add_argument(
        "-o", "--output", type=Path, default=None, help="Output SDF file path"
    )
    args = parser.parse_args()

    try:
        path = smiles_to_3d(args.smiles, args.output)
    except (ValueError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    print(str(path))


if __name__ == "__main__":
    main()
