#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "rdkit",
# ]
# ///
"""Convert molecular notation to 3D JSON using RDKit ETKDGv3."""

import argparse
import json
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

VALID_FORMATS = ("smiles", "inchi", "fasta", "sequence", "helm", "dna", "rna")


def parse_input(input_str: str, fmt: str = "smiles") -> Chem.Mol:
    """Parse molecular notation string into an RDKit Mol object.

    Args:
        input_str: Molecular notation string.
        fmt: Input format. One of: smiles, inchi, fasta, sequence, helm, dna, rna.

    Returns:
        RDKit Mol object.

    Raises:
        ValueError: If format is unsupported or input is invalid.
    """
    parsers = {
        "smiles": lambda s: Chem.MolFromSmiles(s),
        "inchi": lambda s: Chem.MolFromInchi(s),
        "fasta": lambda s: Chem.MolFromFASTA(s),
        "sequence": lambda s: Chem.MolFromSequence(s),
        "helm": lambda s: Chem.MolFromHELM(s),
        "dna": lambda s: Chem.MolFromSequence(s, flavor=6),
        "rna": lambda s: Chem.MolFromSequence(s, flavor=2),
    }
    parser = parsers.get(fmt)
    if parser is None:
        raise ValueError(f"Unsupported format: {fmt}")
    mol = parser(input_str)
    if mol is None:
        raise ValueError(f"Invalid {fmt} input: {input_str}")
    return mol


def mol_to_json(mol: Chem.Mol) -> dict:
    """Serialize an RDKit Mol with 3D coordinates to a JSON-compatible dict.

    Args:
        mol: RDKit Mol object with embedded 3D coordinates and explicit Hs.

    Returns:
        Dict with 'atoms' (list of element/x/y/z) and 'bonds' (list of begin/end/order).
    """
    conf = mol.GetConformer()
    atoms = []
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        atoms.append(
            {
                "element": mol.GetAtomWithIdx(i).GetSymbol(),
                "x": pos.x,
                "y": pos.y,
                "z": pos.z,
            }
        )

    bond_type_map = {
        Chem.BondType.SINGLE: 1.0,
        Chem.BondType.DOUBLE: 2.0,
        Chem.BondType.TRIPLE: 3.0,
        Chem.BondType.AROMATIC: 1.5,
    }
    bonds = []
    for bond in mol.GetBonds():
        bonds.append(
            {
                "begin": bond.GetBeginAtomIdx(),
                "end": bond.GetEndAtomIdx(),
                "order": bond_type_map.get(bond.GetBondType(), 1.0),
            }
        )

    return {"atoms": atoms, "bonds": bonds}


def input_to_json(input_str: str, fmt: str = "smiles") -> dict:
    """Convert molecular notation to a 3D JSON dict using ETKDGv3.

    Args:
        input_str: Molecular notation string.
        fmt: Input format (default: smiles).

    Returns:
        Dict with 'atoms' and 'bonds' keys.

    Raises:
        ValueError: If format is unsupported or input is invalid.
        RuntimeError: If 3D coordinate generation fails.
    """
    mol = parse_input(input_str, fmt)
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        raise RuntimeError(f"Failed to generate 3D coordinates for: {input_str}")

    opt_result = AllChem.MMFFOptimizeMolecule(mol)
    if opt_result == -1:
        import warnings

        warnings.warn(
            f"MMFF force field setup failed for: {input_str}. "
            "3D coordinates are unoptimized.",
            stacklevel=2,
        )

    return mol_to_json(mol)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert molecular notation to 3D JSON"
    )
    parser.add_argument("input", help="Molecular notation string")
    parser.add_argument(
        "-f",
        "--format",
        default="smiles",
        choices=VALID_FORMATS,
        help="Input format (default: smiles)",
    )
    args = parser.parse_args()

    try:
        data = input_to_json(args.input, args.format)
    except (ValueError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {type(e).__name__}: {e}", file=sys.stderr)
        sys.exit(2)

    print(json.dumps(data))


if __name__ == "__main__":
    main()
