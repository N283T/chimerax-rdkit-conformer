#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "rdkit",
# ]
# ///
"""Convert molecular notation to 3D SDF using RDKit ETKDGv3."""

import argparse
import sys
import tempfile
from pathlib import Path

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


def input_to_3d(
    input_str: str, fmt: str = "smiles", output: Path | None = None
) -> Path:
    """Convert molecular notation to a 3D SDF file using ETKDGv3.

    Args:
        input_str: Molecular notation string.
        fmt: Input format (default: smiles).
        output: Output file path. If None, creates a temp file.

    Returns:
        Path to the generated SDF file.
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

    if output is None:
        tmp = tempfile.NamedTemporaryFile(suffix=".sdf", delete=False)
        output = Path(tmp.name)
        tmp.close()

    with Chem.SDWriter(str(output)) as writer:
        writer.write(mol)

    return output


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert molecular notation to 3D SDF")
    parser.add_argument("input", help="Molecular notation string")
    parser.add_argument(
        "-f",
        "--format",
        default="smiles",
        choices=VALID_FORMATS,
        help="Input format (default: smiles)",
    )
    parser.add_argument(
        "-o", "--output", type=Path, default=None, help="Output SDF file path"
    )
    args = parser.parse_args()

    try:
        path = input_to_3d(args.input, args.format, args.output)
    except (ValueError, RuntimeError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error: {type(e).__name__}: {e}", file=sys.stderr)
        sys.exit(2)

    print(str(path))


if __name__ == "__main__":
    main()
