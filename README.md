# ChimeraX-RDKitConformer

Generate 3D conformers from molecular notations using RDKit ETKDGv3, directly in ChimeraX.

Two installable bundles are provided. Both register the same `rdkconf` command — install one OR the other.

| Bundle | How RDKit runs | Pros | Cons |
|--------|---------------|------|------|
| **bundle-direct** | Imported directly into ChimeraX Python | Simple, fast, no extra tools | Adds rdkit to ChimeraX's Python |
| **bundle-uv** | Isolated subprocess via [uv](https://docs.astral.sh/uv/) + [PEP 723](https://peps.python.org/pep-0723/) | Zero environment pollution | Requires uv, slightly slower |

## Requirements

- [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) 1.6+ (1.11+ for `minimize` keyword)

## Installation

### Option A: Direct RDKit (recommended for simplicity)

Install RDKit into ChimeraX's Python, then install the bundle.

**From ChimeraX command line (GUI):**

```
pip install rdkit
cd /path/to/chimerax-rdkit-conformer/bundle-direct
devel install .
```

**From terminal:**

```bash
ChimeraX --nogui --exit --cmd 'pip install rdkit'
cd /path/to/chimerax-rdkit-conformer/bundle-direct
ChimeraX --nogui --exit --cmd 'devel install .'
```

> **Note:** `devel install` must be run from the bundle directory (`cd` then `devel install .`). Passing an absolute path does not work due to a ChimeraX bundle builder limitation.

### Option B: uv subprocess (zero ChimeraX pollution)

Install uv into ChimeraX, then install the bundle.

**From ChimeraX command line (GUI):**

```
pip install uv
cd /path/to/chimerax-rdkit-conformer/bundle-uv
devel install .
```

**From terminal:**

```bash
ChimeraX --nogui --exit --cmd 'pip install uv'
cd /path/to/chimerax-rdkit-conformer/bundle-uv
ChimeraX --nogui --exit --cmd 'devel install .'
```

With this option, RDKit runs in an isolated subprocess managed by uv — ChimeraX's own packages are not affected.

If uv is installed outside ChimeraX, you can point the bundle to it:

```
rdkconf uvPath ~/.local/bin/uv
```

## Usage

```
rdkconf CCO                                        # SMILES (default)
rdkconf c1ccccc1                                   # Benzene
rdkconf "CC(=O)Oc1ccccc1C(=O)O"                   # Aspirin
rdkconf "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"      # InChI (auto-detected)
rdkconf GGG format sequence                        # Peptide sequence
rdkconf ACGT format dna                            # DNA sequence
rdkconf ACGU format rna                            # RNA sequence
rdkconf ">seq\nGGG" format fasta                   # FASTA format
rdkconf "PEPTIDE1{G.G.G}$$$$" format helm          # HELM notation
rdkconf CCO name EtOH                              # Custom residue name
rdkconf CCO hydrogen false                         # Hide hydrogens
rdkconf CCO optimize true                          # RDKit MMFF optimization
rdkconf CCO minimize true                          # ChimeraX minimization (1.11+)
rdkconf CCO optimize true minimize true            # Both: MMFF then minimize
rdkconf c1ccccc1 conformers 5                     # Multiple conformers
rdkconf CCO conformers 10 name EtOH               # Named multi-conformers
```

### Supported Formats

| Format | `format` keyword | Auto-detected |
|--------|-----------------|---------------|
| SMILES | `smiles` (default) | Yes (fallback) |
| InChI | `inchi` | Yes (`InChI=` prefix) |
| Protein sequence | `sequence` | No |
| DNA sequence | `dna` | No |
| RNA sequence | `rna` | No |
| FASTA | `fasta` | No |
| HELM | `helm` | No |

## How It Works

ChimeraX's built-in SMILES support (`open smiles:CCO`) depends on NCI's web service
(requires internet, no control over 3D generation quality). This bundle uses RDKit
ETKDGv3 locally for high-quality 3D coordinate generation without network dependency.
The 3D structure is built directly in memory using ChimeraX's AtomicStructure API — no
intermediate files are written.

## Development

```bash
# Build and install (must cd into the bundle directory first)
cd bundle-direct
ChimeraX --nogui --exit --cmd 'devel install .'

# Run tests (either bundle)
cd bundle-direct && uv run --no-project --with pytest --with rdkit pytest tests/ -v
cd bundle-uv && uv run --no-project --with pytest --with rdkit pytest tests/ -v
```

## License

MIT
