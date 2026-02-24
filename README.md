# ChimeraX-RDKitConformer

Generate 3D conformers from molecular notations using RDKit ETKDGv3, directly in ChimeraX.

Two installable bundles are provided. Both register the same `rdkconf` command — install one OR the other.

| Option | Bundle | How RDKit runs | ChimeraX pollution |
|--------|--------|---------------|-------------------|
| **A** | bundle-direct | Imported directly into ChimeraX Python | rdkit added |
| **B** | bundle-uv | Isolated subprocess via uv (installed in ChimeraX) | uv added |
| **C** | bundle-uv | Isolated subprocess via external [uv](https://docs.astral.sh/uv/) | None |

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

### Option B: uv subprocess (uv installed in ChimeraX)

Install uv into ChimeraX's Python, then install the bundle. RDKit itself runs in an isolated subprocess — only uv is added to ChimeraX.

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

### Option C: uv subprocess (zero ChimeraX pollution)

If uv is already installed on your system, you can skip `pip install uv` entirely — nothing is added to ChimeraX's Python.

**From ChimeraX command line (GUI):**

```
cd /path/to/chimerax-rdkit-conformer/bundle-uv
devel install .
rdkconf uvPath ~/.local/bin/uv
```

**From terminal:**

```bash
cd /path/to/chimerax-rdkit-conformer/bundle-uv
ChimeraX --nogui --exit --cmd 'devel install .'
ChimeraX --nogui --exit --cmd 'rdkconf uvPath ~/.local/bin/uv'
```

The `uvPath` setting is saved persistently — you only need to set it once.

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
