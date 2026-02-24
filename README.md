# ChimeraX-RDKitConformer

Generate 3D conformers from molecular notations using RDKit ETKDGv3, directly in ChimeraX.

## Requirements

- [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) 1.6+ (1.11+ for `minimize` keyword)
- [uv](https://docs.astral.sh/uv/) (manages RDKit dependency automatically)

## Installation

### From source (development)

```bash
ChimeraX --nogui --exit --cmd 'devel install .'
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
rdkconf CCO minimize true                          # ChimeraX AMBER minimization (1.11+)
rdkconf CCO optimize true minimize true            # Both: MMFF then AMBER
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

### uv Setup

This bundle requires [uv](https://docs.astral.sh/uv/) to manage the RDKit subprocess.

**Option A: Install uv into ChimeraX via pip**

Run the following in ChimeraX's command line:

```
pip install uv
```

This places the uv binary in ChimeraX's user directory and is automatically detected by the bundle. The RDKit subprocess runs in an isolated environment managed by uv, so ChimeraX's own packages are not affected. You may need to re-run `pip install uv` after updating ChimeraX.

**Option B: Use an existing uv installation**

If uv is already installed on your system (`which uv` to check), you can point the bundle to it. This is useful when ChimeraX is launched from a GUI where the system PATH may not include user-installed tools:

```
rdkconf uvPath                         # Show current setting and resolved path
rdkconf uvPath ~/.local/bin/uv         # Set path persistently
rdkconf uvPath ""                      # Reset to auto-detect
```

If uv is not installed yet, follow the [installation guide](https://docs.astral.sh/uv/getting-started/installation/).

## How It Works

ChimeraX's built-in SMILES support (`open smiles:CCO`) depends on NCI's web service
(requires internet, no control over 3D generation quality). This bundle uses RDKit
ETKDGv3 locally via `uv run --script` subprocess for high-quality 3D coordinate
generation without network dependency. The 3D structure is built directly in memory
using ChimeraX's AtomicStructure API — no intermediate files are written.

## Development

```bash
ChimeraX --nogui --exit --cmd 'devel build .'    # Build wheel
ChimeraX --nogui --exit --cmd 'devel install .'   # Install to ChimeraX
```

## License

MIT
