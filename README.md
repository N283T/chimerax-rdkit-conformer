# ChimeraX-RDKitConformer

Generate 3D conformers from molecular notations using RDKit ETKDGv3, directly in ChimeraX.

## Requirements

- [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) 1.6+
- [uv](https://docs.astral.sh/uv/) (manages RDKit dependency automatically)

## Installation

### From source (development)

```bash
# Using echidna
echidna install

# Or directly with ChimeraX
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
ETKDGv3 locally via `uv run --script` subprocess for high-quality 3D coordinate
generation without network dependency. The 3D structure is built directly in memory
using ChimeraX's AtomicStructure API — no intermediate files are written.

## Development

```bash
echidna build                          # Build wheel
echidna install                        # Install to ChimeraX
echidna run                            # Build, install, and launch ChimeraX
echidna run --script scripts/smoke.cxc # Run with smoke test
```

## License

MIT
