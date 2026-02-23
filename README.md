# ChimeraX-RDKitSMILES

Generate 3D molecules from SMILES strings using RDKit ETKDGv3, directly in ChimeraX.

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
rdksmiles CCO                          # Ethanol
rdksmiles c1ccccc1                     # Benzene
rdksmiles "CC(=O)Oc1ccccc1C(=O)O"     # Aspirin
rdksmiles CCO name EtOH               # Custom residue name
rdksmiles CCO hydrogen false           # Hide hydrogens
rdksmiles CCO output ~/ethanol.sdf     # Save SDF file
```

## How It Works

ChimeraX's built-in SMILES support (`open smiles:CCO`) depends on NCI's web service
(requires internet, no control over 3D generation quality). This bundle uses RDKit
ETKDGv3 locally via `uv run --script` subprocess for high-quality 3D coordinate
generation without network dependency.

## Development

```bash
echidna build                          # Build wheel
echidna install                        # Install to ChimeraX
echidna run                            # Build, install, and launch ChimeraX
echidna run --script scripts/smoke.cxc # Run with smoke test
```

## License

MIT
