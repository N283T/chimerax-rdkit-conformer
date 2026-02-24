# Multi-Format Input Support Design

Date: 2026-02-23

## Goal

Extend `rdkconf` to accept InChI, FASTA, Sequence, and HELM inputs in addition to SMILES.

## Supported Formats

| Format | RDKit Function | Detection | `format` keyword |
|--------|---------------|-----------|-----------------|
| SMILES | `MolFromSmiles()` | Default (fallback) | `smiles` |
| InChI | `MolFromInchi()` | Auto: starts with `InChI=` | `inchi` |
| FASTA | `MolFromFASTA()` | Explicit only | `fasta` |
| Sequence | `MolFromSequence()` | Explicit only | `sequence` |
| HELM | `MolFromHELM()` | Explicit only | `helm` |

## UX

```
rdkconf CCO                                       # SMILES (default)
rdkconf "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3"     # Auto-detect InChI
rdkconf AGCTTGA format dna                         # DNA sequence
rdkconf PEPTIDE1{A.G.C}$$ format helm              # HELM notation
```

## Architecture

Single command `rdkconf` with optional `format` keyword. No new commands.

### Auto-detect Flow

```
input received
  -> format keyword provided? -> yes -> use that format
  -> no
  -> starts with "InChI=" ? -> yes -> inchi
  -> no
  -> default -> smiles
```

### Changes

1. **`smiles_to_sdf.py` -> `input_to_sdf.py`** (rename)
   - Add `--format` CLI flag
   - Add `parse_input()` dispatcher function
   - Add `rdkit.Chem.inchi.MolFromInchi` import

2. **`cmd.py`**
   - Add `format` keyword argument (`StringArg`)
   - Add auto-detect logic (`InChI=` prefix check)
   - Pass `--format` to subprocess

3. **`pyproject.toml`**
   - Update `package-data` script name

4. **`__init__.py`** - No changes needed

### Error Handling

- Invalid format keyword -> `UserError` with list of valid formats
- Parse failure -> `UserError` with format-specific message
- HELM limitations -> error from RDKit propagated as-is

## Decisions

- Auto-detect only for InChI (unambiguous `InChI=` prefix)
- FASTA, Sequence, HELM require explicit `format` keyword
- HELM included despite limitations (standard monomers only)
- Single command approach (no separate `rdkinchi`, `rdkseq` commands)
