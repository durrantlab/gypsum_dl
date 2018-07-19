Protonation
===========

The Protonation Module takes smiles as either input line stings or from a .smi 
file and returns the protonated smiles, given a specific pH range.

Default pH range is 6.4 to 8.4, considered biologically relevant pH.

Example Run:
```bash
python protonate.py --smiles_file sample_molecules.smi
```