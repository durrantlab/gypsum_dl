Dimorphite-DL 1.0
=================

What is it?
-----------

Dimorphite-DL adds hydrogen atoms to molecular representations, as appropriate
for a user-specified pH range. It is a fast, accurate, accessible, and modular
open-source program for enumerating small-molecule ionization states.

Users can provide SMILES strings from the command line or via an .smi file.

Licensing
---------

Protonation is released under the Apache 2.0 license. See LICENCE.txt for
details.

Usage
-----

```
usage: dimorphite_dl.py [-h] [--min_ph MIN] [--max_ph MAX]
                        [--pka_precision PRE] [--smiles SMI]
                        [--smiles_file FILE] [--output_file FILE]
                        [--label_states] [--test]

Dimorphite 1.0: Creates models of appropriately protonated small moleucles.
Apache 2.0 License. Copyright 2018 Jacob D. Durrant.

optional arguments:
  -h, --help           show this help message and exit
  --min_ph MIN         minimum pH to consider (default: 6.4)
  --max_ph MAX         maximum pH to consider (default: 8.4)
  --pka_precision PRE  pKa precision factor (number of standard devations,
                       default: 1.0)
  --smiles SMI         SMILES string to protonate
  --smiles_file FILE   file that contains SMILES strings to protonate
  --output_file FILE   output file to write protonated SMILES (optional)
  --label_states       label protonated SMILES with target state (i.e.,
                       "DEPROTONATED", "PROTONATED", or "BOTH").
  --test               run unit tests (for debugging)
```

The default pH range is 6.4 to 8.4, considered biologically relevant pH.

Examples
--------

```
  python dimorphite_dl.py --smiles_file sample_molecules.smi
  python dimorphite_dl.py --smiles "CCC(=O)O" --min_ph -3.0 --max_ph -2.0
  python dimorphite_dl.py --smiles "CCCN" --min_ph -3.0 --max_ph -2.0 --output_file output.smi
  python dimorphite_dl.py --smiles_file sample_molecules.smi --pka_precision 2.0 --label_states
  python dimorphite_dl.py --test
```

Authors and Contacts
--------------------

See the `CONTRIBUTORS.md` file for a full list of contributors. Please contact
Jacob Durrant (durrantj@pitt.edu) with any questions.
