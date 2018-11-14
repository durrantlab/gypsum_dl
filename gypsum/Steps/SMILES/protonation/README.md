Protonation 1.0
===============

What is it?
-----------

Protonation 1.0 adds hydrogen atoms to molecular representations, as
appropriate for a user-specified pH range.

Users can provide SMILES strings from the command line or via an .smi file.

Licensing
---------

Protonation is released under the Apache 2.0 license. See LICENCE.txt for 
details.

Usage
-----

```
usage: protonate.py [-h] [--min_ph MIN] [--max_ph MAX] [--st_dev STD]
                    [--smiles SMILES] [--smiles_file SMILES_FILE]
                    [--output_file OUTPUT_FILE] [--label_states]

Protonates small molecules.

optional arguments:
  -h, --help            show this help message and exit
  --min_ph MIN          Minimum pH to consider.
  --max_ph MAX          Maximum pH to consider.
  --st_dev STD          Standard devation range (number of standard devs).
  --smiles SMILES       SMILE string to protonate.
  --smiles_file SMILES_FILE
                        File which contains SMILES strings to protonate.
  --output_file OUTPUT_FILE
                        File to write protonated SMILES. (Optional)
  --label_states        Label protonated SMILES with target state 
                        ("DEPROTONATED", "PROTONATED", or "BOTH").
```

The default pH range is 6.4 to 8.4, considered biologically relevant pH.

Examples
--------

```
python protonate.py --smiles_file sample_molecules.smi
```

```
python protonate.py --smiles "CCC(=O)O" --min_ph -3 --max_ph -2
```

```
python protonate.py --smiles "CCCN" --min_ph -3 --max_ph -2 --output_file output.smi
```

```
python protonate.py --smiles_file sample_molecules.smi --st_dev 2.0 --label_states
```

Authors and Contacts
--------------------

See the `CONTRIBUTORS.md` file for a full list of contributors. Please contact
Jacob Durrant (durrantj@pitt.edu) with any questions.