Changes
=======

1.2.3
-----

* Updated protonation of nitrogen, oxygen, and sulfur atoms to be compatible
  with the latest version of RDKit, which broke backwards compatibility.
* Added "silent" option to suppress all output.
* Added code to suppress unnecessary RDKit warnings.
* Updated copyright to 2020.

1.2.2
-----

* Added a new parameter to limit the number of variants per compound
  (`--max_variants`). The default is 128.

1.2.1
-----

* Corrected a bug that rarely misprotonated/deprotonated compounds with
  multiple ionization sites (e.g., producing a carbanion).

1.2
---

* Corrected a bug that led Dimorphite-DL to sometimes produce output molecules
  that are non-physical.
* Corrected a bug that gave incorrect protonation states for rare molecules
  (aromatic rings with nitrogens that are protonated when electrically
  neutral, e.g. pyridin-4(1H)-one).
* `run_with_mol_list()` now preserves non-string properties.
* `run_with_mol_list()` throws a warning if it cannot process a molecule,
  rather than terminating the program with an error.

1.1
---

* Dimorphite-DL now distinguishes between indoles/pyrroles and
  Aromatic_nitrogen_protonated.
* It is now possible to call Dimorphite-DL from another Python script, in
  addition to the command line. See the `README.md` file for instructions.

1.0
---

The original version described in:

Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019) Dimorphite-DL: An
open-source program for enumerating the ionization states of drug-like small
molecules. J Cheminform 11:14. doi:10.1186/s13321-019-0336-9.
