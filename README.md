# Gypsum

Gypsum is a python program that takes SMILE representations of small molecules and generates 3D models as SDF files. It searches various states of protonation, tautomerization, chirality, and double bonds to find viable, low-energy conformations.

## Getting Started

To run Gypsum, acquire a copy of this repository, either by git clone of this repository or by downloading and unzipping. Install required dependencies via your favorite python package manager, and the program should be good to go.

### Prerequisites

There a few python libraries that gypsum depends on. We suggest using Anaconda to manage packages. The following command will install the necessary packages

```
conda install -c rdkit rdkit numpy scipy
```

## Example

To run a demo of gypsum type the following command in the main gypsum folder:

```
python ./run_gypsum.py -j sample_molecules.json
```

## License

TODO: Add License

## Acknowledgments

TODO: Compile list of authors and contributors
