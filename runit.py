import sys
from gypsum.Start import ConfGenerator

# # Test chiral. One chiral center specified, one not.
# ConfGenerator(['I[C@@](C)(F)C(I)(F)C'], 5, 9)

# # Test double-bond enumeration. Should be two.
# ConfGenerator(['FC(I)=C(Cl)C'], 5, 9)

# # Test Tautomers
# ConfGenerator(['C1=CC=CNC(=O)1'], 5, 9)

# # Test Tautomers from including bad tautomers. This also tests protonation at various states.
# ConfGenerator(['N[C@@H](CC1=CNC=N1)C(O)=O'], 5, 9)

# # Test non-aromatic rings. Should generate multiple ring confomers.
#ConfGenerator(['C1(C=CC=C2)=C2CCCC1'], 5, 9)

# Test all
#ConfGenerator(['C1(C=CC=C2)=C2CCCC1', 'N[C@@H](CC1=CNC=N1)C(O)=O', 'C1=CC=CNC(=O)1', 'FC(I)=C(Cl)C', 'I[C@@](C)(F)C(I)(F)C'], 5, 9)
#ConfGenerator("samples.smi", 5, 9)

#ConfGenerator("sample_parameters.json")
ConfGenerator(sys.argv[1])
