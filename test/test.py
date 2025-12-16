from pypdbio import *
fetch('1aki')
data = PdbReader('1aki.pdb').read()
PdbWriter('1aki_copy.pdb').write(data)
chain = Chain()
residue = Residue('GLY')
atom1 = Atom('ATOM', 'CA', [1, 0, 1], 3.0, 'N')
atom2 = Atom('ATOM', 'N', [1, 2, 4], 3.5, 'N')
residue.add_atom(atom1)
residue.add_atom(atom2)
chain.add_residue(
    residue
)
PdbWriter('test.pdb').write(chain)