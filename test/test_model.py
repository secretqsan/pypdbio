import unittest
from pypdbio import Model, Chain, Residue, Atom, PdbData

def make_atom(name="CA", coord=None):
    if coord is None:
        coord = [1.0, 2.0, 3.0]
    return Atom(
        name=name,
        coord=coord,
        temp_factor=10.0,
        element="C",
        charge=""
    )


class TestAtomCoordinates(unittest.TestCase):
    def test_standalone_atoms_coord(self):
        atom1 = make_atom(name="N", coord=[7.0, 8.0, 9.0])

        self.assertEqual(atom1.coord, [7.0, 8.0, 9.0])
        atom1.coord = [1.0, 2.0, 3.0]
        self.assertEqual(atom1.coord, [1.0, 2.0, 3.0])

    def test_residue_atoms_coord(self):
        residue = Residue("ALA")
        atom1 = make_atom(name="N", coord=[1.1, 2.2, 3.3])
        residue.add_atom(atom1)

        self.assertEqual(atom1.coord, [1.1, 2.2, 3.3])
        atom1.coord = [4.4, 5.5, 6.6]
        self.assertEqual(atom1.coord, [4.4, 5.5, 6.6])
        residue[0].coord = [7.7, 8.8, 9.9]
        self.assertEqual(residue[0].coord, [7.7, 8.8, 9.9])

    def test_chain_atoms_coord(self):
        chain = Chain()
        residue = Residue("GLY")
        atom1 = make_atom(name="N", coord=[9.0, 8.0, 7.0])
        residue.add_atom(atom1)
        chain.add_residue(residue)

        self.assertEqual(atom1.coord, [9.0, 8.0, 7.0])

    def test_model_atoms_coord(self):
        model = Model()
        chain = Chain()
        residue = Residue("SER")
        atom1 = make_atom(name="N", coord=[2.0, 2.0, 2.0])
        residue.add_atom(atom1)
        chain.add_residue(residue)
        model.add_chain(chain)

        self.assertEqual(atom1.coord, [2.0, 2.0, 2.0])
    
    def test_pdb_data_atoms_coord(self):
        pdb_data = PdbData()
        model = Model()
        chain = Chain()
        residue = Residue("SER")
        atom1 = make_atom(name="N", coord=[2.0, 2.0, 2.0])
        residue.add_atom(atom1)
        chain.add_residue(residue)
        model.add_chain(chain)
        pdb_data.add_model(model)

        self.assertEqual(atom1.coord, [2.0, 2.0, 2.0])

if __name__ == "__main__":
    unittest.main()
