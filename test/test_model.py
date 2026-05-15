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

class TestIndexing(unittest.TestCase):
    def test_model_indexing(self):
        pdb_data = PdbData()
        model_1 = Model()
        model_2 = Model()
        pdb_data.add_model(model_1)
        pdb_data.add_model(model_2)
        for index, model in enumerate(pdb_data.models):
            self.assertEqual(model.id, index + 1)

    def test_chain_indexing(self):
        model = Model()
        model.add_chain(Chain())
        model.add_chain(Chain())
        model.add_chain(Chain())
        for index, chain in enumerate(model.chains):
            self.assertEqual(chain.id, chr(ord('A') + index))
        model[1].id = "X"
        self.assertEqual(model[0].id, "A")
        self.assertEqual(model[1].id, "X")
        self.assertEqual(model[2].id, "B")

    def test_residue_indexing(self):
        chain = Chain()
        residue_1 = Residue("ALA")
        residue_2 = Residue("GLY")
        chain.add_residue(residue_1)
        chain.add_residue(residue_2)
        for index, residue in enumerate(chain.residues):
            self.assertEqual(residue.id, index + 1)
        chain[1].icode = "A"
        for index, residue in enumerate(chain.residues):
            self.assertEqual(residue.id, 1)

    def test_atom_indexing(self):
        residue = Residue("ALA")
        atom_1 = make_atom(name="N", coord=[1.1, 2.2, 3.3])
        atom_2 = make_atom(name="CA", coord=[4.4, 5.5, 6.6])
        residue.add_atom(atom_1)
        residue.add_atom(atom_2)
        for index, atom in enumerate(residue.atoms):
            self.assertEqual(atom.id, index + 1)
        residue[1].alt_loc = "A"
        for index, atom in enumerate(residue.atoms):
            self.assertEqual(atom.id, 1)

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
