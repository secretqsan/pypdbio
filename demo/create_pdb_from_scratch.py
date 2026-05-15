from pypdbio import Atom, Chain, Model, PdbData, PdbWriter, Residue, set_unit

def main():
    set_unit("angstrom")

    pdb = PdbData()
    pdb.meta.title = "Synthetic triglycine, 2 models, icode + altLoc demo"
    pdb.meta.remark = {
        "0": "Built programmatically by demo/create_pdb_from_scratch.py"
    }
    for _ in range(2):
        m = Model()
        chain = Chain()
        for j in range(3):
            residue = Residue("GLY")
            x0 = j * 3.8
            y0 = 0.0
            z0 = 0.0
            residue.add_atom(Atom("N", (x0 + 0.87, y0 + 1.20, z0), element="N"))
            residue.add_atom(Atom("CA", (x0, y0, z0), element="C", occupancy=0.5))
            residue.add_atom(Atom("CA", (x0, y0, z0), element="C", occupancy=0.5, alt_loc="A"))
            residue.add_atom(Atom("C", (x0 + 1.52, y0, z0), element="C"))
            residue.add_atom(Atom("O", (x0 + 1.93, y0 + 1.09, z0), element="O"))
            chain.add_residue(residue)
        chain.residues[1].icode = "A"
        chain.residues[-1].end_of_chain = True
        m.add_chain(chain)
        pdb.add_model(m)
    PdbWriter("output/from_scratch.pdb").write(pdb)


if __name__ == "__main__":
    main()
