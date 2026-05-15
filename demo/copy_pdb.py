from pypdbio import fetch, PdbReader, PdbWriter

fetch("1aki", path="output/1aki.pdb")
reader = PdbReader("output/1aki.pdb")
pdb_data = reader.read()
writer = PdbWriter("output/1aki_copy.pdb")
writer.write(pdb_data)
