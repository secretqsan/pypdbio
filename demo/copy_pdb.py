from pypdbio import fetch, PdbReader, PdbWriter

#fetch("1aki", path="1aki.pdb")
reader = PdbReader("pdb_all_feature.pdb")
pdb_data = reader.read()
writer = PdbWriter("pdb_test_copy.pdb")
writer.write(pdb_data)
