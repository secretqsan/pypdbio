import requests
from functools import singledispatchmethod
def fetch(pdb_id, path = ''):
    pdb_filename = pdb_id + '.pdb'
    url = "https://files.rcsb.org/download/" + pdb_filename
    data = requests.get(url)

    if data.status_code != 200:
        raise FileNotFoundError(f"PDB ID {pdb_id} dont exist.")

    if path == '':
        path = pdb_filename

    with open(path, 'wb') as f:
        f.write(data.content)

class PdbMetaData:
    def __init__(self):
        self.classification = ""
        self.date = ""
        self.pdb_id = ""
        self.author = {}
        self.title = {}
        self.remark = {}

class PdbData:
    def __init__(self):
        self.meta = PdbMetaData()
        self.index = 0
        self.models = []
        self.connections = {}

    def __iter__(self):
        return self.models

class Model:
    def __init__(self):
        self.chains = []

    def add_chain(self, chain):
        self.chains.append(chain)

class Chain:
    def __init__(self):
        self.residues = []

    def add_residue(self, residue):
        self.residues.append(residue)

class Residue:
    def __init__(self, name):
        self.name = name
        self.atoms = []

    def add_atom(self, atom):
        self.atoms.append(atom)

class Atom:
    def __init__(self, typ="", name="", coord=(0.0, 0.0, 0.0), temp_factor=0.0, element="", charge=0):
        self.typ = typ
        self.name = name
        self.coord = coord
        self.temp_factor = temp_factor
        self.element = element
        self.charge = charge

class PdbReader:
    def __init__(self, path):
        """ Initializes the PdbReader with the given file path.
        :param path: Path to the PDB file to be read"""
        self.pdb_file = open(path, 'r')
        self.__current_model = -1
        self.__current_chain_id = ''
        self.__current_residue_id = -1

    @staticmethod
    def __parse_header(line):
        classification = line[10:50].strip()
        date = line[50:59].strip()
        pdb_id = line[62:66].strip()
        return classification, date, pdb_id

    @staticmethod
    def __parse_title(line):
        title_id = line[8:10].strip()
        if title_id == '':
            title_id = '0'
        title = line[10:80].strip()
        return title_id, title

    @staticmethod
    def __parse_author(line):
        author_id = line[8:10].strip()
        if author_id == '':
            author_id = '0'
        author_list = line[10:80].strip()
        return author_id, author_list

    @staticmethod
    def __parse_remark(line):
        remark_id = line[7:10].strip()
        if remark_id == '':
            remark_id = '0'
        remark_text = line[11:80].strip()
        return remark_id, remark_text

    @staticmethod
    def __parse_atom(line):
        typ = line[0:6].strip()
        atom_no = line[6:11].strip()
        atom_name = line[12:16].strip()
        res_name = line[17:20].strip()
        chain_id = line[21]
        residue_id = line[22:26].strip()
        coord_x = line[30:38].strip()
        coord_y = line[38:46].strip()
        coord_z = line[46:54].strip()
        temp_factor = line[60:66].strip()
        element = line[76:78].strip()
        charge = line[78:80].strip()
        return (typ, int(atom_no), atom_name, res_name, chain_id, int(residue_id), float(coord_x) / 10,
                     float(coord_y) / 10,
                     float(coord_z) / 10,
                     float(temp_factor),
                     element,
                     0 if charge == '' else int(charge)
        )

    def __parse_connect(self, line):
        atom1 = int(line[6:11].strip())
        bonded_atoms = []
        for i in range(11, 31, 5):
            bonded_atom_str = line[i:i + 5].strip()
            if bonded_atom_str != '':
                bonded_atom_index =  int(bonded_atom_str)
                if bonded_atom_index > atom1:
                    bonded_atoms.append(int(bonded_atom_str))
        return atom1, bonded_atoms

    def read(self):
        """ Reads the PDB file and parses it into a PdbData structure."""
        pdb_data = PdbData()
        for line in self.pdb_file:
            line = line[:-1]
            if line.startswith("HEADER"):
                classification, date, pdb_id = self.__parse_header(line)
                pdb_data.meta.date = date
                pdb_data.meta.classification = classification
                pdb_data.meta.pdb_id = pdb_id
            elif line.startswith("TITLE"):
                title_id, title = self.__parse_title(line)
                if title_id in pdb_data.meta.title:
                    pdb_data.meta.title[title_id] += '\n' + title
                else:
                    pdb_data.meta.title[title_id] = title
            elif line.startswith("AUTHOR"):
                author_id, author_list = self.__parse_author(line)
                if author_id in pdb_data.meta.title:
                    pdb_data.meta.title[author_id] += '\n' + author_list
                else:
                    pdb_data.meta.title[author_id] = author_list
            elif line.startswith("REMARK"):
                remark_id, remark_text = self.__parse_remark(line)
                if remark_id in pdb_data.meta.remark:
                    pdb_data.meta.remark[remark_id] += '\n' + remark_text
                else:
                    pdb_data.meta.remark[remark_id] = remark_text
            elif line.startswith('CONECT'):
                atom1, bonded_atoms = self.__parse_connect(line)
                pdb_data.connections[atom1] = bonded_atoms
            elif line.startswith('MODEL'):
                self.__current_model += 1
                pdb_data.models.append(Model())
                self.__current_chain_id = ''
                self.__current_residue_id = -1

            elif line.startswith('ATOM') or line.startswith('HETATM'):
                if self.__current_model == -1:
                    self.__current_model += 1
                    pdb_data.models.append(Model())
                (typ, atom_no, atom_name, res_name, chain_id, residue_id, coord_x, coord_y, coord_z,
                 temp_factor, element, charge) = self.__parse_atom(line)
                if chain_id != self.__current_chain_id:
                    self.__current_chain_id = chain_id
                    pdb_data.models[self.__current_model].add_chain(Chain())
                    self.__current_residue_id = -1
                if residue_id != self.__current_residue_id:
                    self.__current_residue_id = residue_id
                    pdb_data.models[self.__current_model].chains[-1].add_residue(Residue(res_name))
                pdb_data.models[-1].chains[-1].residues[-1].add_atom(Atom(
                    typ=typ,
                    name=atom_name,
                    coord=(coord_x, coord_y, coord_z),
                    temp_factor=temp_factor,
                    element=element,
                    charge=charge
                ))
        return pdb_data

class PdbWriter:
    def __init__(self, path):
        """ Initializes the PdbWriter with the given file path.
        :param path: Path to the PDB file to be written"""
        self.pdb_file = open(path, 'w')
        self.last_atom_info = {}

    @singledispatchmethod
    def write(self, data):
        """ Writes the given data to the PDB file.
        :param data: Data to be written (PdbData, Model, or Chain)"""
        pass

    @write.register(PdbData)
    def _(self, data):
        self.__write_pdb_data(data)

    @write.register(Model)
    def _(self, data):
        pdb_data = PdbData()
        pdb_data.models.append(data)
        meta = PdbMetaData()
        meta.title = {'0': 'pdb'}
        meta.remark = {'0': 'Generated by pypdbio'}
        pdb_data.meta = meta
        self.write(pdb_data)

    @write.register(Chain)
    def _(self, data):
        model = Model()
        model.add_chain(data)
        self.write(model)

    def __write_pdb_data(self, pdb_data):
        meta = pdb_data.meta
        self.__write_header(meta.classification, meta.date, meta.pdb_id)
        self.__write_title(meta.title)
        self.__write_author(meta.author)
        self.__write_remark(meta.remark)
        chain_id = 'A'
        for model in pdb_data.models:
            if len(pdb_data.models) > 1:
                self.pdb_file.write(f'MODEL     {model.index + 1:>4}\n')

            residue_id = 1
            atom_id = 1
            for chain in model.chains:
                for residue in chain.residues:
                    for atom in residue.atoms:

                        atom_info = {
                            "type": atom.typ,
                            "atom_no": atom_id,
                            "atom_name": atom.name,
                            "residue_name": residue.name,
                            "chain_id": chain_id,
                            "residue_id": residue_id,
                            "coord_x": atom.coord[0],
                            "coord_y": atom.coord[1],
                            "coord_z": atom.coord[2],
                            "temp_factor": atom.temp_factor,
                            "element": atom.element,
                            "charge": atom.charge
                        }
                        if self.last_atom_info.get('type') == 'ATOM' and atom_info["type"] == 'HETATM':
                            self.__write_ter()
                        self.__write_atom(atom_info)
                        atom_id += 1
                    residue_id += 1
                if self.last_atom_info.get('type') == 'ATOM':
                    self.__write_ter()
                chain_id = chr(ord(chain_id) + 1)
            if len(pdb_data.models) > 1:
                self.pdb_file.write('ENDMDL\n')
        self.__write_connect(pdb_data.connections)
        self.pdb_file.write('END\n')

    def __write_header(self, classification, date, pdb_id):
        line = f'HEADER    {classification:<40}{date:>9}   {pdb_id:>4}\n'
        self.pdb_file.write(line)

    def __write_title(self, title_dict):
        for title_id, title in sorted(title_dict.items()):
            lines = title.split('\n')
            for line in lines:
                self.pdb_file.write(f'TITLE   {"" if title_id=='0' else title_id:>2}{line:<70}\n')

    def __write_author(self, author_dict):
        for author_id, author_list in sorted(author_dict.items()):
            lines = author_list.split('\n')
            for line in lines:
                self.pdb_file.write(f'AUTHOR  {"" if author_id=='0' else author_id:>2}{line:<70}\n')

    def __write_remark(self, remark_dict):
        for remark_id, remark_text in sorted(remark_dict.items()):
            lines = remark_text.split('\n')
            for line in lines:
                self.pdb_file.write(f'REMARK {"" if remark_id=='0' else remark_id:>3}{line:<70}\n')

    def __write_atom(self, atom_info):
        line = f'{atom_info["type"]:<6}{atom_info["atom_no"]:>5} {atom_info["atom_name"]:<4} '
        line += f'{atom_info["residue_name"]:>3} {atom_info["chain_id"]}{atom_info['residue_id']:>4}    '
        line += f'{atom_info["coord_x"] * 10:8.3f}{atom_info["coord_y"] * 10:8.3f}{atom_info["coord_z"] * 10:8.3f}'
        line += f'{1.0:6.2f}{atom_info["temp_factor"]:6.2f}          {atom_info["element"]:>2}'
        line += f'{''if atom_info["charge"] == 0 else atom_info["charge"]:>2}\n'
        self.pdb_file.write(line)
        self.last_atom_info = atom_info

    def __write_connect(self, connections):
        connections_with_dup = {}
        for atom1, bonded_atoms in connections.items():
            for atom2 in bonded_atoms:
                if atom1 not in connections_with_dup:
                    connections_with_dup[atom1] = []
                connections_with_dup[atom1].append(atom2)
                if atom2 not in connections_with_dup:
                    connections_with_dup[atom2] = []
                connections_with_dup[atom2].append(atom1)
        for atom1 in sorted(connections_with_dup.keys()):
            bonded_atoms = sorted(set(connections_with_dup[atom1]))
            line = f'CONECT{atom1:>5}'
            for atom2 in bonded_atoms:
                line += f'{atom2:>5}'
            line += '\n'
            self.pdb_file.write(line)

    def __write_ter(self):
        line = f'TER   {self.last_atom_info["atom_no"]:>5}      '
        line += f'{self.last_atom_info["residue_name"]:>3} {self.last_atom_info["chain_id"]}'
        line += f'{self.last_atom_info['residue_id']:>4} \n'
        self.pdb_file.write(line)