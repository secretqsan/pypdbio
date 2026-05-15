# pylint: disable=missing-module-docstring
# pylint: disable=missing-function-docstring
# pylint: disable=protected-access
from math import pi
from dataclasses import dataclass, field

from .utils import index_of_chain_id, chain_id_of_index


@dataclass
class ObsoleteInfo:
    """
    PDB obsolete information class.
    """
    replace_date: str = ""
    new_entry_id: str = ""


@dataclass
class JournalInfo:
    """
    PDB journal information class.
    """
    author: list[str] = field(default_factory=list)
    title: str = ""
    editor: list[str] = field(default_factory=list)
    journal: str = ""
    volume: str = ""
    pages: str = ""
    year: str = ""
    publisher: str = ""
    issn: str = ""
    essn: str = ""
    pmid: str = ""
    doi: str = ""


@dataclass
class RevisionInfo:
    """
    PDB revision information class.
    """
    date: str = ""
    modifications: list[str] = field(default_factory=list)


@dataclass
class ReplaceInfo:
    """
    PDB replace information class.
    """
    date: str = ""
    ids: list[str] = field(default_factory=list)


@dataclass
class PdbMetaData:
    """
    PDB meta data class.
    """
    classification: str = ""
    date: str = ""
    pdb_id: str = ""
    author: list[str] = field(default_factory=list)
    split: list[str] = field(default_factory=list)
    obsolete: ObsoleteInfo = None
    caveat: str = ""
    compounds: list[dict] = field(default_factory=list)
    keywords: list[str] = field(default_factory=list)
    experiment: list[str] = field(default_factory=list)
    model_type: str = ""
    title = ""
    revisions: list[RevisionInfo] = field(default_factory=list)
    replace: ReplaceInfo = None
    journal: JournalInfo = None
    remark: dict = field(default_factory=dict)
    _compound_text: str = ""
    _source_text: str = ""


@dataclass
class CrystalInfo:
    """
    PDB Crystal information class.
    """
    space_group: str = "P 1"
    fake_crystallographic: bool = False
    cell_angles: list[float] = field(default_factory=list)
    cell_lengths: list[float] = field(default_factory=list)
    z: int = 1
    origin_matrix: list[list[float]] = field(default_factory=list)
    scale_matrix: list[list[float]] = field(default_factory=list)
    ncs_matrix: list[list[list[float]]] = field(default_factory=list)


@dataclass
class ConnectivityInfo:
    """Connectivity information class."""
    ss_bond: list[str] = field(default_factory=list)
    link: list[str] = field(default_factory=list)
    cis_peptide: list[str] = field(default_factory=list)
    connections: dict = field(default_factory=dict)


@dataclass
class SecondaryStructureInfo:
    """Secondary structure information class."""
    helix: dict = field(default_factory=dict)
    sheet: dict = field(default_factory=dict)
# secondary structure class


@dataclass
class SheetStrand:
    """Sheet class."""
    init_chain_id: str = ""
    init_seq_num: int = 0
    init_icode: str = ""

    end_chain_id: str = ""
    end_res_seq: int = 0
    end_icode: str = ""

    sense: int = 0

    cur_atom: str = ""
    cur_chain_id: str = ""
    cur_res_seq: int = None
    cur_icode: str = ""

    prev_atom: str = ""
    prev_chain_id: str = ""
    prev_res_seq: int = None
    prev_icode: str = ""


@dataclass
class Helix:
    """Helix class."""
    chain_id: str = ""
    init_seq_num: int = 0
    init_icode: str = ""
    end_seq_num: int = 0
    end_icode: str = ""
    helix_class: int = 0
    comment: str = ""


@dataclass
class Site:
    """Site class."""
    seq_num: int = 0
    chain_id: str = ""
    icode: str = ""


@dataclass
class Heterogen:
    """Heterogen class."""
    name: str = ""
    comment: str = ""
    _alias_text: str = ""
    formula: str = ""
    alias: list[str] = field(default_factory=list)


@dataclass
class SequenceDBInfo:
    """Sequence database information class."""
    init_seq_num: int = 0
    init_icode: str = ""
    end_seq_num: int = 0
    end_icode: str = ""
    database: str = ""
    db_accession: str = ""
    db_id_code: str = ""
    db_init_seq_num: int = 0
    db_end_seq_num: int = 0
    db_init_icode: str = ""
    db_end_icode: str = ""


@dataclass
class SequenceDifferenceInfo:
    """Sequence difference class."""
    alt_res_name: str = ""
    seq_num: int = 0
    icode: str = ""
    database: str = ""
    db_accession: str = ""
    db_res_num: int = 0
    db_res_name: str = ""
    comment: str = ""


@dataclass
class ResidueModificationInfo:
    """Residue modification class."""
    seq_num: int = 0
    std_res: str = ""
    comment: str = ""
    icode: str = ""


@dataclass
class SequenceInfo:
    """Sequence information class."""
    sequence_db: SequenceDBInfo = None
    sequence: list[str] = field(default_factory=list)
    sequence_differences: list[SequenceDifferenceInfo] = field(default_factory=list)
    residue_modifications: list[ResidueModificationInfo] = field(default_factory=list)


@dataclass
class NcsMatrix:
    """NCS matrix information class."""
    matrix: list[list[float]] = field(default_factory=list)
    given: bool = False


@dataclass
class SsBond:
    """SSBOND information class."""
    chain_id_1: str = ""
    seq_num_1: int = 0
    icode_1: str = ""
    chain_id_2: str = ""
    seq_num_2: int = 0
    icode_2: str = ""
    symmetry_operation_1: str = "1555"
    symmetry_operation_2: str = "1555"
    distance: float = 0.0


@dataclass
class Link:
    """LINK information class."""
    name_1: str = ""
    alt_loc_1: str = ""
    chain_id_1: str = ""
    seq_num_1: int = 0
    icode_1: str = ""
    name_2: str = ""
    alt_loc_2: str = ""
    chain_id_2: str = ""
    seq_num_2: int = 0
    icode_2: str = ""
    symmetry_operation_1: str = "1555"
    symmetry_operation_2: str = "1555"
    distance: float = 0.0


@dataclass
class CisPeptide:
    """Cis-peptide information class."""
    chain_id_1: str = ""
    seq_num_1: int = 0
    icode_1: str = ""
    chain_id_2: str = ""
    seq_num_2: int = 0
    icode_2: str = ""
    num_model: int = 0
    measure: float = 0.0


class IterableParentBase:
    """Class with a child list."""

    def __init__(self):
        self._children = []

    def __len__(self):
        return len(self._children)

    def __getitem__(self, index):
        return self._children[index]

    def __setitem__(self, index, value):
        self._children[index] = value
        value._parent = self

    def __delitem__(self, index):
        self._children[index]._parent = None
        del self._children[index]

    def __iter__(self):
        return iter(self._children)

    def add(self, value):
        self._children.append(value)
        value._parent = self

    def _find(self, value, should_increase_id_fn):
        index = 0
        for child in self._children:
            if should_increase_id_fn(child, value):
                index += 1
            if id(child) == id(value):
                break
        return index


class IterableChildBase:
    """Class with a parent."""
    _default_id = 1

    @staticmethod
    def _should_increase_id(a, b):
        return True

    def __init__(self):
        self._parent = None

    @property
    def id(self):
        if self._parent is not None:
            return self._parent._find(self, self._should_increase_id)
        else:
            return self._default_id


class PdbData(IterableParentBase):
    """
    PDB data class.
    """

    def __init__(self):
        super().__init__()
        self._tmp = {
            "primary": {},
        }
        self.secondary_structure = SecondaryStructureInfo()
        self.heterogen = {}
        self.meta = PdbMetaData()
        self.crystallographic = CrystalInfo()
        self.models = self._children
        self.index = 0
        self.connectivity = ConnectivityInfo()
        self.sites = {}
        self._validation_info = {}

    def add_model(self, model):
        self.add(model)


class Model(IterableParentBase, IterableChildBase):
    """Model class."""

    def __init__(self):
        super().__init__()
        self.chains = self._children
        self._occupied_chain_ids = []

    def add_chain(self, chain):
        self.add(chain)

    def __getitem__(self, index):
        if isinstance(index, str):
            self._occupied_chain_ids.clear()
            for chain in self.chains:
                if chain.id != "":
                    self._occupied_chain_ids.append(chain.id)
                    if chain.id == index:
                        return chain
            return self.chains[index_of_chain_id(index, self._occupied_chain_ids)]
        return super().__getitem__(index)

    def _find(self, value, should_increase_id_fn):
        target_index = -1
        occupied_chain_ids = []
        for chain in self.chains:
            if chain._malloc_id != "":
                occupied_chain_ids.append(chain._malloc_id)
            else:
                target_index += 1
            if id(chain) == id(value):
                break
        return chain_id_of_index(target_index, occupied_chain_ids)


class Chain(IterableParentBase, IterableChildBase):
    """Chain class."""
    _default_id = "A"

    def __init__(self):
        super().__init__()
        self.sequence_info = None
        self.residues = self._children
        self._malloc_id = ""

    def add_residue(self, residue):
        self.add(residue)

    def __getitem__(self, index):
        if isinstance(index, str):
            icode = index[-1].strip()
            seq_num = int(index[:-1])
            current_res_num = 0
            for residue in self.residues:
                if residue.icode == "":
                    current_res_num += 1
                if current_res_num == seq_num and residue.icode == icode:
                    return residue
            raise IndexError(f"Residue {index} not found")
        elif isinstance(index, int):
            return super().__getitem__(index)
        else:
            raise TypeError(f"Invalid index type: {type(index)}")

    @property
    def id(self):
        if self._malloc_id != "":
            return self._malloc_id
        else:
            return super().id

    @id.setter
    def id(self, value):
        self._malloc_id = value


class Residue(IterableParentBase, IterableChildBase):
    """Residue class."""
    _default_id = 1

    @staticmethod
    def _should_increase_id(a, b):
        return a.icode == ""

    def __init__(self, name):
        super().__init__()
        self.het = False
        self.solvent = False
        self.name = name
        self.atoms = self._children
        self.icode = ""
        self.end_of_chain = False

    def add_atom(self, atom):
        self.add(atom)


class Atom(IterableChildBase):
    """Atom class."""
    @staticmethod
    def _should_increase_id(a, b):
        return a.alt_loc == ""

    def __init__(
        self,
        name="",
        coord=(0.0, 0.0, 0.0),
        temp_factor=0.0,
        element="",
        charge=0,
        occupancy=1.0,
        alt_loc="",
    ):
        super().__init__()
        self.name = name
        self.coord = coord
        self.__temp_factor = temp_factor
        self.occupancy = occupancy
        self.alt_loc = alt_loc
        self.element = element
        self.charge = charge

    @property
    def temp_factor(self):
        if isinstance(self.__temp_factor, list):
            return 8 * pi ** 2 * sum(self.__temp_factor[0:3]) / 3
        else:
            return self.__temp_factor

    @property
    def anisotropic_temp_factor(self):
        if isinstance(self.__temp_factor, list):
            return self.__temp_factor
        else:
            return None

    @temp_factor.setter
    def temp_factor(self, value):
        self.__temp_factor = value
