"""Unit tests for the PDB Writer.
"""
# pylint: disable=missing-function-docstring
# pylint: disable=missing-class-docstring
# pylint: disable=unused-wildcard-import
# pylint: disable=wildcard-import
import unittest
import os
from contextlib import contextmanager
from pypdbio import set_unit, PdbWriter
from pypdbio.writer_helper import *
from pypdbio.models import *
set_unit("A")

TEMP_PDB_PATH = "test_pdb.pdb"

@contextmanager
def write_temp_pdb(pdb_data, start_line=None, end_line=None):
    writer = PdbWriter(TEMP_PDB_PATH)
    writer.write(pdb_data)
    with open(TEMP_PDB_PATH, 'r', encoding='utf-8') as f:
        content = f.readlines()
    try:
        yield "".join(content[start_line:end_line])
    finally:
        os.remove(TEMP_PDB_PATH)



class TestPdbHelperFunctions(unittest.TestCase):
    def test_split_text_field(self):
        self.assertEqual(
            split_text("This is a test string", width=10, prefix_space=True),
            ["This is a", " test", " string"],
        )

    def test_line_template_field(self):
        self.assertEqual(
            gen_line_from_template(
                "{field1:>8}{field2:>8}",
                {"field1": {"value": "123"}, "field2": {"value": "456"}},
            )[0],
            "     123     456",
        )


class TestPdbHeaderWriter(unittest.TestCase):
    def test_header_generator(self):
        self.assertEqual(
            gen_header(
                "TRANSFERASE/TRANSFERASE INHIBITOR", "17-SEP-04", "1XH6"
            ),
            ["HEADER    TRANSFERASE/TRANSFERASE INHIBITOR       17-SEP-04   1XH6"],
        )

    def test_write_header(self):
        pdb_data = PdbData()
        pdb_data.meta.classification = "HYDROLASE"
        pdb_data.meta.date = "19-MAY-97"
        pdb_data.meta.pdb_id = "1AKI"
        with write_temp_pdb(pdb_data, 0, 1) as content:
            self.assertEqual(
                content,
                "HEADER    HYDROLASE                               19-MAY-97   1AKI              \n"
            )

    def test_title_generator(self):
        self.assertEqual(
            gen_title(
                "RHIZOPUSPEPSIN COMPLEXED WITH REDUCED PEPTIDE INHIBITOR",
            ),
            ["TITLE     RHIZOPUSPEPSIN COMPLEXED WITH REDUCED PEPTIDE INHIBITOR               "],
        )
        self.assertEqual(
            gen_title(
                "THE STRUCTURE OF THE ORTHORHOMBIC FORM OF HEN EGG-WHITE LYSOZYME AT "
                "1.5 ANGSTROMS RESOLUTION",
            ),
            [
                "TITLE     THE STRUCTURE OF THE ORTHORHOMBIC FORM OF HEN EGG-WHITE LYSOZYME AT   ",
                "TITLE    2 1.5 ANGSTROMS RESOLUTION                                             ",
            ],
        )

    def test_write_title(self):
        pdb_data = PdbData()
        pdb_data.meta.title = (
            "THE STRUCTURE OF THE ORTHORHOMBIC FORM OF HEN EGG-WHITE "
            "LYSOZYME AT 1.5 ANGSTROMS RESOLUTION"
        )
        with write_temp_pdb(pdb_data, 0, 2) as content:
            self.assertEqual(
                content,
                "TITLE     THE STRUCTURE OF THE ORTHORHOMBIC FORM OF HEN EGG-WHITE LYSOZYME AT   \n" + \
                "TITLE    2 1.5 ANGSTROMS RESOLUTION                                             \n"
            )

    def test_obslte_generator(self):
        self.assertEqual(
            gen_obslte("31-JAN-94", "2MBP", "1MBP"),
            ["OBSLTE     31-JAN-94 1MBP      2MBP"],
        )

    def test_write_obslte(self):
        pdb_data = PdbData()
        pdb_data.meta.pdb_id = "1MBP"
        pdb_data.meta.obsolete = ObsoleteInfo(
            replace_date="31-JAN-94", new_entry_id="2MBP")
        with write_temp_pdb(pdb_data, 1, 2) as content:
            self.assertEqual(
                content,
                "OBSLTE     31-JAN-94 1MBP      2MBP                                             \n"
            )

    def test_split_generator(self):
        self.assertEqual(
            gen_split(
                [
                    "1VOQ", "1VOR", "1VOS", "1VOU", "1VOV",
                    "1VOW", "1VOX", "1VOY", "1VP0", "1VOZ", "1VPA",
                ],
            ),
            [
                "SPLIT      1VOQ 1VOR 1VOS 1VOU 1VOV 1VOW 1VOX 1VOY 1VP0 1VOZ 1VPA               ",
            ],
        )

    def test_write_split(self):
        pdb_data = PdbData()
        pdb_data.meta.split = [
            "1VOQ", "1VOR", "1VOS", "1VOU", "1VOV",
            "1VOW", "1VOX", "1VOY", "1VP0", "1VOZ", "1VPA",
        ]
        with write_temp_pdb(pdb_data, 0, 1) as content:
            self.assertEqual(
                content,
                "SPLIT      1VOQ 1VOR 1VOS 1VOU 1VOV 1VOW 1VOX 1VOY 1VP0 1VOZ 1VPA               \n"
            )

    def test_caveat_generator(self):
        self.assertEqual(
            gen_caveat(
                "2WMW",
                ["CHAIN A IS MISSING RESIDUES AND CHAIN B IS MISSING RESIDUES"],
            ),
            [
                "CAVEAT     2WMW    CHAIN A IS MISSING RESIDUES AND CHAIN B IS MISSING RESIDUES ",
            ],
        )

    def test_write_caveat(self):
        pdb_data = PdbData()
        pdb_data.meta.pdb_id = "2WMW"
        pdb_data.meta.caveat = "CHAIN A IS MISSING RESIDUES"
        with write_temp_pdb(pdb_data, 1, 2) as content:
            self.assertEqual(
                content,
                "CAVEAT     2WMW    CHAIN A IS MISSING RESIDUES                                  \n"
            )

    def test_compound_generator(self):
        compound_texts = [
            "MOL_ID: 1;",
            "MOLECULE: HEMOGLOBIN ALPHA CHAIN;",
            "CHAIN: A,  C;",
            "MOL_ID: 2;",
            "MOLECULE: HEMOGLOBIN BETA CHAIN",
        ]
        self.assertEqual(
            gen_compnd(compound_texts),
            [
                "COMPND    MOL_ID: 1;                                                            ",
                "COMPND   2 MOLECULE: HEMOGLOBIN ALPHA CHAIN;                                    ",
                "COMPND   3 CHAIN: A,  C;                                                        ",
                "COMPND   4 MOL_ID: 2;                                                           ",
                "COMPND   5 MOLECULE: HEMOGLOBIN BETA CHAIN                                      ",
            ],
        )

    def test_source_generator(self):
        source_texts = [
            "MOL_ID: 1;",
            "ORGANISM_SCIENTIFIC: AVIAN SARCOMA VIRUS;",
            "MOL_ID: 2;",
            "ORGANISM_SCIENTIFIC: AVIAN SARCOMA VIRUS",
        ]
        self.assertEqual(
            gen_source(source_texts),
            [
                "SOURCE    MOL_ID: 1;                                                           ",
                "SOURCE   2 ORGANISM_SCIENTIFIC: AVIAN SARCOMA VIRUS;                           ",
                "SOURCE   3 MOL_ID: 2;                                                          ",
                "SOURCE   4 ORGANISM_SCIENTIFIC: AVIAN SARCOMA VIRUS                            ",
            ],
        )

    def test_write_compound_and_source(self):
        pdb_data = PdbData()
        pdb_data.meta.compounds = [
            {
                "molecule": "HEMOGLOBIN ALPHA CHAIN",
                "chain": "A,  C",
                "source": {"organism_scientific": "AVIAN SARCOMA VIRUS"},
            },
            {
                "molecule": "HEMOGLOBIN BETA CHAIN",
                "source": {"organism_scientific": "AVIAN SARCOMA VIRUS"},
            },
        ]
        with write_temp_pdb(pdb_data, 0, 9) as content:
            self.assertEqual(
                content,
                "COMPND    MOL_ID: 1;                                                            \n" + \
                "COMPND   2 MOLECULE: HEMOGLOBIN ALPHA CHAIN;                                    \n" + \
                "COMPND   3 CHAIN: A,  C;                                                        \n" + \
                "COMPND   4 MOL_ID: 2;                                                           \n" + \
                "COMPND   5 MOLECULE: HEMOGLOBIN BETA CHAIN                                      \n" + \
                "SOURCE    MOL_ID: 1;                                                            \n" + \
                "SOURCE   2 ORGANISM_SCIENTIFIC: AVIAN SARCOMA VIRUS;                            \n" + \
                "SOURCE   3 MOL_ID: 2;                                                           \n" + \
                "SOURCE   4 ORGANISM_SCIENTIFIC: AVIAN SARCOMA VIRUS                             \n"
            )

    def test_keywords_generator(self):
        self.assertEqual(
            gen_keywds(
                ["LYASE", "TRICARBOXYLIC ACID CYCLE",
                    "MITOCHONDRION", "OXIDATIVE", "METABOLISM"],
            ),
            [
                "KEYWDS    LYASE, TRICARBOXYLIC ACID CYCLE, MITOCHONDRION, OXIDATIVE, METABOLISM"
            ],
        )

    def test_write_keywords(self):
        pdb_data = PdbData()
        pdb_data.meta.keywords = [
            "LYASE", "TRICARBOXYLIC ACID CYCLE",
            "MITOCHONDRION", "OXIDATIVE", "METABOLISM",
        ]
        with write_temp_pdb(pdb_data, 0, 1) as content:
            self.assertEqual(
                content,
                "KEYWDS    LYASE, TRICARBOXYLIC ACID CYCLE, MITOCHONDRION, OXIDATIVE, METABOLISM \n"
            )

    def test_nummdl_generator(self):
        self.assertEqual(
            gen_nummdl(20),
            ["NUMMDL    20  "],
        )

    def test_write_nummdl(self):
        pdb_data = PdbData()
        pdb_data.add_model(Model())
        pdb_data.add_model(Model())
        with write_temp_pdb(pdb_data, 0, 1) as content:
            self.assertEqual(
                content,
                "NUMMDL    2                                                                     \n"
            )

    def test_expdta_generator(self):
        self.assertEqual(
            gen_expdta(
                ["NEUTRON DIFFRACTION", "X-RAY DIFFRACTION", "SOLUTION NMR"]),
            [
                "EXPDTA    NEUTRON DIFFRACTION; X-RAY DIFFRACTION; SOLUTION NMR                 ",
            ],
        )

    def test_write_expdta(self):
        pdb_data = PdbData()
        pdb_data.meta.experiment = [
            "NEUTRON DIFFRACTION", "X-RAY DIFFRACTION", "SOLUTION NMR",
        ]
        with write_temp_pdb(pdb_data, 0, 1) as content:
            self.assertEqual(
                content,
                "EXPDTA    NEUTRON DIFFRACTION; X-RAY DIFFRACTION; SOLUTION NMR                  \n"
            )

    def test_mdltyp_generator(self):
        self.assertEqual(
            gen_mdltyp(
                "CA ATOMS ONLY, CHAIN A, B, C, D, E, F, G, H, I, J, K ; P ATOMS ONLY, CHAIN X, Y, Z",
            ),
            [
                "MDLTYP    CA ATOMS ONLY, CHAIN A, B, C, D, E, F, G, H, I, J, K ; P ATOMS ONLY,  ",
                "MDLTYP   2 CHAIN X, Y, Z                                                        ",
            ],
        )

    def test_write_mdltyp(self):
        pdb_data = PdbData()
        pdb_data.meta.model_type = (
            "CA ATOMS ONLY, CHAIN A, B, C, D, E, F, G, H, I, J, K ; P ATOMS ONLY, CHAIN X, Y, Z"
        )
        with write_temp_pdb(pdb_data, 0, 2) as content:
            self.assertEqual(
                content,
                "MDLTYP    CA ATOMS ONLY, CHAIN A, B, C, D, E, F, G, H, I, J, K ; P ATOMS ONLY,  \n" + \
                "MDLTYP   2 CHAIN X, Y, Z                                                        \n"
            )

    def test_author_generator(self):
        self.assertEqual(
            gen_author(
                [
                    "M.B.BERRY", "B.MEADOR",
                    "T.BILDERBACK", "P.LIANG",
                    "M.GLASER", "G.N.PHILLIPS JR.",
                    "T.L.ST. STEVENS",
                ],
            ),
            [
                "AUTHOR    M.B.BERRY, B.MEADOR, T.BILDERBACK, P.LIANG, M.GLASER, G.N.PHILLIPS   ",
                "AUTHOR   2 JR., T.L.ST. STEVENS                                                ",
            ],
        )

    def test_write_author(self):
        pdb_data = PdbData()
        pdb_data.meta.author = [
            "M.B.BERRY", "B.MEADOR",
            "T.BILDERBACK", "P.LIANG",
            "M.GLASER", "G.N.PHILLIPS JR.",
            "T.L.ST. STEVENS",
        ]
        with write_temp_pdb(pdb_data, 0, 2) as content:
            self.assertEqual(
                content,
                "AUTHOR    M.B.BERRY, B.MEADOR, T.BILDERBACK, P.LIANG, M.GLASER, G.N.PHILLIPS    \n" + \
                "AUTHOR   2 JR., T.L.ST. STEVENS                                                 \n"
            )

    def test_revdat_generator(self):
        self.assertEqual(
            gen_revdat(2, "11-MAR-08", "2ABC", 1, ["JRNL", "VERSN"]),
            [
                "REVDAT   2   11-MAR-08 2ABC    1       JRNL   VERSN  ",
            ]
        )

    def test_write_revdat(self):
        pdb_data = PdbData()
        pdb_data.meta.pdb_id = "2ABC"
        pdb_data.meta.revisions = [
            RevisionInfo(date="09-DEC-03", modifications=[]),
            RevisionInfo(date="11-MAR-08", modifications=["JRNL", "VERSN"]),
        ]
        with write_temp_pdb(pdb_data, 1, 3) as content:
            self.assertEqual(
                content,
                "REVDAT   2   11-MAR-08 2ABC    1       JRNL   VERSN                             \n" + \
                "REVDAT   1   09-DEC-03 2ABC    0                                                \n"
            )

    def test_jrnl_auth_generator(self):
        self.assertEqual(
            gen_jrnl_auth(["G.FERMI", "M.F.PERUTZ", "B.SHAANAN", "R.FOURME"]),
            [
                "JRNL        AUTH   G.FERMI, M.F.PERUTZ, B.SHAANAN, R.FOURME                    ",
            ]
        )

    def test_jrnl_titl_generator(self):
        self.assertEqual(
            gen_jrnl_titl("THE CRYSTAL STRUCTURE OF HUMAN DEOXYHAEMOGLOBIN AT 1.74 A RESOLUTION"),
            [
                "JRNL        TITL   THE CRYSTAL STRUCTURE OF HUMAN DEOXYHAEMOGLOBIN AT 1.74 A   ",
                "JRNL        TITL 2 RESOLUTION                                                  ",
            ]
        )

    def test_jrnl_editor_generator(self):
        self.assertEqual(
            gen_jrnl_editor(["G.FERMI", "M.F.PERUTZ", "B.SHAANAN", "R.FOURME"]),
            [
                "JRNL        EDIT   G.FERMI, M.F.PERUTZ, B.SHAANAN, R.FOURME                    ",
            ]
        )

    def test_jrnl_ref_to_be_published_generator(self):
        self.assertEqual(
            gen_jrnl_ref_to_be_published(),
            ["JRNL        REF    TO BE PUBLISHED"],
        )

    def test_jrnl_ref_generator(self):
        self.assertEqual(
            gen_jrnl_ref("J.MOL.BIOL.", "175", "159", "1984"),
            [
                "JRNL        REF    J.MOL.BIOL.                   V. 175   159 1984",
            ]
        )

    def test_jrnl_publisher_generator(self):
        self.assertEqual(
            gen_jrnl_publisher("WILEY-VCH VERLAG WEINHEIM"),
            [
                "JRNL        PUBL   WILEY-VCH VERLAG WEINHEIM                          ",
            ]
        )

    def test_jrnl_refn_generator(self):
        self.assertEqual(
            gen_jrnl_refn("ISSN", "1234-5678"),
            [
                "JRNL        REFN                   ISSN 1234-5678                ",
            ]
        )

    def test_jrnl_pmid_generator(self):
        self.assertEqual(
            gen_jrnl_pmid("12345678"),
            [
                "JRNL        PMID   12345678                                                    ",
            ]
        )

    def test_jrnl_doi_generator(self):
        self.assertEqual(
            gen_jrnl_doi("10.1000/12345678"),
            [
                "JRNL        DOI    10.1000/12345678                                            ",
            ]
        )

    def test_write_journal(self):
        pdb_data = PdbData()
        pdb_data.meta.journal = JournalInfo(
            author=["G.FERMI", "M.F.PERUTZ", "B.SHAANAN", "R.FOURME"],
            title="THE CRYSTAL STRUCTURE OF HUMAN DEOXYHAEMOGLOBIN AT 1.74 A RESOLUTION",
            journal="J.MOL.BIOL.",
            volume="175",
            pages="159",
            year="1984",
            publisher="LONDON : ACADEMIC PRESS",
            issn="0022-2836",
            pmid="6726807",
            doi="10.1016/0022-2836(84)90472-8",
        )
        with write_temp_pdb(pdb_data, 0, 8) as content:
            self.assertEqual(
                content,
                "JRNL        AUTH   G.FERMI, M.F.PERUTZ, B.SHAANAN, R.FOURME                     \n" + \
                "JRNL        TITL   THE CRYSTAL STRUCTURE OF HUMAN DEOXYHAEMOGLOBIN AT 1.74 A    \n" + \
                "JRNL        TITL 2 RESOLUTION                                                   \n" + \
                "JRNL        REF    J.MOL.BIOL.                   V. 175   159 1984              \n" + \
                "JRNL        PUBL   LONDON : ACADEMIC PRESS                                      \n" + \
                "JRNL        REFN                   ISSN 0022-2836                               \n" + \
                "JRNL        PMID   6726807                                                      \n" + \
                "JRNL        DOI    10.1016/0022-2836(84)90472-8                                 \n"
            )

    def test_sprsde_generator(self):
        self.assertEqual(
            gen_sprsde("27-FEB-95", "1GDJ", ["1LH4", "2LH4", "3LH4"]),
            [
                "SPRSDE     27-FEB-95 1GDJ      1LH4 2LH4 3LH4 ",
            ]
        )

    def test_write_sprsde(self):
        pdb_data = PdbData()
        pdb_data.meta.pdb_id = "1GDJ"
        pdb_data.meta.replace = ReplaceInfo(
            date="27-FEB-95", ids=["1LH4", "2LH4", "3LH4"])
        with write_temp_pdb(pdb_data, 1, 2) as content:
            self.assertEqual(
                content,
                "SPRSDE     27-FEB-95 1GDJ      1LH4 2LH4 3LH4                                   \n"
            )

    def test_remark_generator(self):
        self.assertEqual(
            gen_remark("0", "THIS IS A COMMENT"),
            ["REMARK   0 THIS IS A COMMENT                                                   "],
        )

    def test_write_remark(self):
        pdb_data = PdbData()
        pdb_data.meta.remark = {"0": "THIS IS A COMMENT"}
        with write_temp_pdb(pdb_data, 0, 1) as content:
            self.assertEqual(
                content,
                "REMARK   0 THIS IS A COMMENT                                                    \n"
            )


class TestPdbPrimaryStructureWriter(unittest.TestCase):

    def test_dbref_generator(self):
        self.assertEqual(
            gen_dbref(
                "2JHQ", "A", 1, "", 226, "", "UNP", "Q9KPK8", "UNG_VIBCH", 1, "", 226, "",
            ),
            [
                "DBREF  2JHQ A    1   226  UNP    Q9KPK8   UNG_VIBCH        1    226 "
            ],
        )
        self.assertEqual(
            gen_dbref(
                "1ABC", "A", 61, "", 322, "", "UNIMES", "UPI000148A153", "46197919", 1534489, "",1537377, ""
            ),
            [
                "DBREF1 1ABC A   61   322  UNIMES               UPI000148A153       ",
                "DBREF2 1ABC A     46197919                      1534489     1537377"
            ],
        )

    def test_write_dbref(self):
        pdb_data = PdbData()
        pdb_data.meta.pdb_id = "2JHQ"
        chain = Chain()
        chain.sequence_info = SequenceInfo(
            sequence_db=SequenceDBInfo(
                init_seq_num=1,
                init_icode="",
                end_seq_num=226,
                end_icode="",
                database="UNP",
                db_accession="Q9KPK8",
                db_id_code="UNG_VIBCH",
                db_init_seq_num=1,
                db_init_icode="",
                db_end_seq_num=226,
                db_end_icode="",
            ),
        )
        pdb_data.add_model(Model())
        pdb_data.models[0].add_chain(chain)
        with write_temp_pdb(pdb_data, 1, 2) as content:
            self.assertEqual(
                content,
                "DBREF  2JHQ A    1   226  UNP    Q9KPK8   UNG_VIBCH        1    226             \n"
            )

    def test_seqadv_generator(self):
        self.assertEqual(
            gen_seqadv(
                "3ABC", "MET", "A", -1, "", "UNP", "P10725", "", "", "EXPRESSION TAG",
            ),
            [
                "SEQADV 3ABC MET A   -1  UNP  P10725              EXPRESSION TAG       ",
            ],
        )

    def test_seqres_generator(self):
        residues = [
            "GLY", "ILE", "VAL", "GLU", "GLN", "CYS", "CYS", "THR", "SER",
            "ILE", "CYS", "SER", "LEU",
            "TYR", "GLN", "LEU", "GLU", "ASN", "TYR", "CYS", "ASN",
        ]
        self.assertEqual(
            gen_seqres("A", residues),
            [
                "SEQRES   1 A   21  GLY ILE VAL GLU GLN CYS CYS THR SER ILE CYS SER LEU ",
                "SEQRES   2 A   21  TYR GLN LEU GLU ASN TYR CYS ASN                     ",
            ],
        )

    def test_write_seqres(self):
        pdb_data = PdbData()
        pdb_data.meta.pdb_id = "1ABC"
        chain = Chain()
        chain.sequence_info = SequenceInfo(
            sequence=[
                "GLY", "ILE", "VAL", "GLU", "GLN", "CYS", "CYS", "THR", "SER",
                "ILE", "CYS", "SER", "LEU",
                "TYR", "GLN", "LEU", "GLU", "ASN", "TYR", "CYS", "ASN",
            ],
        )
        pdb_data.add_model(Model())
        pdb_data.models[0].add_chain(chain)
        with write_temp_pdb(pdb_data, 1, 3) as content:
            self.assertEqual(
                content,
                "SEQRES   1 A   21  GLY ILE VAL GLU GLN CYS CYS THR SER ILE CYS SER LEU          \n" + \
                "SEQRES   2 A   21  TYR GLN LEU GLU ASN TYR CYS ASN                              \n"
            )

    def test_modres_generator(self):
        self.assertEqual(
            gen_modres(
                "2R0L", "ASN", "A", 74, "", "ASN", "GLYCOSYLATION SITE"
            ),
            ["MODRES 2R0L ASN A   74  ASN  GLYCOSYLATION SITE                       "],
        )


class TestPdbHeterogenWriter(unittest.TestCase):

    def test_het_generator(self):
        self.assertEqual(
            gen_het("TRS", "B", 975, "", 8, ""),
            ["HET    TRS  B 975       8                                             "],
        )

    def test_hetnam_generator(self):
        self.assertEqual(
            gen_hetnam(
                "SAD", "BETA-METHYLENE SELENAZOLE-4-CARBOXAMIDE ADENINE"),
            ["HETNAM     SAD BETA-METHYLENE SELENAZOLE-4-CARBOXAMIDE ADENINE        "],
        )

    def test_hetsyn_generator(self):
        self.assertEqual(
            gen_hetsyn("TRS", "TRIS BUFFER;"),
            ["HETSYN     TRS TRIS BUFFER;                                           "],
        )

    def test_formul_generator(self):
        self.assertEqual(
            gen_formul(3, "MG", "2(MG 2+)"),
            ["FORMUL   3   MG    2(MG 2+)                                           "],
        )
        self.assertEqual(
            gen_formul(3, "HOH", "2(H 2O)"),
            ["FORMUL   3  HOH   *2(H 2O)                                            "],
        )

    def test_write_heterogen_records(self):
        pdb_data = PdbData()
        pdb_data.add_model(Model())
        pdb_data.models[0].add_chain(Chain())
        residue = Residue("MG")
        residue.het = True
        residue.add_atom(Atom("MG"))
        pdb_data.models[0].chains[0].add_residue(residue)
        pdb_data.heterogen["MG"] = Heterogen(
            name="MAGNESIUM",
            comment="",
            alias=["MAGNESIUM ION"],
            formula="2(MG 2+)",
        )
        with write_temp_pdb(pdb_data, 0, 4) as content:
            self.assertEqual(
                content,
                "HET    MG   A   1       1                                                       \n" + \
                "HETNAM     MG  MAGNESIUM                                                        \n" + \
                "HETSYN     MG  MAGNESIUM ION                                                    \n" + \
                "FORMUL   1   MG    2(MG 2+)                                                     \n"
            )


class TestPdbSecondaryStructureWriter(unittest.TestCase):

    def test_helix_generator(self):
        self.assertEqual(
            gen_helix(
                1, "HA", "GLY", "A", 86, "", "GLY", "A", 94, "", 1, "", 9,
            ),
            [
                "HELIX    1  HA GLY A   86  GLY A   94  1                                   9",
            ],
        )

    def test_write_helix(self):
        pdb_data = PdbData()
        pdb_data.secondary_structure = SecondaryStructureInfo()
        pdb_data.secondary_structure.helix = {
            "HA": Helix(
                chain_id="A",
                init_seq_num=1,
                init_icode="",
                end_seq_num=2,
                end_icode="",
                helix_class=1,
                comment="",
            ),
        }
        pdb_data.add_model(Model())
        chain = Chain()
        chain.add_residue(Residue("GLY"))
        chain.add_residue(Residue("GLY"))
        pdb_data.models[0].add_chain(chain)
        with write_temp_pdb(pdb_data, 0, 1) as content:
            self.assertEqual(
                content,
                "HELIX    1  HA GLY A    1  GLY A    2  1                                   2    \n"
            )

    def test_sheet_generator(self):
        self.assertEqual(
            gen_sheet(
                1, "A", 5, "THR", "A", 107, "", "ARG", "A", 110, "", 0,
                "", "", "", None, "", "", "", "", None, "",
            ), 
            [
                "SHEET    1   A 5 THR A 107  ARG A 110  0 ",
            ],
        )
        self.assertEqual(
            gen_sheet(
                2, "A", 2, "ILE", "A", 96, "", "THR", "A", 99, "", -1,
                "N", "LYS", "A", 98, "", "O", "THR", "A", 107, "",
            ),
            [
                "SHEET    2   A 2 ILE A  96  THR A  99 -1  N  LYS A  98   O  THR A 107 ",
            ],
        )

    def test_write_sheet(self):
        pdb_data = PdbData()
        pdb_data.secondary_structure = SecondaryStructureInfo()
        pdb_data.secondary_structure.sheet = {
            "A": [
                SheetStrand(
                    init_chain_id="A",
                    init_seq_num=1,
                    init_icode="",
                    end_chain_id="A",
                    end_res_seq=2,
                    end_icode="",
                    sense=0,
                    cur_atom="",
                    cur_chain_id="",
                    cur_res_seq=None,
                    cur_icode="",
                    prev_atom="",
                    prev_chain_id="",
                    prev_res_seq=None,
                    prev_icode="",
                ),
                SheetStrand(
                    init_chain_id="A",
                    init_seq_num=1,
                    init_icode="",
                    end_chain_id="A",
                    end_res_seq=2,
                    end_icode="",
                    sense=-1,
                    cur_atom="N",
                    cur_chain_id="A",
                    cur_res_seq=1,
                    cur_icode="",
                    prev_atom="O",
                    prev_chain_id="A",
                    prev_res_seq=2,
                    prev_icode="",
                ),
            ],
        }
        pdb_data.add_model(Model())
        chain = Chain()
        chain.add_residue(Residue("THR"))
        chain.add_residue(Residue("ARG"))
        pdb_data.models[0].add_chain(chain)
        with write_temp_pdb(pdb_data, 0, 2) as content:
            self.assertEqual(
                content,
                "SHEET    1   A 2 THR A   1  ARG A   2  0                                        \n" + \
                "SHEET    2   A 2 THR A   1  ARG A   2 -1  N  THR A   1   O  ARG A   2           \n"
            )


class TestPdbConnectivityWriter(unittest.TestCase):

    def test_ssbond_generator(self):
        self.assertEqual(
            gen_ssbond(1, "A", 26, "", "A", 84, "", "1555", "1555", 2.04),
            [
                "SSBOND   1 CYS A   26    CYS A   84                          1555   1555  2.04",
            ],
        )

    def test_write_ssbond(self):
        pdb_data = PdbData()
        pdb_data.connectivity.ss_bond = [
            SsBond(
                chain_id_1="A",
                seq_num_1=26,
                icode_1="",
                chain_id_2="A",
                seq_num_2=84,
                icode_2="",
                symmetry_operation_1="1555",
                symmetry_operation_2="1555",
                distance=2.04,
            ),
        ]
        with write_temp_pdb(pdb_data, 0, 1) as content:
            self.assertEqual(
                content,
                "SSBOND   1 CYS A   26    CYS A   84                          1555   1555  2.04  \n"
            )

    def test_link_generator(self):
        self.assertEqual(
            gen_link(
                "O", "", "SER", "A", 82, "", "OG", "", "SER", "B", 5, "", "1555", "1555", 2.84,
            ),
            [
                "LINK         O   SER A  82                 OG  SER B   5     1555   1555  2.84",
            ],
        )

    def test_write_link(self):
        pdb_data = PdbData()
        pdb_data.connectivity.link = [
            Link(
                name_1="O",
                alt_loc_1="",
                chain_id_1="A",
                seq_num_1=1,
                icode_1="",
                name_2="OG",
                alt_loc_2="",
                chain_id_2="B",
                seq_num_2=2,
                icode_2="",
                distance=2.84,
            ),
        ]
        pdb_data.add_model(Model())
        chain = Chain()
        chain.add_residue(Residue("SER"))
        chain.add_residue(Residue("SER"))
        pdb_data.models[0].add_chain(chain)
        with write_temp_pdb(pdb_data, 0, 1) as content:
            self.assertEqual(
                content,
                "LINK         O   SER A   1                 OG  SER B   2     1555   1555  2.84  \n"
            )

    def test_cispep_generator(self):
        self.assertEqual(
            gen_cispep(1, "SER", "A", 58, "", "GLY", "A", 59, "", 0, 20.91),
            [
                "CISPEP   1 SER A   58    GLY A   59          0        20.91",
            ],
        )

    def test_write_cispep(self):
        pdb_data = PdbData()
        pdb_data.connectivity.cis_peptide = [
            CisPeptide(
                chain_id_1="A",
                seq_num_1=1,
                icode_1="",
                chain_id_2="A",
                seq_num_2=2,
                icode_2="",
                num_model=0,
                measure=-0.24,
            ),
        ]
        pdb_data.add_model(Model())
        chain = Chain()
        chain.add_residue(Residue("ASP"))
        chain.add_residue(Residue("PRO"))
        pdb_data.models[0].add_chain(chain)
        with write_temp_pdb(pdb_data, 0, 1) as content:
            self.assertEqual(
                content,
                "CISPEP   1 ASP A    1    PRO A    2          0        -0.24                     \n"
            )


class TestPdbMiscFeatureWriter(unittest.TestCase):

    def test_site_generator(self):
        residues = [
            {"res_name": "ASN", "chain_id": "A", "seq_num": 62, "icode": ""},
            {"res_name": "GLY", "chain_id": "A", "seq_num": 63, "icode": ""},
            {"res_name": "HIS", "chain_id": "A", "seq_num": 64, "icode": ""},
            {"res_name": "HOH", "chain_id": "A", "seq_num": 328, "icode": ""},
            {"res_name": "HOH", "chain_id": "A", "seq_num": 634, "icode": ""},
        ]
        self.assertEqual(
            gen_site("AC2", residues),
            [
                "SITE     1 AC2  5 ASN A  62  GLY A  63  HIS A  64  HOH A 328  ",
                "SITE     2 AC2  5 HOH A 634                                   ",
            ],
        )

    def test_write_site(self):
        pdb_data = PdbData()
        pdb_data.add_model(Model())
        chain = Chain()
        chain.add_residue(Residue("ASN"))
        pdb_data.models[0].add_chain(chain)
        pdb_data.sites["AC2"] = [
            Site(chain_id="A", seq_num=1, icode="")
        ]
        with write_temp_pdb(pdb_data, 0, 1) as content:
            self.assertEqual(
                content,
                "SITE     1 AC2  1 ASN A   1                                                     \n"
            )


class TestPdbCrystallographicWriter(unittest.TestCase):

    def test_cryst1_generator(self):
        self.assertEqual(
            gen_cryst1(52.000, 58.600, 61.900, 90.00,
                       90.00, 90.00, "P 21 21 21", 8),
            ["CRYST1   52.000   58.600   61.900  90.00  90.00  90.00 P 21 21 21    8"],
        )

    def test_write_cryst1(self):
        pdb_data = PdbData()
        pdb_data.crystallographic = CrystalInfo(
            cell_lengths=[52.0, 58.6, 61.9],
            cell_angles=[90.0, 90.0, 90.0],
            space_group="P 21 21 21",
            z=8,
        )
        with write_temp_pdb(pdb_data, 0, 1) as content:
            self.assertEqual(
                content,
                "CRYST1   52.000   58.600   61.900  90.00  90.00  90.00 P 21 21 21    8          \n"
            )

    def test_origin_matrix_generator(self):
        self.assertEqual(
            gen_matrix_rows(
                "ORIGX",
                [
                    [0.963457, 0.136613, 0.230424, 16.61000],
                    [-0.158977, 0.983924, 0.081383, 13.72000],
                    [-0.215598, -0.115048, 0.969683, 37.00000],
                ],
            ),
            [
                "ORIGX1      0.963457  0.136613  0.230424       16.61000",
                "ORIGX2     -0.158977  0.983924  0.081383       13.72000",
                "ORIGX3     -0.215598 -0.115048  0.969683       37.00000",
            ],
        )

    def test_write_origin_matrix(self):
        pdb_data = PdbData()
        pdb_data.crystallographic.origin_matrix = [
            [0.963457, 0.136613, 0.230424, 16.61000],
            [-0.158977, 0.983924, 0.081383, 13.72000],
            [-0.215598, -0.115048, 0.969683, 37.00000],
        ]
        with write_temp_pdb(pdb_data, 0, 3) as content:
            self.assertEqual(
                content,
                "ORIGX1      0.963457  0.136613  0.230424       16.61000                         \n" + \
                "ORIGX2     -0.158977  0.983924  0.081383       13.72000                         \n" + \
                "ORIGX3     -0.215598 -0.115048  0.969683       37.00000                         \n"
            )

    def test_scale_matrix_generator(self):
        self.assertEqual(
            gen_matrix_rows(
                "SCALE",
                [
                    [0.019231, 0.0, 0.0, 0.0],
                    [0.0, 0.017065, 0.0, 0.0],
                    [0.0, 0.0, 0.016155, 0.0],
                ],
            ),
            [
                "SCALE1      0.019231  0.000000  0.000000        0.00000",
                "SCALE2      0.000000  0.017065  0.000000        0.00000",
                "SCALE3      0.000000  0.000000  0.016155        0.00000",
            ],
        )

    def test_write_scale_matrix(self):
        pdb_data = PdbData()
        pdb_data.crystallographic.scale_matrix = [
            [0.019231, 0.0, 0.0, 0.0],
            [0.0, 0.017065, 0.0, 0.0],
            [0.0, 0.0, 0.016155, 0.0],
        ]
        with write_temp_pdb(pdb_data, 0, 3) as content:
            self.assertEqual(
                content,
                "SCALE1      0.019231  0.000000  0.000000        0.00000                         \n" + \
                "SCALE2      0.000000  0.017065  0.000000        0.00000                         \n" + \
                "SCALE3      0.000000  0.000000  0.016155        0.00000                         \n"
            )

    def test_ncs_matrix_generator(self):
        matrix=[
            [-1.000000, 0.000000, 0.000000, 0.000000],
            [0.000000, 1.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, -1.000000, 0.000000],
        ]
        self.assertEqual(
            gen_ncs_matrix(1, matrix, True),
            [
                "MTRIX1   1 -1.000000  0.000000  0.000000        0.00000    1",
                "MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1",
                "MTRIX3   1  0.000000  0.000000 -1.000000        0.00000    1",
            ],
        )

    def test_write_ncs_matrix(self):
        pdb_data = PdbData()
        matrix = [
            [-1.000000, 0.000000, 0.000000, 0.000000],
            [0.000000, 1.000000, 0.000000, 0.000000],
            [0.000000, 0.000000, -1.000000, 0.000000],
        ]
        pdb_data.crystallographic.ncs_matrix = [
            NcsMatrix(matrix=matrix, given=True),
            NcsMatrix(
                matrix=[
                    [1.000000, 0.000000, 0.000000, 0.000000],
                    [0.000000, 1.000000, 0.000000, 0.000000],
                    [0.000000, 0.000000, 1.000000, 0.000000],
                ],
                given=False,
            ),
        ]
        with write_temp_pdb(pdb_data, 0, 6) as content:
            self.assertEqual(
                content,
                "MTRIX1   1 -1.000000  0.000000  0.000000        0.00000    1                    \n" + \
                "MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1                    \n" + \
                "MTRIX3   1  0.000000  0.000000 -1.000000        0.00000    1                    \n" + \
                "MTRIX1   2  1.000000  0.000000  0.000000        0.00000    0                    \n" + \
                "MTRIX2   2  0.000000  1.000000  0.000000        0.00000    0                    \n" + \
                "MTRIX3   2  0.000000  0.000000  1.000000        0.00000    0                    \n"
            )


class TestPdbCoordinateWriter(unittest.TestCase):
    def test_atom_generator(self):
        self.assertEqual(
            gen_atom(
                "ATOM", 1, "N", "", "ALA", "A", 1, "", 11.104, 6.134, -6.504, 1.0, 0.0, "N", 0
            ),
            [
                "ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N  ",
            ],
        )

    def test_anisou_generator(self):
        u11 = 0.2406
        u22 = 0.1892
        u33 = 0.1614
        u12 = 0.0198
        u13 = 0.0519
        u23 = -0.0328
        self.assertEqual(
            gen_anisou(
                107, "N", "", "GLY", "A", 13, "", u11, u22, u33, u12, u13, u23, "N", 0
            ), 
            [
                "ANISOU  107  N   GLY A  13     2406   1892   1613    198    519   -328       N  "
            ],
        )

    def test_ter_generator(self):
        self.assertEqual(
            gen_ter(1, "ALA", "A", 1),
            [
                "TER       1      ALA A   1",
            ],
        )

    def test_write_atom_coordinate(self):
        pdb_data = PdbData()
        pdb_data.add_model(Model())
        chain = Chain()
        residue = Residue("MET")
        residue.add_atom(
            Atom(
                name="N",
                coord=(11.104, 13.207, 10.000),
                occupancy=1.0,
                element="N",
                charge=0,
                temp_factor=[0.2406, 0.1892, 0.1613, 0.0198, 0.0519, -0.0328],
            )
        )
        chain.add_residue(residue)
        pdb_data.models[0].add_chain(chain)
        with write_temp_pdb(pdb_data, 0, 2) as content:
            print(content)
            self.assertEqual(
                content,
                "ATOM      1  N   MET A   1      11.104  13.207  10.000  1.00 15.56           N  \n" + \
                "ANISOU    1  N   MET A   1     2406   1892   1613    198    519   -328       N  \n"
            )


class TestPdbConectWriter(unittest.TestCase):

    def test_conect_generator(self):
        self.assertEqual(
            gen_conect(1, [2, 1]),
            ["CONECT    1    2    1"],
        )

    def test_write_conect(self):
        pdb_data = PdbData()
        pdb_data.connectivity.connections = {1: [2, 3]}
        with write_temp_pdb(pdb_data, 0, 3) as content:
            self.assertEqual(
                content,
                "CONECT    1    2    3                                                           \n" + \
                "CONECT    2    1                                                                \n" + \
                "CONECT    3    1                                                                \n"
            )


class TestPdbBookkeepingWriter(unittest.TestCase):

    def test_master_generator(self):
        self.assertEqual(
            gen_master(40, 0, 0, 0, 0, 6, 2930, 2, 0, 29),
            ["MASTER       40    0    0    0    0    0    0    6 2930    2    0   29"],
        )

    def test_write_master(self):
        pdb_data = PdbData()
        pdb_data._validation_info = {
            "num_remark": 40,
            "num_het": 0,
            "num_helix": 0,
            "num_sheet": 0,
            "num_site": 0,
        }
        with write_temp_pdb(pdb_data, 0, 1) as content:
            self.assertEqual(
                content,
                "MASTER       40    0    0    0    0    0    0    0    0    0    0    0          \n"
            )

if __name__ == "__main__":
    unittest.main()
