"""Unit tests for the PDB Writer.
"""
# pylint: disable=missing-function-docstring
# pylint: disable=missing-class-docstring
# pylint: disable=unused-wildcard-import
# pylint: disable=wildcard-import
import unittest

from pypdbio import set_unit
from pypdbio.models import NcsMatrix
from pypdbio.writer_helper import *

set_unit("A")


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
    def test_header_field_generator(self):
        self.assertEqual(
            gen_header("TRANSFERASE/TRANSFERASE INHIBITOR",
                       "17-SEP-04", "1XH6"),
            ["HEADER    TRANSFERASE/TRANSFERASE INHIBITOR       17-SEP-04   1XH6"],
        )

    def test_title_field_generator(self):
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

    def test_obslte_field_generator(self):
        self.assertEqual(
            gen_obslte("31-JAN-94", "2MBP", "1MBP"),
            ["OBSLTE     31-JAN-94 1MBP      2MBP"],
        )

    def test_split_field_generator(self):
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

    def test_caveat_field_generator(self):
        self.assertEqual(
            gen_caveat(
                "2WMW",
                ["CHAIN A IS MISSING RESIDUES AND CHAIN B IS MISSING RESIDUES"],
            ),
            [
                "CAVEAT     2WMW    CHAIN A IS MISSING RESIDUES AND CHAIN B IS MISSING RESIDUES ",
            ],
        )

    def test_compound_field_generator(self):
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

    def test_source_field_generator(self):
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

    def test_keywords_field_generator(self):
        self.assertEqual(
            gen_keywds(
                ["LYASE", "TRICARBOXYLIC ACID CYCLE",
                    "MITOCHONDRION", "OXIDATIVE", "METABOLISM"],
            ),
            [
                "KEYWDS    LYASE, TRICARBOXYLIC ACID CYCLE, MITOCHONDRION, OXIDATIVE, METABOLISM"
            ],
        )

    def test_nummdl_field_generator(self):
        self.assertEqual(
            gen_nummdl(20),
            ["NUMMDL    20  "],
        )

    def test_expdta_field_generator(self):
        self.assertEqual(
            gen_expdta(
                ["NEUTRON DIFFRACTION", "X-RAY DIFFRACTION", "SOLUTION NMR"]),
            [
                "EXPDTA    NEUTRON DIFFRACTION; X-RAY DIFFRACTION; SOLUTION NMR                 ",
            ],
        )

    def test_mdltyp_field_generator(self):
        self.assertEqual(
            gen_mdltyp(
                "CA ATOMS ONLY, CHAIN A, B, C, D, E, F, G, H, I, J, K ; P ATOMS ONLY, CHAIN X, Y, Z",
            ),
            [
                "MDLTYP    CA ATOMS ONLY, CHAIN A, B, C, D, E, F, G, H, I, J, K ; P ATOMS ONLY,  ",
                "MDLTYP   2 CHAIN X, Y, Z                                                        ",
            ],
        )

    def test_author_field_generator(self):
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

    def test_revdat_field_generator(self):
        self.assertEqual(
            gen_revdat(2, "11-MAR-08", "2ABC", 1, ["JRNL", "VERSN"]),
            [
                "REVDAT   2   11-MAR-08 2ABC    1       JRNL   VERSN  ",
            ]
        )

    def test_jrnl_auth_field_generator(self):
        self.assertEqual(
            gen_jrnl_auth(["G.FERMI", "M.F.PERUTZ", "B.SHAANAN", "R.FOURME"]),
            [
                "JRNL        AUTH   G.FERMI, M.F.PERUTZ, B.SHAANAN, R.FOURME                    ",
            ]
        )

    def test_jrnl_titl_field_generator(self):
        self.assertEqual(
            gen_jrnl_titl("THE CRYSTAL STRUCTURE OF HUMAN DEOXYHAEMOGLOBIN AT 1.74 A RESOLUTION"),
            [
                "JRNL        TITL   THE CRYSTAL STRUCTURE OF HUMAN DEOXYHAEMOGLOBIN AT 1.74 A   ",
                "JRNL        TITL 2 RESOLUTION                                                  ",
            ]
        )

    def test_jrnl_editor_field_generator(self):
        self.assertEqual(
            gen_jrnl_editor(["G.FERMI", "M.F.PERUTZ", "B.SHAANAN", "R.FOURME"]),
            [
                "JRNL        EDIT   G.FERMI, M.F.PERUTZ, B.SHAANAN, R.FOURME                    ",
            ]
        )

    def test_jrnl_ref_to_be_published_field_generator(self):
        self.assertEqual(
            gen_jrnl_ref_to_be_published(),
            ["JRNL        REF    TO BE PUBLISHED"],
        )

    def test_jrnl_ref_field_generator(self):
        self.assertEqual(
            gen_jrnl_ref("J.MOL.BIOL.", "175", "159", "1984"),
            [
                "JRNL        REF    J.MOL.BIOL.                   V. 175   159 1984",
            ]
        )

    def test_jrnl_publisher_field_generator(self):
        self.assertEqual(
            gen_jrnl_publisher("WILEY-VCH VERLAG WEINHEIM"),
            [
                "JRNL        PUBL   WILEY-VCH VERLAG WEINHEIM                          ",
            ]
        )

    def test_jrnl_refn_field_generator(self):
        self.assertEqual(
            gen_jrnl_refn("ISSN", "1234-5678"),
            [
                "JRNL        REFN                   ISSN 1234-5678                ",
            ]
        )

    def test_jrnl_pmid_field_generator(self):
        self.assertEqual(
            gen_jrnl_pmid("12345678"),
            [
                "JRNL        PMID   12345678                                                    ",
            ]
        )

    def test_jrnl_doi_field_generator(self):
        self.assertEqual(
            gen_jrnl_doi("10.1000/12345678"),
            [
                "JRNL        DOI    10.1000/12345678                                            ",
            ]
        )

    def test_sprsde_field_generator(self):
        self.assertEqual(
            gen_sprsde("27-FEB-95", "1GDJ", ["1LH4", "2LH4", "3LH4"]),
            [
                "SPRSDE     27-FEB-95 1GDJ      1LH4 2LH4 3LH4 ",
            ]
        )

    def test_remark_field_generator(self):
        self.assertEqual(
            gen_remark("0", "THIS IS A COMMENT"),
            ["REMARK   0 THIS IS A COMMENT                                                   "],
        )


class TestPdbPrimaryStructureWriter(unittest.TestCase):

    def test_dbref_field_generator(self):
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

    def test_seqadv_field_generator(self):
        self.assertEqual(
            gen_seqadv(
                "3ABC", "MET", "A", -1, "", "UNP", "P10725", "", "", "EXPRESSION TAG",
            ),
            [
                "SEQADV 3ABC MET A   -1  UNP  P10725              EXPRESSION TAG       ",
            ],
        )

    def test_seqres_field_generator(self):
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

    def test_modres_field_generator(self):
        self.assertEqual(
            gen_modres(
                "2R0L", "ASN", "A", 74, "", "ASN", "GLYCOSYLATION SITE"
            ),
            ["MODRES 2R0L ASN A   74  ASN  GLYCOSYLATION SITE                       "],
        )


class TestPdbHeterogenWriter(unittest.TestCase):

    def test_het_field_generator(self):
        self.assertEqual(
            gen_het("TRS", "B", 975, "", 8, ""),
            ["HET    TRS  B 975       8"],
        )

    def test_hetnam_field_generator(self):
        self.assertEqual(
            gen_hetnam(
                "SAD", "BETA-METHYLENE SELENAZOLE-4-CARBOXAMIDE ADENINE"),
            ["HETNAM     SAD BETA-METHYLENE SELENAZOLE-4-CARBOXAMIDE ADENINE               "],
        )

    def test_hetsyn_field_generator(self):
        self.assertEqual(
            gen_hetsyn("TRS", "TRIS BUFFER;"),
            ["HETSYN     TRS TRIS BUFFER;                                                   "],
        )

    def test_formul_field_generator(self):
        self.assertEqual(
            gen_formul(3, "MG", 0, "", "2(MG 2+)"),
            ["FORMUL   3   MG    2(MG 2+)"],
        )

    def test_het_integration_field_generator(self):
        self.assertEqual(
            gen_het("MG", "B", 975, "", 8, "")
            + gen_hetnam("MG", "MAGNESIUM")
            + gen_hetsyn("MG", "MAGNESIUM ION;")
            + gen_formul(3, "MG", 0, "", "2(MG 2+)"),
            [
                "HET    MG   B 975       8",
                "HETNAM     MG  MAGNESIUM                                                      ",
                "HETSYN     MG  MAGNESIUM ION;                                                 ",
                "FORMUL   3   MG    2(MG 2+)",
            ],
        )


class TestPdbSecondaryStructureWriter(unittest.TestCase):

    def test_helix_field_generator(self):
        self.assertEqual(
            gen_helix(
                1, "HA", "GLY", "A", 86, "", "GLY", "A", 94, "", 1, "", 9,
            ),
            [
                "HELIX    1  HA GLY A   86  GLY A   94  1                                   9 ",
            ],
        )

    def test_sheet_field_generator(self):
        self.assertEqual(
            gen_sheet(
                1, "A", 2, "THR", "A", 107, "", "ARG", "A", 110, "", 0,
                "", "", "", None, "", "", "", "", None, "",
            )
            + gen_sheet(
                2, "A", 2, "ILE", "A", 96, "", "THR", "A", 99, "", -1,
                "N", "LYS", "A", 98, "", "O", "THR", "A", 107, "",
            ),
            [
                "SHEET    1   A 2 THR A 107  ARG A 110  0",
                "SHEET    2   A 2 ILE A  96  THR A  99 -1  N  LYS A  98   O  THR A 107",
            ],
        )


class TestPdbConnectivityWriter(unittest.TestCase):

    def test_ssbond_field_generator(self):
        self.assertEqual(
            gen_ssbond(1, "A", 26, "", "A", 84, "", "1555", "1555", 2.04),
            [
                "SSBOND   1 CYS A   26    CYS A   84                          1555   1555  2.04",
            ],
        )

    def test_link_field_generator(self):
        self.assertEqual(
            gen_link(
                "O", "", "SER", "A", 82, "", "OG", "", "SER", "B", 5, "", "1555", "1555", 2.84,
            ),
            [
                "LINK         O   SER A  82                 OG  SER B   5     1555   1555  2.84",
            ],
        )

    def test_cispep_field_generator(self):
        self.assertEqual(
            gen_cispep(1, "ASP", "A", 118, "", "PRO", "A", 119, "", 0, -0.24),
            [
                "CISPEP   1 ASP A 118     PRO A 119          0        -0.24",
            ],
        )


class TestPdbMiscFeatureWriter(unittest.TestCase):

    def test_site_field_generator(self):
        residues = [
            {"res_name": "HIS", "chain_id": "A", "seq_num": 94, "icode": ""},
            {"res_name": "HIS", "chain_id": "A", "seq_num": 96, "icode": ""},
            {"res_name": "HIS", "chain_id": "A", "seq_num": 119, "icode": ""},
        ]
        self.assertEqual(
            gen_site("AC1", residues),
            [
                "SITE   1 AC1  3 HIS A  94 HIS A  96 HIS A 119                                 ",
            ],
        )

    def test_site_integration_field_generator(self):
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
                "SITE   1 AC2  5 ASN A  62 GLY A  63 HIS A  64 HOH A 328                     ",
                "SITE   2 AC2  5 HOH A 634                                                     ",
            ],
        )


class TestPdbCrystallographicWriter(unittest.TestCase):

    def test_cryst1_field_generator(self):
        self.assertEqual(
            gen_cryst1(52.000, 58.600, 61.900, 90.00,
                       90.00, 90.00, "P 21 21 21", 8),
            ["CRYST1   52.000   58.600   61.900  90.00  90.00  90.00 P 21 21 21    8"],
        )

    def test_origin_matrix_field_generator(self):
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
                "ORIGX1      0.963457  0.136613  0.230424      16.61000",
                "ORIGX2     -0.158977  0.983924  0.081383      13.72000",
                "ORIGX3     -0.215598 -0.115048  0.969683      37.00000",
            ],
        )

    def test_scale_matrix_field_generator(self):
        self.assertEqual(
            gen_scale_matrix_rows(
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

    def test_ncs_matrix_field_generator(self):
        ncs = [
            NcsMatrix(
                matrix=[
                    [-1.000000, 0.000000, 0.000000, 0.000000],
                    [0.000000, 1.000000, 0.000000, 0.000000],
                    [0.000000, 0.000000, -1.000000, 0.000000],
                ],
                given=True,
            ),
            NcsMatrix(
                matrix=[
                    [1.000000, 0.000000, 0.000000, 0.000000],
                    [0.000000, 1.000000, 0.000000, 0.000000],
                    [0.000000, 0.000000, 1.000000, 0.000000],
                ],
                given=False,
            ),
        ]
        self.assertEqual(
            gen_ncs_matrix(ncs),
            [
                "MTRIX1   1 -1.000000  0.000000  0.000000        0.00000    1",
                "MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1",
                "MTRIX3   1  0.000000  0.000000 -1.000000        0.00000    1",
                "MTRIX1   2  1.000000  0.000000  0.000000        0.00000    0",
                "MTRIX2   2  0.000000  1.000000  0.000000        0.00000    0",
                "MTRIX3   2  0.000000  0.000000  1.000000        0.00000    0",
            ],
        )


class TestPdbCoordinateWriter(unittest.TestCase):
    """对应 TestPdbCoordinateReader；尚无 gen_atom / gen_anisou，占位以便与 reader 样例对齐。"""

    @unittest.skip("尚无 gen_atom：输入/输出见 TestPdbCoordinateReader.test_atom_parser")
    def test_atom_field(self):
        pass

    @unittest.skip("尚无 gen_anisou：输入/输出见 TestPdbCoordinateReader.test_anisou_parser")
    def test_anisou_field(self):
        pass


class TestPdbConectWriter(unittest.TestCase):

    def test_conect_field(self):
        self.assertEqual(
            gen_conect({1: [2], 2: [1]}),
            ["CONECT    1    2", "CONECT    2    1"],
        )


class TestPdbBookkeepingWriter(unittest.TestCase):

    def test_master_field(self):
        self.assertEqual(
            gen_master({
                "num_remark": 40,
                "num_het": 0,
                "num_helix": 0,
                "num_sheet": 0,
                "num_site": 0,
                "num_xform": 6,
                "num_coord": 2930,
                "num_ter": 2,
                "num_conect": 0,
                "num_seq": 29,
            }),
            ["MASTER       40    0    0    0    0    0    0    6 2930    2    0   29"],
        )


if __name__ == "__main__":
    unittest.main()
