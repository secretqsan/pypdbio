# pylint: disable=wildcard-import
# pylint: disable=protected-access
# pylint: disable=missing-function-docstring
# pylint: disable=missing-class-docstring
# pylint: disable=unused-wildcard-import
"""
Unit tests for the PDB Reader.
"""
import os
import unittest
from contextlib import contextmanager
from pypdbio import PdbReader, set_unit
from pypdbio.models import *
from pypdbio.reader_helper import *


TEMP_PDB_PATH = "test_pdb.pdb"
set_unit("A")

@contextmanager
def write_temp_pdb(content):
    with open(TEMP_PDB_PATH, 'w', encoding='utf-8') as f:
        f.write(content)
    try:
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        yield pdb_data
    finally:
        os.remove(TEMP_PDB_PATH)


class TestPdbReaderHelper(unittest.TestCase):
    def test_continuous_field_parser(self):
        def parser1(line):
            return line
        result1 = "123"
        result1 = parse_continuous_field("123", result1, parser1)
        self.assertEqual(result1, "123 123")

        def parser2(line):
            return [1, line]
        result2 = [1, "123"]
        result2 = parse_continuous_field(
            "123", result2, parser2, place=1)
        self.assertEqual(result2, [1, "123 123"])

class TestPdbHeaderReader(unittest.TestCase):

    def test_header_parser(self):
        self.assertEqual(
            parse_header(
                "HEADER    TRANSFERASE/TRANSFERASE INHIBITOR       17-SEP-04   1XH6              "
            ),
            ("TRANSFERASE/TRANSFERASE INHIBITOR", "17-SEP-04", "1XH6"),
        )

    def test_header(self):
        with write_temp_pdb(
            "HEADER    HYDROLASE                               19-MAY-97   1AKI              "
        ) as pdb_data:
            self.assertEqual(pdb_data.meta.classification, "HYDROLASE")
            self.assertEqual(pdb_data.meta.date, "19-MAY-97")
            self.assertEqual(pdb_data.meta.pdb_id, "1AKI")

    def test_remark_parser(self):
        self.assertEqual(
            parse_remark(
                "REMARK  0  THIS IS A COMMENT"
            ),
            ["0", "THIS IS A COMMENT"],
        )

    def test_remark(self):
        with write_temp_pdb(
            "REMARK  0  THIS IS A COMMENT\n"
        ) as pdb_data:
            self.assertEqual(pdb_data.meta.remark, {
                "0": "THIS IS A COMMENT",
            })
    def test_title_parser(self):
        self.assertEqual(
            parse_title(
                "TITLE     RHIZOPUSPEPSIN COMPLEXED WITH REDUCED PEPTIDE INHIBITOR           \n"
            ),
            "RHIZOPUSPEPSIN COMPLEXED WITH REDUCED PEPTIDE INHIBITOR",
        )

    def test_title(self):
        with write_temp_pdb(
            "TITLE     THE STRUCTURE OF THE ORTHORHOMBIC FORM OF HEN EGG-WHITE LYSOZYME AT   \n"
            "TITLE    2 1.5 ANGSTROMS RESOLUTION                                             \n"
        ) as pdb_data:
            self.assertEqual(
                pdb_data.meta.title,
                "THE STRUCTURE OF THE ORTHORHOMBIC FORM OF HEN EGG-WHITE LYSOZYME AT 1.5 ANGSTROMS RESOLUTION"
            )

    def test_obslte_parser(self):
        self.assertEqual(
            parse_obslte(
                "OBSLTE     31-JAN-94 1MBP      2MBP"
            ),
            ("31-JAN-94", "2MBP"),
        )

    def test_obslte(self):
        with write_temp_pdb(
            "OBSLTE     31-JAN-94 1MBP      2MBP\n"
        ) as pdb_data:
            self.assertEqual(
                pdb_data.meta.obsolete,
                ObsoleteInfo(replace_date="31-JAN-94", new_entry_id="2MBP"),
            )

    def test_split_parser(self):
        self.assertEqual(
            parse_split(
                "SPLIT      1VOQ 1VOR 1VOS 1VOU 1VOV 1VOW 1VOX 1VOY 1VP0 1VOZ"
            ),
            ["1VOQ", "1VOR", "1VOS", "1VOU", "1VOV",
                "1VOW", "1VOX", "1VOY", "1VP0", "1VOZ"],
        )

    def test_split(self):
        with write_temp_pdb(
            "SPLIT      1VOQ 1VOR 1VOS 1VOU 1VOV 1VOW 1VOX 1VOY 1VP0 1VOZ\n"
            "SPLIT      1VPA\n"
        ) as pdb_data:
            self.assertEqual(
                pdb_data.meta.split,
                ["1VOQ", "1VOR", "1VOS", "1VOU", "1VOV", "1VOW",
                    "1VOX", "1VOY", "1VP0", "1VOZ", "1VPA"]
            )

    def test_caveat_parser(self):
        self.assertEqual(
            parse_caveat(
                "CAVEAT    2WMW     CHAIN A IS MISSING RESIDUES"
            ),
            "CHAIN A IS MISSING RESIDUES",
        )

    def test_caveat(self):
        with write_temp_pdb(
            "CAVEAT    2WMW     CHAIN A IS MISSING RESIDUES\n"
        ) as pdb_data:
            self.assertEqual(
                pdb_data.meta.caveat,
                "CHAIN A IS MISSING RESIDUES"
            )

    def test_compond_parser(self):
        self.assertEqual(
            parse_compond(
                "COMPND    MOL_ID:  1; MOLECULE:  HEMOGLOBIN ALPHA CHAIN; CHAIN: A,  C"
            ),
            "MOL_ID:  1; MOLECULE:  HEMOGLOBIN ALPHA CHAIN; CHAIN: A,  C",
        )

    def test_source_parser(self):
        self.assertEqual(
            parse_source(
                "SOURCE    MOL_ID: 1; ORGANISM_SCIENTIFIC: AVIAN SARCOMA VIRUS;"
            ),
            "MOL_ID: 1; ORGANISM_SCIENTIFIC: AVIAN SARCOMA VIRUS;",
        )

    def test_compound(self):
        with write_temp_pdb(
            "COMPND    MOL_ID:  1;\n"
            "COMPND   2 MOLECULE:  HEMOGLOBIN ALPHA CHAIN;\n"
            "COMPND   3 CHAIN: A,  C;\n"
            "COMPND   7 MOL_ID:  2;\n"
            "COMPND   8 MOLECULE:  HEMOGLOBIN BETA CHAIN;\n"
            "SOURCE    MOL_ID: 1;\n"
            "SOURCE   2 ORGANISM_SCIENTIFIC: AVIAN SARCOMA VIRUS;\n"
            "SOURCE   3 MOL_ID: 2;\n"
            "SOURCE   4 ORGANISM_SCIENTIFIC: AVIAN SARCOMA VIRUS;\n"
        ) as pdb_data:
            self.assertEqual(
            pdb_data.meta.compounds,
                [
                    {
                        "molecule": "HEMOGLOBIN ALPHA CHAIN",
                        "chain": "A,  C",
                        "source": {
                            "organism_scientific": "AVIAN SARCOMA VIRUS",
                        },
                    },
                    {
                        "molecule": "HEMOGLOBIN BETA CHAIN",
                        "source": {
                            "organism_scientific": "AVIAN SARCOMA VIRUS",
                        },
                    },
                ]
            )

    def test_compound_warning(self):
        with open(TEMP_PDB_PATH, 'w', encoding='utf-8') as f:
            f.write(
                "COMPND    MOL_ID:  1;\n"
                "COMPND   2 MOL_ID:  2\n"
                "SOURCE    MOL_ID: 1\n"
            )
        reader = PdbReader(TEMP_PDB_PATH)
        with self.assertWarns(Warning):
            reader.read()
        os.remove(TEMP_PDB_PATH)

    def test_keywords_parser(self):
        self.assertEqual(
            parse_keywords(
                "KEYWDS    LYASE, TRICARBOXYLIC ACID CYCLE, MITOCHONDRION"
            ),
            ["LYASE", "TRICARBOXYLIC ACID CYCLE", "MITOCHONDRION"],
        )

    def test_keywords(self):
        with open(TEMP_PDB_PATH, 'w', encoding='utf-8') as f:
            f.write(
                "KEYWDS    LYASE,  TRICARBOXYLIC ACID CYCLE, MITOCHONDRION, OXIDATIVE\n"
                "KEYWDS   2 METABOLISM\n"
            )
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(
            pdb_data.meta.keywords,
            ["LYASE", "TRICARBOXYLIC ACID CYCLE",
                "MITOCHONDRION", "OXIDATIVE", "METABOLISM"]
        )
        os.remove(TEMP_PDB_PATH)

    def test_nummdl_parser(self):
        self.assertEqual(
            parse_nummdl("NUMMDL    20"),
            20,
        )

    def test_expdta_parser(self):
        self.assertEqual(
            parse_expdta(
                "EXPDTA    NEUTRON DIFFRACTION; X-RAY DIFFRACTION"
            ),
            ["NEUTRON DIFFRACTION", "X-RAY DIFFRACTION"],
        )

    def test_expdta(self):
        with open(TEMP_PDB_PATH, 'w', encoding='utf-8') as f:
            f.write(
                "EXPDTA    NEUTRON DIFFRACTION; X-RAY DIFFRACTION\n"
                "EXPDTA    SOLUTION NMR\n"
            )
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(
            pdb_data.meta.experiment,
            ["NEUTRON DIFFRACTION", "X-RAY DIFFRACTION", "SOLUTION NMR"]
        )
        os.remove(TEMP_PDB_PATH)

    def test_mdltyp_parser(self):
        self.assertEqual(
            parse_mdltyp("MDLTYP   2 CHAIN X, Y, Z"),
            "CHAIN X, Y, Z"
        )

    def test_mdltyp(self):
        with open(TEMP_PDB_PATH, 'w', encoding='utf-8') as f:
            f.write(
                "MDLTYP    CA ATOMS ONLY, CHAIN A, B, C, D, E, F, G, H, I, J, K ; P ATOMS ONLY,\n"
                "MDLTYP   2 CHAIN X, Y, Z"
            )
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(
            pdb_data.meta.model_type,
            "CA ATOMS ONLY, CHAIN A, B, C, D, E, F, G, H, I, J, K ; P ATOMS ONLY, CHAIN X, Y, Z"
        )
        os.remove(TEMP_PDB_PATH)

    def test_author_parser(self):
        self.assertEqual(
            parse_author(
                "AUTHOR    FIRST AUTHOR, LAST AUTHOR"
            ),
            ["FIRST AUTHOR", "LAST AUTHOR"],
        )

    def test_author(self):
        with open(TEMP_PDB_PATH, 'w', encoding='utf-8') as f:
            f.write(
                "AUTHOR    M.B.BERRY,B.MEADOR,T.BILDERBACK,P.LIANG,M.GLASER,\n"
                "AUTHOR   2 G.N.PHILLIPS JR.,T.L.ST. STEVENS\n"
            )
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(
            pdb_data.meta.author,
            [
                "M.B.BERRY", "B.MEADOR",
                "T.BILDERBACK", "P.LIANG",
                "M.GLASER", "G.N.PHILLIPS JR.",
                "T.L.ST. STEVENS"
            ]
        )
        os.remove(TEMP_PDB_PATH)

    def test_jrnl_auth_parser(self):
        self.assertEqual(
            parse_jrnl_auth(
                "JRNL        AUTH   G.FERMI,M.F.PERUTZ,B.SHAANAN,R.FOURME"
            ),
            ["G.FERMI", "M.F.PERUTZ", "B.SHAANAN", "R.FOURME"],
        )

    def test_jrnl_titl_parser(self):
        self.assertEqual(
            parse_jrnl_titl(
                "JRNL        TITL   THE CRYSTAL STRUCTURE OF HUMAN DEOXYHAEMOGLOBIN AT"
            ),
            "THE CRYSTAL STRUCTURE OF HUMAN DEOXYHAEMOGLOBIN AT",
        )

    def test_jrnl_edit_parser(self):
        self.assertEqual(
            parse_jrnl_edit(
                "JRNL        EDIT   J.DOE,A.SMITH"
            ),
            "J.DOE,A.SMITH",
        )

    def test_jrnl_ref_parser(self):
        self.assertEqual(
            parse_jrnl_ref(
                "JRNL        REF    J.MOL.BIOL.                   V. 175   159 1984"
            ),
            ("J.MOL.BIOL.", "175", "159", "1984"),
        )

    def test_jrnl_publ_parser(self):
        self.assertEqual(
            parse_jrnl_publ(
                "JRNL        PUBL   LONDON : ACADEMIC PRESS"
            ),
            "LONDON : ACADEMIC PRESS",
        )

    def test_jrnl_refn_parser(self):
        self.assertEqual(
            parse_jrnl_refn(
                "JRNL        REFN                   ISSN 0022-2836"
            ),
            ("ISSN", "0022-2836"),
        )

    def test_jrnl_pmid_parser(self):
        self.assertEqual(
            parse_jrnl_pmid(
                "JRNL        PMID   6726807"
            ),
            "6726807",
        )

    def test_jrnl_doi_parser(self):
        self.assertEqual(
            parse_jrnl_doi(
                "JRNL        DOI    10.1016/0022-2836(84)90472-8"
            ),
            "10.1016/0022-2836(84)90472-8",
        )

    def test_sprsde_parser(self):
        self.assertEqual(
            parse_sprsde(
                "SPRSDE     27-FEB-95 1GDJ      1LH4 2LH4"
            ),
            ("27-FEB-95", ["1LH4", "2LH4"]),
        )

    def test_jrnl(self):
        with open(TEMP_PDB_PATH, 'w', encoding='utf-8') as f:
            f.write(
                "JRNL        AUTH   G.FERMI,M.F.PERUTZ,B.SHAANAN,R.FOURME\n"
                "JRNL        TITL   THE CRYSTAL STRUCTURE OF HUMAN DEOXYHAEMOGLOBIN AT\n"
                "JRNL        TITL 2 1.74 A RESOLUTION\n"
                "JRNL        REF    J.MOL.BIOL.                   V. 175   159 1984\n"
                "JRNL        PUBL   LONDON : ACADEMIC PRESS\n"
                "JRNL        REFN                   ISSN 0022-2836\n"
                "JRNL        PMID   6726807\n"
                "JRNL        DOI    10.1016/0022-2836(84)90472-8\n"
            )
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(
            pdb_data.meta.journal.author,
            ["G.FERMI", "M.F.PERUTZ", "B.SHAANAN", "R.FOURME"]
        )
        self.assertEqual(
            pdb_data.meta.journal.title,
            "THE CRYSTAL STRUCTURE OF HUMAN DEOXYHAEMOGLOBIN AT 1.74 A RESOLUTION"
        )
        self.assertEqual(
            pdb_data.meta.journal.journal,
            "J.MOL.BIOL."
        )
        self.assertEqual(
            pdb_data.meta.journal.volume,
            "175"
        )
        self.assertEqual(
            pdb_data.meta.journal.pages,
            "159"
        )
        self.assertEqual(
            pdb_data.meta.journal.year,
            "1984"
        )
        self.assertEqual(
            pdb_data.meta.journal.publisher,
            "LONDON : ACADEMIC PRESS"
        )
        self.assertEqual(
            pdb_data.meta.journal.issn,
            "0022-2836"
        )
        self.assertEqual(
            pdb_data.meta.journal.pmid,
            "6726807"
        )
        self.assertEqual(
            pdb_data.meta.journal.doi,
            "10.1016/0022-2836(84)90472-8"
        )
        os.remove(TEMP_PDB_PATH)

    def test_revdat_parser(self):
        self.assertEqual(
            parse_revdat(
                "REVDAT   2   15-OCT-99 1ABC    1       REMARK"),
            ("15-OCT-99", ["REMARK"]),
        )

    def test_revdat(self):
        with open(TEMP_PDB_PATH, 'w', encoding='utf-8') as f:
            f.write(
                "REVDAT   2   11-MAR-08 2ABC    1       JRNL   VERSN\n"
                "REVDAT   1   09-DEC-03 2ABC    0"
            )
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(
            pdb_data.meta.revisions,
            [
                RevisionInfo(
                    date="09-DEC-03",
                    modifications=[],
                ),
                RevisionInfo(
                    date="11-MAR-08",
                    modifications=["JRNL", "VERSN"],
                ),
            ]
        )
        os.remove(TEMP_PDB_PATH)

    def test_sprsde(self):
        with open(TEMP_PDB_PATH, 'w', encoding='utf-8') as f:
            f.write(
                "SPRSDE     27-FEB-95 1GDJ      1LH4 2LH4\n"
                "SPRSDE                         3LH4\n"
            )
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(
            pdb_data.meta.replace,
            ReplaceInfo(
                date="27-FEB-95",
                ids=["1LH4", "2LH4", "3LH4"],
            )
        )
        os.remove(TEMP_PDB_PATH)


class TestPdbPrimaryStructureReader(unittest.TestCase):
    def test_dbref_parser(self):
        self.assertEqual(
            parse_dbref(
                "DBREF  2JHQ A    1   226  UNP    Q9KPK8   UNG_VIBCH        1    226 "
            ),
            ["A", 1, "", 226, "", "UNP", "Q9KPK8", "UNG_VIBCH", 1, "", 226, ""],
        )

    def test_dbref1_parser(self):
        line = (
            "DBREF1 1ABC A   61   322  UNIMES               UPI000148A153                   "
        )
        self.assertEqual(
            parse_dbref1(line),
            ["A", 61, "", 322, "", "UNIMES", "UPI000148A153"],
        )

    def test_dbref2_parser(self):
        line = (
            "DBREF2 1ABC A     46197919                      1534489     1537377"
        )
        self.assertEqual(
            parse_dbref2(line),
            ["A", "46197919", 1534489, 1537377],
        )

    def test_seqadv_parser(self):
        self.assertEqual(
            parse_seqadv(
                "SEQADV 3ABC MET A   -1  UNP  P10725              EXPRESSION TAG"
            ),
            ["MET", "A", -1, "", "UNP", "P10725", "", None, "EXPRESSION TAG"],
        )

    def test_seqres_parser(self):
        self.assertEqual(
            parse_seqres(
                "SEQRES   1 A   21  GLY ILE VAL GLU GLN CYS CYS THR SER ILE CYS SER LEU\n"
            ),
            ["A", 21, ["GLY", "ILE", "VAL", "GLU", "GLN", "CYS",
                       "CYS", "THR", "SER", "ILE", "CYS", "SER", "LEU"]],
        )

    def test_modres_parser(self):
        self.assertEqual(
            parse_modres(
                "MODRES 2R0L ASN A   74  ASN  GLYCOSYLATION SITE "
            ),
            ["ASN", "A", 74, "", "ASN", "GLYCOSYLATION SITE"],
        )

    def test_dbref(self):
        with open(TEMP_PDB_PATH, "w", encoding="utf-8") as f:
            f.write(
                "DBREF  2JHQ A    1   226  UNP    Q9KPK8   UNG_VIBCH        1    226 \n")
            f.write(
                "ATOM      1  N   LYS A   1      35.365  22.342 -11.980  1.00 22.28           N  \n")
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(
            pdb_data._tmp["primary"]["A"].sequence_db,
            SequenceDBInfo(
                init_seq_num=1,
                end_seq_num=226,
                database="UNP",
                db_accession="Q9KPK8",
                db_id_code="UNG_VIBCH",
                db_init_seq_num=1,
                db_end_seq_num=226
            )
        )
        os.remove(TEMP_PDB_PATH)

    def test_split_dbref(self):
        with open(TEMP_PDB_PATH, "w", encoding="utf-8") as f:
            f.write(
                "DBREF1 1ABC A   61   322  UNIMES               UPI000148A153                   \n")
            f.write(
                "DBREF2 1ABC A     46197919                      1534489     1537377\n")
            f.write(
                "ATOM      1  N   LYS A   1      35.365  22.342 -11.980  1.00 22.28           N  \n")
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(
            pdb_data._tmp["primary"]["A"].sequence_db,
            SequenceDBInfo(
                init_seq_num=61,
                end_seq_num=322,
                database="UNIMES",
                db_accession="46197919",
                db_id_code="UPI000148A153",
                db_init_seq_num=1534489,
                db_end_seq_num=1537377
            )
        )
        os.remove(TEMP_PDB_PATH)

    def test_seqadv(self):
        with open(TEMP_PDB_PATH, "w", encoding="utf-8") as f:
            f.write(
                "SEQADV 3ABC MET A   -1  UNP  P10725              EXPRESSION TAG\n")
            f.write(
                "ATOM      1  N   LYS A   1      35.365  22.342 -11.980  1.00 22.28           N  \n")
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(
            pdb_data._tmp["primary"]["A"].sequence_differences,
            [SequenceDifferenceInfo(
                seq_num=-1,
                alt_res_name="MET",
                database="UNP",
                db_accession="P10725",
                db_res_num=None,
                db_res_name="",
                comment="EXPRESSION TAG",
            )],
        )
        os.remove(TEMP_PDB_PATH)

    def test_modres(self):
        with open(TEMP_PDB_PATH, "w", encoding="utf-8") as f:
            f.write("MODRES 2R0L ASN A   74  ASN  GLYCOSYLATION SITE \n")
            f.write(
                "ATOM      1  N   LYS A   1      35.365  22.342 -11.980  1.00 22.28           N  \n")
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(
            pdb_data._tmp["primary"]["A"].residue_modifications,
            [ResidueModificationInfo(
                seq_num=74,
                std_res="ASN",
                comment="GLYCOSYLATION SITE",
            )]
        )
        os.remove(TEMP_PDB_PATH)

    def test_seqres(self):
        with open(TEMP_PDB_PATH, "w", encoding="utf-8") as f:
            f.write(
                "SEQRES   1 A   21  GLY ILE VAL GLU GLN CYS CYS THR SER ILE CYS SER LEU\n")
            f.write("SEQRES   2 A   21  TYR GLN LEU GLU ASN TYR CYS ASN\n")
            f.write(
                "ATOM      1  N   LYS A   1      35.365  22.342 -11.980  1.00 22.28           N  \n")
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(
            pdb_data._tmp["primary"]["A"].sequence,
            [
                "GLY", "ILE", "VAL", "GLU", "GLN", "CYS", "CYS", "THR", "SER",
                "ILE", "CYS", "SER", "LEU",
                "TYR", "GLN", "LEU", "GLU", "ASN", "TYR", "CYS", "ASN",
            ],
        )
        os.remove(TEMP_PDB_PATH)


class TestPdbHeterogenSectionReader(unittest.TestCase):
    def test_het_parser(self):
        self.assertEqual(
            parse_het(
                "HET    TRS  B 975       8"
            ),
            ["TRS", ""],
        )

    def test_hetnam_parser(self):
        self.assertEqual(
            parse_hetnam(
                "HETNAM     SAD BETA-METHYLENE SELENAZOLE-4-CARBOXAMIDE ADENINE"
            ),
            ["SAD", "BETA-METHYLENE SELENAZOLE-4-CARBOXAMIDE ADENINE"],
        )

    def test_hetsyn_parser(self):
        self.assertEqual(
            parse_hetsyn(
                "HETSYN     TRS TRIS BUFFER; "
            ),
            ["TRS", "TRIS BUFFER;"],
        )

    def test_formul_parser(self):
        self.assertEqual(
            parse_formul(
                "FORMUL   3   MG    2(MG 2+)"
            ),
            ["MG", "2(MG 2+)"],
        )

    def test_het(self):
        with open(TEMP_PDB_PATH, "w", encoding="utf-8") as f:
            f.write(
                "HET    MG   B 975       8\n"
                "HETNAM     MG  MAGNESIUM \n"
                "HETSYN     MG  MAGNESIUM ION; \n"
                "FORMUL   3   MG    2(MG 2+)\n"
            )
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(pdb_data.heterogen["MG"], Heterogen(
            name="MAGNESIUM",
            comment="",
            alias=["MAGNESIUM ION"],
            formula="2(MG 2+)",
        ))
        os.remove(TEMP_PDB_PATH)


class TestPdbSecondaryStructureReader(unittest.TestCase):
    def test_helix_parser(self):
        self.assertEqual(
            parse_helix(
                "HELIX    1  HA GLY A   86  GLY A   94  1                                   9 "
            ),
            ["HA", 'A', 86, '', 94, '', 1, ''],
        )

    def test_helix(self):
        with open(TEMP_PDB_PATH, "w", encoding="utf-8") as f:
            f.write(
                "HELIX    1  HA GLY A   86  GLY A   94  1                                   9 "
            )
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(
            pdb_data.secondary_structure.helix["HA"],
            Helix(
                chain_id="A",
                init_seq_num=86,
                end_seq_num=94,
                helix_class=1,
                comment=""
            )
        )
        os.remove(TEMP_PDB_PATH)

    def test_sheet_parser(self):
        self.assertEqual(
            parse_sheet(
                "SHEET    1   A 5 THR A 107  ARG A 110  0"
            ),
            ['A', 'A', 107, '', 0, '', '', None, '', '', '', None, '', 'A', 110, ''])

        self.assertEqual(
            parse_sheet(
                "SHEET    2   A 5 ILE A  96  THR A  99 -1  N  LYS A  98   O  THR A 107"
            ),
            ['A', 'A', 96, '', -1, 'N', 'A', 98,
                '', 'O', 'A', 107, '', 'A', 99, ''],
        )

    def test_sheet(self):
        with open(TEMP_PDB_PATH, "w", encoding="utf-8") as f:
            f.write(
                "SHEET    1   A 2 THR A 107  ARG A 110  0\n"
                "SHEET    2   A 2 ILE A  96  THR A  99 -1  N  LYS A  98   O  THR A 107"
            )
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(str(pdb_data.secondary_structure.sheet["A"][0]), str(
            SheetStrand(
                init_chain_id="A",
                init_seq_num=107,
                end_chain_id="A",
                end_res_seq=110,
                sense=0
            )
        ))


class TestPdbMiscFeatureReader(unittest.TestCase):
    def test_site_parser(self):
        self.assertEqual(
            parse_site(
                "SITE     1 AC1  3 HIS A  94  HIS A  96  HIS A 119                    "
            ),
            ["AC1", [["A", 94, ""], ["A", 96, ""], ["A", 119, ""]]],
        )

    def test_site(self):
        with open(TEMP_PDB_PATH, "w", encoding="utf-8") as f:
            f.write(
                "SITE     1 AC2  5 ASN A  62  GLY A  63  HIS A  64  HOH A 328                    \n"
                "SITE     2 AC2  5 HOH A 634                     "
            )
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(pdb_data.sites["AC2"][0], Site(
            seq_num=62, chain_id="A", icode=""))
        self.assertEqual(pdb_data.sites["AC2"][1], Site(
            seq_num=63, chain_id="A", icode=""))
        self.assertEqual(pdb_data.sites["AC2"][2], Site(
            seq_num=64, chain_id="A", icode=""))
        self.assertEqual(pdb_data.sites["AC2"][3], Site(
            seq_num=328, chain_id="A", icode=""))
        self.assertEqual(pdb_data.sites["AC2"][4], Site(
            seq_num=634, chain_id="A", icode=""))
        os.remove(TEMP_PDB_PATH)


class TestPdbCrystallographicReader(unittest.TestCase):
    def test_cryst1_parser(self):
        self.assertEqual(
            parse_cryst1(
                "CRYST1   52.000   58.600   61.900  90.00  90.00  90.00 P 21 21 21    8 \n"
            ),
            [52.000, 58.600, 61.900, 90.00, 90.00, 90.00, "P 21 21 21", 8],
        )

    def test_cryst1(self):
        with open(TEMP_PDB_PATH, "w", encoding="utf-8") as f:
            f.write(
                "CRYST1   52.000   58.600   61.900  90.00  90.00  90.00 P 21 21 21    8 \n"
            )
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(pdb_data.crystallographic.cell_lengths, [
                         52.000, 58.600, 61.900])
        self.assertEqual(pdb_data.crystallographic.cell_angles, [
                         90.00, 90.00, 90.00])
        self.assertEqual(pdb_data.crystallographic.space_group, "P 21 21 21")
        self.assertEqual(pdb_data.crystallographic.z, 8)
        os.remove(TEMP_PDB_PATH)

    def test_matrix_line_parser(self):
        self.assertEqual(
            parse_matrix_line(
                "ORIGX1      0.963457  0.136613  0.230424      16.61000\n"
            ),
            [0.963457, 0.136613, 0.230424, 16.61000],
        )

    def test_origin_matrix(self):
        with open(TEMP_PDB_PATH, "w", encoding="utf-8") as f:
            f.write(
                'ORIGX1      0.963457  0.136613  0.230424       16.61000\n'
                'ORIGX2     -0.158977  0.983924  0.081383       13.72000\n'
                'ORIGX3     -0.215598 -0.115048  0.969683       37.00000\n'
            )
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(pdb_data.crystallographic.origin_matrix, [
            [0.963457, 0.136613, 0.230424, 16.61000],
            [-0.158977, 0.983924, 0.081383, 13.72000],
            [-0.215598, -0.115048, 0.969683, 37.00000],
        ])
        os.remove(TEMP_PDB_PATH)

    def test_scale_matrix(self):
        with open(TEMP_PDB_PATH, "w", encoding="utf-8") as f:
            f.write(
                "SCALE1      0.019231  0.000000  0.000000        0.00000\n"
                "SCALE2      0.000000  0.017065  0.000000        0.00000\n"
                "SCALE3      0.000000  0.000000  0.016155        0.00000\n"
            )
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(pdb_data.crystallographic.scale_matrix, [
            [0.019231, 0.0, 0.0, 0.0],
            [0.0, 0.017065, 0.0, 0.0],
            [0.0, 0.0, 0.016155, 0.0],
        ])
        os.remove(TEMP_PDB_PATH)

    def test_ncs_matrix_line_parser(self):
        self.assertEqual(
            parse_ncs_matrix_line(
                "MTRIX1   1 -1.000000  0.000000  0.000000        0.00000    1"
            ),
            [-1.000000, 0.000000, 0.000000, 0.000000, 1],
        )

    def test_ncs_matrix(self):
        with open(TEMP_PDB_PATH, "w", encoding="utf-8") as f:
            f.write(
                "MTRIX1   1 -1.000000  0.000000  0.000000        0.00000    1\n"
                "MTRIX2   1  0.000000  1.000000  0.000000        0.00000    1\n"
                "MTRIX3   1  0.000000  0.000000 -1.000000        0.00000    1\n"
                "MTRIX1   2  1.000000  0.000000  0.000000        0.00000    0\n"
                "MTRIX2   2  0.000000  1.000000  0.000000        0.00000    0\n"
                "MTRIX3   2  0.000000  0.000000  1.000000        0.00000    0\n"
            )
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(
            pdb_data.crystallographic.ncs_matrix,
            [
                NcsMatrix(
                    matrix=[[-1.000000, 0.000000, 0.000000, 0.000000],
                            [0.000000, 1.000000, 0.000000, 0.000000],
                            [0.000000, 0.000000, -1.000000, 0.000000]],
                    given=True
                ),
                NcsMatrix(
                    matrix=[[1.000000, 0.000000, 0.000000, 0.000000],
                            [0.000000, 1.000000, 0.000000, 0.000000],
                            [0.000000, 0.000000, 1.000000, 0.000000]],
                    given=False
                )
            ]
        )
        os.remove(TEMP_PDB_PATH)


class TestPdbConnectivityAnnotationReader(unittest.TestCase):
    def test_ssbond_parser(self):
        line = "SSBOND   1 CYS A   26    CYS A   84                          1555   1555  2.04"
        self.assertEqual(
            parse_ssbond(line.ljust(80)),
            ["A", 26, "", "A", 84, "", "1555", "1555", 2.04],
        )
    
    def test_ssbond(self):
        with open(TEMP_PDB_PATH, "w", encoding="utf-8") as f:
            f.write(
                "SSBOND   1 CYS A   26    CYS A   84                          1555   1555  2.04\n"
            )
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(
            pdb_data.connectivity.ss_bond[0],
            SsBond(
                chain_id_1="A",
                seq_num_1=26,
                chain_id_2="A",
                seq_num_2=84,
                symmetry_operation_1="1555",
                symmetry_operation_2="1555",
                distance=2.04
            ))
        os.remove(TEMP_PDB_PATH)

    def test_link_parser(self):
        line = "LINK         O   SER A  82                 OG  SER B   5     1555   1555  2.84"
        self.assertEqual(
            parse_link(line.ljust(80)),
            ["O", "", "SER", "A", 82, "", "OG", "",
                "SER", "B", 5, "", "1555", "1555", 2.84],
        )

    def test_link(self):
        with open(TEMP_PDB_PATH, "w", encoding="utf-8") as f:
            f.write(
                "LINK         O   SER A  82                 OG  SER B   5     1555   1555  2.84\n"
            )
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(
            pdb_data.connectivity.link[0],
            Link(
                name_1="O",
                alt_loc_1="",
                chain_id_1="A",
                seq_num_1=82,
                icode_1="",
                name_2="OG",
                alt_loc_2="",
                chain_id_2="B",
                seq_num_2=5,
                icode_2="",
                distance=2.84,
            ))
        os.remove(TEMP_PDB_PATH)

    def test_cispep_parser(self):
        line = "CISPEP   1 ASP A 118     PRO A 119          0        -0.24"
        self.assertEqual(
            parse_cispep(line.ljust(80)),
            ["A", 118, "", "A", 119, "", 0, -0.24],
        )

    def test_cispep(self):
        with open(TEMP_PDB_PATH, "w", encoding="utf-8") as f:
            f.write(
                "CISPEP   1 ASP A 118     PRO A 119          0        -0.24\n"
            )
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(
            pdb_data.connectivity.cis_peptide[0],
            CisPeptide(
                chain_id_1="A", seq_num_1=118, icode_1="",
                chain_id_2="A", seq_num_2=119, icode_2="",
                num_model=0, measure=-0.24
            ))
        os.remove(TEMP_PDB_PATH)

class TestPdbCoordinateReader(unittest.TestCase):
    def test_atom_parser(self):
        line = "ATOM      1  N   MET A   1      11.104  13.207  10.000  1.00 20.00           N  "
        self.assertEqual(
            parse_atom(line),
            ["ATOM", 1, "N", "", "MET", "A", 1, "", 11.104, 13.207, 10.0, 1.0, 20.0, "N", 0],
        )

    def test_anisou_parser(self):
        line = "ANISOU  107  N   GLY A  13     2406   1892   1614    198    519   -328       N  "
        self.assertEqual(
            parse_anisou(line),
            [
                0.2406,
                0.1892,
                0.1614,
                0.0198,
                0.0519,
                -0.0328,
            ],
        )

    def test_connect_parser(self):
        line = "CONECT 1179  746 1184 1195 1203"
        self.assertEqual(parse_connect(line), (1179, [1184, 1195, 1203]))


class TestPdbBookkeepingReader(unittest.TestCase):
    def test_master_parser(self):
        line = "MASTER       40    0    0    0    0    0    0    6 2930    2    0   29   "
        self.assertEqual(
            parse_master(line.ljust(80)),
            [40, 0, 0, 0, 0, 6, 2930, 2, 0, 29],
        )

    def test_master(self):
        with open(TEMP_PDB_PATH, "w", encoding="utf-8") as f:
            f.write(
                "MASTER       40    0    0    0    0    0    0    6 2930    2    0   29  "
            )
        reader = PdbReader(TEMP_PDB_PATH)
        pdb_data = reader.read()
        self.assertEqual(pdb_data._validation_info, {
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
        })
        os.remove(TEMP_PDB_PATH)


if __name__ == "__main__":
    unittest.main()
