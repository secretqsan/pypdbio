# pylint: disable=missing-module-docstring
# pylint: disable=protected-access
# pylint: disable=unused-wildcard-import
# pylint: disable=wildcard-import
import warnings
from .models import *
from .reader_helper import *
from .utils import next_chain_id

class PdbReader:
    """
    PDB reader class.
    """

    def __init__(self, path):
        """ Initializes the PdbReader with the given file path.
        :param path: Path to the PDB file to be read"""
        self.pdb_path = path
        self.__validation_info = {
            "num_remark": 0,
            "num_het": 0,
            "num_helix": 0,
            "num_sheet": 0,
            "num_site": 0,
            "num_xform": 0,
            "num_coord": 0,
            "num_ter": 0,
            "num_conect": 0,
            "num_seq": 0,
        }
        self.__n_model_expected = -1
        self.__current_model = -1
        self.__current_chain_id = ''
        self.__current_chain_index = -1
        self.__occupied_chain_ids = []
        self.__current_residue_id = -1

    #post processing functions
    def __check_status(self, pdb_data):
        if pdb_data.meta.obsolete is not None:
            obslte_info = pdb_data.meta.obsolete
            warnings.warn(
                f"It seems that the PDB file you use is outdated at {obslte_info.replace_date}."
                f"The new entry ID is {obslte_info.new_entry_id}."
            )
        if pdb_data.meta.split is not None:
            warnings.warn(
                "It seems that the PDB file you use is incomplete."
                f" The other parts can be found by id {', '.join(pdb_data.meta.split)}."
            )
        if len(pdb_data.meta.caveat) > 0:
            warnings.warn(
                f"Something wrong is found in the PDB file:\n {pdb_data.meta.caveat}"
            )
        if self.__n_model_expected != -1 and len(pdb_data.models) != self.__n_model_expected:
            warnings.warn(
                "The number of models in the PDB file does not match the expected number."
            )
        if pdb_data._validation_info != {}:
            def check_field_ok(field_name):
                field_name_str = f"num_{field_name}"
                if pdb_data._validation_info[field_name_str] != self.__validation_info[field_name_str]:
                    warnings.warn(
                        f"The number of {field_name} in the PDB file does not match the expected number."
                    )
            target_fields = [
                "remark", "het", "helix",
                "sheet", "site", "xform",
                "coord", "ter", "conect",
                "seq"
            ]
            for field_name in target_fields:
                check_field_ok(field_name)

    def __post_process(self, pdb_data):
        primary = pdb_data._tmp["primary"]
        for sequence_id, data in primary.items():
            sequence = pdb_data.models[0][sequence_id]
            sequence.sequence_info = data
            sequence_length = len(sequence)
            if sequence.sequence_info.sequence_differences is not None:
                for entry in sequence.sequence_info.sequence_differences:
                    if entry.seq_num < sequence_length and entry.seq_num > 1:
                        entry.alt_res_name = ""
        primary = None
        meta = pdb_data.meta
        if meta._compound_text != "":
            compound_entries = parse_kv_entries(meta._compound_text)
            for entry in compound_entries:
                key, value = entry
                if "mol_id" == key:
                    safe_append(meta, "compounds", {})
                else:
                    meta.compounds[-1][key] = value
        if meta._source_text != "":
            source_entries = parse_kv_entries(meta._source_text)
            n_source = sum(1 for entry in source_entries if entry[0] == "mol_id")
            if len(meta.compounds) != n_source:
                warnings.warn("The number of compounds and sources are not the same.")
            else:
                current_compound_index = -1
                for entry in source_entries:
                    key, value = entry
                    if key == "mol_id" :
                        current_compound_index += 1
                        meta.compounds[current_compound_index]["source"] = {}
                    else:
                        meta.compounds[current_compound_index]["source"][key] = value
        meta._source_text = ""
        meta._compound_text = ""
        if meta.revisions is not None:
            meta.revisions.reverse()
        if pdb_data.heterogen is not None:
            for _, het in pdb_data.heterogen.items():
                het.alias = [
                    alias.strip()
                    for alias
                    in het._alias_text.split(";")
                    if alias.strip() != ""
                ]
                het._alias_text = ""

    # section parsing functions
    def __parse_title_section(self, line, pdb_data):
        meta = pdb_data.meta
        if line.startswith("HEADER"):
            classification, date, pdb_id = parse_header(line)
            meta.date = date
            meta.classification = classification
            meta.pdb_id = pdb_id
        elif line.startswith("TITLE"):
            meta.title = parse_continuous_field(
                line, meta.title, parse_title)
        elif line.startswith("OBSLTE"):
            rep_date, new_entry_id = parse_obslte(line)
            meta.obsolete = ObsoleteInfo(
                replace_date=rep_date,
                new_entry_id=new_entry_id
            )
        elif line.startswith("SPLIT"):
            meta.split = parse_continuous_field(
                line,
                meta.split,
                parse_split
            )
        elif line.startswith("CAVEAT"):
            meta.caveat = parse_continuous_field(
                line,
                meta.caveat,
                parse_caveat
            )
        elif line.startswith("COMPND"):
            meta._compound_text = parse_continuous_field(
                line, meta._compound_text, parse_compond
            )
        elif line.startswith("SOURCE"):
            meta._source_text = parse_continuous_field(
                line, meta._source_text, parse_source
            )
        elif line.startswith("KEYWDS"):
            line_keywords = parse_keywords(line)
            if line[8:10].strip() != '' and len(line_keywords) > 0:
                line_keywords[0] = line_keywords[0].split(maxsplit=1)[-1]
            safe_extend(meta, "keywords", [keyword for keyword in line_keywords if keyword != ''])
        elif line.startswith("EXPDTA"):
            meta.experiment = parse_continuous_field(
                line, meta.experiment, parse_expdta
            )
        elif line.startswith("NUMMDL"):
            self.__n_model_expected = parse_nummdl(line)
        elif line.startswith("MDLTYP"):
            meta.model_type = parse_continuous_field(
                line, meta.model_type, parse_mdltyp
            )
        elif line.startswith("AUTHOR"):
            meta.author = parse_continuous_field(
                line, meta.author, parse_author
            )
        elif line.startswith("REVDAT"):
            date, modifications = parse_revdat(line)
            safe_append(meta, "revisions", RevisionInfo(
                date=date,
                modifications=modifications
            ))
        elif line.startswith("SPRSDE"):
            if meta.replace is None:
                meta.replace = ReplaceInfo()
            meta.replace.date, meta.replace.ids = parse_continuous_field(
                line,
                [
                    meta.replace.date,
                    meta.replace.ids
                ],
                parse_sprsde,
                place=1
            )
        elif line.startswith("JRNL"):
            if meta.journal is None:
                meta.journal = JournalInfo()
            jrnl_type, _ = parse_jrnl(line)
            if jrnl_type == "AUTH":
                meta.journal.author = parse_continuous_field(
                    line, meta.journal.author, parse_jrnl_auth
                )
            elif jrnl_type == "TITL":
                meta.journal.title = parse_continuous_field(
                    line, meta.journal.title, parse_jrnl_titl
                )
            elif jrnl_type == "EDIT":
                meta.journal.editor = parse_continuous_field(
                    line, meta.journal.editor, parse_jrnl_edit
                )
            elif jrnl_type == "REF":
                j = meta.journal
                j.journal, j.volume, j.pages, j.year = parse_continuous_field(
                    line,
                    [
                        j.journal, j.volume, j.pages, j.year,
                    ],
                    parse_jrnl_ref,
                    place=0
                )
            elif jrnl_type == "PUBL":
                meta.journal.publisher = parse_continuous_field(
                    line, meta.journal.publisher, parse_jrnl_publ
                )
            elif jrnl_type == "REFN":
                refn_type, refn_value = parse_jrnl_refn(line)
                if refn_type != "" and refn_value != "":
                    if refn_type == "ISSN":
                        meta.journal.issn = refn_value
                    elif refn_type == "ESSN":
                        meta.journal.essn = refn_value
            elif jrnl_type == "PMID":
                meta.journal.pmid = parse_continuous_field(
                    line, meta.journal.pmid, parse_jrnl_pmid
                )
            elif jrnl_type == "DOI":
                meta.journal.doi = parse_continuous_field(
                    line, meta.journal.doi, parse_jrnl_doi, separator=""
                )
        elif line.startswith("REMARK"):
            self.__validation_info["num_remark"] += 1
            comment_id, comment_text = parse_remark(line)
            if pdb_data.meta.remark is None:
                pdb_data.meta.remark = {}
            if comment_id not in pdb_data.meta.remark:
                meta.remark[comment_id] = comment_text
            else:
                meta.remark[comment_id] += f"\n{comment_text}"

    def __parse_primary_structure_section(self, line, pdb_data):
        primary_structure = pdb_data._tmp["primary"]
        def check_chain_id(chain_id):
            if chain_id not in primary_structure:
                primary_structure[chain_id] = SequenceInfo()
        if line.startswith("DBREF "):
            chain_id, seq_begin, insert_begin, seq_end, insert_end,\
            database, db_accession, db_id_code, db_seq_begin,\
            db_begin_icode, db_seq_end, db_end_icode = parse_dbref(line)
            check_chain_id(chain_id)
            primary_structure[chain_id].sequence_db = SequenceDBInfo(
                init_seq_num=seq_begin,
                init_icode=insert_begin,
                end_seq_num=seq_end,
                end_icode=insert_end,
                database=database,
                db_accession=db_accession,
                db_id_code=db_id_code,
                db_init_seq_num=db_seq_begin,
                db_init_icode=db_begin_icode,
                db_end_seq_num=db_seq_end,
                db_end_icode=db_end_icode,
            )
        elif line.startswith("DBREF1"):
            chain_id, seq_begin, insert_begin, seq_end, insert_end, database, db_id_code = parse_dbref1(line)
            check_chain_id(chain_id)
            primary_structure[chain_id].sequence_db = SequenceDBInfo(
                init_seq_num=seq_begin,
                end_seq_num=seq_end,
                database=database,
                db_accession="",
                db_init_seq_num=0,
                db_end_seq_num=0,
                db_id_code=db_id_code,
                init_icode=insert_begin,
                end_icode=insert_end,
            )
        elif line.startswith("DBREF2"):
            chain_id, db_accession, db_seq_begin, db_seq_end = parse_dbref2(line)
            check_chain_id(chain_id)
            dbref = primary_structure[chain_id].sequence_db
            dbref.db_accession = db_accession
            dbref.db_init_seq_num = db_seq_begin
            dbref.db_end_seq_num = db_seq_end
        elif line.startswith("SEQADV"):
            res_name, chain_id, seq_num, icode, database, db_accession, db_res, db_seq, conflict = parse_seqadv(line)
            check_chain_id(chain_id)
            safe_append(primary_structure[chain_id], "sequence_differences", SequenceDifferenceInfo(
                seq_num=seq_num,
                database=database,
                alt_res_name=res_name,
                db_accession=db_accession,
                db_res_num=db_seq,
                db_res_name=db_res,
                comment=conflict,
                icode=icode,
            ))
        elif line.startswith("SEQRES"):
            self.__validation_info["num_seq"] += 1
            chain_id, _, sequence = parse_seqres(line)
            check_chain_id(chain_id)
            safe_extend(primary_structure[chain_id], "sequence", sequence)
        elif line.startswith("MODRES"):
            _, chain_id, seq_num, icode, std_res, comment = parse_modres(line)
            check_chain_id(chain_id)
            safe_append(primary_structure[chain_id], "residue_modifications", ResidueModificationInfo(
                seq_num=seq_num,
                icode=icode,
                std_res=std_res,
                comment=comment,
            ))

    def __parse_coordinate_section(self, line, pdb_data):
        if line.startswith('MODEL'):
            self.__current_model += 1
            pdb_data.models.append(Model())
            self.__current_chain_id = ''
            self.__current_residue_id = -1

        elif line.startswith('ATOM') or line.startswith('HETATM'):
            self.__validation_info["num_coord"] += 1
            if self.__current_model == -1:
                self.__current_model += 1
                pdb_data.models.append(Model())
            typ, _, atom_name, alt_loc,\
                res_name, chain_id, residue_id, icode,\
                coord_x, coord_y, coord_z,\
                occupancy, temp_factor, element, charge = parse_atom(line)
            if chain_id != self.__current_chain_id:
                new_chain = Chain()
                target_chain_id = next_chain_id(
                    self.__current_chain_id, self.__occupied_chain_ids
                )
                if chain_id != target_chain_id:
                    new_chain.id = target_chain_id
                    self.__occupied_chain_ids.append(chain_id)
                self.__current_chain_id = chain_id
                pdb_data.models[self.__current_model].add_chain(new_chain)
                self.__current_residue_id = -1
            if residue_id != self.__current_residue_id:
                self.__current_residue_id = residue_id
                pdb_data.models[self.__current_model].chains[-1].add_residue(Residue(res_name))
                mol_type = "het" if typ == "HETATM" else "std"
                pdb_data.models[self.__current_model].chains[-1].residues[-1].het = mol_type
                pdb_data.models[self.__current_model].chains[-1].residues[-1].icode = icode
            pdb_data.models[-1].chains[-1].residues[-1].add_atom(Atom(
                name=atom_name,
                coord=(coord_x, coord_y, coord_z),
                temp_factor=temp_factor,
                occupancy=occupancy,
                alt_loc=alt_loc,
                element=element,
                charge=charge
            ))
        elif line.startswith("TER   "):
            self.__validation_info["num_ter"] += 1
            pdb_data.models[-1].chains[-1].residues[-1].is_ter = True
        elif line.startswith("ANISOU"):
            u_list = parse_anisou(line)
            pdb_data.models[-1].chains[-1].residues[-1].atoms[-1].temp_factor = u_list

    def __parse_heterogen_section(self, line, pdb_data):
        heterogen = pdb_data.heterogen
        def check_het_id(het_id):
            if het_id not in heterogen:
                heterogen[het_id] = Heterogen()
        if line.startswith("HET   "):
            self.__validation_info["num_het"] += 1
            het_id, text = parse_het(line)
            check_het_id(het_id)
            heterogen[het_id].comment = text
        elif line.startswith("HETNAM"):
            het_id, text = parse_hetnam(line)
            check_het_id(het_id)
            heterogen[het_id].name += text
        elif line.startswith("HETSYN"):
            het_id, text = parse_hetsyn(line)
            check_het_id(het_id)
            heterogen[het_id]._alias_text += text
        elif line.startswith("FORMUL"):
            het_id, text = parse_formul(line)
            check_het_id(het_id)
            heterogen[het_id].formula += text

    def __parse_secondary_structure_section(self, line, pdb_data):
        if line.startswith("HELIX "):
            self.__validation_info["num_helix"] += 1
            if pdb_data.secondary_structure.helix is None :
                pdb_data.secondary_structure.helix = {}
            helix = pdb_data.secondary_structure.helix
            helix_id, chain_id, init_seq_num, init_icode,\
                end_res_num, end_icode, helix_class, comment = parse_helix(line)
            helix[helix_id] = Helix(
                chain_id=chain_id,
                init_seq_num=init_seq_num,
                end_seq_num=end_res_num,
                helix_class=helix_class,
                comment=comment,
                init_icode=init_icode,
                end_icode=end_icode,
            )
        elif line.startswith("SHEET "):
            self.__validation_info["num_sheet"] += 1
            if pdb_data.secondary_structure.sheet is None :
                pdb_data.secondary_structure.sheet = {}
            sheet = pdb_data.secondary_structure.sheet
            sheet_id, init_chain_id, init_seq_num, init_icode,\
            sense, cur_atom, cur_chain_id, cur_res_seq, cur_icode,\
            prev_atom, prev_chain_id, prev_res_seq, prev_icode,\
            end_chain_id, end_res_seq, end_icode = parse_sheet(line)
            if sheet_id not in sheet:
                sheet[sheet_id] = []
            sheet[sheet_id].append(SheetStrand(
                init_chain_id=init_chain_id,
                init_seq_num=init_seq_num,
                init_icode=init_icode,
                sense=sense,
                cur_atom=cur_atom,
                cur_chain_id=cur_chain_id,
                cur_res_seq=cur_res_seq,
                cur_icode=cur_icode,
                prev_atom=prev_atom,
                prev_chain_id=prev_chain_id,
                prev_res_seq=prev_res_seq,
                prev_icode=prev_icode,
                end_chain_id=end_chain_id,
                end_res_seq=end_res_seq,
                end_icode=end_icode,
            ))

    def __parse_connectivity_section(self, line, pdb_data):
        connectivity = pdb_data.connectivity
        if line.startswith("SSBOND"):
            chain_id, seq_num, icode, chain_id2, seq_num2, icode2, sym1, sym2, distance = parse_ssbond(line)
            safe_append(connectivity, "ss_bond", SsBond(
                chain_id_1=chain_id,
                seq_num_1=seq_num,
                icode_1=icode,
                chain_id_2=chain_id2,
                seq_num_2=seq_num2,
                icode_2=icode2,
                symmetry_operation_1=sym1,
                symmetry_operation_2=sym2,
                distance=distance,
            ))
        elif line.startswith("LINK  "):
            name1, alt_loc1, _, chain_id1, res_seq1, icode1,\
            name2, alt_loc2, _, chain_id2, res_seq2, icode2,\
            sym1, sym2, length = parse_link(line)
            safe_append(connectivity, "link", Link(
                name_1=name1,
                alt_loc_1=alt_loc1,
                chain_id_1=chain_id1,
                seq_num_1=res_seq1,
                icode_1=icode1,
                name_2=name2,
                alt_loc_2=alt_loc2,
                chain_id_2=chain_id2,
                seq_num_2=res_seq2,
                icode_2=icode2,
                distance=length,
            ))
        elif line.startswith("CISPEP"):
            chain_id1, seq_num1, icode1,\
            chain_id2, seq_num2, icode2,\
            num_model, measure = parse_cispep(line)
            safe_append(connectivity, "cis_peptide", CisPeptide(
                chain_id_1=chain_id1,
                seq_num_1=seq_num1,
                icode_1=icode1,
                chain_id_2=chain_id2,
                seq_num_2=seq_num2,
                icode_2=icode2,
                num_model=num_model,
                measure=measure,
            ))
        elif line.startswith('CONECT'):
            self.__validation_info["num_conect"] += 1
            atom1, bonded_atoms = parse_connect(line)
            if connectivity.connections is None:
                connectivity.connections = {}
            connectivity.connections[atom1] = bonded_atoms

    def __parse_misc_feature_section(self, line, pdb_data):
        if line.startswith("SITE  "):
            self.__validation_info["num_site"] += 1
            site_id, residues = parse_site(line)
            if site_id not in pdb_data.sites:
                pdb_data.sites[site_id] = []
            for residue in residues:
                pdb_data.sites[site_id].append(Site(
                    seq_num=residue[1],
                    chain_id=residue[0],
                    icode=residue[2],
                ))

    def __parse_crystallographic_section(self, line, pdb_data):
        crystallographic = pdb_data.crystallographic
        if line.startswith("CRYST1"):
            a, b, c, alpha, beta, gamma, s_group, z  = parse_cryst1(line)
            crystallographic.cell_lengths = [a, b, c]
            crystallographic.cell_angles = [alpha, beta, gamma]
            crystallographic.space_group = s_group
            crystallographic.z = z
        elif line.startswith("ORIGX"):
            self.__validation_info["num_xform"] += 1
            o1, o2, o3, t = parse_matrix_line(line)
            safe_append(crystallographic, "origin_matrix", [o1, o2, o3, t])
        elif line.startswith("SCALE"):
            self.__validation_info["num_xform"] += 1
            o1, o2, o3, t = parse_matrix_line(line)
            safe_append(crystallographic, "scale_matrix", [o1, o2, o3, t])
        elif line.startswith("MTRIX"):
            self.__validation_info["num_xform"] += 1
            o1, o2, o3, t, i_given = parse_ncs_matrix_line(line)
            if crystallographic.ncs_matrix is None or len(crystallographic.ncs_matrix[-1].matrix) == 3:
                safe_append(crystallographic, "ncs_matrix", NcsMatrix(
                    matrix=[],
                    given=i_given == 1
                ))
            safe_append(crystallographic.ncs_matrix[-1], "matrix", [o1, o2, o3, t])

    def __parse_bookkeeping_section(self, line, pdb_data):
        if line.startswith("MASTER"):
            num_remark, num_het, num_helix, num_sheet, num_site, num_xform,\
                num_coord, num_ter, num_conect, num_seq = parse_master(line)
            pdb_data._validation_info = {
                "num_remark": num_remark,
                "num_het": num_het,
                "num_helix": num_helix,
                "num_sheet": num_sheet,
                "num_site": num_site,
                "num_xform": num_xform,
                "num_coord": num_coord,
                "num_ter": num_ter,
                "num_conect": num_conect,
                "num_seq": num_seq,
            }

    def read(self):
        """ Reads the PDB file and parses it into a PdbData structure."""
        pdb_data = PdbData()
        with open(self.pdb_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.rstrip('\n')
                self.__parse_title_section(line, pdb_data)
                self.__parse_primary_structure_section(line, pdb_data)
                self.__parse_heterogen_section(line, pdb_data)
                self.__parse_secondary_structure_section(line, pdb_data)
                self.__parse_misc_feature_section(line, pdb_data)
                self.__parse_crystallographic_section(line, pdb_data)
                self.__parse_coordinate_section(line, pdb_data)
                self.__parse_connectivity_section(line, pdb_data)
                self.__parse_bookkeeping_section(line, pdb_data)

        self.__post_process(pdb_data)
        self.__check_status(pdb_data)
        return pdb_data
