# pylint: disable=missing-module-docstring
# pylint: disable=wildcard-import
# pylint: disable=unused-wildcard-import
# pylint: disable=protected-access
from functools import singledispatchmethod
from .writer_helper import *
from .models import (
    PdbData,
    PdbMetaData,
    Model,
    Chain,
)
from .utils import chain_id_of_index
from .unit import unit_config

class PdbWriter:
    """
    PDB writer class.
    """

    def __init__(self, path):
        """ Initializes the PdbWriter with the given file path.
        :param path: Path to the PDB file to be written"""
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

    def __write_header_section(self, pdb_data, f):
        meta = pdb_data.meta
        if meta.classification or meta.date or meta.pdb_id:
            f.write(warp_lines(
                gen_header(meta.classification, meta.date, meta.pdb_id)
            ))
        if meta.obsolete:
            f.write(
                warp_lines(
                    gen_obslte(
                        meta.obsolete.replace_date,
                        meta.obsolete.new_entry_id,
                        meta.pdb_id,
                    )
                )
            )
        if meta.title:
            f.write(warp_lines(gen_title(meta.title)))
        if meta.split:
            f.write(warp_lines(gen_split(meta.split)))
        if meta.caveat:
            f.write(warp_lines(gen_caveat(meta.pdb_id, meta.caveat)))
        if meta.compounds:
            compound_texts = []
            source_texts = []
            for index, comp in enumerate(meta.compounds):
                compound_texts.append(f"MOL_ID: {index + 1};")
                source_texts.append(f"MOL_ID: {index + 1};")
                for key, value in comp.items():
                    if key == "source":
                        for src_key, src_value in value.items():
                            source_texts.append(
                                f"{src_key.upper()}: {src_value};")
                    else:
                        compound_texts.append(f"{key.upper()}: {value};")
            compound_texts[-1] = compound_texts[-1].rstrip(";")
            source_texts[-1] = source_texts[-1].rstrip(";")
            f.write(warp_lines(gen_compnd(compound_texts)))
            f.write(warp_lines(gen_source(source_texts)))
        if meta.keywords:
            f.write(warp_lines(gen_keywds(meta.keywords)))
        if meta.experiment:
            f.write(warp_lines(gen_expdta(meta.experiment)))
        if len(pdb_data.models) > 1:
            f.write(warp_lines(gen_nummdl(len(pdb_data.models))))
        if meta.model_type:
            f.write(warp_lines(gen_mdltyp(meta.model_type)))
        if meta.author:
            f.write(warp_lines(gen_author(meta.author)))
        if meta.revisions:
            revision_sorted = meta.revisions[::-1]
            for index, revision in enumerate(revision_sorted):
                f.write(warp_lines(gen_revdat(
                    len(meta.revisions) - index,
                    revision.date,
                    meta.pdb_id,
                    1 if len(meta.revisions) - index > 1 else 0,
                    revision.modifications,
                )))
        if meta.replace:
            f.write(warp_lines(gen_sprsde(
                meta.replace.date,
                meta.pdb_id,
                meta.replace.ids,
            )))
        if meta.journal:
            if meta.journal.author and len(meta.journal.author) > 0:
                f.write(warp_lines(gen_jrnl_auth(meta.journal.author)))
            if meta.journal.title:
                f.write(warp_lines(gen_jrnl_titl(meta.journal.title)))
            if meta.journal.editor:
                f.write(warp_lines(gen_jrnl_editor(meta.journal.editor)))
            if meta.journal.journal == "":
                f.write(warp_lines(gen_jrnl_ref_to_be_published()))
            else:
                f.write(warp_lines(gen_jrnl_ref(
                    meta.journal.journal,
                    meta.journal.volume,
                    meta.journal.pages,
                    meta.journal.year
                )))
            if meta.journal.publisher:
                f.write(warp_lines(gen_jrnl_publisher(meta.journal.publisher)))
            if meta.journal.issn != "":
                f.write(warp_lines(gen_jrnl_refn("ISSN", meta.journal.issn)))
            if meta.journal.essn != "":
                f.write(warp_lines(gen_jrnl_refn("ESSN", meta.journal.essn)))
            if meta.journal.pmid != "":
                f.write(warp_lines(gen_jrnl_pmid(meta.journal.pmid)))
            if meta.journal.doi != "":
                f.write(warp_lines(gen_jrnl_doi(meta.journal.doi)))
        if meta.remark:
            for remark_id, remark_text in sorted(meta.remark.items(), key=lambda x: int(x[0])):
                text = gen_remark(remark_id, remark_text)
                self.__validation_info["num_remark"] += len(text)
                f.write(warp_lines(text))

    def __write_primary_structure_section(self, pdb_data, f):
        if not pdb_data.models:
            return
        pdb_id = pdb_data.meta.pdb_id
        model = pdb_data.models[0]
        occupied_chain_ids = []
        dbref_dict = {}
        seqadv_dict = {}
        seq_dict = {}
        modres_dict = {}
        for index, chain in enumerate(model.chains):
            if chain.id == "":
                chain_id = chain_id_of_index(index, occupied_chain_ids)
            else:
                chain_id = chain.id
                occupied_chain_ids.append(chain_id)
            if chain.sequence_info is not None:
                dbref_dict[chain_id] = chain.sequence_info.sequence_db
                seqadv_dict[chain_id] = chain.sequence_info.sequence_differences
                seq_dict[chain_id] = chain.sequence_info.sequence
                modres_dict[chain_id] = chain.sequence_info.residue_modifications

        if dbref_dict != {}:
            for chain_id, dbref in dbref_dict.items():
                if dbref is None:
                    continue
                f.write(warp_lines(
                    gen_dbref(
                        pdb_id, chain_id,
                        dbref.init_seq_num, dbref.init_icode,
                        dbref.end_seq_num, dbref.end_icode,
                        dbref.database, dbref.db_accession, dbref.db_id_code,
                        dbref.db_init_seq_num, dbref.db_init_icode,
                        dbref.db_end_seq_num, dbref.db_end_icode
                    )
                ))

        if seqadv_dict != {}:
            for chain_id, seqadv in seqadv_dict.items():
                if seqadv is None:
                    continue
                for entry in seqadv:
                    try:
                        res_name = pdb_data[0][chain_id][f"{entry.seq_num}{entry.icode:>1}"].name
                    except IndexError:
                        res_name = entry.alt_res_name
                    f.write(warp_lines(
                        gen_seqadv(
                            pdb_id, res_name, chain_id, entry.seq_num, entry.icode,
                            entry.database, entry.db_accession, entry.db_res_name,
                            entry.db_res_num, entry.comment
                        )
                    ))

        if seq_dict != {}:
            for chain_id, seq in seq_dict.items():
                if seq is not None:
                    sequence = seq
                else:
                    sequence = []
                    for residue in model[chain_id]:
                        sequence.append(residue.name)
                        if residue.end_of_chain:
                            break
                text = gen_seqres(chain_id, sequence)
                self.__validation_info["num_seq"] += len(text)
                f.write(warp_lines(text))

        if modres_dict != {}:
            for chain_id, modres in modres_dict.items():
                if modres is None:
                    continue
                for entry in modres:
                    res_name = model[chain_id][f"{entry.seq_num}{entry.icode:>1}"].name
                    f.write(warp_lines(
                        gen_modres(
                            pdb_id,
                            res_name, chain_id, entry.seq_num, entry.icode,
                            entry.std_res,
                            entry.comment
                        )
                    ))

    def __write_heterogen_section(self, pdb_data, f):
        if len(pdb_data.models) > 0:
            for chain in pdb_data.models[0]:
                for residue in chain.residues:
                    if residue.het:
                        het_id = residue.name
                        if het_id == "HOH" or residue.solvent:
                            continue
                        het = pdb_data.heterogen.get(het_id, None)
                        if het is None:
                            continue
                        text = gen_het(
                            het_id,
                            chain.id,
                            residue.id,
                            residue.icode,
                            len(residue),
                            het.comment,
                        )
                        self.__validation_info["num_het"] += len(text)
                        f.write(warp_lines(text))
        for het_id, het in pdb_data.heterogen.items():
            if getattr(het, "name", ""):
                f.write(warp_lines(gen_hetnam(het_id, het.name)))
            alias_src = getattr(het, "_alias_text", None) or getattr(
                het, "alias", None)
            if alias_src:
                if isinstance(alias_src, list):
                    hetsyn_text = "; ".join(alias_src)
                else:
                    hetsyn_text = str(alias_src)
                f.write(warp_lines(gen_hetsyn(het_id, hetsyn_text)))
            if getattr(het, "formula", ""):
                comp_num = getattr(het, "formula_comp_num", 1)
                f.write(warp_lines(gen_formul(comp_num, het_id, het.formula)))

    def __write_secondary_structure_section(self, pdb_data, f):
        secondary = pdb_data.secondary_structure
        if secondary.helix:
            model = pdb_data.models[0]
            for serial, (hid, hx) in enumerate(secondary.helix.items(), start=1):
                length = abs(hx.end_seq_num - hx.init_seq_num) + 1
                text = gen_helix(
                    serial,
                    hid,
                    model[hx.chain_id][f"{hx.init_seq_num}{hx.init_icode:1}"].name,
                    hx.chain_id,
                    hx.init_seq_num,
                    hx.init_icode,
                    model[hx.chain_id][f"{hx.end_seq_num}{hx.end_icode:1}"].name,
                    hx.chain_id,
                    hx.end_seq_num,
                    hx.end_icode,
                    hx.helix_class,
                    hx.comment,
                    length,
                )
                self.__validation_info["num_helix"] += len(text)
                f.write(warp_lines(text))
        if secondary.sheet:
            model = pdb_data.models[0]
            for sid, strands in secondary.sheet.items():
                n = len(strands)
                for i, strand in enumerate(strands, start=1):
                    cur_rs = (
                        strand.cur_res_seq
                        if strand.cur_res_seq is not None else 0)
                    prev_rs = (
                        strand.prev_res_seq
                        if strand.prev_res_seq is not None else 0)
                    text = gen_sheet(
                        i,
                        sid,
                        n % 100,
                        model[strand.init_chain_id][f"{strand.init_seq_num}{strand.init_icode:1}"].name,
                        strand.init_chain_id,
                        strand.init_seq_num,
                        strand.init_icode,
                        model[strand.end_chain_id][f"{strand.end_res_seq}{strand.end_icode:1}"].name,
                        strand.end_chain_id,
                        strand.end_res_seq,
                        strand.end_icode,
                        strand.sense,
                        strand.cur_atom,
                        model[strand.cur_chain_id][f"{cur_rs}{strand.cur_icode:1}"].name if cur_rs != 0 else "",
                        strand.cur_chain_id,
                        cur_rs,
                        strand.cur_icode,
                        strand.prev_atom,
                        model[strand.prev_chain_id][f"{prev_rs}{strand.prev_icode:1}"].name if prev_rs != 0 else "",
                        strand.prev_chain_id,
                        prev_rs,
                        strand.prev_icode,
                    )
                    self.__validation_info["num_sheet"] += len(text)
                    f.write(warp_lines(text))

    def __write_connectivity_section(self, pdb_data, f):
        conn = pdb_data.connectivity
        for i, ssbond in enumerate(conn.ss_bond, start=1):
            if ssbond.distance == 0:
                residue1 = pdb_data.models[0][ssbond.chain_id_1][f"{ssbond.seq_num_1}{ssbond.icode_1:1}"]
                residue2 = pdb_data.models[0][ssbond.chain_id_2][f"{ssbond.seq_num_2}{ssbond.icode_2:1}"]
                coord_atom1 = residue1.atoms['SG '].coord
                coord_atom2 = residue2.atoms['SG '].coord
                distance = ((coord_atom1[0] - coord_atom2[0]) ** 2 + \
                    (coord_atom1[1] - coord_atom2[1]) ** 2 + \
                    (coord_atom1[2] - coord_atom2[2])**2) ** 0.5
            else:
                distance = ssbond.distance
            distance = distance / unit_config.conversion_factor
            f.write(warp_lines(gen_ssbond(
                i,
                ssbond.chain_id_1,
                ssbond.seq_num_1,
                ssbond.icode_1,
                ssbond.chain_id_2,
                ssbond.seq_num_2,
                ssbond.icode_2,
                ssbond.symmetry_operation_1,
                ssbond.symmetry_operation_2,
                ssbond.distance,
            )))
        for link in conn.link:
            if link.distance == 0:
                residue1 = pdb_data.models[0][link.chain_id_1][f"{link.seq_num_1}{link.icode_1:1}"]
                residue2 = pdb_data.models[0][link.chain_id_2][f"{link.seq_num_2}{link.icode_2:1}"]
                coord_atom1 = residue1.atoms['OG '].coord
                coord_atom2 = residue2.atoms['OG '].coord
                distance = ((coord_atom1[0] - coord_atom2[0]) ** 2 + \
                    (coord_atom1[1] - coord_atom2[1]) ** 2 + \
                    (coord_atom1[2] - coord_atom2[2])**2) ** 0.5
            else:
                distance = link.distance
            distance = distance / unit_config.conversion_factor
            f.write(warp_lines(gen_link(
                link.name_1,
                link.alt_loc_1,
                pdb_data.models[0][link.chain_id_1][f"{link.seq_num_1}{link.icode_1:1}"].name,
                link.chain_id_1,
                link.seq_num_1,
                link.icode_1,
                link.name_2,
                link.alt_loc_2,
                pdb_data.models[0][link.chain_id_2][f"{link.seq_num_2}{link.icode_2:1}"].name,
                link.chain_id_2,
                link.seq_num_2,
                link.icode_2,
                link.symmetry_operation_1,
                link.symmetry_operation_2,
                link.distance,
            )))
        for i, cispep in enumerate(conn.cis_peptide, start=1):
            f.write(warp_lines(gen_cispep(
                i,
                pdb_data.models[0][cispep.chain_id_1][f"{cispep.seq_num_1}{cispep.icode_1:1}"].name,
                cispep.chain_id_1,
                cispep.seq_num_1,
                cispep.icode_1,
                pdb_data.models[0][cispep.chain_id_2][f"{cispep.seq_num_2}{cispep.icode_2:1}"].name,
                cispep.chain_id_2,
                cispep.seq_num_2,
                cispep.icode_2,
                cispep.num_model,
                cispep.measure,
            )))

    def __write_misc_feature_section(self, pdb_data, f):
        for site_id, site_list in pdb_data.sites.items():
            residues = []
            for s in site_list:
                residues.append({
                    "res_name": pdb_data.models[0][s.chain_id][f"{s.seq_num}{s.icode:1}"].name,
                    "chain_id": s.chain_id,
                    "seq_num": s.seq_num,
                    "icode": s.icode,
                })
            text = gen_site(site_id, residues)
            self.__validation_info["num_site"] += len(text)
            f.write(warp_lines(text))

    def __write_crystallographic_section(self, pdb_data, f):
        crystal = pdb_data.crystallographic
        if crystal.fake_crystallographic:
            f.write(warp_lines(gen_cryst1(
                1,
                1,
                1,
                90,
                90,
                90,
                "P 1",
                1,
            )))
        elif (
            crystal.cell_lengths and len(crystal.cell_lengths) == 3
            and crystal.cell_angles and len(crystal.cell_angles) == 3
        ):
            text = gen_cryst1(
                crystal.cell_lengths[0] / unit_config.conversion_factor,
                crystal.cell_lengths[1] / unit_config.conversion_factor,
                crystal.cell_lengths[2] / unit_config.conversion_factor,
                crystal.cell_angles[0],
                crystal.cell_angles[1],
                crystal.cell_angles[2],
                crystal.space_group,
                crystal.z,
            )
            f.write(warp_lines(text))
        if crystal.origin_matrix:
            converted_origin_matrix = []
            for row in crystal.origin_matrix:
                converted_origin_matrix.append([
                    row[0],
                    row[1],
                    row[2],
                    row[3] / unit_config.conversion_factor,
                ])
            text = gen_matrix_rows("ORIGX", converted_origin_matrix)
            self.__validation_info["num_xform"] += len(text)
            f.write(warp_lines(text))
        if crystal.scale_matrix:
            converted_scale_matrix = []
            for row in crystal.scale_matrix:
                converted_scale_matrix.append([
                    row[0] * unit_config.conversion_factor,
                    row[1] * unit_config.conversion_factor,
                    row[2] * unit_config.conversion_factor,
                    row[3],
                ])
            text = gen_matrix_rows("SCALE", converted_scale_matrix)
            self.__validation_info["num_xform"] += len(text)
            f.write(warp_lines(text))
        if crystal.ncs_matrix:
            for i, ncs in enumerate(crystal.ncs_matrix, start=1):
                converted_matrix = []
                for row in ncs.matrix:
                    converted_matrix.append([
                        row[0],
                        row[1],
                        row[2],
                        row[3] / unit_config.conversion_factor,
                    ])
                text = gen_ncs_matrix(
                    i, converted_matrix, 1 if ncs.given else 0
                )
                self.__validation_info["num_xform"] += len(text)
                f.write(warp_lines(text))

    def __write_conect_section(self, pdb_data, f):
        if pdb_data.connectivity.connections is None:
            return
        connections_with_dup = {}
        for atom1, bonded_atoms in pdb_data.connectivity.connections.items():
            for atom2 in bonded_atoms:
                if atom1 not in connections_with_dup:
                    connections_with_dup[atom1] = []
                connections_with_dup[atom1].append(atom2)
                if atom2 not in connections_with_dup:
                    connections_with_dup[atom2] = []
                connections_with_dup[atom2].append(atom1)
        for atom1 in sorted(connections_with_dup.keys()):
            bonded_atoms = sorted(set(connections_with_dup[atom1]))
            text = gen_conect(atom1, bonded_atoms)
            self.__validation_info["num_conect"] += len(text)
            f.write(warp_lines(text))

    @singledispatchmethod
    def write(self, data):
        """ Writes the given data to the PDB file.
        :param data: Data to be written (PdbData, Model, or Chain)"""
        raise NotImplementedError(
            f"This method is not implemented for the given data type: {type(data)}.")

    @write.register(Model)
    def _(self, data):
        pdb_data = PdbData()
        pdb_data.add_model(data)
        meta = PdbMetaData()
        meta.title = "pdb"
        meta.remark = {"0": "Generated by pypdbio"}
        pdb_data.meta = meta
        self.write(pdb_data)

    @write.register(Chain)
    def _(self, data):
        model = Model()
        model.add_chain(data)
        self.write(model)

    def __write_coordinate_section(self, pdb_data, f):
        for model_id, model in enumerate(pdb_data.models):
            if len(pdb_data.models) > 1:
                f.write(warp_lines([f"MODEL     {model_id + 1:>4}"]))
            for chain in model.chains:
                chain_id = chain.id
                atom_id = 1
                for residue in chain.residues:
                    residue_id = residue.id
                    atom_typ = "HETATM" if residue.het else "ATOM"
                    for atom in residue.atoms:
                        f.write(warp_lines(
                            gen_atom(
                                atom_typ,
                                atom_id,
                                atom.name,
                                atom.alt_loc,
                                residue.name,
                                chain_id,
                                residue_id,
                                residue.icode,
                                atom.coord[0] / unit_config.conversion_factor,
                                atom.coord[1] / unit_config.conversion_factor,
                                atom.coord[2] / unit_config.conversion_factor,
                                atom.occupancy,
                                atom.temp_factor / unit_config.conversion_factor ** 2,
                                atom.element,
                                atom.charge,
                            )))
                        if atom.anisotropic_temp_factor:
                            f.write(warp_lines(gen_anisou(
                                atom_id,
                                atom.name,
                                atom.alt_loc,
                                residue.name,
                                chain_id,
                                residue_id,
                                residue.icode,
                                atom.anisotropic_temp_factor[0] / unit_config.conversion_factor ** 2,
                                atom.anisotropic_temp_factor[1] / unit_config.conversion_factor ** 2,
                                atom.anisotropic_temp_factor[2] / unit_config.conversion_factor ** 2,
                                atom.anisotropic_temp_factor[3] / unit_config.conversion_factor ** 2,
                                atom.anisotropic_temp_factor[4] / unit_config.conversion_factor ** 2,
                                atom.anisotropic_temp_factor[5] / unit_config.conversion_factor ** 2,
                                atom.element,
                                atom.charge,
                            )))
                        atom_id += 1
                        self.__validation_info["num_coord"] += 1
                    if residue.end_of_chain:
                        f.write(warp_lines(
                            gen_ter(
                                atom_id,
                                residue.name,
                                chain_id,
                                residue_id,
                                residue.icode,
                            )
                        ))
                        atom_id += 1
                        self.__validation_info["num_ter"] += 1
            if len(pdb_data.models) > 1:
                f.write(warp_lines(["ENDMDL"]))

    def __write_bookkeeping_section(self, f):
        f.write(warp_lines(gen_master(
            self.__validation_info["num_remark"],
            self.__validation_info["num_het"],
            self.__validation_info["num_helix"],
            self.__validation_info["num_sheet"],
            self.__validation_info["num_site"],
            self.__validation_info["num_xform"],
            self.__validation_info["num_coord"],
            self.__validation_info["num_ter"],
            self.__validation_info["num_conect"],
            self.__validation_info["num_seq"],
        )))
        f.write(warp_lines(["END"]))

    @write.register(PdbData)
    def _(self, data):
        with open(self.pdb_path, "w", encoding="utf-8") as f:
            self.__write_header_section(data, f)
            self.__write_primary_structure_section(data, f)
            self.__write_heterogen_section(data, f)
            self.__write_secondary_structure_section(data, f)
            self.__write_connectivity_section(data, f)
            self.__write_misc_feature_section(data, f)
            self.__write_crystallographic_section(data, f)
            self.__write_coordinate_section(data, f)
            self.__write_conect_section(data, f)
            self.__write_bookkeeping_section(f)
