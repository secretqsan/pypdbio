# pylint: disable=missing-module-docstring
# pylint: disable=wildcard-import
# pylint: disable=unused-wildcard-import
from functools import singledispatchmethod
from .writer_helper import *
from .models import (
    PdbData,
    PdbMetaData,
    Model,
    Chain,
    Site,
)
from .utils import chain_id_of_index


class PdbWriter:
    """
    PDB writer class.
    """

    def __init__(self, path):
        """ Initializes the PdbWriter with the given file path.
        :param path: Path to the PDB file to be written"""
        self.pdb_path = path

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
                meta.replace.ids or [],
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
                f.write(warp_lines(gen_remark(remark_id, remark_text)))

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
                f.write(warp_lines(gen_seqres(chain_id, sequence)))

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
                        het = pdb_data.heterogen[het_id]
                        f.write(warp_lines(gen_het(
                            het_id,
                            chain.id,
                            residue.id,
                            residue.icode,
                            len(residue),
                            het.comment,
                        )))
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
                f.write(warp_lines(gen_helix(
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
                )))
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
                    f.write(warp_lines(gen_sheet(
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
                    )))

    def __write_connectivity_section(self, pdb_data, f):
        conn = pdb_data.connectivity
        for i, ssbond in enumerate(conn.ss_bond or [], start=1):
            sym1 = (ssbond.symmetry_operation_1 or "1555")[:6].rjust(6)
            sym2 = (ssbond.symmetry_operation_2 or "1555")[:6].rjust(6)
            f.write(warp_lines(gen_ssbond(
                i,
                ssbond.chain_id_1,
                ssbond.seq_num_1,
                ssbond.icode_1,
                ssbond.chain_id_2,
                ssbond.seq_num_2,
                ssbond.icode_2,
                sym1,
                sym2,
                ssbond.distance,
            )))
        for link in conn.link or []:
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
        for i, cispep in enumerate(conn.cis_peptide or [], start=1):
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
            f.write(warp_lines(gen_site(site_id, residues)))

    def __write_crystallographic_section(self, pdb_data, f):
        crystal = pdb_data.crystallographic
        if (
            crystal.cell_lengths and len(crystal.cell_lengths) == 3
            and crystal.cell_angles and len(crystal.cell_angles) == 3
        ):
            cryst_lines = gen_cryst1(
                crystal.cell_lengths[0],
                crystal.cell_lengths[1],
                crystal.cell_lengths[2],
                crystal.cell_angles[0],
                crystal.cell_angles[1],
                crystal.cell_angles[2],
                crystal.space_group,
                crystal.z,
            )
        else:
            cryst_lines = []
        if cryst_lines:
            f.write(warp_lines(cryst_lines))
        if crystal.origin_matrix:
            f.write(warp_lines(gen_matrix_rows("ORIGX", crystal.origin_matrix)))
        if crystal.scale_matrix:
            f.write(warp_lines(gen_matrix_rows("SCALE", crystal.scale_matrix)))
        if crystal.ncs_matrix:
            for i, ncs in enumerate(crystal.ncs_matrix, start=1):
                f.write(warp_lines(gen_ncs_matrix(
                    i, ncs.matrix, 1 if ncs.given else 0)))

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
            f.write(warp_lines(gen_conect(atom1, bonded_atoms)))

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
        last_atom_info = {}
        for model_id, model in enumerate(pdb_data.models):
            if len(pdb_data.models) > 1:
                f.write(f"MODEL     {model_id + 1:>4}\n")
            residue_id = 1
            atom_id = 1
            for chain in model.chains:
                chain_id = getattr(chain, "id", None) or "A"
                for residue in chain.residues:
                    atom_typ = "HETATM" if residue.het else "ATOM"
                    for atom in residue.atoms:
                        atom_info = {
                            "type": atom_typ,
                            "atom_no": atom_id,
                            "atom_name": atom.name,
                            "alt_loc": atom.alt_loc,
                            "residue_name": residue.name,
                            "chain_id": chain_id,
                            "residue_id": residue_id,
                            "icode": residue.icode,
                            "coord_x": atom.coord[0],
                            "coord_y": atom.coord[1],
                            "coord_z": atom.coord[2],
                            "temp_factor": atom.temp_factor,
                            "element": atom.element,
                            "charge": atom.charge,
                            "occupancy": atom.occupancy,
                        }
                        if last_atom_info.get('type') == 'ATOM' and atom_info["type"] == 'HETATM':
                            f.write(warp_lines(gen_ter(
                                last_atom_info["atom_no"],
                                last_atom_info["residue_name"],
                                last_atom_info["chain_id"],
                                last_atom_info["residue_id"],
                            )))
                        f.write(warp_lines(
                            gen_atom(
                                atom_info["type"],
                                atom_info["atom_no"],
                                atom_info["atom_name"],
                                atom_info["alt_loc"],
                                atom_info["residue_name"],
                                atom_info["chain_id"],
                                atom_info["residue_id"],
                                atom_info["icode"],
                                atom_info["coord_x"],
                                atom_info["coord_y"],
                                atom_info["coord_z"],
                                atom_info["occupancy"],
                                atom_info["temp_factor"],
                                atom_info["element"],
                                atom_info["charge"],
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
                                atom.anisotropic_temp_factor[0],
                                atom.anisotropic_temp_factor[1],
                                atom.anisotropic_temp_factor[2],
                                atom.anisotropic_temp_factor[3],
                                atom.anisotropic_temp_factor[4],
                                atom.anisotropic_temp_factor[5],
                                atom.element,
                                atom.charge,
                            )))
                        last_atom_info = atom_info
                        atom_id += 1
                    residue_id += 1

                if last_atom_info.get("type") == "ATOM":
                    f.write(warp_lines(
                        gen_ter(
                            last_atom_info["atom_no"],
                            last_atom_info["residue_name"],
                            last_atom_info["chain_id"],
                            last_atom_info["residue_id"],
                        )
                    ))
            if len(pdb_data.models) > 1:
                f.write("ENDMDL\n")

    def __write_bookkeeping_section(self, pdb_data, f):
        vi = getattr(pdb_data, "_validation_info", None) or {}
        f.write(warp_lines(gen_master(
            vi.get("num_remark", 0),
            vi.get("num_het", 0),
            vi.get("num_helix", 0),
            vi.get("num_sheet", 0),
            vi.get("num_site", 0),
            vi.get("num_xform", 0),
            vi.get("num_coord", 0),
            vi.get("num_ter", 0),
            vi.get("num_conect", 0),
            vi.get("num_seq", 0),
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
            self.__write_bookkeeping_section(data, f)
