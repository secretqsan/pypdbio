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
from .unit import unit_config
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
                        meta.obsolete.new_entry_id,
                        meta.obsolete.old_entry_id,
                        meta.obsolete.pdb_id
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
            f.write(warp_lines(gen_sprsde(meta.replace, meta.pdb_id)))
        if meta.journal:
            if len(meta.journal.author) > 0:
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
        for het_id, het in pdb_data.heterogen.items():
            het_rec = {
                "het_id": het_id,
                "chain_id": getattr(het, "chain_id", " ") or " ",
                "seq_num": getattr(het, "seq_num", None),
                "icode": getattr(het, "icode", "") or "",
                "num_het_atoms": getattr(het, "num_het_atoms", None),
                "text": getattr(het, "comment", "") or "",
            }
            f.write(warp_lines(gen_het(
                het_rec["het_id"],
                het_rec["chain_id"],
                het_rec["seq_num"],
                het_rec["icode"],
                het_rec["num_het_atoms"],
                het_rec["text"],
            )))
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
                formul_rec = {
                    "comp_num": getattr(het, "formula_comp_num", 1),
                    "het_id": het_id,
                    "continuation": getattr(het, "formula_continuation", 0),
                    "asterisk": getattr(het, "formula_asterisk", "") or " ",
                    "text": het.formula,
                }
                f.write(warp_lines(gen_formul(
                    formul_rec["comp_num"],
                    formul_rec["het_id"],
                    formul_rec["continuation"],
                    formul_rec["asterisk"],
                    formul_rec["text"],
                )))

    def __write_secondary_structure_section(self, pdb_data, f):
        secondary = pdb_data.secondary_structure
        if secondary.helix:
            for serial, (hid, hx) in enumerate(secondary.helix.items(), start=1):
                d = helix_to_writer_dict(hid, hx, serial)
                f.write(warp_lines(gen_helix(
                    d["ser_num"],
                    d["helix_id"],
                    d["init_res_name"],
                    d["init_chain_id"],
                    d["init_seq_num"],
                    d["init_icode"],
                    d["end_res_name"],
                    d["end_chain_id"],
                    d["end_seq_num"],
                    d["end_icode"],
                    d["helix_class"],
                    d["comment"],
                    d["length"],
                )))
        if secondary.sheet:
            for sid, strands in secondary.sheet.items():
                n = len(strands)
                for i, strand in enumerate(strands, start=1):
                    d = sheet_strand_to_writer_dict(sid, strand, i, n)
                    f.write(warp_lines(gen_sheet(
                        d["strand"],
                        d["sheet_id"],
                        d["num_strands"],
                        d["init_res_name"],
                        d["init_chain_id"],
                        d["init_seq_num"],
                        d["init_icode"],
                        d["end_res_name"],
                        d["end_chain_id"],
                        d["end_seq_num"],
                        d["end_icode"],
                        d["sense"],
                        d["cur_atom"],
                        d["cur_res_name"],
                        d["cur_chain_id"],
                        d["cur_res_seq"],
                        d["cur_icode"],
                        d["prev_atom"],
                        d["prev_res_name"],
                        d["prev_chain_id"],
                        d["prev_res_seq"],
                        d["prev_icode"],
                    )))

    def __write_connectivity_section(self, pdb_data, f):
        conn = pdb_data.connectivity
        for i, ssbond in enumerate(conn.ss_bond or [], start=1):
            d = ssbond_to_writer_dict(i, ssbond)
            f.write(warp_lines(gen_ssbond(
                d["ser_num"],
                d["chain_id1"],
                d["seq_num1"],
                d["icode1"],
                d["chain_id2"],
                d["seq_num2"],
                d["icode2"],
                d["sym1"],
                d["sym2"],
                d["length"],
            )))
        for link in conn.link or []:
            d = link_to_writer_dict(link)
            f.write(warp_lines(gen_link(
                d["name1"],
                d["alt_loc1"],
                d["res_name1"],
                d["chain_id1"],
                d["res_seq1"],
                d["icode1"],
                d["name2"],
                d["alt_loc2"],
                d["res_name2"],
                d["chain_id2"],
                d["res_seq2"],
                d["icode2"],
                d["sym1"],
                d["sym2"],
                d["length"],
            )))
        for i, cispep in enumerate(conn.cis_peptide or [], start=1):
            d = cispep_to_writer_dict(i, cispep)
            f.write(warp_lines(gen_cispep(
                d["ser_num"],
                d["pep1"],
                d["chain_id1"],
                d["seq_num1"],
                d["icode1"],
                d["pep2"],
                d["chain_id2"],
                d["seq_num2"],
                d["icode2"],
                d["mod_num"],
                d["measure"],
            )))

    def __write_misc_feature_section(self, pdb_data, f):
        for site_id, site_list in pdb_data.sites.items():
            residues = []
            for s in site_list:
                if isinstance(s, Site):
                    residues.append({
                        "res_name": "   ",
                        "chain_id": s.chain_id or " ",
                        "seq_num": s.seq_num,
                        "icode": (s.icode or " ")[:1],
                    })
                else:
                    residues.append({
                        "res_name": s.get("res_name", "   "),
                        "chain_id": s.get("chain_id", " "),
                        "seq_num": s.get("seq_num"),
                        "icode": s.get("icode", " "),
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
            f.write(warp_lines(gen_scale_matrix_rows(
                "SCALE", crystal.scale_matrix)))
        if crystal.ncs_matrix:
            f.write(warp_lines(gen_ncs_matrix(crystal.ncs_matrix)))

    def __write_conect_section(self, pdb_data, f):
        conn = pdb_data.connectivity
        if conn.connections:
            f.write(warp_lines(gen_conect(conn.connections)))

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

    def _gen_ter_line(self, last_atom_info):
        line = f'TER   {last_atom_info["atom_no"]:>5}      '
        line += f'{last_atom_info["residue_name"]:>3} {last_atom_info["chain_id"]}'
        line += f'{last_atom_info["residue_id"]:>4} \n'
        return line

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
                    atom_typ = "HETATM" if getattr(
                        residue, "het", "") == "het" else "ATOM"
                    for atom in residue.atoms:
                        atom_info = {
                            "type": atom_typ,
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
                            "charge": atom.charge,
                            "occupancy": getattr(atom, "occupancy", 1.0),
                        }
                        if last_atom_info.get('type') == 'ATOM' and atom_info["type"] == 'HETATM':
                            f.write(self._gen_ter_line(last_atom_info))
                        f.write(self._gen_atom_line(atom_info))
                        last_atom_info = atom_info
                        atom_id += 1
                    residue_id += 1
                if last_atom_info.get("type") == "ATOM":
                    f.write(self._gen_ter_line(last_atom_info))
            if len(pdb_data.models) > 1:
                f.write("ENDMDL\n")

    def __write_bookkeeping_section(self, pdb_data, f):
        f.write(warp_lines(gen_master(getattr(pdb_data, "_validation_info", None))))
        f.write(warp_lines(["END"]))

    @write.register(PdbData)
    def _(self, data):
        with open(self.pdb_path, "w", encoding="utf-8") as f:
            self.__write_header_section(data, f)
            self.__write_primary_structure_section(data, f)
            # self.__write_heterogen_section(data, f)
            # self.__write_secondary_structure_section(data, f)
            # self.__write_connectivity_section(data, f)
            # self.__write_misc_feature_section(data, f)
            # self.__write_crystallographic_section(data, f)
            # self.__write_coordinate_section(data, f)
            # self.__write_conect_section(data, f)
            # self.__write_bookkeeping_section(data, f)
