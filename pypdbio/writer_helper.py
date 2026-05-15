# pylint: disable=missing-function-docstring
# pylint: disable=missing-module-docstring
# pylint: disable=wildcard-import
# pylint: disable=unused-wildcard-import
import warnings
from .models import *

# helper functions


def gen_fixed_length_list(lst, n):
    return [lst[i: i + n] for i in range(0, len(lst), n)]


def gen_fixed_length_list_str(lst, width, n):
    target_list = gen_fixed_length_list(lst, n)
    result = [
        "".join(f"{item.strip():<{width}}" for item in item_chunk)
        for item_chunk
        in target_list
    ]
    return result if result else [""]


def split_text(text, width=80, prefix_space=True, split_words=False):
    first_line = True
    lines = []
    if isinstance(text, list):
        return [
            f"{' ' if index != 0 and prefix_space else ''}{line}"
            for index, line
            in enumerate(text)
        ]
    if split_words:
        def max_length():
            if not first_line and prefix_space:
                return width - 1
            else:
                return width
        while len(text) > max_length():
            lines.append(
                f"{' ' if prefix_space and not first_line else ''}{text[:max_length()]}"
            )
            text = text[max_length():]
            first_line = False
        lines.append(text)
    else:
        current_line = ""
        words = text.split()
        while len(words) > 0:
            word = words.pop(0)
            len_word = len(word)
            len_current_line = len(current_line)
            if len_word > width:
                lines.append(current_line)
                lines.append(word)
                first_line = False
                current_line = ""
                warnings.warn(
                    f"The word '{word}' in '{text}' is too long to be split")
            elif len_current_line + len_word + 1 > width:
                lines.append(current_line)
                first_line = False
                current_line = f"{' ' if prefix_space else ''}{word}"
            else:
                current_line += f"{' ' if current_line else ''}{word}"
        if current_line:
            lines.append(current_line)
    return lines


def warp_lines(lines):
    return "".join([f"{line:<80}\n" for line in lines])


def gen_line_from_template(template, fields_template):
    """
    Generates a line from a template and a dictionary of fields.
    :param template: The template to use.
    :param fields_template: The dictionary of fields to use.
    :return: The generated line.
    """
    split_field = ""
    splitted_values = []
    for field, value in fields_template.items():
        if value.get("split", False):
            split_field = field
            splitted_values = split_text(
                value["value"],
                **value.get("kwargs", {})
            )
    if split_field == "":
        splitted_values = [""]  # make it a list with one element
    results = []
    for i, sub_str in enumerate(splitted_values):
        fields = {}
        for field, value in fields_template.items():
            if value.get("no_repetition", False) and i > 0:
                fields[field] = ""
            elif value.get("continuous_id", False):
                fields[field] = "" if i == 0 and value.get(
                    "hide_first", True) else i + 1
            elif field == split_field:
                fields[field] = sub_str
            else:
                fields[field] = value["value"]
        results.append(template.format(**fields))
    return results

# title section generation functions


def gen_header(classification, date, pdb_id):
    return gen_line_from_template(
        "HEADER    {classification:<40}{date:>9}   {pdb_id:>4}",
        {
            "classification": {"value": classification},
            "date": {"value": date},
            "pdb_id": {"value": pdb_id},
        }
    )


def gen_obslte(replace_date, new_entry_id, pdb_id):
    return gen_line_from_template(
        "OBSLTE     {replace_date:>9} {pdb_id:>4}      {new_entry_id:>4}",
        {
            "replace_date": {"value": replace_date},
            "pdb_id": {"value": pdb_id},
            "new_entry_id": {"value": new_entry_id},
        }
    )


def gen_title(title):
    return gen_line_from_template(
        "TITLE   {continuation:>2}{title:<70}",
        {
            "continuation": {"continuous_id": True},
            "title": {
                "value": title,
                "split": True,
                "kwargs": {
                    "width": 70,
                }
            },
        },
    )


def gen_split(split_ids):
    text = gen_fixed_length_list_str(split_ids, 5, 14)
    return gen_line_from_template(
        "SPLIT   {continuation:>2} {split_ids:<69}",
        {
            "continuation": {"continuous_id": True},
            "split_ids": {"value": text, "split": True},
        },
    )


def gen_caveat(pdb_id, caveat):
    return gen_line_from_template(
        "CAVEAT  {continuation:>2} {pdb_id:<4}    {caveat:<60}",
        {
            "continuation": {"continuous_id": True},
            "pdb_id": {"value": pdb_id},
            "caveat": {
                "value": caveat,
                "split": True,
                "kwargs": {
                    "width": 60,
                },
            },
        },
    )


def gen_compnd(texts):
    return gen_line_from_template(
        "COMPND {continuation:>3}{compound:<70}",
        {
            "continuation": {"continuous_id": True},
            "compound": {"value": texts, "split": True},
        },
    )


def gen_source(texts):
    return gen_line_from_template(
        "SOURCE {continuation:>3}{source:<69}",
        {
            "continuation": {"continuous_id": True},
            "source": {"value": texts, "split": True},
        },
    )


def gen_keywds(keywords):
    return gen_line_from_template(
        "KEYWDS  {continuation:>2}{keywords:<69}",
        {
            "continuation": {"continuous_id": True},
            "keywords": {
                "value": ", ".join(keywords),
                "split": True,
                "kwargs": {"width": 69},
            },
        },
    )


def gen_expdta(experiments):
    return gen_line_from_template(
        "EXPDTA  {continuation:>2}{experiments:<69}",
        {
            "continuation": {"continuous_id": True},
            "experiments": {
                "value": "; ".join(experiments),
                "split": True,
                "kwargs": {"width": 69},
            },
        },
    )


def gen_nummdl(n_models):
    return [f"NUMMDL    {n_models:<4}"]


def gen_mdltyp(model_type):
    text = model_type.strip()
    return gen_line_from_template(
        "MDLTYP  {continuation:>2}{mdltyp:<70}",
        {
            "continuation": {"continuous_id": True},
            "mdltyp": {
                "value": text,
                "split": True,
                "kwargs": {"width": 70},
            },
        },
    )


def gen_author(author):
    return gen_line_from_template(
        "AUTHOR  {continuation:>2}{author:<69}",
        {
            "continuation": {"continuous_id": True},
            "author": {
                "value": ", ".join(author),
                "split": True,
                "kwargs": {"width": 69},
            },
        },
    )


def gen_revdat(serial, date, pdb_id, mod_type, mods):
    return gen_line_from_template(
        "REVDAT {serial:>3}{continuation:>2} {date:>9} {pdb_id:>4}    {mod_type}       {mods}",
        {
            "serial": {"value": serial},
            "continuation": {"continuous_id": True},
            "date": {"value": date},
            "pdb_id": {"value": pdb_id},
            "mod_type": {"value": mod_type},
            "mods": {"value": gen_fixed_length_list_str(mods, 7, 4)[0]},
        },
    )


def gen_sprsde(date, pdb_id, replace_pdb_id_list):
    ids_str_list = gen_fixed_length_list_str(replace_pdb_id_list, 5, 9)
    return gen_line_from_template(
        "SPRSDE  {continuation:>2} {date:>9} {pdb_id:>4}      {ids}",
        {
            "continuation": {"continuous_id": True},
            "date": {"value": date, "no_repetition": True},
            "pdb_id": {"value": pdb_id.strip(), "no_repetition": True},
            "ids": {
                "value": ids_str_list,
                "split": True,
                "kwargs": {
                    "prefix_space": False,
                },
            },
        },
    )


def gen_jrnl_auth(author):
    return gen_line_from_template(
        "JRNL        AUTH{continuation:>2} {auth:<60}",
        {
            "continuation": {"continuous_id": True},
            "auth": {
                "value": ", ".join(author),
                "split": True,
                "kwargs": {
                    "width": 60,
                    "prefix_space": False,
                },
            },
        },
    )


def gen_jrnl_titl(title):
    return gen_line_from_template(
        "JRNL        TITL{continuation:>2} {titl:<60}",
        {
            "continuation": {"continuous_id": True},
            "titl": {
                "value": title,
                "split": True,
                "kwargs": {
                    "width": 60,
                    "prefix_space": False,
                },
            },
        },
    )


def gen_jrnl_editor(editor):
    return gen_line_from_template(
        "JRNL        EDIT{continuation:>2} {editor:<60}",
        {
            "continuation": {"continuous_id": True},
            "editor": {
                "value": ", ".join(editor),
                "split": True,
                "kwargs": {
                    "width": 60,
                    "prefix_space": False,
                },
            },
        },
    )


def gen_jrnl_ref_to_be_published():
    return ["JRNL        REF    TO BE PUBLISHED"]


def gen_jrnl_ref(pub_name, volume, pages, year):
    return gen_line_from_template(
        "JRNL        REF {continuation:>2} {pub_name:<28}  {v:<2}{volume:>4} {pages:>5} {year:<4}",
        {
            "continuation": {"continuous_id": True},
            "pub_name": {
                "value": pub_name,
                "split": True,
                "kwargs": {
                    "width": 28,
                    "prefix_space": False,
                },
            },
            "v": {
                "value": "V." if volume else "",
                "repetition": False,
            },
            "volume": {
                "value": volume,
                "repetition": False,
            },
            "pages": {
                "value": pages,
                "repetition": False,
            },
            "year": {
                "value": year,
                "repetition": False,
            },
        }
    )


def gen_jrnl_publisher(publisher):
    return gen_line_from_template(
        "JRNL        PUBL{continuation:>2} {publisher:<51}",
        {
            "continuation": {"continuous_id": True},
            "publisher": {
                "value": publisher,
                "split": True,
                "kwargs": {
                    "width": 51,
                    "prefix_space": False,
                },
            },
        },
    )


def gen_jrnl_refn(target, value):
    return gen_line_from_template(
        "JRNL        REFN                   {target:<4} {value:<25}",
        {
            "target": {"value": target},
            "value": {
                "value": value,
            },
        },
    )


def gen_jrnl_pmid(pmid):
    return gen_line_from_template(
        "JRNL        PMID{continuation:>2} {pmid:<60}",
        {
            "continuation": {"continuous_id": True},
            "pmid": {
                "value": pmid.strip(),
                "split": True,
                "kwargs": {
                    "width": 60,
                    "prefix_space": False,
                    "split_words": True,
                },
            },
        },
    )


def gen_jrnl_doi(doi):
    return gen_line_from_template(
        "JRNL        DOI {continuation:>2} {doi:<60}",
        {
            "continuation": {"continuous_id": True},
            "doi": {
                "value": doi.strip(),
                "split": True,
                "kwargs": {
                    "width": 60,
                    "prefix_space": False,
                    "split_words": True,
                },
            },
        },
    )


def gen_remark(remark_id, remark_text):
    out = []
    for line in remark_text.split("\n"):
        new_line = " " + line
        if len(new_line) > 69:
            warnings.warn(f"The remark text is too long: {line}")
        out.extend(gen_line_from_template(
            "REMARK {remark_id:>3}{remark:<69}",
            {
                "remark_id": {"value": remark_id},
                "remark": {"value": new_line},
            },
        ))
    return out


def gen_dbref(
    id_code, chain_id, seq_begin, insert_begin, seq_end, insert_end,
    database, db_accession, db_id_code, db_seq_begin, db_ins_begin,
    db_seq_end, db_ins_end
):
    if database == "GB" or database == "UNIMES":
        line1 = f"DBREF1 {id_code:>4} {chain_id} {seq_begin:>4}{insert_begin:>1} " + \
            f"{seq_end:>4}{insert_end:>1} {database:<6}               {db_accession:<20}"
        line2 = f"DBREF2 {id_code:>4} {chain_id}     {db_id_code:<22}     " + \
            f"{db_seq_begin:>10}  {db_seq_end:>10}"
        return [line1, line2]
    return [
        f"DBREF  {id_code:>4} {chain_id} {seq_begin:>4}{insert_begin:>1} " +
        f"{seq_end:>4}{insert_end:>1} {database:<6} {db_accession:<8} " +
        f"{db_id_code:<12} {db_seq_begin:>5}{db_ins_begin:>1} " +
        f"{db_seq_end:>5}{db_ins_end:>1}"
    ]


def gen_seqadv(
    id_code, res_name, chain_id, seq_num, icode,
    database, db_accession, db_res, db_seq, conflict
):
    return [
        f"SEQADV {id_code:>4} {res_name:<3} {chain_id:1} {seq_num:>4}{icode:>1} " +
        f"{database:<4} {db_accession:<9} {db_res:<3} " +
        f"{" " if db_seq is None else db_seq:>5} {conflict:<21}"
    ]


def gen_seqres(chain_id, residues):
    seq_str = gen_fixed_length_list_str(residues, 4, 13)
    return gen_line_from_template(
        "SEQRES {index:>3} {chain_id:1} {num_res:>4}  {residues:<52}",
        {
            "index": {"continuous_id": True, "hide_first": False},
            "chain_id": {"value": chain_id},
            "num_res": {"value": len(residues)},
            "residues": {
                "value": seq_str,
                "split": True,
                "kwargs": {"prefix_space": False}
            }
        },
    )


def gen_modres(id_code, res_name, chain_id, seq_num, icode, std_res, comment):
    return [
        f"MODRES {id_code:>4} {res_name:<3} {chain_id:1} {seq_num:>4}{icode:1} " +
        f"{std_res:<3}  {comment:<41}"
    ]


def gen_het(het_id, chain_id, seq_num, icode, num_het_atoms, text):
    return [
        f"HET    {het_id:<3}  {chain_id:1}{seq_num:>4}{icode:1}  {num_het_atoms:>5}     {text:<40}"
    ]


def gen_hetnam(het_id, text):
    return gen_line_from_template(
        "HETNAM  {continuation:>2} {het_id:<3} {text:<55}",
        {
            "continuation": {"continuous_id": True},
            "het_id": {"value": het_id, "no_repetition": True},
            "text": {"value": text, "split": True, "kwargs": {"width": 55}},
        },
    )


def gen_hetsyn(het_id, text):
    return gen_line_from_template(
        "HETSYN  {continuation:>2} {het_id:<3} {text:<55}",
        {
            "continuation": {"continuous_id": True},
            "het_id": {"value": het_id, "no_repetition": True},
            "text": {"value": text, "split": True, "kwargs": {"width": 55}},
        },
    )


def gen_formul(comp_num, het_id, text):
    return gen_line_from_template(
        "FORMUL  {comp_num:>2}  {het_id:>3} {continuation:>2}{asterisk:1}{text:<51}",
        {
            "comp_num": {"value": comp_num},
            "het_id": {"value": het_id},
            "continuation": {"continuous_id": True},
            "asterisk": {"value": "*" if het_id == "HOH" else ""},
            "text": {"value": text, "split": True, "kwargs": {"width": 51}},
        },
    )


def gen_helix(
    ser_num, helix_id, init_res_name, init_chain_id, init_seq_num, init_icode,
    end_res_name, end_chain_id, end_seq_num, end_icode, helix_class, comment, length
):
    return [
        f"HELIX  {ser_num:>3} {helix_id:>3} " +
        f"{init_res_name:<3} {init_chain_id:1} {init_seq_num:>4}{init_icode:1} " +
        f"{end_res_name:<3} {end_chain_id:1} {end_seq_num:>4}{end_icode:1}" +
        f"{helix_class:>2}{comment:<30} {length:>5}"
    ]


def gen_atom_name(atom_name):
    idx = atom_name[0]
    if idx.isdigit():
        return atom_name
    return f" {atom_name}"


def gen_sheet(
    strand, sheet_id, num_strands, init_res_name, init_chain_id, init_seq_num,
    init_icode, end_res_name, end_chain_id, end_seq_num, end_icode, sense,
    cur_atom, cur_res_name, cur_chain_id, cur_res_seq, cur_icode,
    prev_atom, prev_res_name, prev_chain_id, prev_res_seq, prev_icode
):
    prefix = f"SHEET  {strand:>3} {sheet_id:>3}{num_strands:>2} " + \
        f"{init_res_name:<3} {init_chain_id:1}{init_seq_num:>4}{init_icode:1} " + \
        f"{end_res_name:<3} {end_chain_id:1}{end_seq_num:>4}{end_icode:1}" + \
        f"{sense:>2} "
    if sense == 0:
        registration = ""
    else:
        registration = f"{gen_atom_name(cur_atom):<4}{cur_res_name:>3} " + \
            f"{cur_chain_id:1}{cur_res_seq:>4}{cur_icode:1} " + \
            f"{gen_atom_name(prev_atom):<4}{prev_res_name:>3} " + \
            f"{prev_chain_id:1}{prev_res_seq:>4}{prev_icode:1}"
    return [
        prefix + registration
    ]


def gen_ssbond(
    ser_num, chain_id1, seq_num1, icode1, chain_id2, seq_num2, icode2,
    sym1, sym2, length
):
    return [
        f"SSBOND {ser_num:>3} CYS {chain_id1:1} {seq_num1:>4}{icode1:1}   " +
        f"CYS {chain_id2:1} {seq_num2:>4}{icode2:1}                       " +
        f"{sym1:>6} {sym2:>6} {length:5.2f}"
    ]


def gen_link(
    name1, alt_loc1, res_name1, chain_id1, res_seq1, icode1, name2, alt_loc2,
    res_name2, chain_id2, res_seq2, icode2, sym1, sym2, length
):
    return [
        f"LINK        {gen_atom_name(name1):<4}{alt_loc1:1}" +
        f"{res_name1:>3} {chain_id1:1}{res_seq1:>4}{icode1:1}               " +
        f"{gen_atom_name(name2):<4}{alt_loc2:1}" +
        f"{res_name2:>3} {chain_id2:1}{res_seq2:>4}{icode2:1}  " +
        f"{sym1:>6} {sym2:>6} {length:5.2f}"
    ]


def gen_cispep(
    ser_num, pep1, chain_id1, seq_num1, icode1,
    pep2, chain_id2, seq_num2, icode2, mod_num, measure
):
    return [
        f"CISPEP {ser_num:>3} {pep1:<3} {chain_id1:1} {seq_num1:>4}{icode1:1}   " +
        f"{pep2:<3} {chain_id2:1} {seq_num2:>4}{icode2:1}       {mod_num:>3}       {measure:6.2f}"
    ]


def gen_site(site_id, residues):
    texts = [
        f"{residue.get('res_name', ''):>3} {residue.get('chain_id', ' '):1}{residue.get('seq_num'):>4} {residue.get('icode', ' '):1} "
        for residue
        in residues
    ]
    texts = gen_fixed_length_list_str(texts, 11, 4)
    return gen_line_from_template(
        "SITE   {index:>3} {site_id:<3} {num_res:>2} {texts:<44}",
        {
            "index": {"continuous_id": True, "hide_first": False},
            "site_id": {"value": site_id},
            "num_res": {"value": len(residues)},
            "texts": {"value": texts, "split": True, "kwargs": {"prefix_space": False}},
        },
    )


def gen_cryst1(a, b, c, alpha, beta, gamma, s_group, z):
    return [
        f"CRYST1{a:>9.3f}{b:>9.3f}{c:>9.3f}{alpha:>7.2f}{beta:>7.2f}{gamma:>7.2f} " + \
        f"{s_group:<11}{z:>4}"
    ]


def gen_matrix_rows(prefix, matrix):
    rows = [
        f"{row[0]:10.6f}{row[1]:10.6f}{row[2]:10.6f}     {row[3]:10.5f}"
        for row in matrix
    ]
    return gen_line_from_template(
        "{prefix}{n}    {rows:<44}",
        {
            "prefix": {"value": prefix},
            "n": {"continuous_id": True, "hide_first": False},
            "rows": {"value": rows, "split": True, "kwargs": {"prefix_space": False}},
        },
    )


def gen_ncs_matrix(serial, ncs_matrix, igiven):
    rows = [
        f"{row[0]:10.6f}{row[1]:10.6f}{row[2]:10.6f}     {row[3]:10.5f}    {igiven:1}"
        for row in ncs_matrix
    ]
    return gen_line_from_template(
        "MTRIX{n} {serial:>3}{rows:<44}",
        {
            "n": {"continuous_id": True, "hide_first": False},
            "serial": {"value": serial},
            "rows": {"value": rows, "split": True, "kwargs": {"prefix_space": False}},
        },
    )


def gen_master(
    num_remark, num_het, num_helix, num_sheet,
    num_site, num_xform, num_coord, num_ter,
    num_conect, num_seq
):
    return [
        f"MASTER    {num_remark:>5}    0{num_het:>5}{num_helix:>5}" +
        f"{num_sheet:>5}    0{num_site:>5}{num_xform:>5}" +
        f"{num_coord:>5}{num_ter:>5}{num_conect:>5}{num_seq:>5}"
    ]


def gen_conect(atom1, bonded_atoms):
    line = f"CONECT{atom1:>5}"
    for atom2 in bonded_atoms:
        line += f"{atom2:>5}"
    return [line]

def gen_charge(charge):
    if charge == 1:
        return "1+"
    elif charge == -1:
        return "1-"
    else:
        return ""

def gen_atom(
    atom_type, atom_no, atom_name, alt_loc, residue_name,
    chain_id, residue_id, icode, coord_x, coord_y, coord_z,
    occupancy, temp_factor, element, charge
):
    line = f'{atom_type:<6}{atom_no:>5} {gen_atom_name(atom_name):<4}'
    line += f'{alt_loc:1}{residue_name:>3} {chain_id}'
    cx = coord_x
    cy = coord_y
    cz = coord_z
    line += f"{residue_id:>4}{icode:1}   {cx:8.3f}"
    line += f"{cy:8.3f}{cz:8.3f}"
    occ = occupancy
    line += f"{occ:6.2f}{temp_factor:6.2f}"
    line += f"          {element:>2}"
    charge_str = gen_charge(charge)
    line += f"{charge_str:>2}"
    return [line]

def gen_anisou(
    atom_no,atom_name, alt_loc, residue_name, chain_id, residue_id, icode, u11, u22, u33, u12, u13, u23, element, charge
):
    return [
        f"ANISOU{atom_no:>5} {gen_atom_name(atom_name):<4}{alt_loc:1}" + \
        f"{residue_name:>3} {chain_id}{residue_id:>4}{icode:1} " + \
        f"{int(u11 * 10000):>7}{int(u22 * 10000):>7}{int(u33 * 10000):>7}" + \
        f"{int(u12 * 10000):>7}{int(u13 * 10000):>7}{int(u23 * 10000):>7}      " + \
        f"{element:>2}{gen_charge(charge):>2}"
    ]

def gen_ter(atom_no, residue_name, chain_id, residue_id, icode):
    return [f'TER   {atom_no:>5}      {residue_name:>3} {chain_id}{residue_id:>4}{icode:1}']
