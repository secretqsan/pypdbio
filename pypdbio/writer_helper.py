# pylint: disable=missing-function-docstring
# pylint: disable=missing-module-docstring
# pylint: disable=wildcard-import
# pylint: disable=unused-wildcard-import
import warnings
from .models import *

# helper functions


def gen_fixed_length_list(lst, width, n):
    target_list = [lst[i: i + n] for i in range(0, len(lst), n)]
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
                fields[field] = "" if i == 0 and value.get("hide_first", True) else i + 1
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
    text = gen_fixed_length_list(split_ids, 5, 14)
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
            "mods": {"value": gen_fixed_length_list(mods, 7, 4)[0]},
        },
    )


def gen_sprsde(date, pdb_id, replace_pdb_id_list):
    ids_str_list = gen_fixed_length_list(replace_pdb_id_list, 5, 9)
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
        f"DBREF  {id_code:>4} {chain_id} {seq_begin:>4}{insert_begin:>1} " + \
        f"{seq_end:>4}{insert_end:>1} {database:<6} {db_accession:<8} " + \
        f"{db_id_code:<12} {db_seq_begin:>5}{db_ins_begin:>1} " + \
        f"{db_seq_end:>5}{db_ins_end:>1}"
    ]

def gen_seqadv(
    id_code, res_name, chain_id, seq_num, icode,
    database, db_accession, db_res, db_seq, conflict
):
    return [
        f"SEQADV {id_code:>4} {res_name:<3} {chain_id:1} {seq_num:>4}{icode:>1} " + \
        f"{database:<4} {db_accession:<9} {db_res:<3} " + \
        f"{" " if db_seq is None else db_seq:>5} {conflict:<21}"
    ]


def gen_seqres(chain_id, residues):
    seq_str = gen_fixed_length_list(residues, 4, 13)
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
        f"MODRES {id_code:>4} {res_name:<3} {chain_id:1} {seq_num:>4}{icode:1} " + \
        f"{std_res:<3}  {comment:<41}"
    ]


def gen_het(het_id, chain_id, seq_num, icode, num_het_atoms, text):
    return [
        f"HET   {het_id:>3} {chain_id:1} {seq_num:>4}{icode:1} {num_het_atoms:>5} {text:<40}"
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


def gen_formul(comp_num, het_id, continuation, asterisk, text):
    return [
        f"FORMUL {comp_num:>2} {het_id:>3} {continuation:>2}{asterisk:1} {text:<51}"
    ]


def gen_helix(
    ser_num, helix_id, init_res_name, init_chain_id, init_seq_num, init_icode,
    end_res_name, end_chain_id, end_seq_num, end_icode, helix_class, comment, length
):
    return [
        (
            "HELIX "
            f"{ser_num:>3}"
            f" {helix_id:<3}"
            f" {init_res_name:<3}"
            f" {init_chain_id:1}"
            f"{init_seq_num:>4}{init_icode:1}"
            f" {end_res_name:<3}"
            f" {end_chain_id:1}"
            f"{end_seq_num:>4}{end_icode:1}"
            f"{helix_class:>2}"
            f"{comment:<30}"
            f" {length:>5}"
        )
    ]


def gen_sheet(
    strand, sheet_id, num_strands, init_res_name, init_chain_id, init_seq_num,
    init_icode, end_res_name, end_chain_id, end_seq_num, end_icode, sense,
    cur_atom, cur_res_name, cur_chain_id, cur_res_seq, cur_icode,
    prev_atom, prev_res_name, prev_chain_id, prev_res_seq, prev_icode
):
    return [
        f"SHEET {strand:>3} {sheet_id:<3} {num_strands:>2} {init_res_name:<3} {init_chain_id:1} {init_seq_num:>4}{init_icode:1} {end_res_name:<3} {end_chain_id:1} {end_seq_num:>4}{end_icode:1} {sense:>2} {cur_atom:<4} {cur_res_name:>3} {cur_chain_id:1} {cur_res_seq:>4}{cur_icode:1} {prev_atom:<4} {prev_res_name:>3} {prev_chain_id:1} {prev_res_seq:>4}{prev_icode:1}"
    ]


def helix_to_writer_dict(helix_id: str, hx: Helix, serial: int) -> dict:
    length = abs((hx.end_seq_num or 0) - (hx.init_seq_num or 0)) + 1
    return {
        "ser_num": serial,
        "helix_id": helix_id[:3].ljust(3),
        "init_res_name": "   ",
        "init_chain_id": hx.chain_id or " ",
        "init_seq_num": hx.init_seq_num,
        "init_icode": (hx.init_icode or " ")[:1],
        "end_res_name": "   ",
        "end_chain_id": hx.chain_id or " ",
        "end_seq_num": hx.end_seq_num,
        "end_icode": (hx.end_icode or " ")[:1],
        "helix_class": hx.helix_class,
        "comment": (hx.comment or "")[:30],
        "length": length,
    }


def sheet_strand_to_writer_dict(sheet_id: str, strand: SheetStrand,
                                strand_idx: int, num_strands: int) -> dict:
    def rs(res_seq):
        return res_seq if res_seq is not None else 0

    return {
        "strand": strand_idx,
        "sheet_id": sheet_id[:3].ljust(3),
        "num_strands": num_strands,
        "init_res_name": "   ",
        "init_chain_id": strand.init_chain_id or " ",
        "init_seq_num": strand.init_seq_num,
        "init_icode": (strand.init_icode or " ")[:1],
        "end_res_name": "   ",
        "end_chain_id": strand.end_chain_id or " ",
        "end_seq_num": strand.end_res_seq,
        "end_icode": (strand.end_icode or " ")[:1],
        "sense": strand.sense,
        "cur_atom": (strand.cur_atom or " ")[:4].ljust(4),
        "cur_res_name": "   ",
        "cur_chain_id": strand.cur_chain_id or " ",
        "cur_res_seq": rs(strand.cur_res_seq),
        "cur_icode": (strand.cur_icode or " ")[:1],
        "prev_atom": (strand.prev_atom or " ")[:4].ljust(4),
        "prev_res_name": "   ",
        "prev_chain_id": strand.prev_chain_id or " ",
        "prev_res_seq": rs(strand.prev_res_seq),
        "prev_icode": (strand.prev_icode or " ")[:1],
    }


def gen_ssbond(
    ser_num, chain_id1, seq_num1, icode1, chain_id2, seq_num2, icode2,
    sym1, sym2, length
):
    length_str = f"{length:5.2f}" if length is not None else " " * 5
    return [
        (
            "SSBOND "
            f"{fmt_int(ser_num, 3)} CYS"
            f" {chain_id1:1}"
            f"{fmt_int(seq_num1, 4)}"
            f"{icode1:1}"
            f"   CYS"
            f" {chain_id2:1}"
            f"{fmt_int(seq_num2, 4)}"
            f"{icode2:1}"
            f"                        {sym1:>6}"
            f" {sym2:>6} {length_str}"
        )
    ]


def gen_link(
    name1, alt_loc1, res_name1, chain_id1, res_seq1, icode1, name2, alt_loc2,
    res_name2, chain_id2, res_seq2, icode2, sym1, sym2, length
):
    length_str = f"{length:5.2f}" if length is not None else " " * 5
    return [
        (
            "LINK        "
            f"{name1:<4}"
            f"{alt_loc1:1}"
            f"{res_name1:>3}"
            f" {chain_id1:1}"
            f"{fmt_int(res_seq1, 4)}"
            f"{icode1:1}"
            f"               "
            f"{name2:<4}"
            f"{alt_loc2:1}"
            f"{res_name2:>3}"
            f" {chain_id2:1}"
            f"{fmt_int(res_seq2, 4)}"
            f"{icode2:1}"
            f"  {sym1:>6}"
            f" {sym2:>6} {length_str}"
        )
    ]


def gen_cispep(
    ser_num, pep1, chain_id1, seq_num1, icode1,
    pep2, chain_id2, seq_num2, icode2, mod_num, measure
):
    measure_str = f"{measure:6.2f}" if measure is not None else " " * 6
    return [
        (
            "CISPEP "
            f"{fmt_int(ser_num, 3)}"
            f" {pep1:<3}"
            f" {chain_id1:1}"
            f"{fmt_int(seq_num1, 4)}"
            f"{icode1:1}"
            f"   {pep2:<3}"
            f" {chain_id2:1}"
            f"{fmt_int(seq_num2, 4)}"
            f"{icode2:1}"
            f"       {fmt_int(mod_num, 3)}"
            f"      {measure_str}"
        )
    ]


def ssbond_to_writer_dict(serial: int, sb: SsBond) -> dict:
    return {
        "ser_num": serial,
        "chain_id1": sb.chain_id_1 or " ",
        "seq_num1": sb.seq_num_1,
        "icode1": (sb.icode_1 or " ")[:1],
        "chain_id2": sb.chain_id_2 or " ",
        "seq_num2": sb.seq_num_2,
        "icode2": (sb.icode_2 or " ")[:1],
        "sym1": sb.symmetry_operation_1 or "1555",
        "sym2": sb.symmetry_operation_2 or "1555",
        "length": sb.distance,
    }


def link_to_writer_dict(link: Link) -> dict:
    return {
        "name1": (link.name_1 or "")[:4].ljust(4),
        "alt_loc1": (link.alt_loc_1 or " ")[:1],
        "res_name1": "   ",
        "chain_id1": link.chain_id_1 or " ",
        "res_seq1": link.seq_num_1,
        "icode1": (link.icode_1 or " ")[:1],
        "name2": (link.name_2 or "")[:4].ljust(4),
        "alt_loc2": (link.alt_loc_2 or " ")[:1],
        "res_name2": "   ",
        "chain_id2": link.chain_id_2 or " ",
        "res_seq2": link.seq_num_2,
        "icode2": (link.icode_2 or " ")[:1],
        "sym1": "1555",
        "sym2": "1555",
        "length": link.distance,
    }


def cispep_to_writer_dict(serial: int, cp: CisPeptide) -> dict:
    return {
        "ser_num": serial,
        "pep1": "   ",
        "chain_id1": cp.chain_id_1 or " ",
        "seq_num1": cp.seq_num_1,
        "icode1": (cp.icode_1 or " ")[:1],
        "pep2": "   ",
        "chain_id2": cp.chain_id_2 or " ",
        "seq_num2": cp.seq_num_2,
        "icode2": (cp.icode_2 or " ")[:1],
        "mod_num": cp.num_model,
        "measure": cp.measure,
    }


def gen_site(site_id, residues):
    num_res = len(residues)
    chunks = [residues[i:i + 4] for i in range(0, len(residues), 4)] or [[]]
    out = []
    for index, chunk in enumerate(chunks, start=1):
        line = f"SITE   {index:>3} {site_id:<3}{num_res:>2}"
        for residue in chunk:
            line += (
                f" {residue.get('res_name', ''):>3}"
                f" {residue.get('chain_id', ' '):1}"
                f"{fmt_int(residue.get('seq_num'), 4)}"
                f"{residue.get('icode', ' '):1}"
            )
        out.append(line)
    return out


def sites_to_site_records(site_id: str, sites: list) -> list:
    """Build SITE lines from reader Site list (residue names omitted -> blank)."""
    residues = []
    for s in sites:
        if isinstance(s, Site):
            residues.append({
                "res_name": "   ",
                "chain_id": s.chain_id or " ",
                "seq_num": s.seq_num,
                "icode": (s.icode or " ")[:1],
            })
        else:
            residues.append(s)
    num_res = len(residues)
    chunks = [residues[i:i + 4] for i in range(0, len(residues), 4)] or [[]]
    out = []
    for index, chunk in enumerate(chunks, start=1):
        line = f"SITE   {index:>3} {site_id:<3}{num_res:>2}"
        for residue in chunk:
            line += (
                f" {residue.get('res_name', '   '):>3}"
                f" {residue.get('chain_id', ' '):1}"
                f"{fmt_int(residue.get('seq_num'), 4)}"
                f"{residue.get('icode', ' '):1}"
            )
        out.append(line)
    return out


def gen_cryst1(a, b, c, alpha, beta, gamma, s_group, z):
    return gen_line_from_template(
        "CRYST1{a:>9.3f}{b:>9.3f}{c:>9.3f}{alpha:>7.2f}{beta:>7.2f}{gamma:>7.2f} {s_group:<11}{z:>4}",
        {
            "a": {"value": a},
            "b": {"value": b},
            "c": {"value": c},
            "alpha": {"value": alpha},
            "beta": {"value": beta},
            "gamma": {"value": gamma},
            "s_group": {"value": s_group},
            "z": {"value": z}
        }
    )


def gen_matrix_rows(prefix, rows):
    out = []
    for n, row in enumerate(rows, start=1):
        out.extend(gen_line_from_template(
            "{prefix}{n}    {o1:10.6f}{o2:10.6f}{o3:10.6f}     {t_field:10.5f}",
            {
                "prefix": {"value": prefix},
                "n": {"value": n},
                "o1": {"value": row[0]},
                "o2": {"value": row[1]},
                "o3": {"value": row[2]},
                "t_field": {"value": row[3]}
            }
        ))
    return out


def gen_scale_matrix_rows(prefix, rows):
    out = []
    for n, row in enumerate(rows, start=1):
        out.extend(gen_line_from_template(
            "{prefix}{n}    {o1:10.6f}{o2:10.6f}{o3:10.6f}     {t_field:10.5f}",
            {
                "prefix": {"value": prefix},
                "n": {"value": n},
                "o1": {"value": row[0]},
                "o2": {"value": row[1]},
                "o3": {"value": row[2]},
                "t_field": {"value": row[3]}
            }
        ))
    return out


def gen_ncs_matrix(ncs_matrix):
    out = []
    for serial, block in enumerate(ncs_matrix, start=1):
        if isinstance(block, NcsMatrix):
            rows = block.matrix or []
            given = block.given
        else:
            rows = block.get("matrix", [])
            if "given" in block:
                given = block["given"]
            elif "calculated" in block:
                given = block["calculated"]
            else:
                given = True
        igiven = "1" if given else ""
        for n, row in enumerate(rows, start=1):
            out.extend(gen_line_from_template(
                "MTRIX{n}{serial:>3}{o1:10.6f}{o2:10.6f}{o3:10.6f}     {t_field:10.5f}    {igiven:1}",
                {
                    "n": {"value": n},
                    "serial": {"value": serial},
                    "o1": {"value": row[0]},
                    "o2": {"value": row[1]},
                    "o3": {"value": row[2]},
                    "t_field": {"value": row[3]},
                    "igiven": {"value": igiven}
                }
            ))
    return out


def gen_master(validation_info):
    vi = validation_info or {}

    def gv(key, default=0):
        return vi.get(key, default)

    return [
        (
            "MASTER    "
            f"{gv('num_remark'):>5}"
            f"{gv('num_het'):>5}"
            f"{gv('num_helix'):>5}"
            f"{gv('num_sheet'):>5}"
            f"{gv('num_site'):>5}"
            f"{gv('num_xform'):>5}"
            f"{gv('num_coord'):>5}"
            f"{gv('num_ter'):>5}"
            f"{gv('num_conect'):>5}"
            f"{gv('num_seq'):>5}"
            f"{0:>5}"
            f"{0:>5}"
        )
    ]


def gen_conect(connections):
    connections_with_dup = {}
    for atom1, bonded_atoms in connections.items():
        for atom2 in bonded_atoms:
            if atom1 not in connections_with_dup:
                connections_with_dup[atom1] = []
            connections_with_dup[atom1].append(atom2)
            if atom2 not in connections_with_dup:
                connections_with_dup[atom2] = []
            connections_with_dup[atom2].append(atom1)
    out = []
    for atom1 in sorted(connections_with_dup.keys()):
        bonded_atoms = sorted(set(connections_with_dup[atom1]))
        line = f"CONECT{atom1:>5}"
        for atom2 in bonded_atoms:
            line += f"{atom2:>5}"
        out.append(line)
    return out


def _gen_atom_line(atom_info):
    line = f'{atom_info["type"]:<6}{atom_info["atom_no"]:>5} {atom_info["atom_name"]:<4} '
    line += f'{atom_info["residue_name"]:>3} {atom_info["chain_id"]}'
    cx = atom_info["coord_x"] * unit_config.conversion_factor
    cy = atom_info["coord_y"] * unit_config.conversion_factor
    cz = atom_info["coord_z"] * unit_config.conversion_factor
    line += f"{atom_info['residue_id']:>4}    {cx:8.3f}"
    line += f"{cy:8.3f}{cz:8.3f}"
    occ = atom_info.get("occupancy", 1.0)
    line += f"{occ:6.2f}{atom_info['temp_factor']:6.2f}"
    line += f"          {atom_info['element']:>2}"
    if atom_info["charge"] == 1:
        charge_str = "1+"
    elif atom_info["charge"] == -1:
        charge_str = "1-"
    else:
        charge_str = ""
    line += f"{charge_str:>2}\n"
    return line
