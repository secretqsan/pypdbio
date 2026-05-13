# pylint: disable=missing-function-docstring
"""field parser functions"""
# Helper functions
def safe_append(obj, attribute, item):
    if getattr(obj, attribute) is None:
        setattr(obj, attribute, [])
    getattr(obj, attribute).append(item)


def safe_extend(obj, attribute, item):
    if getattr(obj, attribute) is None:
        setattr(obj, attribute, [])
    getattr(obj, attribute).extend(item)


def parse_int_value(value):
    value = value.strip()
    return int(value) if value else None


def parse_kv_entries(merged_text):
    entries = []
    for item in merged_text.split(';'):
        item = item.strip()
        if item == '':
            continue
        if ':' in item:
            key, value = item.split(':', maxsplit=1)
            key = key.strip().lower()
            value = value.strip()
            entries.append([key, value])
        else:
            entries.append([item, None])
    return entries


def parse_continuous_field(line, field, parser=None, parsed=None, separator=' ', place=-1):
    """
    Parses a continuous field of the PDB file.
    :param line: The line to parse.
    :param field: The field to parse.
    :param parser: The parser to use.
    :param separator: The separator to use.
    :param place: The place to insert the value.
    :return: The parsed value.
    """
    # get the parsed value
    if parser is not None:
        value = parser(line)
    else:
        value = parsed
    # uniform the shape of the input
    if place == -1:
        field = [field]
        value = [value]
    for i in range(len(field)):  # pylint: disable=consider-using-enumerate
        if not field[i]:
            field[i] = value[i]
        elif i == place or place == -1:
            if isinstance(value[i], str):
                field[place] += f"{separator}{value[i]}"
            elif isinstance(value[i], list):
                field[place].extend(value[i])
    return field if place != -1 else field[0]


def parse_float_value(value):
    value = value.strip()
    return float(value) if value else None

# Meta data related functions


def parse_header(line):
    classification = line[10:50].strip()
    date = line[50:59].strip()
    pdb_id = line[62:66].strip()
    return classification, date, pdb_id


def parse_obslte(line):
    obslte_date = line[11:20].strip()
    new_entry_id = line[31:35].strip()
    return obslte_date, new_entry_id


def parse_title(line):
    title = line[10:80].strip()
    return title


def parse_split(line):
    split_list = line[11:80].strip().split()
    return split_list


def parse_caveat(line):
    caveat = line[19:79].strip()
    return caveat


def parse_author(line):
    author_list = line[10:80].strip().split(',')
    return [author.strip() for author in author_list if author.strip() != '']


def parse_remark(line):
    remark_id = line[7:10].strip()
    remark_text = line[11:80].strip()
    return [remark_id, remark_text]


def parse_keywords(line):
    keywords = line[10:79].strip().split(',')
    keywords = [keyword.strip() for keyword in keywords]
    return keywords


def parse_compond(line):
    return line[10:80].strip()


def parse_source(line):
    return line[10:79].strip()


def parse_expdta(line):
    experiments = line[10:79].strip().split(';')
    return [experiment.strip() for experiment in experiments]


def parse_mdltyp(line):
    return line[10:80].strip()


def parse_nummdl(line):
    return parse_int_value(line[10:14])


def parse_jrnl(line):
    jrnl_type = line[12:16].strip()
    jrnl_text = line[19:79].strip()
    return jrnl_type, jrnl_text


def parse_jrnl_auth(line):
    author_list = line[19:79].strip().split(',')
    return [author.strip() for author in author_list if author.strip() != '']


def parse_jrnl_titl(line):
    return line[19:79].strip()


def parse_jrnl_edit(line):
    return line[19:79].strip()


def parse_jrnl_ref(line):
    published = line[19:34].strip() != "TO BE PUBLISHED"
    if not published:
        return "", "", "", ""
    pub_name = line[19:47].strip()
    volume = line[51:55].strip()
    page = line[56:61].strip()
    year = line[62:66].strip()
    return pub_name, volume, page, year


def parse_jrnl_publ(line):
    return line[19:70].strip()


def parse_jrnl_refn(line):
    refn_type = line[35:39].strip()
    refn_value = line[40:65].strip()
    return refn_type, refn_value


def parse_jrnl_pmid(line):
    return line[19:79].strip()


def parse_jrnl_doi(line):
    return line[19:79].strip()


def parse_sprsde(line):
    sprsde_date = line[11:20].strip()
    superseded_ids = []
    for start in range(31, 76, 5):
        sid_code = line[start:start + 4].strip()
        if sid_code != "":
            superseded_ids.append(sid_code)
    return sprsde_date, superseded_ids


def parse_revdat(line):
    date = line[13:22].strip()
    modifications = []
    for start in [39, 46, 53, 60]:
        modification = line[start:start + 6].strip()
        if modification != "":
            modifications.append(modification)
    return date, modifications

# Primary structure related functions


def parse_dbref(line):
    chain_id = line[12:13].strip()
    seq_begin = parse_int_value(line[14:18])
    insert_begin = line[18:19].strip()
    seq_end = parse_int_value(line[20:24])
    insert_end = line[24:25].strip()
    database = line[26:32].strip()
    db_accession = line[33:41].strip()
    db_id_code = line[42:54].strip()
    db_seq_begin = parse_int_value(line[55:60])
    db_ins_begin = line[60:61].strip()
    db_seq_end = parse_int_value(line[62:67])
    db_ins_end = line[67:68].strip()
    return [
        chain_id, seq_begin, insert_begin, seq_end, insert_end,
        database, db_accession, db_id_code, db_seq_begin,
        db_ins_begin, db_seq_end, db_ins_end
    ]


def parse_dbref1(line):
    chain_id = line[12:13].strip()
    seq_begin = parse_int_value(line[14:18])
    insert_begin = line[18:19].strip()
    seq_end = parse_int_value(line[20:24])
    insert_end = line[24:25].strip()
    database = line[26:32].strip()
    db_id_code = line[47:67].strip()
    return [
        chain_id, seq_begin, insert_begin, seq_end, insert_end,
        database, db_id_code
    ]


def parse_dbref2(line):
    chain_id = line[12:13].strip()
    db_accession = line[18:40].strip()
    db_seq_begin = parse_int_value(line[45:55])
    db_seq_end = parse_int_value(line[57:67])
    return [
        chain_id, db_accession, db_seq_begin, db_seq_end
    ]


def parse_seqadv(line):
    res_name = line[12:15].strip()
    chain_id = line[16:17].strip()
    seq_num = parse_int_value(line[18:22])
    icode = line[22:23].strip()
    database = line[24:28].strip()
    db_accession = line[29:38].strip()
    db_res = line[39:42].strip()
    db_seq = parse_int_value(line[43:48])
    conflict = line[49:70].strip()
    return [res_name, chain_id, seq_num, icode, database, db_accession, db_res, db_seq, conflict]


def parse_seqres(line):
    chain_id = line[11:12]
    num_res = parse_int_value(line[13:17])
    residues = []
    for start in range(19, 70, 4):
        name = line[start:start + 3].strip()
        if name:
            residues.append(name)
    return [chain_id, num_res, residues]


def parse_modres(line):
    res_name = line[12:15].strip()
    chain_id = line[16:17].strip()
    seq_num = parse_int_value(line[18:22])
    icode = line[22:23].strip()
    std_res = line[24:27].strip()
    comment = line[29:70].strip()
    return [res_name, chain_id, seq_num, icode, std_res, comment]

# Heterogen related functions


def parse_het(line):
    het_id = line[7:10].strip()
    text = line[30:70].strip()
    return [het_id, text]


def parse_hetnam(line):
    het_id = line[11:14].strip()
    het_name = line[15:70].strip()
    return [het_id, het_name]


def parse_hetsyn(line):
    het_id = line[11:14].strip()
    het_syn = line[15:70].strip()
    return [het_id, het_syn]


def parse_formul(line):
    het_id = line[12:15].strip()
    text = line[19:70].strip()
    return [het_id, text]

# Secondary structure related functions


def parse_helix(line):
    helix_id = line[11:14].strip()
    chain_id = line[19:20]
    init_seq_num = parse_int_value(line[21:25])
    init_icode = line[25:26].strip()
    end_res_num = parse_int_value(line[33:37])
    end_icode = line[37:38].strip()
    helix_class = parse_int_value(line[38:40])
    comment = line[40:70].strip()
    return [
        helix_id,
        chain_id,
        init_seq_num,
        init_icode,
        end_res_num,
        end_icode,
        helix_class,
        comment
    ]


def parse_sheet(line):
    sheet_id = line[11:14].strip()
    init_chain_id = line[21:22]
    init_seq_num = parse_int_value(line[22:26])
    init_icode = line[26:27].strip()
    end_chain_id = line[32:33]
    end_res_seq = parse_int_value(line[33:37])
    end_icode = line[37:38].strip()
    sense = parse_int_value(line[38:40])
    # registration information
    cur_atom = line[41:45].strip()
    cur_chain_id = line[49:50]
    cur_res_seq = parse_int_value(line[50:54])
    cur_icode = line[54:55].strip()
    prev_atom = line[56:60].strip()
    prev_chain_id = line[64:65]
    prev_res_seq = parse_int_value(line[65:69])
    prev_icode = line[69:70].strip()
    return [
        sheet_id,
        init_chain_id,
        init_seq_num,
        init_icode,
        sense,
        cur_atom,
        cur_chain_id,
        cur_res_seq,
        cur_icode,
        prev_atom,
        prev_chain_id,
        prev_res_seq,
        prev_icode,
        end_chain_id,
        end_res_seq,
        end_icode,
    ]

# Connectivity annotation related functions


def parse_ssbond(line):
    chain_id1 = line[15:16]
    seq_num1 = parse_int_value(line[17:21])
    icode1 = line[21:22].strip()
    chain_id2 = line[29:30]
    seq_num2 = parse_int_value(line[31:35])
    icode2 = line[35:36].strip()
    sym1 = line[59:65].strip()
    sym2 = line[66:72].strip()
    length = parse_float_value(line[73:78])
    return [
        chain_id1, seq_num1, icode1, chain_id2, seq_num2, icode2, sym1, sym2, length
    ]


def parse_link(line):
    name1 = line[12:16].strip()
    alt_loc1 = line[16:17].strip()
    res_name1 = line[17:20].strip()
    chain_id1 = line[21:22]
    res_seq1 = parse_int_value(line[22:26])
    icode1 = line[26:27].strip()
    name2 = line[42:46].strip()
    alt_loc2 = line[46:47].strip()
    res_name2 = line[47:50].strip()
    chain_id2 = line[51:52]
    res_seq2 = parse_int_value(line[52:56])
    icode2 = line[56:57].strip()
    sym1 = line[59:65].strip()
    sym2 = line[66:72].strip()
    length = parse_float_value(line[73:78])
    return [
        name1,
        alt_loc1,
        res_name1,
        chain_id1,
        res_seq1,
        icode1, name2,
        alt_loc2,
        res_name2,
        chain_id2,
        res_seq2,
        icode2,
        sym1,
        sym2,
        length,
    ]


def parse_cispep(line):
    chain_id1 = line[15:16]
    seq_num1 = parse_int_value(line[17:21])
    icode1 = line[21:22].strip()
    chain_id2 = line[29:30]
    seq_num2 = parse_int_value(line[31:35])
    icode2 = line[35:36].strip()
    mod_num = parse_int_value(line[43:46])
    measure = parse_float_value(line[53:59])
    return [
        chain_id1,
        seq_num1,
        icode1,
        chain_id2,
        seq_num2,
        icode2,
        mod_num,
        measure,
    ]
# Misc features related functions


def parse_site(line):
    site_id = line[11:14].strip()
    residues = []
    for start in [18, 29, 40, 51]:
        chain_id = line[start + 4:start + 5].strip()
        seq_num = parse_int_value(line[start + 5:start + 9])
        icode = line[start + 9:start + 10].strip()
        if seq_num is not None:
            residues.append([chain_id, seq_num, icode])
    return [site_id, residues]
# Crystallographic related functions


def parse_cryst1(line):
    a = float(line[6:15].strip())
    b = float(line[15:24].strip())
    c = float(line[24:33].strip())
    alpha = float(line[33:40].strip())
    beta = float(line[40:47].strip())
    gamma = float(line[47:54].strip())
    s_group = line[55:66].strip()
    z = int(line[66:70].strip())
    return [a, b, c, alpha, beta, gamma, s_group, z]


def parse_matrix_line(line):
    o1 = parse_float_value(line[10:20])
    o2 = parse_float_value(line[20:30])
    o3 = parse_float_value(line[30:40])
    t = parse_float_value(line[45:55])
    return [o1, o2, o3, t]


def parse_ncs_matrix_line(line):
    o1 = parse_float_value(line[10:20])
    o2 = parse_float_value(line[20:30])
    o3 = parse_float_value(line[30:40])
    t = parse_float_value(line[45:55])
    i_given = parse_int_value(line[59:60])
    return [o1, o2, o3, t, i_given]

# Coordinate related functions


def parse_atom(line):
    typ = line[0:6].strip()
    atom_no = line[6:11].strip()
    atom_name = line[12:16].strip()
    alt_loc = line[16:17].strip()
    res_name = line[17:20].strip()
    chain_id = line[21]
    residue_id = line[22:26].strip()
    icode = line[26:27].strip()
    coord_x = line[30:38].strip()
    coord_y = line[38:46].strip()
    coord_z = line[46:54].strip()
    occupancy = line[54:60].strip()
    temp_factor = line[60:66].strip()
    element = line[76:78].strip()
    charge = line[78:80].strip()
    return [
        typ,
        int(atom_no),
        atom_name,
        alt_loc,
        res_name,
        chain_id,
        int(residue_id),
        icode,
        float(coord_x),
        float(coord_y),
        float(coord_z),
        float(occupancy),
        float(temp_factor),
        element,
        0 if charge == '' else 1 if charge == '1+' else -1
    ]


def parse_anisou(line):
    u11 = parse_int_value(line[28:35]) / 10000
    u22 = parse_int_value(line[35:42]) / 10000
    u33 = parse_int_value(line[42:49]) / 10000
    u12 = parse_int_value(line[49:56]) / 10000
    u13 = parse_int_value(line[56:63]) / 10000
    u23 = parse_int_value(line[63:70]) / 10000
    return [
        u11,
        u22,
        u33,
        u12,
        u13,
        u23,
    ]

# Connectivity related functions


def parse_connect(line):
    atom1 = int(line[6:11].strip())
    bonded_atoms = []
    for i in range(11, 31, 5):
        bonded_atom_str = line[i:i + 5].strip()
        if bonded_atom_str != '':
            bonded_atom_index = int(bonded_atom_str)
            if bonded_atom_index > atom1:
                bonded_atoms.append(int(bonded_atom_str))
    return atom1, bonded_atoms

# Bookkeeping related functions


def parse_master(line):
    num_remark = parse_int_value(line[10:15])
    num_het = parse_int_value(line[20:25])
    num_helix = parse_int_value(line[25:30])
    num_sheet = parse_int_value(line[30:35])
    num_site = parse_int_value(line[40:45])
    num_xform = parse_int_value(line[45:50])
    num_coord = parse_int_value(line[50:55])
    num_ter = parse_int_value(line[55:60])
    num_conect = parse_int_value(line[60:65])
    num_seq = parse_int_value(line[65:70])
    return [
        num_remark,
        num_het,
        num_helix,
        num_sheet,
        num_site,
        num_xform,
        num_coord,
        num_ter,
        num_conect,
        num_seq,
    ]
