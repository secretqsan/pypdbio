# pypdbio API 参考

本文档描述 `pypdbio` 包面向用户的 **公共 API** 与常用数据模型。与「快速上手」相关的简介见仓库根目录的 [README.md](../README.md)。

---

## 目录

1. [包入口与导出](#包入口与导出)
2. [导出数据类速览](#导出数据类速览)
3. [坐标单位](#坐标单位)
4. [下载：`fetch`](#下载：`fetch`)
5. [读取：`PdbReader`](#读取：`pdbreader`)
6. [写入：`PdbWriter`](#写入：`pdbwriter`)
7. [数据模型概览](#数据模型概览)
8. [层级与索引](#层级与索引)
9. [贡献与反馈](#贡献与反馈)

---

## 包入口与导出

以下名称由包根 `__all__` 声明，推荐 `from pypdbio import ...` 直接导入（与 `import pypdbio` 后 `pypdbio.<name>` 一致）。

### I/O 与全局配置


| 名称          | 说明                    |
| ----------- | --------------------- |
| `PdbReader` | PDB 文本文件解析器           |
| `PdbWriter` | PDB 文本文件写出器           |
| `fetch`     | 从 RCSB 按 ID 下载 `.pdb` |
| `set_unit`  | 设置全局长度单位（影响读写换算）      |


### 坐标层级


| 名称        | 说明      |
| --------- | ------- |
| `PdbData` | 完整条目根对象 |
| `Model`   | 单套坐标模型  |
| `Chain`   | 链       |
| `Residue` | 残基      |
| `Atom`    | 原子      |


### 数据类


| 名称                        | 说明                                               |
| ------------------------- | ------------------------------------------------ |
| `PdbMetaData`             | HEADER 区主元数据容器                                   |
| `ObsoleteInfo`            | `OBSLTE`：替换日期、新条目 ID                             |
| `JournalInfo`             | `JRNL`：作者、题名、期刊、DOI 等                            |
| `RevisionInfo`            | `REVDAT`：修订日期与说明列表                               |
| `ReplaceInfo`             | `SPRSDE`：取代关系                                    |
| `CrystalInfo`             | 晶胞、`CRYST1`、`ORIGX`/`SCALE`、占位晶胞等                |
| `NcsMatrix`               | 非晶体学对称矩阵（`igiven` 标志 + `matrix`）                 |
| `ConnectivityInfo`        | 连接信息，包含`SSBOND`/`LINK`/`CISPEP` 列表与 `CONECT` 字典等 |
| `SsBond`                  | 二硫键结构化字段                                         |
| `Link`                    | 化学键 `LINK` 结构化字段                                 |
| `CisPeptide`              | `CISPEP` （顺式氨基酸）结构化字段                            |
| `SecondaryStructureInfo`  | 二级结构容器（`helix`、`sheet` 字典）                       |
| `Helix`                   | 单条螺旋几何/分类信息                                      |
| `SheetStrand`             | β 片层中一条链股                                        |
| `Site`                    | `SITE` 中一处残基引用                                   |
| `Heterogen`               | `HET`/`FORMUL` 等异质组分描述                           |
| `SequenceInfo`            | `DBREF`/`SEQRES`/`SEQADV`/`MODRES` 等集合           |
| `SequenceDBInfo`          | `DBREF` 数据库区间与 accession                         |
| `SequenceDifferenceInfo`  | `SEQADV` 差异项                                     |
| `ResidueModificationInfo` | `MODRES` 修饰项                                     |


每一个字段可能有几种表示方式，即列表、字典和结构化对象，具体区分标准为：如果字段依赖ID，则为以ID为键的字典，例如pdb文件的注释；如果存在多个并列的描述，则为列表，例如pdb文件的修订记录，否则为结构化对象，例如pdb文件的期刊信息。也可以通过文本编辑器的Intellisense信息获取相关内容。

每一字段的详细含义可参考 [PDB 规范](https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html)。缺失的字段通常可以由程序自动计算。

---

## 坐标单位

库在 **读入与写出** 时，在「内存中的数值」与「PDB 文件中的埃（Å）」之间做换算。

### 默认与环境变量

- 进程启动时读取环境变量 `**PYPDBIO_UNIT`**。
- 未设置时默认为 `**nm`（纳米）**：与文件交互时，坐标与部分晶体学长度会按约定换算（与 `unit_config.conversion_factor` 一致）。
- 设为 `angstrom` 或 `A` 时，内部数值与 PDB 中 Å 一致（因子为 `1.0`）。若 `unit` 非法，会发出 `warnings` 并回退为因子 `1.0` 的行为（与「明确 Å」接近，但建议始终传入合法值）。

### `set_unit(unit)`

```python
from pypdbio import set_unit

set_unit("nm")         # 与默认一致（若未改环境变量）
set_unit("angstrom")   # 或 set_unit("A")
```

- 在 **同一进程** 内为全局配置；多线程场景需注意共享状态，请勿在运行过程中修改单位设置，否则可能因发生重复换算而导致错误。
- 请注意各变换矩阵的单位换算，依赖默认行为可能会出现错误。

---

## 下载：`fetch`

```python
def fetch(pdb_id: str, path: str | None = None) -> None
```


| 参数       | 说明                                 |
| -------- | ---------------------------------- |
| `pdb_id` | 四位 PDB ID；内部会转为 **大写**。            |
| `path`   | 保存路径。为 `None` 时保存为 `<PDB_ID>.pdb`。 |


- 请求地址：`https://files.rcsb.org/download/<ID>.pdb`。
- HTTP 状态非 200 时抛出 `FileNotFoundError`，提示该 ID 不存在。
- 以 **二进制** 写入文件，与 RCSB 返回内容一致。

---

## 读取：`PdbReader`

### 构造

```python
PdbReader(path: str)
```

- `path`：PDB 文件路径。
- 文件以 `utf-8` 打开；若文件含非 UTF-8 字节，可能需在系统侧转换编码。

### `read() -> PdbData`

```python
reader = PdbReader("structure.pdb")
pdb_data = reader.read()
```

按行扫描文件，依次解析 HEADER、一级结构、异质、二级结构、杂项、晶体学、坐标、连通性、簿记等记录，经后处理与校验后返回 `PdbData`。

### 读入后行为与告警

`read()` 末尾会调用内部校验，在以下情况可能触发警告（不中断解析）：

- 条目已废弃（`OBSLTE`）、拆分（`SPLIT`）、存在 `CAVEAT`；
- `NUMMDL` 与实际模型数不一致；
- 文件中 `MASTER` 等与内部计数不一致等

---

## 写入：`PdbWriter`

### 构造

```python
PdbWriter(path: str)
```

- `path`：输出 PDB 路径。
- 写出时使用 `utf-8` 编码。

### `write(data)`

支持三种类型：


| 参数类型      | 行为                                                                                                           |
| --------- | ------------------------------------------------------------------------------------------------------------ |
| `PdbData` | 完整写出：HEADER、一级结构、异质、二级结构、连通性、SITE、晶体学、坐标、`CONECT`、`MASTER`、`END` 等（按数据中是否具备相应字段决定段落）。                        |
| `Model`   | 包装为新的 `PdbData`：含单个模型；`meta.title` 设为 `"pdb"`，`meta.remark["0"]` 为 `"Generated by pypdbio"`，再按 `PdbData` 写出。 |
| `Chain`   | 再包一层：单模型、单链，再同上。                                                                                             |


---

## 数据模型概览

### `PdbData`


| 属性 / 成员               | 说明                                                                                     |
| --------------------- | -------------------------------------------------------------------------------------- |
| `models`              | **模型列表**。                                                                              |
| `meta`                | HEADER 信息。                                                                             |
| `crystallographic`    | 晶体学信息，包含选 `ORIGX`/`SCALE` 等。                                                           |
| `connectivity`        | 键连信息，包含`SSBOND`/`LINK`/`CISPEP`/`CONECT` 等相关部分。注意，在pdb文件中多肽相关的键连关系与异质部分的键连关系分别存储在不同字段。 |
| `secondary_structure` | 二级结构信息（螺旋、片层等）。                                                                        |
| `heterogen`           | 异质组分字典（读入时填充）。                                                                         |
| `sites`               | SITE 特征相关。                                                                             |
| `add_model(model)`    | 添加 `Model`。                                                                            |


内部另有 `_tmp`、`_validation_info` 等，供解析过程使用，一般应用代码不应修改。

### `Model`


| 属性 / 方法            | 说明   |
| ------------------ | ---- |
| `chains`           | 链列表。 |
| `add_chain(chain)` | 添加链。 |


### `Chain`


| 属性 / 方法                | 说明                                      |
| ---------------------- | --------------------------------------- |
| `residues`             | 残基列表。                                   |
| `id`                   | 链 ID（**属性**）。若未显式设置，则按在模型中的顺序推导默认字母 ID。 |
| `sequence_info`        | 可选 `SequenceInfo`（DBREF、SEQRES 等）。      |
| `add_residue(residue)` | 添加残基。                                   |


### `Residue`

```python
Residue(name: str)
```


| 属性               | 说明                      |
| ---------------- | ----------------------- |
| `name`           | 残基名（如 `"ALA"`、`"HOH"`）。 |
| `atoms`          | 原子列表。                   |
| `icode`          | 插入码（PDB 单列）；空字符串表示无插入码。 |
| `het`            | 是否为异质组分。                |
| `solvent`        | 溶剂等标记（会影响异质部分的写出）。      |
| `end_of_chain`   | 是否为链的最后一个残基。            |
| `add_atom(atom)` | 添加原子。                   |


`id` 由在链中的位置及 `icode` 规则在父类中推导，在没有插入码的情况下不建议使用，会大幅延长迭代速度。

### `Atom`

```python
Atom(
    name: str = "",
    coord: tuple = (0.0, 0.0, 0.0),
    temp_factor: float = 0.0,
    element: str = "",
    charge: int = 0,
    occupancy: float = 1.0,
    alt_loc: str = "",
)
```


| 属性            | 说明                                                                        |
| ------------- | ------------------------------------------------------------------------- |
| `name`        | 原子名（如 `" CA "` ，无左右空格）。                                                   |
| `coord`       | `(x, y, z)`，单位受 `set_unit` 影响。                                            |
| `temp_factor` | 各向同性 B 因子。在包含`anisotropic_temp_factor` 时，通过`anisotropic_temp_factor`计算得到。 |
| `element`     | 元素符号。                                                                     |
| `charge`      | 形式电荷（写出列）。                                                                |
| `occupancy`   | 占有率。                                                                      |
| `alt_loc`     | 交替构象 ID。                                                                  |


---

## 层级与索引

### 遍历

```python
for model in pdb_data.models:
    for chain in model.chains:
        for residue in chain.residues:
            for atom in residue.atoms:
                ...
```

### `PdbData[model_index]`

- **整数** `i`：第 `i` 个模型（与 `pdb_data.models[i]` 一致）。

### `Model[index]`

- **字符串**（如 `"A"`）：按链 ID 取链。已显式声明 ID 的链优先匹配；否则按未声明链的顺序与字母占位规则解析（与 `pypdbio.utils` 中链 ID 推算一致）。
- **整数** `i`：第 `i` 条链。

### `Chain[index]`

- **整数** `i`：第 `i` 个残基（0 起）。
- **字符串**：按 **「标准残基序号 + 插入码」** 查找：解析方式为 `icode = index[-1].strip()`，`seq_num = int(index[:-1])`，再在链中累计 `icode == ""` 的残基以得到「序号」，与 `icode` 联合匹配。
  - **无插入码** 时，字符串最后一位应为 **空白**，例如 `"12 "` 表示第 12 号标准位、插入码为空（注意与 `"12"` 区分，后者会把 `"2"` 误当作插入码）。
  - **有插入码** 时形如 `"12A"`。

### 链式访问示例

读入后常见写法（与 writer 内逻辑一致）：

```python
res = pdb_data.models[0]["A"]["23 "]   # 链 A，标准序号 23，无插入码
# 或
res = pdb_data.models[0][0][22]       # 按残基列表下标（若与文件顺序一致）
```

---

## 贡献与反馈

若你发现与 PDB 规范边缘情况相关的问题，建议附带 **最小 PDB 样例** 与 **期望行为** 提交 issue 或 PR。

---

