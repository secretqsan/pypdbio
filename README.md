# pypdbio

`pypdbio` 是一个用于 **PDB 文件下载、读取与写入** 的轻量 Python 包。

## 安装

```bash
pip install pypdbio
```

## 快速开始

```python
from pypdbio import fetch, PdbReader, PdbWriter

# 1) 下载 PDB 文件（示例：1aki）
fetch("1aki", path="1aki.pdb")

# 2) 读取 PDB 文件
reader = PdbReader("1aki.pdb")
pdb_data = reader.read()

# 3) 写出 PDB 文件
writer = PdbWriter("1aki_copy.pdb")
writer.write(pdb_data)
```

## 主要能力

- 通过 PDB ID 从 RCSB 下载结构文件
- 将 PDB 文件解析为可遍历的数据结构（模型、链、残基、原子）
- 将数据结构写回标准 PDB 文本文件

## 许可证

MIT
