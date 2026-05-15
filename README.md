# pypdbio

`pypdbio` 是一个用于 **PDB 文件下载、读取与写入** 的轻量 Python 包，兼容PDB标准的所有特征。

## 安装

```bash
pip install pypdbio
```

依赖：`requests`（用于 `fetch`）。

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

## API 参考

包根 `from pypdbio import ...` 已导出解析/写出类型及 HEADER、连通性、二级结构、序列等相关数据类（以 `pypdbio.__all__` 为准）。完整说明见 [docs/reference.md](docs/reference.md)。

## 许可证

MIT