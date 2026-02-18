# smiles2group

SMILES から [GenIce3](https://github.com/genice-dev/GenIce3) の **Group プラグイン**用の原子配置（sites, labels, bonds）を生成するツールです。  
[RDKit](https://www.rdkit.org/) で 3D 座標を生成し、GenIce3 の内部単位（nm）で出力します。

## インストール

```bash
pip install smiles2group
```

## 使い方

### コマンドライン

```bash
# [*] で結合点を指定（推奨）。アンカーは [*] の隣の原子、[*] とアンカーの水素1つは出力から除かれる
smiles2group "[*]CCCC" --name butyl -o butyl.py

# [*] なしの場合はアンカー原子の 0-based インデックスを指定
smiles2group "C" 0 --name methyl -o methyl.py
```

### Python API

```python
from smiles2group import smiles_to_group_data

data = smiles_to_group_data("[*]CCCC", name="butyl")
# data["sites"], data["labels"], data["bonds"], data["name"]
```

### GenIce3 で使う

生成した `.py` を `~/.genice/group/` に置くか、`genice3_group` の entry point で提供すると、  
`--spot_cation "[51=N --group 12=butyl]"` のように指定できます。GenIce3 は別途インストールしてください。

## 規約（GenIce3 Group）

- **アンカー原子**：置換イオンに直結する原子。原点。
- **除く水素**：アンカーに結合した水素を1つ除く（イオン結合で置き換えるため）。
- **z 軸**：残基の重心が z 軸方向になるように回転。

## 開発・GitHub への公開

リポジトリは `genice-dev/smiles2group` を想定しています。

1. GitHub で [genice-dev/smiles2group](https://github.com/genice-dev/smiles2group) を新規作成（空のリポジトリ、README なしで可）。
2. 初回プッシュ:
   ```bash
   cd /path/to/smiles2group
   git push -u origin main
   ```

## ライセンス

MIT
