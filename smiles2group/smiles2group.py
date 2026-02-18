"""
SMILES から GenIce3 Group 互換の sites / labels / bonds を生成するコア処理。

RDKit で 3D 座標を生成し、アンカー原子を原点に、残基の重心が z 軸方向になるように回転する。
アンカーに結合した水素を1つ除き、イオン結合で置き換わる価数にする。
"""

from __future__ import annotations

import numpy as np
from typing import Any

try:
    from rdkit import Chem
    from rdkit.Chem import rdDistGeom
    _HAS_RDKIT = True
except ImportError:
    _HAS_RDKIT = False


def _rotation_to_z_positive(v: np.ndarray) -> np.ndarray:
    """ベクトル v を z 軸正方向 (0,0,1) に重ねる回転行列（行ベクトル適用: new = old @ R.T）。"""
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v)
    if n < 1e-10:
        return np.eye(3)
    v = v / n
    if abs(v[2]) >= 1.0 - 1e-9:
        x = np.array([1.0, 0.0, 0.0])
    else:
        x = np.array([0.0, 0.0, 1.0])
    y = np.cross(v, x)
    yn = np.linalg.norm(y)
    if yn < 1e-10:
        x = np.array([1.0, 0.0, 0.0])
        y = np.cross(v, x)
        yn = np.linalg.norm(y)
    y = y / yn
    x = np.cross(y, v)
    x = x / np.linalg.norm(x)
    R = np.array([x, y, v])
    return R


def _find_dummy_and_anchor(mol) -> tuple[int | None, int | None]:
    """[*] ダミー原子と、その隣（アンカー）のインデックスを返す。"""
    for i in range(mol.GetNumAtoms()):
        if mol.GetAtomWithIdx(i).GetAtomicNum() == 0:
            neighbors = [a.GetIdx() for a in mol.GetAtomWithIdx(i).GetNeighbors()]
            if neighbors:
                return (i, neighbors[0])
            return (i, None)
    return (None, None)


def smiles_to_group_data(
    smiles: str,
    anchor_atom_index: int | None = None,
    z_neighbor_index: int | None = None,
    name: str | None = None,
    embed_seed: int = 0,
) -> dict[str, Any]:
    """
    SMILES から GenIce3 Group 互換の辞書を生成する。

    規約:
      - アンカー原子の位置が原点になる。
      - アンカーに結合した水素を1つ除く（イオン結合で置き換わるため）。
      - 残基の重心が z 軸方向になるように座標変換する（除く水素の向きにはよらない）。

    Parameters
    ----------
    smiles : str
        SMILES 文字列。結合点を [*] で書ける（例: "[*]CCCC" でブチル基）。
    anchor_atom_index : int, optional
        置換イオンに直結する原子の 0-based インデックス。[*] がない場合に使用。
    z_neighbor_index : int, optional
        [*] がない場合に、除く水素（アンカーの隣接 H）のインデックス。
    name : str, optional
        Group の名前。
    embed_seed : int
        RDKit 3D 埋め込みの乱数シード。

    Returns
    -------
    dict
        "sites": (N,3) ndarray（単位は nm）、"labels", "bonds", "name"
    """
    if not _HAS_RDKIT:
        raise RuntimeError("RDKit が必要です: pip install rdkit")

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"無効な SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    n_atoms = mol.GetNumAtoms()

    dummy_idx, anchor_from_dummy = _find_dummy_and_anchor(mol)
    if dummy_idx is not None and anchor_from_dummy is not None:
        anchor_atom_index = anchor_from_dummy
        drop_dummy = True
    else:
        drop_dummy = False
        if anchor_atom_index is None:
            raise ValueError(
                "SMILES に [*] が含まれていない場合は、アンカー原子のインデックスを指定してください"
            )

    if anchor_atom_index < 0 or anchor_atom_index >= n_atoms:
        raise ValueError(f"anchor_atom_index は 0..{n_atoms - 1} の範囲で指定してください")

    params = rdDistGeom.ETKDGv3()
    params.randomSeed = embed_seed
    cid = rdDistGeom.EmbedMolecule(mol, params)
    if cid != 0:
        raise RuntimeError("3D 埋め込みに失敗しました。SMILES やシードを変更してみてください。")

    conf = mol.GetConformer()
    ANGSTROM_TO_NM = 0.1
    sites = np.array(conf.GetPositions(), dtype=float) * ANGSTROM_TO_NM

    bonds: list[tuple[int, int]] = []
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bonds.append((i, j))

    old_anchor_idx = anchor_atom_index
    anchor_atom = mol.GetAtomWithIdx(old_anchor_idx)
    anchor_h_neighbors = [
        a.GetIdx()
        for a in anchor_atom.GetNeighbors()
        if a.GetSymbol() == "H" and (dummy_idx is None or a.GetIdx() != dummy_idx)
    ]
    if not anchor_h_neighbors:
        raise ValueError(
            "アンカー原子に結合した水素がありません。置換点となる水素が必要です。"
        )
    if dummy_idx is not None:
        anchor_pos = sites[anchor_atom_index]
        dummy_dir = sites[dummy_idx] - anchor_pos
        if np.linalg.norm(dummy_dir) > 1e-10:
            remove_h_idx = max(
                anchor_h_neighbors,
                key=lambda h: np.dot(sites[h] - anchor_pos, dummy_dir),
            )
        else:
            remove_h_idx = anchor_h_neighbors[0]
    else:
        remove_h_idx = (
            z_neighbor_index
            if z_neighbor_index is not None and z_neighbor_index in anchor_h_neighbors
            else anchor_h_neighbors[0]
        )

    drop_indices = {dummy_idx} if dummy_idx is not None else set()
    drop_indices.add(remove_h_idx)
    real_indices = [i for i in range(n_atoms) if i not in drop_indices]
    old_to_new: dict[int, int] = {idx: k for k, idx in enumerate(real_indices)}
    sites = sites[real_indices]
    labels = [mol.GetAtomWithIdx(i).GetSymbol() for i in real_indices]
    bonds = [
        (old_to_new[i], old_to_new[j])
        for i, j in bonds
        if i not in drop_indices and j not in drop_indices
    ]
    n_atoms = len(real_indices)
    anchor_atom_index = old_to_new[anchor_atom_index]

    sites = sites - sites[anchor_atom_index]

    centroid = np.mean(sites, axis=0)
    z_dir_norm = np.linalg.norm(centroid)
    if z_dir_norm >= 1e-10:
        z_dir = centroid / z_dir_norm
        R = _rotation_to_z_positive(z_dir)
        sites = sites @ R.T

    group_name = name if name else f"SMILES_{smiles[:20].replace(' ', '_')}"

    return {
        "sites": sites,
        "labels": labels,
        "bonds": bonds,
        "name": group_name,
    }
