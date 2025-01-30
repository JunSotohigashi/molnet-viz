"""
EQ,TS 可視化プログラム
pandas,pyvis,rdkit
"""

import networkx as nx
import json
import pandas as pd
import pyvis
import rdkit.Chem
import rdkit.Chem.rdmolfiles

# 読み込むファイルのパス
EQ_PATH = "manu_EQ_list.log"
TS_PATH = "manu_TS_list.log"
BOND_DB_PATH = "bond_length_database.json"
# 保存先のパス
SAVE_CSV_PATH = "manu.csv"
SAVE_GRAPH_PATH = "manu.html"

# データベースの結合長に対して、何倍まで許容するか
BOND_TOLERANCE = 1.1


class EQ:
    def __init__(self, eq_block_str: str) -> None:
        block_split = eq_block_str.split()
        self.name = block_split[block_split.index("Geometry") + 2] + block_split[block_split.index("Geometry") + 3][:-1]
        self.symmetry = block_split[block_split.index("SYMMETRY") + 2]
        self.n_atoms = (block_split.index("Energy") - (block_split.index("SYMMETRY") + 3)) // 4
        # Energy
        self.energy = float(block_split[block_split.index("Energy") + 2])
        self.energy_0 = float(block_split[block_split.index("Energy") + 3][1:])
        self.energy_1 = float(block_split[block_split.index("Energy") + 5][:-1])
        # Spin, ZPVE
        self.spin = float(block_split[block_split.index("Spin(**2)") + 2])
        self.zpve = float(block_split[block_split.index("ZPVE") + 2])
        # 分子情報読み込み、Molオブジェクト作成
        xyz_lines = [str(self.n_atoms), str(self.name)]
        for i in range(block_split.index("SYMMETRY") + 3, block_split.index("Energy"), 4):
            xyz_lines.append(f"{block_split[i]} {block_split[i+1]} {block_split[i+2]} {block_split[i+3]}")
        xyz_block = "\n".join(xyz_lines)
        self.mol = rdkit.Chem.RWMol(rdkit.Chem.MolFromXYZBlock(xyz_block))
        # 固有値
        self.nmode = int(block_split[block_split.index("nmode") + 2])
        self.eigenvalues = []
        for i in range(block_split.index("nmode") + 3, block_split.index("nmode") + 3 + self.nmode):
            self.eigenvalues.append(float(block_split[i]))
        # 結合を推測する
        self.determine_bonds()

    def determine_bonds(self) -> None:
        """
        原子間の距離から結合を推定する
        """
        bond_db = None
        with open(BOND_DB_PATH) as f:
            bond_db = json.load(f)
        bond_db.sort(key=lambda x: x["length"])
        distance_matrix = rdkit.Chem.rdmolops.Get3DDistanceMatrix(self.mol)
        for x in range(self.n_atoms):
            for y in range(x + 1, self.n_atoms):
                pair = set([self.mol.GetAtoms()[x].GetSymbol(), self.mol.GetAtoms()[y].GetSymbol()])
                for bond in bond_db:
                    if pair == set(bond["elements"]) and distance_matrix[x, y] <= bond["length"] * BOND_TOLERANCE:
                        # print(pair, distance_matrix[x,y], distance_matrix[x,y]/bond["length"])
                        if bond["bond_order"] == 1:
                            self.mol.AddBond(x, y, rdkit.Chem.BondType.SINGLE)
                        elif bond["bond_order"] == 2:
                            self.mol.AddBond(x, y, rdkit.Chem.BondType.DOUBLE)
                        elif bond["bond_order"] == 3:
                            self.mol.AddBond(x, y, rdkit.Chem.BondType.TRIPLE)
                        break

    def get_smiles(self) -> str:
        """
        分子構造をSMILES記法で出力する
        """
        return rdkit.Chem.MolToSmiles(self.mol)

    def __str__(self) -> str:
        return self.name


class TS(EQ):
    def __init__(self, ts_block_str: str, eq_list: list) -> None:
        super().__init__(ts_block_str)
        block_split = ts_block_str.split()
        self.con_from = None
        self.con_to = None
        con_from_name = block_split[block_split.index("CONNECTION") + 2]
        con_to_name = block_split[block_split.index("CONNECTION") + 4]
        eq_name_list = [eq.name for eq in eq_list]
        if con_from_name != "DC":
            self.con_from = eq_list[eq_name_list.index("EQ" + con_from_name)]
        if con_to_name != "DC":
            self.con_to = eq_list[eq_name_list.index("EQ" + con_to_name)]

    def has_dc_connection(self) -> bool:
        if self.con_from is None or self.con_to is None:
            return True
        return False


def select_color(u: EQ, v: EQ) -> str:
    """
    2つの状態のエネルギー差からエッジの色を決定する
    """
    delta_E = abs(u.energy - v.energy) * 627.51
    color = ""
    if delta_E < 50:
        color = "blue"
    elif delta_E < 100:
        color = "green"
    else:
        color = "#FF0000"
    return color


if __name__ == "__main__":
    # ファイルの読み込み
    eq_list = []
    ts_list = []
    with open(EQ_PATH, mode="r") as f:
        eq_blocks = f.read().split("\n\n")[1:]
        eq_list = [EQ(eq_block) for eq_block in eq_blocks]
    with open(TS_PATH, mode="r") as f:
        ts_blocks = f.read().split("\n\n")[1:]
        ts_list = [TS(ts_block, eq_list) for ts_block in ts_blocks]

    print(rdkit.Chem.rdmolfiles.MolToMolBlock(eq_list[56].mol))
    print(rdkit.Chem.rdmolfiles.MolToSmiles(eq_list[56].mol))

    # DCが含まれるTSを削除
    ts_list = [ts for ts in ts_list if ts.con_from is not None and ts.con_to is not None]
    # TSのfrom,toに含まれるEQのみを抽出
    eq_list = list(set([ts.con_from for ts in ts_list] + [ts.con_to for ts in ts_list]))

    # CSVに出力
    eq_df = pd.DataFrame(
        {
            "name": [eq.name for eq in eq_list],
            "energy": [eq.energy for eq in eq_list],
            "smiles": [eq.get_smiles() for eq in eq_list],
        }
    )
    ts_df = pd.DataFrame(
        {
            "name": [ts.name for ts in ts_list],
            "energy": [ts.energy for ts in ts_list],
            "smiles": [ts.get_smiles() for ts in ts_list],
        }
    )
    eq_ts_df = pd.concat([eq_df, ts_df], ignore_index=True)
    eq_ts_df.to_csv(SAVE_CSV_PATH)


    # グラフの構築
    G = nx.Graph()
    for eq in eq_list:
        G.add_node(eq.name, shape="box")
    for ts in ts_list:
        # TSをノードとして追加
        G.add_node(ts.name, color="gray")
        # EQからTSへのエッジ
        G.add_edge(ts.con_from.name, ts.name, color=select_color(ts.con_from, ts))
        # TSからEQへのエッジ
        G.add_edge(ts.con_to.name, ts.name, color=select_color(ts.con_to, ts))

    # グラフの可視化
    visG = pyvis.network.Network()
    visG.from_nx(G)
    visG.show_buttons(True)
    visG.inherit_edge_colors(False)
    visG.show(SAVE_GRAPH_PATH, notebook=False)
