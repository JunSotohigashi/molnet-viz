"""
EQ,TS 可視化プログラム
pandas,pyvis,rdkit
"""

import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.colorbar
import matplotlib.colors
import pandas as pd
import pyvis

from EQ_TS_loader import load, EQ

# 読み込むファイルのパス
EQ_PATH = "./data/manu_EQ_list.log"
TS_PATH = "./data/manu_TS_list.log"
# 保存先のパス
SAVE_CSV_PATH = "./out/manu.csv"
SAVE_GRAPH_PATH = "./out/manu.html"


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


def arrange_nodes_by_smiles(G):
    """networkx のグラフ G を、smiles ごとにクラスタ化してレイアウトを適用する"""

    # クラスタリング（smiles ごとにノードをグループ化）
    clusters = {}
    for node in G.nodes:
        smiles = G.nodes[node].get("smiles", "unknown")
        clusters.setdefault(smiles, []).append(node)

    # クラスタごとに円形レイアウトを適用
    layout = {}
    cluster_positions = {}

    for smiles, cluster_nodes in clusters.items():
        subgraph = G.subgraph(cluster_nodes)
        circular_pos = nx.circular_layout(subgraph, scale=50)  # 円形配置

        # クラスタの仮の中心座標（後でspring_layoutで調整）
        cluster_positions[smiles] = (0, 0)

        for node, (x, y) in circular_pos.items():
            layout[node] = (x, y)

    # クラスタ同士を配置するためのspring_layoutを適用
    cluster_graph = nx.Graph()
    for smiles in clusters.keys():
        cluster_graph.add_node(smiles)

    cluster_pos = nx.spring_layout(
        cluster_graph, scale=1000
    )  # クラスタ全体を適度に分散

    # 各クラスタの中心座標を設定
    for smiles, (cx, cy) in cluster_pos.items():
        cluster_positions[smiles] = (cx, cy)

    # 各ノードの座標を調整
    for smiles, cluster_nodes in clusters.items():
        cx, cy = cluster_positions[smiles]  # クラスタ中心
        for node in cluster_nodes:
            x, y = layout[node]  # 円形内の位置
            layout[node] = (cx + x, cy + y)  # 全体の位置に調整

    # networkx のグラフに座標を設定
    nx.set_node_attributes(
        G,
        {
            node: {"x": layout[node][0], "y": layout[node][1]}
            for node in G.nodes()
        },
    )

    return G


if __name__ == "__main__":

    # ファイルの読み込み
    eq_list, ts_list = load(EQ_PATH, TS_PATH)

    # 水素の結合先原子の変化を抽出
    print("Counting moved H atoms...")
    # 注目する原子のインデックス C,H,H,C,Nの順なので 水素原子は1,2
    target_atoms = [1, 2]
    ts_atom_moved = []
    for ts in ts_list:
        atom_moved = []
        # for target_atom in target_atoms:
        #     # from側の
        #     from_eq = ts.con_from.mol.GetAtoms()[target_atom]
        #     from_eq_bond = None
        #     if from_eq.GetBonds():
        #         if from_eq.GetBonds()[0].GetBeginAtomIdx() == target_atom:
        #             from_eq_bond = from_eq.GetBonds()[0].GetEndAtomIdx()
        #         else:
        #             from_eq_bond = from_eq.GetBonds()[0].GetBeginAtomIdx()
        #     to_eq = ts.con_to.mol.GetAtoms()[target_atom]
        #     to_eq_bond = None
        #     if to_eq.GetBonds():
        #         if to_eq.GetBonds()[0].GetBeginAtomIdx() == target_atom:
        #             to_eq_bond = to_eq.GetBonds()[0].GetEndAtomIdx()
        #         else:
        #             to_eq_bond = to_eq.GetBonds()[0].GetBeginAtomIdx()
        #     atom_moved.append(from_eq_bond != to_eq_bond)
        ts_atom_moved.append(atom_moved.count(True))

    # CSVに出力
    print("Saving data as CSV...")
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
            "from_EQ": [ts.con_from.name for ts in ts_list],
            "from_EQ_smiles": [ts.con_from.get_smiles() for ts in ts_list],
            "to_EQ": [ts.con_to.name for ts in ts_list],
            "to_EQ_smiles": [ts.con_to.get_smiles() for ts in ts_list],
            "H_moved": ts_atom_moved,
        }
    )
    eq_ts_df = pd.concat([eq_df, ts_df], ignore_index=True)
    eq_ts_df.to_csv(SAVE_CSV_PATH)
    print(f"CSV file has saved to {SAVE_CSV_PATH}")

    # グラフの構築
    print("Visualizing...")
    G = nx.Graph()

    # カラーマップを作製
    color_map = plt.get_cmap("turbo")
    color_norm = matplotlib.colors.Normalize(
        vmin=min(min(eq_df["energy"]), min(ts_df["energy"])),
        vmax=max(max(eq_df["energy"]), max(ts_df["energy"])),
    )
    fig, ax = plt.subplots(figsize=(1, 5))
    matplotlib.colorbar.Colorbar(
        mappable=matplotlib.cm.ScalarMappable(color_norm, color_map),
        ax=ax,
        orientation="vertical",
    ).set_label("energy")
    plt.savefig("./out/colorbar.svg", format="svg", bbox_inches="tight")

    # 指定した色・形状でノードとエッジを追加
    for eq in eq_list:
        G.add_node(
            eq.name,
            shape="square",
            color=matplotlib.colors.rgb2hex(color_map(color_norm(eq.energy))),
            size=20,
            font="25px arial black",
            smiles=eq.get_smiles(),
        )
    for ts in ts_list:
        # # TSをノードとして追加
        # G.add_node(
        #     ts.name,
        #     color=matplotlib.colors.rgb2hex(color_map(color_norm(ts.energy))),
        #     size=20,
        #     font="25px arial black",
        #     smiles=ts.get_smiles(),
        # )
        # # EQからTSへのエッジ
        # G.add_edge(ts.con_from.name, ts.name, color="black")
        # # TSからEQへのエッジ
        # G.add_edge(ts.con_to.name, ts.name, color="black")
        # EQ間のエッジ
        G.add_edge(
            ts.con_from.name,
            ts.con_to.name,
        )

    G = arrange_nodes_by_smiles(G)

    # すべてのノードの配置が完了後、カラーマップをグラフに追加
    G.add_node(
        "colorbar", shape="image", image="./colorbar.svg", size=120, label=""
    )

    # グラフの可視化
    visG = pyvis.network.Network(height="1080px", width="1920px")
    visG.from_nx(G)
    # visG.show_buttons(True)
    # visG.inherit_edge_colors(False)
    visG.set_options(
        """
        var options = {
            "physics": {
                "enabled": false
            }
        }
        """
    )
    visG.show(SAVE_GRAPH_PATH, notebook=False)

    print("Finished!")
