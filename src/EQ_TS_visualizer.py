"""
EQ,TS 可視化プログラム
"""

import os
from datetime import datetime

import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.colorbar
import matplotlib.colors
import pandas as pd
import pyvis

from EQ_TS_loader import load, convert_to_df

# 読み込むファイルのパス
EQ_PATH = "./data/manu_EQ_list.log"
TS_PATH = "./data/manu_TS_list.log"
# 保存先のパス
SAVE_DIR = f"./out/{datetime.now().strftime('%Y%m%d_%H%M%S')}"
# TSを表示するかどうか
SHOW_TS = False


class ColorMapper:
    def __init__(self, target_data: list):
        self.color_map = plt.get_cmap("turbo")
        self.color_norm = matplotlib.colors.Normalize(
            vmin=min(target_data), vmax=max(target_data)
        )

    def get_color(self, value: float) -> str:
        return matplotlib.colors.rgb2hex(self.color_map(self.color_norm(value)))

    def save_colorbar(self, save_path: str, label: str = ""):
        fig, ax = plt.subplots(figsize=(1, 5))
        matplotlib.colorbar.Colorbar(
            mappable=matplotlib.cm.ScalarMappable(self.color_norm, self.color_map),
            ax=ax,
            orientation="vertical",
        ).set_label(label)
        fig.savefig(
            save_path,
            format="svg",
            bbox_inches="tight",
        )


def arrange_nodes_by_smiles(G: nx.Graph) -> nx.Graph:
    """networkx のグラフ G を、smiles ごとにクラスタ化してレイアウトを適用する"""

    # クラスタリング（smiles ごとにノードをグループ化）
    clusters = {}
    for node in G.nodes:
        smiles = G.nodes[node]["smiles"]
        if smiles not in clusters:
            clusters[smiles] = {
                "nodes": [node],
                "atoms": G.nodes[node]["atoms"],
            }
        else:
            clusters[smiles]["nodes"].append(node)

    # クラスタごとに円形レイアウトを適用
    layout = {}
    cluster_positions = {}

    for smiles, cluster in clusters.items():
        subgraph = G.subgraph(cluster["nodes"])
        circular_pos = nx.circular_layout(
            subgraph,
            scale=(len(cluster["nodes"]) - 1) ** 0.5 * 40.0,
        )

        # クラスタの仮の中心座標（後でspring_layoutで調整）
        cluster_positions[smiles] = (0, 0)

        for node, (x, y) in circular_pos.items():
            layout[node] = (x, y)

    # クラスタ同士を配置するためのspring_layoutを適用
    cluster_graph = nx.Graph()
    for smiles, cluster in clusters.items():
        cluster_graph.add_node(smiles)
        # 同じ原子数のクラスタ同士をエッジで結ぶ
        for other_smiles, other_cluster in clusters.items():
            if smiles != other_smiles and cluster["atoms"] == other_cluster["atoms"]:
                cluster_graph.add_edge(smiles, other_smiles, weight=300.0)
            elif smiles != other_smiles:
                cluster_graph.add_edge(smiles, other_smiles, weight=0.1)

    cluster_pos = nx.kamada_kawai_layout(
        cluster_graph,
        scale=(cluster_graph.number_of_nodes() - 1) ** 0.5 * 150.0,
    )  # クラスタ全体を適度に分散

    # 各クラスタの中心座標を設定
    for smiles, (cx, cy) in cluster_pos.items():
        cluster_positions[smiles] = (cx, cy)

    # 各ノードの座標を調整
    for smiles, cluster in clusters.items():
        cx, cy = cluster_positions[smiles]  # クラスタ中心
        for node in cluster["nodes"]:
            x, y = layout[node]  # 円形内の位置
            layout[node] = (cx + x, cy + y)  # 全体の位置に調整

    # networkx のグラフに座標を設定
    nx.set_node_attributes(
        G,
        {node: {"x": layout[node][0], "y": layout[node][1]} for node in G.nodes()},
    )
    return G


if __name__ == "__main__":
    # 保存先のディレクトリを作成
    os.makedirs(SAVE_DIR, exist_ok=True)

    # ファイルの読み込み
    eq_list, ts_list = load(EQ_PATH, TS_PATH)

    # カラーマップを作製
    color_map = ColorMapper(
        [eq.energy for eq in eq_list] + [ts.energy for ts in ts_list]
    )
    color_map.save_colorbar(SAVE_DIR + "/colorbar.svg", "Energy")

    # グラフの構築
    G = nx.Graph()
    # 指定した色・形状でノードとエッジを追加
    for eq in eq_list:
        G.add_node(
            eq.name,
            shape="square",
            color=color_map.get_color(eq.energy),
            size=16,
            font="25px arial black",
            title=f"{eq.name}\n{eq.energy}\n{eq.smiles}\n{str(eq.atoms_in_fragments)}",
            smiles=eq.smiles,
            atoms=eq.atoms_in_fragments,
        )
    if SHOW_TS:
        for ts in ts_list:
            # TSをノードとして追加
            G.add_node(
                ts.name,
                color=color_map.get_color(ts.energy),
                size=16,
                font="25px arial black",
                title=f"{ts.name}\n{ts.energy}\n{ts.smiles}\n{str(ts.atoms_in_fragments)}",
                smiles=ts.smiles,
                atoms=ts.atoms_in_fragments,
            )
            # EQからTSへのエッジ
            G.add_edge(ts.con_from.name, ts.name, color="black")
            # TSからEQへのエッジ
            G.add_edge(ts.con_to.name, ts.name, color="black")
    else:
        for ts in ts_list:
            # EQ間のエッジ
            G.add_edge(
                ts.con_from.name,
                ts.con_to.name,
            )

    G = arrange_nodes_by_smiles(G)

    # すべてのノードの配置が完了後、カラーマップをグラフに追加
    G.add_node(
        "colorbar",
        shape="image",
        image="./colorbar.svg",
        size=120,
        label=" ",
        x=1300,
        y=0,
    )

    # グラフの可視化
    visG = pyvis.network.Network(
        height="100vh",
        width="100vw",
    )
    visG.from_nx(G)
    visG.set_options(
        """
        var options = {
            "configure": {
                "enabled": false,
                "filter": true
            },
            "edges": {
                "color": {
                    "inherit": false
                },
                "smooth": {
                    "enabled": false
                }
            },
            "physics": {
                "enabled": false
            }
        }
        """
    )
    visG.show(SAVE_DIR + "/out.html", notebook=False)

    # CSVとして保存
    eq_df, ts_df = convert_to_df(eq_list, ts_list)
    pd.concat([eq_df, ts_df], ignore_index=True).to_csv(SAVE_DIR + "/out.csv")

    print("Done!")
