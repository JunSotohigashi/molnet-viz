"""
EQ,TS 可視化プログラム
"""

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
SAVE_CSV_PATH = "./out/manu.csv"
SAVE_GRAPH_PATH = "./out/manu.html"
SHOW_TS = False


def arrange_nodes_by_smiles(G: nx.Graph) -> nx.Graph:
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

    # カラーマップを作製
    color_map = plt.get_cmap("turbo")
    color_norm = matplotlib.colors.Normalize(
        vmin=min(
            min([eq.energy for eq in eq_list]),
            min([ts.energy for ts in ts_list]),
        ),
        vmax=max(
            max([eq.energy for eq in eq_list]),
            max([ts.energy for ts in ts_list]),
        ),
    )
    fig, ax = plt.subplots(figsize=(1, 5))
    matplotlib.colorbar.Colorbar(
        mappable=matplotlib.cm.ScalarMappable(color_norm, color_map),
        ax=ax,
        orientation="vertical",
    ).set_label("energy")
    plt.savefig("./out/colorbar.svg", format="svg", bbox_inches="tight")

    # グラフの構築
    print("Visualizing...")
    G = nx.Graph()
    # 指定した色・形状でノードとエッジを追加
    for eq in eq_list:
        pass
        G.add_node(
            eq.name,
            shape="square",
            color=matplotlib.colors.rgb2hex(color_map(color_norm(eq.energy))),
            size=20,
            font="25px arial black",
            smiles=eq.smiles,
            title=f"{eq.name}\n{eq.smiles}",
        )
    if SHOW_TS:
        for ts in ts_list:
            # TSをノードとして追加
            G.add_node(
                ts.name,
                color=matplotlib.colors.rgb2hex(
                    color_map(color_norm(ts.energy))
                ),
                size=20,
                font="25px arial black",
                smiles=ts.smiles,
                title=f"{eq.name}\n{eq.smiles}",
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
        "colorbar", shape="image", image="./colorbar.svg", size=120, label=""
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
                "enabled": true,
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
    visG.show(SAVE_GRAPH_PATH, notebook=False)

    # CSVとして保存
    eq_df, ts_df = convert_to_df(eq_list, ts_list)
    pd.concat([eq_df, ts_df], ignore_index=True).to_csv(SAVE_CSV_PATH)
    print(f"CSV file has saved to {SAVE_CSV_PATH}")

    print("Finished!")
