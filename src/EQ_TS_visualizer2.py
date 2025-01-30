"""
EQ,TS 可視化プログラム
pandas,pyvis をpipでインストールすること
"""
import networkx as nx
import pandas as pd
import pyvis

EQ_path = "manu_EQ_list.log"
TS_path = "manu_TS_list.log"


def read_EQ_block(raw_text: str) -> pd.DataFrame:
    """
    単一のEQに対応するの文字列のブロックを渡すと、pandas.Dataframeに格納する
    """
    # 値の取り出しのためスペースでsplit
    # 例: ['#', 'Geometry', 'of', 'EQ', '0,', 'SYMMETRY', '=', 'C2v', 'C', '-0.000112000000', '1.198714000000', '0.000000000000', 'H', '-0.000007000000', '1.750081000000', '0.933062000000', 'H', '-0.000007000000', '1.750081000000', '-0.933062000000', 'C', '0.000000000000', '-0.182179000000', '0.000000000000', 'N', '0.000098000000', '-1.371339000000', '0.000000000000', 'Energy', '=', '-132.055704962587', '(-132.055704962587', ':', '0.000000000000)', 'Spin(**2)', '=', '0.770506670235', 'ZPVE', '=', '0.031274079270', 'Normal', 'mode', 'eigenvalues', ':', 'nmode', '=', '9', '0.005453116', '0.006931664', '0.019747598', '0.042820693', '0.043936614', '0.083153588', '0.165421172', '0.384847516', '0.410750967']
    raw_text_split = raw_text.split()

    # 結果を格納するDataframe
    df = pd.DataFrame()

    # EQの識別番号 - "EQ" に続く "0," の末尾カンマを除去して整数化
    df["id"] = [raw_text_split[raw_text_split.index("Geometry") + 3][:-1]]
    # "SYMMETRY", "=", "C2v" なので+2
    df["symmetry"] = raw_text_split[raw_text_split.index("SYMMETRY") + 2]
    # 各元素のデータ取り出し
    for element_id, index in enumerate(range(raw_text_split.index("SYMMETRY") + 3, raw_text_split.index("Energy"), 4)):
        df[f"element{element_id}"] = raw_text_split[index]
        df[f"element{element_id}_0"] = float(raw_text_split[index + 1])
        df[f"element{element_id}_1"] = float(raw_text_split[index + 2])
        df[f"element{element_id}_2"] = float(raw_text_split[index + 3])
    # Energy
    df["energy"] = float(raw_text_split[raw_text_split.index("Energy") + 2])
    df["energy_0"] = float(raw_text_split[raw_text_split.index("Energy") + 3][1:])
    df["energy_1"] = float(raw_text_split[raw_text_split.index("Energy") + 5][:-1])
    # Spin, ZPVE
    df["spin"] = float(raw_text_split[raw_text_split.index("Spin(**2)") + 2])
    df["zpve"] = float(raw_text_split[raw_text_split.index("ZPVE") + 2])
    # Normal mode eigenvalues
    nmode = int(raw_text_split[raw_text_split.index("nmode") + 2])
    df["nmode"] = nmode
    for ev_id, index in enumerate(range(raw_text_split.index("nmode") + 3, raw_text_split.index("nmode") + 3 + nmode)):
        df[f"eigenvalue{ev_id}"] = float(raw_text_split[index])
    return df


def read_TS_block(raw_text: str) -> pd.DataFrame:
    """
    単一のTSに対応するの文字列のブロックを渡すと、pandas.Dataframeに格納する
    EQと同様だが、CONNECTIONの情報を追加する
    """
    df = read_EQ_block(raw_text)
    raw_text_split = raw_text.split()
    connection_from = raw_text_split[raw_text_split.index("CONNECTION") + 2]
    if connection_from == "DC":
        df["from"] = None
    else:
        df["from"] = connection_from
    connection_to = raw_text_split[raw_text_split.index("CONNECTION") + 4]
    if connection_to == "DC":
        df["to"] = None
    else:
        df["to"] = connection_to
    return df


def open_EQ(path) -> pd.DataFrame:
    """
    EQのファイルを開いてDataframeを作成する
    """
    result = pd.DataFrame()
    with open(path, mode="r") as f:
        eq_blocks = f.read().split("\n\n")[1:]
        result = pd.concat([read_EQ_block(eq_block) for eq_block in eq_blocks])
    result["id"] = "EQ"+result["id"]
    result.set_index("id", drop=False, inplace=True)
    return result


def open_TS(path) -> pd.DataFrame:
    """
    TSのファイルを開いてDataframeを作成する
    """
    result = pd.DataFrame()
    with open(path, mode="r") as f:
        ts_blocks = f.read().split("\n\n")[1:]
        result = pd.concat([read_TS_block(ts_block) for ts_block in ts_blocks])
    result["id"] = "TS"+result["id"]
    result["from"] = "EQ"+result["from"]
    result["to"] = "EQ"+result["to"]
    result.set_index("id", drop=False, inplace=True)
    return result


if __name__ == "__main__":
    # ファイルの読み込み
    eq_df = open_EQ(EQ_path)
    print(eq_df)
    ts_df = open_TS(TS_path)
    print(ts_df)

    # おまけ Excel出力 pipでopenpyxlを入れると動く
    eq_df.to_excel(EQ_path+".xlsx")
    ts_df.to_excel(TS_path+".xlsx")

    # グラフ用に要素を抽出
    # from,toの組が重複しているエッジを削除
    edge_df = ts_df.drop_duplicates(subset=["from", "to"])
    # None(TS_listではDC)を含むエッジを削除
    edge_df = edge_df[edge_df.notnull().all(axis=1)]
    # EQ0-TS0-EQ1 を EQ0-TS0,EQ1-TS0 に変換
    edge_df1 = edge_df.loc[:, ["from", "id"]]
    edge_df1.columns = ["from", "to"]
    edge_df2 = edge_df.loc[:, ["to", "id"]]
    edge_df2.columns = ["from", "to"]
    edge_df = pd.concat([edge_df1, edge_df2])
    edge_df.reset_index(drop=True, inplace=True)
    # エネルギー差 算出
    edge_df["E_EQ"] = eq_df.loc[edge_df["from"], ["energy"]].reset_index(drop=True)
    edge_df["E_TS"] = ts_df.loc[edge_df["to"], ["energy"]].reset_index(drop=True)
    # hartree to kcal
    edge_df["delta_E"] = (edge_df["E_TS"] - edge_df["E_EQ"]) * 627.51

    # 色を決定する関数
    def select_color(delta_E) -> str:
        color = ""
        if delta_E < 50:
            color = "blue"
        elif delta_E < 100:
            color = "green"
        else:
            color = "#FF0000"
        return color

    edge_df["color"] = edge_df["delta_E"].apply(select_color)
    edge_df.to_excel("manu_Edge_list.xlsx")
    # print(edge_df)

    node_df1 = eq_df
    node_df1["shape"] = "box"
    node_df2 = ts_df
    node_df2["color"] = "gray"
    node_df = pd.concat([node_df1, node_df2])
    node_df.to_excel("manu_Node_list.xlsx")
    # print(node_df)

    G = nx.from_pandas_edgelist(edge_df, source="from", target="to", edge_attr=True)
    nx.set_node_attributes(G, node_df.set_index("id").to_dict(orient="index"))
    visG = pyvis.network.Network()
    visG.from_nx(G)
    visG.show_buttons(True)
    visG.inherit_edge_colors(False)
    visG.show("pyvis.html", notebook=False)

