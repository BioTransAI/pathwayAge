## notwork data preparsion

import json
import pandas as pd 
import numpy as np
import json
from pyecharts import options as opts
from pyecharts.charts import Graph, Grid
from sklearn.preprocessing import MinMaxScaler
from snapshot_selenium import snapshot as driver
from pyecharts.render import make_snapshot


#  calculate the correlation between the GOs and prepare data for plot

def moduleDataForPlot():
    colmap = {
        6: "Cytokine Production Regulation", 
        24: "Homeostasis Response",
        19: "Transduction Regulation",
        25: "DNA repair",
        10: "Regulation of Cellular Organelle Organization",
        22: "Positive Regulation of Protein Phosphorylation",
        26: "Negative Regulation of mRNA Catabolic Process",
        27: "Regulation of Cellular Process"

    }

    GOs = pd.read_csv("./Result/GOWGCNACluster.csv", index_col = 0)
    GOs["Tag"] = GOs['Cluster'].replace(colmap)


    AgeGapPerGOPath = "./result/AgeAccPerGO.csv"
    AgeGapPerGO = pd.read_csv(AgeGapPerGOPath, index_col=0)
    AgeGapPerGO = AgeGapPerGO[list(GOs.GO.values)]
    heatmap = AgeGapPerGO.corr()

    scaler = MinMaxScaler()
    normalized_data = scaler.fit_transform(heatmap)
    normalizd_heatmap = pd.DataFrame(normalized_data, columns=heatmap.columns, index=heatmap.index)

    upper_triangle_values =  normalizd_heatmap.mask(np.triu(np.ones(heatmap.shape)).astype(bool)).fillna(0)
    nodes = GOs[["GO", "Tag"]].sort_values(by=["Tag"]).rename(columns={"GO": "id", "Tag": "category"})
    nodes["name"] = nodes["id"]

    nodes['category'] = nodes['category'].replace(colmap)

    link = []
    linkedNodes = []
    for index, row in upper_triangle_values.iterrows():
        for column in upper_triangle_values.columns:
            if row[column] < 0 or row[column] > 0.35:
                dict = {}
                # print(f"Row {index}, Column {column}: Value is not equal to 1 ({row[column]})")
                dict["source"] = index 
                dict["target"] = column
                dict["value"] =  row[column]
                linkedNodes.append(index)
                linkedNodes.append(column)
                link.append(dict)

    # remove the not linked nodes
    link_df = pd.DataFrame(link)
    # Count the occurrences of each node in both the 'source' and 'target' columns
    source_count = link_df['source'].value_counts()
    target_count = link_df['target'].value_counts()
    # Combine the counts of 'source' and 'target'
    node_count = source_count.add(target_count, fill_value=0)
    node_count.name = "symbolSize"

    nodes = nodes[nodes.id.isin(linkedNodes)]
    nodes.set_index("id", inplace = True)
    nodes = nodes.join(node_count)

    def categorize(x):
        if x < 4:
            return 1
        elif 4 < x < 8:
            return 5
        elif 8 < x < 14:
            return 10
        elif 14 < x:
            return 15
            
    nodes['symbolSize'] = nodes['symbolSize'].apply(categorize)
    nodes = nodes.to_dict(orient='records')

    result = {}
    result["nodes"] = nodes
    result["links"] = link
    result["categories"] = GOs[["Tag"]].drop_duplicates(keep="first").rename(columns = {"Tag": "name"}).to_dict(orient='records')

    file_name = "./Result/Aging.json"
    print(len(result["links"]))

    with open(file_name, "w") as json_file:
        json.dump(result, json_file)
     
    return "Modules data for plot been saved."


def load_data(file_path):
    with open(file_path, "r", encoding="utf-8") as f:
        data = json.load(f)
        return data["nodes"], data["links"], data["categories"]

nodes1, links1, categories1 = load_data("./Result/Aging.json")
nodes2, links2, categories2 = load_data("./Result/AgingKEGG.json")

def create_network(nodes, links, categories, legend_pos):
    node_category_map = {node['name']: node['category'] for node in nodes}
    all_links = []
    same_category_links = []

    for link in links:
        source_category = node_category_map.get(link['source'])
        target_category = node_category_map.get(link['target'])
        
        all_links.append(link)
        
        if source_category and target_category and source_category == target_category:
            same_category_links.append(link)

    network = (
        Graph(init_opts=opts.InitOpts(width="500px", height="500px", bg_color="white"))
        .add(
            "",
            nodes=nodes,
            links=all_links,
            categories=categories,
            layout="circular",
            is_rotate_label=True,
            linestyle_opts=opts.LineStyleOpts(color="source", curve=0.7, opacity=0.5),
            gravity=10,  
            repulsion=5000,  
        )
        .set_global_opts(
            toolbox_opts=opts.ToolboxOpts(is_show=True),
            legend_opts=opts.LegendOpts(
                orient="vertical",
                pos_left=legend_pos[0],
                pos_top=legend_pos[1],
                textstyle_opts=opts.TextStyleOpts(font_size=20)
            ),
        )
        .set_series_opts(
            node_scale_ratio=0.5,
            label_opts=opts.LabelOpts(position="right", is_show=False, font_size=10),
        )
    )

    for link in same_category_links:
        link['lineStyle'] = {'opacity': 0.5, 'curveness': 0.2}

    return network

network1 = create_network(nodes1, links1, categories1, legend_pos=("20%", "72%"))  # 图例1位置
network2 = create_network(nodes2, links2, categories2, legend_pos=("60%", "72%"))  # 图例2位置

grid = (
    Grid(init_opts=opts.InitOpts(width="2000px", height="900px", bg_color="white"))
    .add(network1, grid_opts=opts.GridOpts(pos_left="0%", pos_right="50%", pos_top="0%", pos_bottom="30%"))  # 设置网络图1的位置和底部空间
    .add(network2, grid_opts=opts.GridOpts(pos_left="50%", pos_right="0%", pos_top="0%", pos_bottom="30%"))  # 设置网络图2的位置和底部空间
)

grid.render("combined_network.html")

make_snapshot(driver, grid.render(), "./Result/AgingModule.png", pixel_ratio=3)

