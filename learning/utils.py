from logging import exception
import torch
from torch_geometric.data import Data


class Mydataset(torch.utils.data.Dataset):
    def __init__(self, dataset):
        self.dataset = dataset

    def __len__(self):
        return len(self.dataset)

    def convert_img2graph(self, img, top_k=5): # TODO(nodeとedgeの決め方)
        mask = img != -1 # マスク準備
        x = img.view(-1, 1).float() # 輝度情報の特徴量化
        x = x[mask.view(-1)]
        coords = torch.stack(torch.meshgrid(torch.arange(img.size(1)), torch.arange(img.size(2))), dim=-1) # xy座標情報の特徴量
        coords = coords.view(-1, 2)
        coords = coords[mask.view(-1)]
        x = torch.cat([x, coords.float()], dim=1) # 特徴量の作成

        distance_matrix = torch.cdist(x, x, p=2) # ユークリッド距離の計算
        _, indices = torch.topk(distance_matrix, k=top_k, largest=False) # 類似度がTOP5のノードを取得
        edge_from = torch.arange(x.size(0)).view(-1, 1).repeat(1, top_k) # エッジインデックスの生成
        edge_to = indices
        edge_index = torch.stack([edge_from.view(-1), edge_to.view(-1)], dim=0)
        edge_index = torch.cat([edge_index, edge_index.flip(0)], dim=1)  # 双方向エッジ
        return Data(x=x, y=None, edge_index=edge_index)

    def __getitem__(self, index):
        img, label = self.dataset[index]
        graph_data = self.convert_img2graph(img, top_k = 5)

        return graph_data, label