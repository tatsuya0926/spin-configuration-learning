import numpy as np
import torch
from decimal import Decimal
from torch_geometric.data import Data


class Mydataset(torch.utils.data.Dataset):
    def __init__(self, dataset):
        self.dataset = dataset

    def __len__(self):
        return len(self.dataset)

    def convert_img2graph(self, img, top_k=5):  # TODO(nodeとedgeの決め方)
        mask = img != -1  # マスク準備
        x = img.view(-1, 1).float()  # 輝度情報の特徴量化
        x = x[mask.view(-1)]
        coords = torch.stack(torch.meshgrid(torch.arange(
            img.size(1)), torch.arange(img.size(2))), dim=-1)  # xy座標情報の特徴量
        coords = coords.view(-1, 2)
        coords = coords[mask.view(-1)]
        x = torch.cat([x, coords.float()], dim=1)  # 特徴量の作成

        distance_matrix = torch.cdist(x, x, p=2)  # ユークリッド距離の計算
        _, indices = torch.topk(distance_matrix, k=top_k,
                                largest=False)  # 類似度がTOP5のノードを取得
        edge_from = torch.arange(
            x.size(0)).view(-1, 1).repeat(1, top_k)  # エッジインデックスの生成
        edge_to = indices
        edge_index = torch.stack([edge_from.view(-1), edge_to.view(-1)], dim=0)
        edge_index = torch.cat(
            [edge_index, edge_index.flip(0)], dim=1)  # 双方向エッジ
        return Data(x=x, y=None, edge_index=edge_index)

    def __getitem__(self, index):
        img, label = self.dataset[index]
        graph_data = self.convert_img2graph(img, top_k=5)

        return graph_data, label


def create_param_list(nconf, t_start, L, model_name, q=None):
    prm_list = []
    t_start = Decimal(str(t_start))
    for i in range(nconf):
        if q == None:
            filename = f"../dataset/{model_name}/L{L}/L{L}T{i}_"
        else:
            filename = f"../dataset/{model_name}/L{L}_q={q}/L{L}T{i}_"
        prm_list.append([float(t_start), filename])
        t_start += Decimal("0.01")
    t_end = t_start

    return prm_list, t_end


def create_train_data_hold_out(prm_list, ndata, T_cr, exclude_T=None):
    train_dataset, valid_dataset, exclude_dataset = [], [], []
    if exclude_T != None: 
        (t_start, t_end) = exclude_T
    for itemp in range(len(prm_list)):
        temp, fname = prm_list[itemp]
        if temp < T_cr:
            label = 0
        else:
            label = 1

        for itrj in range(ndata):
            npsc = np.load(f"{fname}{itrj}.npy")
            if temp >= t_start and temp <= t_end:
                exclude_dataset.append((torch.tensor(npsc, dtype=torch.float32).unsqueeze(0), temp, label))
                if itrj == 99: break
            else:
                if itrj < 150:
                    train_dataset.append((torch.tensor(npsc, dtype=torch.float32).unsqueeze(0), temp, label))
                else:
                    valid_dataset.append((torch.tensor(npsc, dtype=torch.float32).unsqueeze(0), temp, label))
    valid_dataset.extend(exclude_dataset)
    valid_dataset = sorted(valid_dataset, reverse=False, key=lambda x: x[1])

    return train_dataset, valid_dataset

def create_train_data_CV(prm_list, ndata, T_cr_1, T_cr_2=None, exclude_T=None):
    dataset, exclude_dataset = [], []
    if exclude_T != None: 
        (t_start, t_end) = exclude_T
    for itemp in range(len(prm_list)):
        temp, fname = prm_list[itemp]
        if T_cr_2 == None:
            if temp < T_cr_1:
                label = 0
            else:
                label = 1
        else:
            if temp < T_cr_1:
                label = 0
            elif temp > T_cr_1 and temp < T_cr_2:
                label = 1
            else:
                label = 2

        for itrj in range(ndata):
            npsc = np.load(f"{fname}{itrj}.npy")
            if temp >= t_start and temp <= t_end:
                exclude_dataset.append((torch.tensor(npsc, dtype=torch.float32).unsqueeze(0), temp, label))
                if itrj == 99: break
            else:
                dataset.append((torch.tensor(npsc, dtype=torch.float32).unsqueeze(0), temp, label))

    return dataset, exclude_dataset
