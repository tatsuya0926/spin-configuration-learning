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
        img, temp, label = self.dataset[index]
        graph_data = self.convert_img2graph(img, top_k=5)

        return graph_data, temp, label


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


def create_train_data_hold_out(
    prm_list,
    ndata,
    T_cr_1,
    exclude_T,
    total_label,
    split_num,
    Q=None,
    T_cr_2=None,
    normalize=False
):
    train_dataset, valid_dataset, merge_valid_dataset, exclude_dataset = [], [], [], []
    if total_label == 2:
        (t_start1, t_end1) = exclude_T
    elif total_label == 3:
        (t_start1, t_end1, t_start2, t_end2) = exclude_T
    else:
        print("Please set argument:'total_label'")
    for itemp in range(len(prm_list)):
        temp, fname = prm_list[itemp]
        if total_label == 2:
            condition = temp >= t_start1 and temp <= t_end1
            if temp < T_cr_1:
                label = 0
            else:
                label = 1
        elif total_label == 3:
            condition = (temp >= t_start1 and temp <= t_end1) or (
                temp >= t_start2 and temp <= t_end2)
            if temp < T_cr_1:
                label = 0
            elif temp > T_cr_1 and temp < T_cr_2:
                label = 1
            else:
                label = 2
        else:
            print("Please set argument:'total_label'")

        for itrj in range(ndata):
            npsc = np.load(f"{fname}{itrj}.npy")
            if (normalize):
                npsc = npsc / Q
            if condition:
                exclude_dataset.append(
                    (torch.tensor(npsc, dtype=torch.float32).unsqueeze(0), temp, label))
                if itrj == split_num:
                    break
            else:
                if itrj < int(ndata*0.7):
                    train_dataset.append(
                        (torch.tensor(npsc, dtype=torch.float32).unsqueeze(0), temp, label))
                else:
                    valid_dataset.append(
                        (torch.tensor(npsc, dtype=torch.float32).unsqueeze(0), temp, label))
                    merge_valid_dataset.append(
                        (torch.tensor(npsc, dtype=torch.float32).unsqueeze(0), temp, label))
    merge_valid_dataset.extend(exclude_dataset)
    merge_valid_dataset = sorted(
        merge_valid_dataset, reverse=False, key=lambda x: x[1])

    return train_dataset, valid_dataset, merge_valid_dataset


def create_train_data_CV(
    prm_list,
    ndata,
    T_cr_1,
    exclude_T,
    total_label,
    split_num,
    Q=None,
    T_cr_2=None,
    normalize=False
):
    dataset, exclude_dataset = [], []
    if total_label == 2:
        (t_start1, t_end1) = exclude_T
    elif total_label == 3:
        (t_start1, t_end1, t_start2, t_end2) = exclude_T
    else:
        print("Please set argument:'total_label'")
    for itemp in range(len(prm_list)):
        temp, fname = prm_list[itemp]
        if total_label == 2:
            condition = temp >= t_start1 and temp <= t_end1
            if temp < T_cr_1:
                label = 0
            else:
                label = 1
        elif total_label == 3:
            condition = (temp >= t_start1 and temp <= t_end1) or (
                temp >= t_start2 and temp <= t_end2)
            if temp < T_cr_1:
                label = 0
            elif temp > T_cr_1 and temp < T_cr_2:
                label = 1
            else:
                label = 2
        else:
            print("Please set argument:'total_label'")

        for itrj in range(ndata):
            npsc = np.load(f"{fname}{itrj}.npy")
            if (normalize):
                npsc = npsc / Q
            if condition:
                exclude_dataset.append(
                    (torch.tensor(npsc, dtype=torch.float32).unsqueeze(0), temp, label))
                if itrj == split_num:
                    break
            else:
                dataset.append(
                    (torch.tensor(npsc, dtype=torch.float32).unsqueeze(0), temp, label))

    return dataset, exclude_dataset


def inference(total_test_dataset, temps, prediction_test, target_size):
    xs, y1s, y2s, y3s = [], [], [], []
    sum_pred_0, sum_pred_1, sum_pred_2 = 0, 0, 0
    count = 0

    for i in range(total_test_dataset):
        if i == 0:
            sum_pred_0, sum_pred_1, sum_pred_2 = __pred_count(
                sum_pred_0, sum_pred_1, sum_pred_2, prediction_test[i])
            count += 1
            xs.append(temps[i])
        else:
            if temps[i] != temps[i-1]:
                # 格納
                y1s.append(sum_pred_0/count)
                y2s.append(sum_pred_1/count)
                if target_size == 3:
                    y3s.append(sum_pred_2/count)

                sum_pred_0, sum_pred_1, sum_pred_2 = 0, 0, 0
                count = 0
                sum_pred_0, sum_pred_1, sum_pred_2 = __pred_count(
                    sum_pred_0, sum_pred_1, sum_pred_2, prediction_test[i])
                count += 1
                xs.append(temps[i])
            elif i == total_test_dataset-1:
                # 格納
                y1s.append(sum_pred_0/count)
                y2s.append(sum_pred_1/count)
                if target_size == 3:
                    y3s.append(sum_pred_2/count)
            else:
                sum_pred_0, sum_pred_1, sum_pred_2 = __pred_count(
                    sum_pred_0, sum_pred_1, sum_pred_2, prediction_test[i])
                count += 1
    if target_size == 2:
        return np.array(xs), np.array(y1s), np.array(y2s)
    elif target_size == 3:
        return np.array(xs), np.array(y1s), np.array(y2s), np.array(y3s)


def __pred_count(sum_pred_0, sum_pred_1, sum_pred_2, prediction_test):
    if prediction_test == 0:
        sum_pred_0 += 1
    elif prediction_test == 1:
        sum_pred_1 += 1
    else:
        sum_pred_2 += 1
    return sum_pred_0, sum_pred_1, sum_pred_2
