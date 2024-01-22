import torch
import torch.nn as nn
from torch_geometric.nn import GCNConv, global_mean_pool


class CNNClassifier(nn.Module):
    def __init__(self):
        super(CNNClassifier, self).__init__()
        self.conv1 = nn.Conv2d(1, 64, 3, 1)
        self.conv2 = nn.Conv2d(64, 128, 3, 1)
        self.conv3 = nn.Conv2d(128, 256, 3, 1)
        self.dropout1 = nn.Dropout(0.25)
        self.dropout2 = nn.Dropout(0.5)
        self.fc1 = nn.Linear(12*12*256, 32)
        self.fc2 = nn.Linear(32, 2)
        self.relu = nn.ReLU()

    def forward(self, x):
        x = self.relu(self.conv1(x))
        x = nn.MaxPool2d(2)(x)
        x = self.relu(self.conv2(x))
        x = nn.MaxPool2d(2)(x)
        x = self.dropout1(x)
        x = self.relu(self.conv3(x))
        x = torch.flatten(x, 1)
        x = self.relu(self.fc1(x))
        x = self.dropout2(x)
        x = self.fc2(x)

        return x


class GCNClassifier(nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = GCNConv(3, 16)
        self.conv2 = GCNConv(16, 64)
        self.fc1 = nn.Linear(64, 32)
        self.fc2 = nn.Linear(32, 2)
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(p=0.5)

    def edge_reconnect(self, x, k):
        distance_matrix = torch.cdist(x, x, p=2)
        _, indices = torch.topk(distance_matrix, k)
        edge_from = torch.arange(x.shape[0]).view(-1, 1).repeat(1, k)
        edge_to = indices
        edge_index = torch.stack([edge_from.view(-1), edge_to.view(-1)], dim=0)
        return edge_index

    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        batch = data.batch

        x = self.relu(self.conv1(x, edge_index))
        # edge_index = self.edge_reconnect(x, k = 5)
        x = self.relu(self.conv2(x, edge_index))
        # edge_index = self.edge_reconnect(x, k = 5)
        x = global_mean_pool(x, batch)
        x = self.fc1(x)
        x = self.fc2(x)

        return x
