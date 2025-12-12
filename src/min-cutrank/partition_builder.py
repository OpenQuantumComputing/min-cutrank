import random
from matrix_tools import create_zero_matrix
from graph_partition import GraphPartition


def set_edge(adjacency_matrix : list[list[int]], n_from : int, n_to : int) -> None:

    if n_from != n_to:
        adjacency_matrix[n_from][n_to] = 1
        adjacency_matrix[n_to][n_from] = 1


def grid_graph(rows : int, columns : int) -> list[list[int]]:

    adj_matrix = create_zero_matrix(rows * columns, rows * columns)
    for r in range(rows):
        for c in range(columns):
            pos = c + r * columns
            if r > 0:
                set_edge(adj_matrix, pos, pos - columns)
            if c > 0:
                set_edge(adj_matrix, pos, pos - 1)
    return adj_matrix


def random_graph(nodes : int, edge_probability : float) -> list[list[int]]:

    adj_mat = create_zero_matrix(nodes, nodes)
    for i in range(nodes - 1):
        for j in range(i + 1, nodes):
            if random.random() < edge_probability:
                set_edge(adj_mat, i, j)
    return adj_mat


def random_partition(adjacency_matrix : list[list[int]], portion : float) -> GraphPartition:

    nmb_nodes = len(adjacency_matrix)
    nmb_part1 = round(nmb_nodes * portion)
    partition_flags = [True] * nmb_part1 + [False] * (nmb_nodes - nmb_part1)
    random.shuffle(partition_flags)
    return GraphPartition(adjacency_matrix, partition_flags)


def random_partition_on_random_graph(nodes : int, edge_probability : float, portion : float) -> GraphPartition:

    return random_partition(random_graph(nodes, edge_probability), portion)
