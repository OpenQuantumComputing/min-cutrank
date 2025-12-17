import random
from min_cutrank.matrix_tools import create_zero_matrix
from min_cutrank.graph_partition import GraphPartition


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


def random_partition(adjacency_matrix : list[list[int]], portion1 : float, portion2 : float = None) -> GraphPartition:

    portion1and2 = 1 if portion2 == None else portion1 + portion2

    nmb_nodes = len(adjacency_matrix)
    nodes = [n for n in range(nmb_nodes)]
    random.shuffle(nodes)
    
    nmb_part1 = round(nmb_nodes * portion1)
    nmb_part1and2 = round(nmb_nodes * portion1and2)
    part1 = nodes[0:nmb_part1]
    part2 = nodes[nmb_part1:nmb_part1and2]
    return GraphPartition(adjacency_matrix, part1, part2)


def random_partition_on_random_graph(nodes : int, edge_probability : float, portion : float) -> GraphPartition:

    return random_partition(random_graph(nodes, edge_probability), portion)
