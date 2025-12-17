
import random
from min_cutrank.matrix_tools import create_zero_matrix

class Graph:

    """
    A representation of a simple graph, defined by an adjacency matrix.
    """

    nmb_nodes : int
    """The number of nodes in the graph."""

    nodes : list[int]
    """List of all nodes, i.e. the range from 0 to nmb_nodes exclusive."""

    adjacencies : list[list[int]]
    """The adjacency matrix of the graph. Should be a square symmetric matrix with a row and column for each graph node, 
    0 on main diagonal, 1 in position (i,j) if (i,j) is an edge in the graph, 0 if not."""

    def __init__(self, adjacencies : list[list[int]]):
        self.adjacencies = adjacencies
        self.nmb_nodes = len(adjacencies)
        self.nodes = list(range(self.nmb_nodes))


    def grid_graph(rows : int, columns : int) -> 'Graph':

        adj_matrix = create_zero_matrix(rows * columns, rows * columns)
        for r in range(rows):
            for c in range(columns):
                pos = c + r * columns
                if r > 0:
                    set_edge(adj_matrix, pos, pos - columns)
                if c > 0:
                    set_edge(adj_matrix, pos, pos - 1)
        return Graph(adj_matrix)


    def random_graph(nodes : int, edge_probability : float) -> 'Graph':

        adj_mat = create_zero_matrix(nodes, nodes)
        for i in range(nodes - 1):
            for j in range(i + 1, nodes):
                if random.random() < edge_probability:
                    set_edge(adj_mat, i, j)
        return Graph(adj_mat)


def set_edge(adjacency_matrix : list[list[int]], n_from : int, n_to : int) -> None:

    if n_from != n_to:
        adjacency_matrix[n_from][n_to] = 1
        adjacency_matrix[n_to][n_from] = 1

