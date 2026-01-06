import random
from min_cutrank.graph import Graph
from min_cutrank.matrix_tools import create_zero_matrix
from min_cutrank.graph_partition import GraphPartition


def random_partition(graph: Graph, portion1 : float, portion2 : float = None) -> GraphPartition:

    portion1and2 = 1 if portion2 == None else portion1 + portion2

    nmb_nodes = graph.nmb_nodes
    nodes = [n for n in range(nmb_nodes)]
    random.shuffle(nodes)
    
    nmb_part1 = round(nmb_nodes * portion1)
    nmb_part1and2 = round(nmb_nodes * portion1and2)
    part1 = nodes[0:nmb_part1]
    part2 = nodes[nmb_part1:nmb_part1and2]
    return GraphPartition(graph, [part1, part2])


def random_partition_on_random_graph(nodes : int, edge_probability : float, portion : float) -> GraphPartition:

    return random_partition(Graph.random_graph(nodes, edge_probability), portion)
