import random
from min_cutrank.graph import Graph
from min_cutrank.matrix_tools import create_zero_matrix
from min_cutrank.graph_partition import GraphPartition


def random_bipartition(graph: Graph, portion1 : float, portion2 : float = None) -> GraphPartition:

    portions = [portion1, 1 - portion1] if portion2 == None \
        else [portion1, portion2]
    return random_partition(graph, portions)

def random_partition(graph: Graph, portions : list[float]) -> GraphPartition:

    nmb_nodes = graph.nmb_nodes
    nodes = [n for n in range(nmb_nodes)]
    random.shuffle(nodes)
    
    parts = []
    taken_portion = 0.0
    
    for portion in portions:
        taken = round(nmb_nodes * taken_portion)
        next_taken = round(nmb_nodes * (taken_portion + portion))
        parts.append(nodes[taken:next_taken])
        taken_portion += portion

    if taken_portion > 1.0:
        raise Exception("Sum of portions exceeds 1.0")

    return GraphPartition(graph, parts)


def random_partition_on_random_graph(nodes : int, edge_probability : float, portion : float) -> GraphPartition:

    return random_bipartition(Graph.random_graph(nodes, edge_probability), portion)
