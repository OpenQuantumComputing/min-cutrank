from min_cutrank.graph import Graph
from min_cutrank.matrix_tools import create_zero_matrix
from min_cutrank.sub_matrix import SubMatrix

class GraphPartition:

    """
    A partition of the nodes of a graph into a collection of two or more subsets
    """

    graph : Graph
    """The graph being partitioned."""
    
    cut_rank : int
    """The cut-rank of the partition."""

    subsets : list[list[int]]
    """The subsets of nodes in the partition."""
    
    subset_index: list[int]
    """For each node in the graph, the index of the subset it belongs to, or -1 if it is not in any subset."""

    matrices: list[list[SubMatrix]]
    """matrices[i][j] where i < j is the SubMatrix whose rows are subsets[i] and whose columns are subsets[j]."""

    matrix_list: list[SubMatrix]
    """A flat list of all SubMatrix objects in this partition"""

    buffer : list[list[int]]
    """A square nmb_nodes x nmb_nodes used for caching intermediate calculations when updating the variables after the partition has been changed."""

    def __init__(self, graph: Graph, subsets : list[list[int]]):
        if subsets is None or len(subsets) < 2:
            raise Exception("GraphPartition must have at least two subsets")
        if len(set(node for subset in subsets for node in subset)) < sum(len(subset) for subset in subsets):
            raise Exception("Subsets in GraphPartition must be disjoint")

        self.graph = graph
        self.subsets = [subset[:] for subset in subsets]

        self.subset_index = [-1] * graph.nmb_nodes
        for i, subset in enumerate(self.subsets):
            for node in subset:
                self.subset_index[node] = i

        self.matrices = [[SubMatrix(self, s, t) if i < j else None for j, t in enumerate(self.subsets)] for i, s in enumerate(self.subsets)]
        self.matrix_list = [m for i, line in enumerate(self.matrices) for m in line[i+1:]]
        self.cut_rank = sum(m.rank for m in self.matrix_list)
        self.buffer = self._empty_matrix()

    def clone(self) -> 'GraphPartition':
        """Creates a deep copy of this GraphPartition."""
        return GraphPartition(self.graph, self.subsets)

    def copy(self, other: 'GraphPartition') -> None:
        """Copy all partition state from other into self."""
        if self.graph != other.graph:
            raise Exception("Cannot copy partition from different graph")
        
        self.cut_rank = other.cut_rank
        self.subset_index[:] = other.subset_index[:]

        for j in range(len(self.subsets)):
            self.subsets[j][:] = other.subsets[j][:]
            for i in range(0, j):
                self.matrices[i][j].copy(other.matrices[i][j])

    def fromFlags(graph: Graph, subset1_flags : list[bool], subset2_flags : list[bool] = None) -> 'GraphPartition':
        """Creates a GraphPartition with two subsets, from flags for which nodes belong to the subsets.
        If the subset2_flags is not given, it is assumed to be the negation of the subset1_flags.
        """
        subset2_flags = subset2_flags if subset2_flags is not None else [not flag for flag in subset1_flags]

        subset1 = [n for n in graph.nodes if subset1_flags[n]]
        subset2 = [n for n in graph.nodes if subset2_flags[n]]

        return GraphPartition(graph, [subset1, subset2])

    def _empty_matrix(self) -> list[list[int]]:

        return create_zero_matrix(self.graph.nmb_nodes, self.graph.nmb_nodes)

    def rows_and_columns_copy(self) -> tuple[list[int], list[int]]:
        """Returns a copy of the rows and columns of this partition.
        This makes sense only for bipartitions, where the rows are the first subset and the columns the second subset.
        """
        if self.subsets is None or len(self.subsets) != 2:
            raise Exception("Partition does not have exactly two subsets")
        return self.subsets[0][:], self.subsets[1][:]

    def apply_swap(self, node1: int, node2: int) -> None:
        """Applies the swap of the given nodes to this partition.
        
        args:
            - node1: 'int' The first node to be swapped.
            - node2: 'int' The second node to be swapped.
        """
        subset1_index = self.subset_index[node1]
        subset2_index = self.subset_index[node2]

        for i in range(len(self.subsets)):
            if subset1_index != -1:
                if i < subset1_index:
                    self.cut_rank += self.matrices[i][subset1_index].apply_swap(node2, node1, i == subset2_index, True)
                elif subset1_index < i:
                    self.cut_rank += self.matrices[subset1_index][i].apply_swap(node1, node2, True, i == subset2_index)
            if i != subset1_index and subset2_index != -1:
                if i < subset2_index:
                    self.cut_rank += self.matrices[i][subset2_index].apply_swap(node1, node2, False, True)
                elif subset2_index < i:
                    self.cut_rank += self.matrices[subset2_index][i].apply_swap(node2, node1, True, False)

        if subset1_index != -1:
            subset1 = self.subsets[subset1_index]
            node1_idx = subset1.index(node1)
            subset1[node1_idx] = node2
        self.subset_index[node2] = subset1_index

        if subset2_index != -1:
            subset2 = self.subsets[subset2_index]
            node2_idx = subset2.index(node2)
            subset2[node2_idx] = node1
        self.subset_index[node1] = subset2_index
        