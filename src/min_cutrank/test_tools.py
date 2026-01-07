from abc import ABC, abstractmethod
import numpy as np
import time
import random
from min_cutrank.graph import Graph, set_edge
from min_cutrank.graph_partition import GraphPartition, SubMatrix
from min_cutrank.matrix_tools import create_zero_matrix, copy_matrix, rank_matrix_positions, set_common_matrix_value, insert_zero_matrix, add_matrix, add_product_matrix, is_zero_matrix, is_identity_matrix
from min_cutrank.swap_rank_calculator import all_swap_cut_ranks, one_to_many_swap_cut_ranks, single_swap_cut_rank_delta


def parse_int(value: str, default: int) -> int:

    try:
        result = int(value)
    except ValueError:
        result = default
    return result


def parse_float(value: str, default: float) -> float:

    try:
        result = float(value)
    except ValueError:
        result = default
    return result



def parse_bool(value: str, default: bool) -> bool:

    value_up = value.upper()

    if value_up in ("T", "TRUE", "1", "Y", "YES"):
        return True
    elif value_up in ("F", "FALSE", "0", "N", "NO"):
        return False
    else:
        return default


def graph_from_description(description : str) -> list[list[int]]:

    # 'gNxM' for grid with N rows and M colmns, like 'g5x6'
    # 'rN[eP]' for graph with N nodes and random edge probability of P (default 0.5), like 'r20' for graph of 20 nodes with edge probability 0.5, or 'r16P0.3' for graph of 16 nodes with edge probability 0.3

    gr_type = description[0]

    if gr_type == "g":
        idx_x = description.index("x")
        rows = int(description[1 : idx_x])
        cols = int(description[(idx_x + 1) :])
        return Graph.grid_graph(rows, cols)

    elif gr_type == "r":
        idx_e = description.find("e")
        if idx_e >= 0:
            nodes = int(description[1 : idx_e])
            edge_prob = float(description[(idx_e + 1) :])
        else:
            nodes = int(description[1 :])
            edge_prob = 0.5
        return Graph.random_graph(nodes, edge_prob)

    else:
        raise Exception(f"Unknown graph type : {gr_type}")


def temperatures_from_description(description : str) -> np.ndarray[float]:

    # Temperatures given at format 'BeEsS' for S samples beginning at temperature B and ending at temperature E. Example: '1.0e0.1s10' for 10 samples from 1.0 to 0.1
    idx_e = description.index("e")
    idx_s = description.index("s")
    start = float(description[: idx_e])
    end = float(description[(idx_e + 1) : idx_s])
    samples = int(description[(idx_s + 1) :])
    return np.linspace(start, end, samples)



class RankCollector(ABC):

    partition : GraphPartition

    @abstractmethod
    def collect_ranks(self, cut_ranks : list[list[int]], nodes_to_swap1: list[int], nodes_to_swap2: list[int]) -> None:
        pass
    @abstractmethod
    def name(self) -> str:
        return None

def swap(list, item1, item2):
    for i in range(len(list)):
        if list[i] == item1:
            list[i] = item2
        elif list[i] == item2:
            list[i] = item1


class DirectSwapRankCollector(RankCollector):

    buffer : list[list[int]]

    def __init__(self, partition : GraphPartition):
        self.partition = partition
        self.buffer = create_zero_matrix(partition.graph.nmb_nodes, partition.graph.nmb_nodes)

    def collect_ranks(self, cut_ranks : list[list[int]], nodes_to_swap1: list[int], nodes_to_swap2: list[int]) -> None:
        subsets_copy = [s[:] for s in self.partition.subsets]
        for node1 in nodes_to_swap1:
            for node2 in nodes_to_swap2:
                cut_ranks[node1][node2] = 0
                for i, rows in enumerate(subsets_copy):
                    for j in range(i):
                        cols = subsets_copy[j]
                        swap(rows, node1, node2)
                        swap(cols, node1, node2)
                        copy_matrix(self.partition.graph.adjacencies, self.buffer, rows, cols)
                        base_rows, _ = rank_matrix_positions(self.buffer, rows, cols)
                        cut_ranks[node1][node2] += len(base_rows)
                        swap(rows, node1, node2)
                        swap(cols, node1, node2)

    def name(self) -> str:
        return "Gauss-Jordan elimination rank calculation"


class FormulaRankCollector(RankCollector):

    single_ranks : bool

    row_ranks : bool

    def __init__(self, partition : GraphPartition, single_ranks : bool, row_ranks : bool):
        self.partition = partition
        self.single_ranks = single_ranks
        self.row_ranks = row_ranks and not single_ranks
        self.call_count = 0

    def collect_ranks(self, cut_ranks : list[list[int]], nodes_to_swap1: list[int], nodes_to_swap2: list[int]) -> None:
        
        old_rank = self.partition.cut_rank
        
        if self.single_ranks:
            for node1 in nodes_to_swap1:
                for node2 in nodes_to_swap2:
                    cut_ranks[node1][node2] = old_rank + single_swap_cut_rank_delta(self.partition, node1, node2)
        elif self.row_ranks:
            for node1 in nodes_to_swap1:
                one_to_many_swap_cut_ranks(self.partition, node1, nodes_to_swap2, cut_ranks[node1])
        else:
            all_swap_cut_ranks(self.partition, nodes_to_swap1, nodes_to_swap2, cut_ranks)

    def name(self) -> str:
        return "Single ranks by formulas" if self.single_ranks else ("Row ranks by formulas" if self.row_ranks else "All ranks by formulas")


def print_base_matrices(heading : str, partition : GraphPartition):
    print()
    print(heading)
    print_matrix("C^(-1):", partition.base_inverse)
    print_matrix("D = A^YB * C^(-1):", partition.adj_b_inverse)
    print_matrix("E = C^(-1) * A_XB:", partition.b_inverse_adj)
    print_matrix("F = A^YB * C^(-1) * A_XB + A:", partition.adj_b_inv_adj)

def print_matrix(heading : str, matrix : list[list[int]]):
    print(heading)
    for r in matrix:
        row = "   "
        for n in r:
            row += str(n)
        print(row)

class ApplySwapRankCollector(RankCollector):

    backup : GraphPartition

    validate : bool

    buffer_flag : list[bool]

    def __init__(self, partition : GraphPartition, validate : bool):
        self.partition = partition
        self.backup = partition.clone()

        self.validate = validate
        self.buffer_flag = [False] * partition.graph.nmb_nodes

    def collect_ranks(self, cut_ranks : list[list[int]], nodes_to_swap1: list[int], nodes_to_swap2: list[int]) -> None:
        self.backup.copy(self.partition)
        for node1 in nodes_to_swap1:
            for node2 in nodes_to_swap2:
                self.partition.apply_swap(node1, node2)
                cut_ranks[node1][node2] = self.partition.cut_rank
                if self.validate:
                    self._validate_partition()
                self.partition.copy(self.backup)

    def name(self) -> str:
        return "Apply swap with validation" if self.validate else "Apply swap without validation"

    def _validate_partition(self) -> None:

        p = self.partition

        for subset_index, subset in enumerate(p.subsets):
            for node in subset:
                if p.subset_index[node] != subset_index:
                    raise Exception("Subset index inconsistent with subsets")
        if sum(i != -1 for i in p.subset_index) != sum(len(s) for s in p.subsets):
            raise Exception("Subset index has invalid entries")

        for i, subset1 in enumerate(p.subsets):
            for j, subset2 in enumerate(p.subsets):
                matrix = p.matrices[i][j]
                if i < j:
                    if matrix.rows != subset1 or matrix.columns != subset2:
                        raise Exception("Submatrix rows or columns inconsistent with subsets")
                else:
                    if matrix is not None:
                        raise Exception("Unexpected submatrix")

        for list in p.matrices:
            for m in list:
                if m is not None:
                    self._validate_matrix(m)

    def _validate_matrix(self, m: SubMatrix) -> None:

        # Test matrix consistency

        p = self.partition
        nodes = p.graph.nodes
        row_flag = self.buffer_flag
        for n in m.rows:
            row_flag[n] = True
        for n in m.columns:
            if row_flag[n]:
                raise Exception("Node used both as row and column")

        for n in m.base_rows:
            if not m.base_flag[n] or not row_flag[n]:
                raise Exception("Unexpected element in base_rows")
            row_flag[n] = False
        for n in m.free_rows:
            if m.base_flag[n] or not row_flag[n]:
                raise Exception("Unexpected element in free_rows")
            row_flag[n] = False
        if any(row_flag):
            raise Exception("Row not classified as base or free")

        column_flag = self.buffer_flag
        for n in m.columns:
            column_flag[n] = True

        for n in m.base_columns:
            if not m.base_flag[n] or not column_flag[n]:
                raise Exception("Unexpected element in base_columns")
            column_flag[n] = False
        for n in m.free_columns:
            if m.base_flag[n] or not column_flag[n]:
                raise Exception("Unexpected element in free_columns")
            column_flag[n] = False
        if any(column_flag):
            raise Exception("Column not classified as base or free")
                
        if m.rank != len(m.base_rows):
            raise Exception("Number of base rows differs from rank")
        if m.rank != len(m.base_columns):
            raise Exception("Number of base columns differs from rank")
        for n in nodes:
            self.buffer_flag[n] = False

        adjacencies = p.graph.adjacencies
        buffer = p.buffer

        # Test C * C^(-1) = Id
        insert_zero_matrix(buffer, m.base_rows, m.base_rows)
        add_product_matrix(adjacencies, m.base_inverse, buffer, m.base_rows, m.base_columns, m.base_rows)
        if not is_identity_matrix(buffer, m.base_rows):
            raise Exception("Wrong inverse of rank matrix")

        # Test D-matrix in base set
        copy_matrix(m.adj_b_inverse, buffer, nodes, m.base_rows)
        add_product_matrix(adjacencies, m.base_inverse, buffer, nodes, m.base_columns, m.base_rows)
        if not is_zero_matrix(buffer, nodes, m.base_rows):
            raise Exception("Wrong value of A^(YB) * C^(-1)")
        if not is_identity_matrix(m.adj_b_inverse, m.base_rows):
            raise Exception("XB x XB submatrix of A^(YB) * C^(-1) is not identity")

        # Test E-matrix in base set
        copy_matrix(m.b_inverse_adj, buffer, m.base_columns, nodes)
        add_product_matrix(m.base_inverse, adjacencies, buffer, m.base_columns, m.base_rows, nodes)
        if not is_zero_matrix(buffer, nodes, m.base_rows):
            raise Exception("Wrong value of C^(-1) * A_(XB)")
        if not is_identity_matrix(m.b_inverse_adj, m.base_columns):
            raise Exception("YB x YB submatrix of C^(-1) * A_(XB) is not identity")

        # Test F-matrix in base set
        copy_matrix(m.adj_b_inv_adj, buffer, nodes, nodes)
        add_matrix(adjacencies, buffer, nodes, nodes)
        add_product_matrix(m.adj_b_inverse, adjacencies, buffer, nodes, m.base_rows, nodes)
        if not is_zero_matrix(buffer, nodes, nodes):
            raise Exception("Wrong value of A^(YB) * C^(-1) * A_(XB) + A")
        if not is_zero_matrix(m.adj_b_inv_adj, m.base_rows, nodes):
            raise Exception("XB rows of A^(YB) * C^(-1) * A_(XB) + A is not zero")
        if not is_zero_matrix(m.adj_b_inv_adj, nodes, m.base_columns):
            raise Exception("YB columns of A^(YB) * C^(-1) * A_(XB) + A is not zero")

        # Test that all other rows and columns in adjacency matrix between partition sets is generated by C
        copy_matrix(adjacencies, buffer, m.free_rows, m.free_columns)
        insert_zero_matrix(buffer, m.free_rows, m.base_rows)
        add_product_matrix(adjacencies, m.base_inverse, buffer, m.free_rows, m.base_columns, m.base_rows)
        add_product_matrix(buffer, adjacencies, buffer, m.free_rows, m.base_rows, m.free_columns)
        if not is_zero_matrix(buffer, m.free_rows, m.free_columns):
            raise Exception("Not a full rank matrix")

class CutRankCalculatorComparer:

    partition : GraphPartition
    """The partition for which to calculate cut rank deltas."""

    nodes1: list[int]
    """The first subset of nodes to calculate for."""

    nodes2: list[int]
    """The second subset of nodes to calculate for."""

    first_cut_ranks : list[list[int]]
    """The cut-ranks calculated in the first calculation pass."""

    second_cut_ranks : list[list[int]]
    """The cut-ranks calculated in the second calculation pass."""

    first_calculations_name : str

    def __init__(self, partition : GraphPartition):
        self.partition = partition
        self.first_cut_ranks = create_zero_matrix(partition.graph.nmb_nodes, partition.graph.nmb_nodes)
        self.second_cut_ranks = create_zero_matrix(partition.graph.nmb_nodes, partition.graph.nmb_nodes)

    def reset(self, nodes1: list[int], nodes2: list[int]) -> None:
        self.first_calculations_name = None
        self.nodes1 = nodes1
        self.nodes2 = nodes2

    def is_reset(self) -> bool:
        return self.first_calculations_name == None

    def calculate_and_compare(self, collector : RankCollector, log: bool = True) -> None:

        name = collector.name()
        is_first = self.is_reset()
        if is_first:
            set_common_matrix_value(-1, self.first_cut_ranks, self.partition.graph.nodes, self.partition.graph.nodes)
            self.first_calculations_name = name
            start = time.time()
            collector.collect_ranks(self.first_cut_ranks, self.nodes1, self.nodes2)
            end = time.time()
        else:
            set_common_matrix_value(-1, self.second_cut_ranks, self.partition.graph.nodes, self.partition.graph.nodes)
            start = time.time()
            collector.collect_ranks(self.second_cut_ranks, self.nodes1, self.nodes2)
            end = time.time()

        if log:
            print(f"Cut rank method '{name}' executed in {end - start} sec")

        if is_first:
            p_rank = self.partition.cut_rank
            max_rank = sum(min(len(m.rows), len(m.columns)) for m in self.partition.matrix_list)
            # The rank of the matrix containing the row and column can change by up to 2.
            # Each of the other 2*(n-2) matrices that contain either the row or the column can change rank by up to 1.
            max_change = 2 + 2 * (len(self.partition.subsets) - 2)
            for node1 in self.nodes1:
                for node2 in self.nodes2:
                    rank = self.first_cut_ranks[node1][node2]
                    if rank < 0 or rank > max_rank:
                        print(f"Cut-rank for position ({node1},{node2}) is {rank}, outside allowed range of [0,{max_rank}]")
                        raise Exception("Cut-rank outside allowed range")
                    if rank < p_rank - max_change or rank > p_rank + max_change:
                        print(f"Cut-rank for position ({node1},{node2}) is {rank}, too far from current cut-rank {p_rank}")
                        raise Exception("New cut-rank too far from current cut-rank")
                for node2 in self.nodes1:
                    if self.first_cut_ranks[node1][node2] != -1:
                        print(f"Cut-rank for position ({node1},{node2}) is {rank}, should be -1 since {node2} is not a column position")
                        raise Exception("Cut-rank set outside Rows x Columns")
            for node1 in self.nodes2:
                for node2 in self.nodes2:
                    if self.first_cut_ranks[node1][node2] != -1:
                        print(f"Cut-rank for position ({node1},{node2}) is {rank}, should be -1 since {node1} is not a row position")
                        raise Exception("Cut-rank set outside Rows x Columns")
                for node2 in self.nodes1:
                    if self.first_cut_ranks[node1][node2] != -1:
                        print(f"Cut-rank for position ({node1},{node2}) is {rank}, should be -1 since {node1} is not a row position and {node2} is not a column position")
                        raise Exception("Cut-rank set outside Rows x Columns")
        else:
            for i in self.partition.graph.nodes:
                for j in self.partition.graph.nodes:
                    if self.first_cut_ranks[i][j] != self.second_cut_ranks[i][j]:
                        print(f"Cut-rank mismatch for position ({i},{j}):")
                        print(f"   {self.first_calculations_name}: {self.first_cut_ranks[i][j]}")
                        print(f"   {name}: {self.second_cut_ranks[i][j]}")
                        raise Exception("Cut-rank mismatch")


def triangle_example() -> GraphPartition:
    matr = create_zero_matrix(6, 6)
    set_edge(matr, 0, 1)
    set_edge(matr, 0, 2)
    set_edge(matr, 1, 2)
    set_edge(matr, 3, 4)
    set_edge(matr, 3, 5)
    set_edge(matr, 4, 5)
    set_edge(matr, 0, 3)
    partition_flags = [True, True, False, True, False, False]
    return GraphPartition.fromFlags(matr, partition_flags)



def rank_collector_from_name(method_name : str, partition : GraphPartition) -> RankCollector:

    if method_name == "gauss":
        return DirectSwapRankCollector(partition)
    elif method_name == "single":
        return FormulaRankCollector(partition, True, False)
    elif method_name == "row":
        return FormulaRankCollector(partition, False, True)
    elif method_name == "all":
        return FormulaRankCollector(partition, False, False)
    elif method_name == "apply":
        return ApplySwapRankCollector(partition, False)
    elif method_name == "validate":
        return ApplySwapRankCollector(partition, True)
    else:
        raise Exception(f"Unknown cut-rank calculation method: '{method_name}'")


def run_greedy_min_rank(partition : GraphPartition, rank_calculation_methods: list[str]) -> None:

    rank_collectors = [rank_collector_from_name(method_name, partition) for method_name in rank_calculation_methods]
    rank_comparer = CutRankCalculatorComparer(partition)
    local_minimum_found = False
    cut_rank = partition.cut_rank
    iteration = 0


    while not local_minimum_found:

        iteration += 1
        print()
        print(f"Iteration {iteration} starting, current rank = {cut_rank}")

        rows, columns = partition.rows_and_columns_copy()
        rank_comparer.reset(rows, columns)
        for coll in rank_collectors:
            rank_comparer.calculate_and_compare(coll)
        if rank_comparer.is_reset():
            raise Exception("No ranks have been calculated")

        ranks_grouped = [[] for _ in range(5)]
        for row in rows:
            for col in columns:
                new_rank_pos = rank_comparer.first_cut_ranks[row][col] + 2 - cut_rank
                ranks_grouped[new_rank_pos].append((row, col))
        best_rank_pos = -1
        for i, ranks_at_level in enumerate(ranks_grouped):
            if len(ranks_at_level) > 0:
                print(f"Cut-rank = {i + cut_rank - 2} for {len(ranks_at_level)} swaps")
                if best_rank_pos == -1:
                    best_rank_pos = i
        local_minimum_found = best_rank_pos >= 2
        if not local_minimum_found:
            sw_pos = random.randrange(len(ranks_grouped[best_rank_pos]))
            sw_i, sw_j = ranks_grouped[best_rank_pos][sw_pos]
            print(f"Applying swap ({sw_i}, {sw_j})")
            partition.apply_swap(sw_i, sw_j)
            cut_rank += best_rank_pos - 2
            if cut_rank != partition.cut_rank:
                raise Exception("Unexpected rank after swap")

    print()
    print(f"Stopped after {iteration} iterations")
    print(f"Final rank is {cut_rank}")
    print("Final partition:")
    for i, subset in enumerate(partition.subsets):
        print(f"Set {i+1} = {subset}")
