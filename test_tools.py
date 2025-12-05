import numpy as np
import time
import random
from graph_partition import GraphPartition
from partition_builder import set_edge, grid_graph, random_graph
from matrix_tools import create_zero_matrix, copy_matrix, rank_matrix_positions, set_common_matrix_value, insert_zero_matrix, add_matrix, add_product_matrix, is_zero_matrix, is_identity_matrix
from swap_rank_calculator import all_swap_cut_ranks, row_swap_cut_ranks, single_swap_cut_rank


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
        return grid_graph(rows, cols)

    elif gr_type == "r":
        idx_e = description.find("e")
        if idx_e >= 0:
            nodes = int(description[1 : idx_e])
            edge_prob = float(description[(idx_e + 1) :])
        else:
            nodes = int(description[1 :])
            edge_prob = 0.5
        return random_graph(nodes, edge_prob)

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


def clone_partition(partition : GraphPartition) -> GraphPartition:
    return GraphPartition(partition.adjacencies, partition.row_flag)


class RankCollector:

    def collect_ranks(self, cut_ranks : list[list[int]]) -> None:
        pass
    def name(self) -> str:
        return None


class DirectSwapRankCollector(RankCollector):

    partition : GraphPartition

    buffer : list[list[int]]

    def __init__(self, partition : GraphPartition):
        self.partition = partition
        self.buffer = create_zero_matrix(partition.nmb_nodes, partition.nmb_nodes)

    def collect_ranks(self, cut_ranks : list[list[int]]) -> None:
        rows_copy = self.partition.rows[:]
        cols_copy = self.partition.columns[:]
        for i in range(len(self.partition.rows)):
            row = self.partition.rows[i]
            for j in range(len(self.partition.columns)):
                col = self.partition.columns[j]
                rows_copy[i], cols_copy[j] = col, row
                copy_matrix(self.partition.adjacencies, self.buffer, rows_copy, cols_copy)
                base_rows, _ = rank_matrix_positions(self.buffer, rows_copy, cols_copy)
                cut_ranks[row][col] = len(base_rows)
                rows_copy[i], cols_copy[j] = row, col

    def name(self) -> str:
        return "Gauss-Jordan elimination rank calculation"


class FormulaRankCollector(RankCollector):

    partition : GraphPartition

    single_ranks : bool

    row_ranks : bool

    def __init__(self, partition : GraphPartition, single_ranks : bool, row_ranks : bool):
        self.partition = partition
        self.single_ranks = single_ranks
        self.row_ranks = row_ranks and not single_ranks
        row_ranks = [0] * self.partition.nmb_nodes

    def collect_ranks(self, cut_ranks : list[list[int]]) -> None:
        if self.single_ranks:
            for row in self.partition.rows:
                for col in self.partition.columns:
                    cut_ranks[row][col] = single_swap_cut_rank(self.partition, row, col)
        elif self.row_ranks:
            for row in self.partition.rows:
                row_swap_cut_ranks(self.partition, row, cut_ranks[row])
        else:
            all_swap_cut_ranks(self.partition, cut_ranks)

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

    partition : GraphPartition

    validate : bool

    row_flag : list[bool]

    rows : list[int]

    columns : list[int]

    base_flag : list[bool]

    cut_rank : int

    base_rows : list[int]

    base_columns : list[int]

    free_rows : list[int]

    free_columns : list[int]

    base_inverse : list[list[int]]

    adj_b_inverse: list[list[int]]

    b_inverse_adj: list[list[int]]

    adj_b_inv_adj: list[list[int]]

    buffer_flag : list[bool]

    def __init__(self, partition : GraphPartition, validate : bool):
        self.partition = partition
        self.validate = validate
        self.base_inverse = create_zero_matrix(partition.nmb_nodes, partition.nmb_nodes)
        self.adj_b_inverse = create_zero_matrix(partition.nmb_nodes, partition.nmb_nodes)
        self.b_inverse_adj = create_zero_matrix(partition.nmb_nodes, partition.nmb_nodes)
        self.adj_b_inv_adj = create_zero_matrix(partition.nmb_nodes, partition.nmb_nodes)
        self.row_flag = [False] * partition.nmb_nodes
        self.base_flag = [False] * partition.nmb_nodes
        self.rows = [0] * len(partition.rows)
        self.columns = [0] * len(partition.columns)
        self.buffer_flag = [False] * partition.nmb_nodes

    def collect_ranks(self, cut_ranks : list[list[int]]) -> None:
        self._backup()
        for row in self.rows:
            for col in self.columns:
                self.partition.apply_swap(row, col)
                cut_ranks[row][col] = self.partition.cut_rank
                if self.validate:
                    self._validate_partition()
                self._restore()

    def name(self) -> str:
        return "Apply swap with validation" if self.validate else "Apply swap without validation"

    def _validate_partition(self) -> None:

        # Test partition sets: Every node occurs exactly once in exactly one of the four partition sets
        # Nodes have row and base flag set according to the set they are in
        # Equally many base rows and base columns
        p = self.partition
        buffer = p.buffer
        for n in p.base_rows:
            if not p.base_flag[n] or not p.row_flag[n]:
                raise Exception("Unexpected element in base_rows")
            if self.buffer_flag[n]:
                raise Exception("Node used several times in partition sets")
            self.buffer_flag[n] = True
        for n in p.base_columns:
            if not p.base_flag[n] or p.row_flag[n]:
                raise Exception("Unexpected element in base_columns")
            if self.buffer_flag[n]:
                raise Exception("Node used several times in partition sets")
            self.buffer_flag[n] = True
        for n in p.free_rows:
            if p.base_flag[n] or not p.row_flag[n]:
                raise Exception("Unexpected element in free_rows")
            if self.buffer_flag[n]:
                raise Exception("Node used several times in partition sets")
            self.buffer_flag[n] = True
        for n in p.free_columns:
            if p.base_flag[n] or p.row_flag[n]:
                raise Exception("Unexpected element in free_columns")
            if self.buffer_flag[n]:
                raise Exception("Node used several times in partition sets")
            self.buffer_flag[n] = True
        if not all(self.buffer_flag):
            raise Exception("Not all nodes found in partition sets")
        if p.cut_rank != len(p.base_rows):
            raise Exception("Number of base rows differs from cut-rank")
        if p.cut_rank != len(p.base_columns):
            raise Exception("Number of base columns differs from cut-rank")
        for n in p.nodes:
            self.buffer_flag[n] = False

        # Test C * C^(-1) = Id
        insert_zero_matrix(buffer, p.base_rows, p.base_rows)
        add_product_matrix(p.adjacencies, p.base_inverse, buffer, p.base_rows, p.base_columns, p.base_rows)
        if not is_identity_matrix(buffer, p.base_rows):
            raise Exception("Wrong inverse of rank matrix")

        # Test D-matrix in base set
        copy_matrix(p.adj_b_inverse, buffer, p.nodes, p.base_rows)
        add_product_matrix(p.adjacencies, p.base_inverse, buffer, p.nodes, p.base_columns, p.base_rows)
        if not is_zero_matrix(buffer, p.nodes, p.base_rows):
            raise Exception("Wrong value of A^(YB) * C^(-1)")
        if not is_identity_matrix(p.adj_b_inverse, p.base_rows):
            raise Exception("XB x XB submatrix of A^(YB) * C^(-1) is not identity")

        # Test E-matrix in base set
        copy_matrix(p.b_inverse_adj, buffer, p.base_columns, p.nodes)
        add_product_matrix(p.base_inverse, p.adjacencies, buffer, p.base_columns, p.base_rows, p.nodes)
        if not is_zero_matrix(buffer, p.nodes, p.base_rows):
            raise Exception("Wrong value of C^(-1) * A_(XB)")
        if not is_identity_matrix(p.b_inverse_adj, p.base_columns):
            raise Exception("YB x YB submatrix of C^(-1) * A_(XB) is not identity")

        # Test F-matrix in base set
        copy_matrix(p.adj_b_inv_adj, buffer, p.nodes, p.nodes)
        add_matrix(p.adjacencies, buffer, p.nodes, p.nodes)
        add_product_matrix(p.adj_b_inverse, p.adjacencies, buffer, p.nodes, p.base_rows, p.nodes)
        if not is_zero_matrix(buffer, p.nodes, p.nodes):
            raise Exception("Wrong value of A^(YB) * C^(-1) * A_(XB) + A")
        if not is_zero_matrix(p.adj_b_inv_adj, p.base_rows, p.nodes):
            raise Exception("XB rows of A^(YB) * C^(-1) * A_(XB) + A is not zero")
        if not is_zero_matrix(p.adj_b_inv_adj, p.nodes, p.base_columns):
            raise Exception("YB columns of A^(YB) * C^(-1) * A_(XB) + A is not zero")

        # Test that all other rows and columns in adjacency matrix between partition sets is generated by C
        copy_matrix(p.adjacencies, buffer, p.free_rows, p.free_columns)
        insert_zero_matrix(buffer, p.free_rows, p.base_rows)
        add_product_matrix(p.adjacencies, p.base_inverse, buffer, p.free_rows, p.base_columns, p.base_rows)
        add_product_matrix(buffer, p.adjacencies, buffer, p.free_rows, p.base_rows, p.free_columns)
        if not is_zero_matrix(buffer, p.free_rows, p.free_columns):
            raise Exception("Not a full rank matrix")

    def _copy_list(self, from_l : list, to_l : list) -> None:
        for i in range(len(from_l)):
            to_l[i] = from_l[i]

    def _backup(self) -> None:
        copy_matrix(self.partition.base_inverse, self.base_inverse, self.partition.nodes, self.partition.nodes)
        copy_matrix(self.partition.adj_b_inverse, self.adj_b_inverse, self.partition.nodes, self.partition.nodes)
        copy_matrix(self.partition.b_inverse_adj, self.b_inverse_adj, self.partition.nodes, self.partition.nodes)
        copy_matrix(self.partition.adj_b_inv_adj, self.adj_b_inv_adj, self.partition.nodes, self.partition.nodes)
        self._copy_list(self.partition.row_flag, self.row_flag)
        self._copy_list(self.partition.base_flag, self.base_flag)
        self._copy_list(self.partition.rows, self.rows)
        self._copy_list(self.partition.columns, self.columns)
        self.cut_rank = self.partition.cut_rank
        self.base_rows = self.partition.base_rows
        self.base_columns = self.partition.base_columns
        self.free_rows = self.partition.free_rows
        self.free_columns = self.partition.free_columns

    def _restore(self) -> None:
        copy_matrix(self.base_inverse, self.partition.base_inverse, self.partition.nodes, self.partition.nodes)
        copy_matrix(self.adj_b_inverse, self.partition.adj_b_inverse, self.partition.nodes, self.partition.nodes)
        copy_matrix(self.b_inverse_adj, self.partition.b_inverse_adj, self.partition.nodes, self.partition.nodes)
        copy_matrix(self.adj_b_inv_adj, self.partition.adj_b_inv_adj, self.partition.nodes, self.partition.nodes)
        self._copy_list(self.row_flag, self.partition.row_flag)
        self._copy_list(self.base_flag, self.partition.base_flag)
        self._copy_list(self.rows, self.partition.rows)
        self._copy_list(self.columns, self.partition.columns)
        self.partition.cut_rank = self.cut_rank
        self.partition.base_rows = self.base_rows
        self.partition.base_columns = self.base_columns
        self.partition.free_rows = self.free_rows
        self.partition.free_columns = self.free_columns

class CutRankCalculatorComparer:

    partition : GraphPartition

    first_cut_ranks : list[list[int]]

    second_cut_ranks : list[list[int]]

    first_calculations_name : str

    def __init__(self, partition : GraphPartition):
        self.partition = partition
        self.first_cut_ranks = create_zero_matrix(partition.nmb_nodes, partition.nmb_nodes)
        self.second_cut_ranks = create_zero_matrix(partition.nmb_nodes, partition.nmb_nodes)
        self.reset()

    def reset(self) -> None:
        self.first_calculations_name = None

    def is_reset(self) -> bool:
        return self.first_calculations_name == None

    def calculate_and_compare(self, collector : RankCollector) -> None:

        name = collector.name()
        is_first = self.is_reset()
        if is_first:
            set_common_matrix_value(-1, self.first_cut_ranks, self.partition.nodes, self.partition.nodes)
            self.first_calculations_name = name
            start = time.time()
            collector.collect_ranks(self.first_cut_ranks)
            end = time.time()
        else:
            set_common_matrix_value(-1, self.second_cut_ranks, self.partition.nodes, self.partition.nodes)
            start = time.time()
            collector.collect_ranks(self.second_cut_ranks)
            end = time.time()

        print(f"Cut rank method '{name}' executed in {end - start} sec")

        if is_first:
            p_rank = self.partition.cut_rank
            max_rank = min(len(self.partition.rows), len(self.partition.columns))
            for row in self.partition.rows:
                for col in self.partition.columns:
                    rank = self.first_cut_ranks[row][col]
                    if rank < 0 or rank > max_rank:
                        print(f"Cut-rank for position ({row},{col}) is {rank}, outside allowed range of [0,{max_rank}]")
                        raise Exception("Cut-rank outside allowed range")
                    if rank < p_rank - 2 or rank > p_rank + 2:
                        print(f"Cut-rank for position ({row},{col}) is {rank}, too far from current cut-rank {p_rank}")
                        raise Exception("New cut-rank too far from current cut-rank")
                for col in self.partition.rows:
                    if self.first_cut_ranks[row][col] != -1:
                        print(f"Cut-rank for position ({row},{col}) is {rank}, should be -1 since {col} is not a column position")
                        raise Exception("Cut-rank set outside Rows x Columns")
            for row in self.partition.columns:
                for col in self.partition.columns:
                    if self.first_cut_ranks[row][col] != -1:
                        print(f"Cut-rank for position ({row},{col}) is {rank}, should be -1 since {row} is not a row position")
                        raise Exception("Cut-rank set outside Rows x Columns")
                for col in self.partition.rows:
                    if self.first_cut_ranks[row][col] != -1:
                        print(f"Cut-rank for position ({row},{col}) is {rank}, should be -1 since {row} is not a row position and {col} is not a column position")
                        raise Exception("Cut-rank set outside Rows x Columns")
        else:
            for i in self.partition.nodes:
                for j in self.partition.nodes:
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
    return GraphPartition(matr, partition_flags)



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

        rank_comparer.reset()
        for coll in rank_collectors:
            rank_comparer.calculate_and_compare(coll)
        if rank_comparer.is_reset():
            raise Exception("No ranks have been calculated")

        ranks_grouped = [[] for _ in range(5)]
        for row in partition.rows:
            for col in partition.columns:
                new_rank_pos = rank_comparer.first_cut_ranks[row][col] + 2 - cut_rank
                ranks_grouped[new_rank_pos].append((row, col))
        best_rank_pos = -1
        for i in range(5):
            if len(ranks_grouped[i]) > 0:
                print(f"Cut-rank = {i + cut_rank - 2} for {len(ranks_grouped[i])} swaps")
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
    print(f"Set 1 = {partition.rows}")
    print(f"Set 2 = {partition.columns}")
