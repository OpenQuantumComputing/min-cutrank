import unittest
import sys
import random
from min_cutrank.partition_builder import random_partition
from min_cutrank.graph_partition import GraphPartition
from min_cutrank.test_tools import CutRankCalculatorComparer, graph_from_description, rank_collector_from_name

class TestCutRankCalculation(unittest.TestCase):
    """Tests for verifying that different cut-rank calculation methods give the same results."""

    def test_rank_calculation_for_swap_matches_gauss(self):

        for method in ["single", "row", "all", "validate"]:
            for partial in [False, True]:
                with self.subTest(method=method, partial=partial):
            
                    graph_partition = self.createPartition(partial_partition=partial)
                    rank_calculation_methods = ["gauss", method]
                    
                    iterations = 100 if method != "validate" else 10
                    verify_random_swaps(graph_partition, rank_calculation_methods, iterations)
    
    
    def test_rank_calculation_for_replace_matches_gauss(self):

        for method in ["single", "row", "all", "validate"]:
            for target in ["row", "column"]:
                with self.subTest(method=method, target=target):

#                    random.seed(12345)
                    graph_partition = self.createPartition(partial_partition=True)
                    rank_calculation_methods = ["gauss", method]
                    
                    iterations = 100 if method != "validate" else 10
                    verify_random_replaces(graph_partition, rank_calculation_methods, iterations, target)
    
    
    def createPartition(self, size: int = 20, partial_partition: bool = False) -> GraphPartition:
        graph_setup = f"r{size}"
        graph_adj_matrix = graph_from_description(graph_setup)
        graph_partition = random_partition(graph_adj_matrix, 0.3, 0.3) if partial_partition \
            else random_partition(graph_adj_matrix, 0.5)
        return graph_partition



def verify_random_swaps(partition : GraphPartition, rank_calculation_methods: list[str], iterations: int) -> None:
    """Runs random swaps on the partition and verifies that the cut-rank calculated by different methods match."""

    rank_collectors = [rank_collector_from_name(method_name, partition) for method_name in rank_calculation_methods]
    rank_comparer = CutRankCalculatorComparer(partition)
    cut_rank = partition.cut_rank
    report(partition, "Initial")

    for _ in range(iterations):

        rows, cols = partition.rows_and_columns_copy()
        rank_comparer.reset(rows, cols)
        for coll in rank_collectors:
            rank_comparer.calculate_and_compare(coll, log=False)
        if rank_comparer.is_reset():
            raise Exception("No ranks have been calculated")

        row = random.choice(rows)
        col = random.choice(cols)

        print(f"Swapping row {row} and column {col}")
        partition.apply_swap(row, col)
        cut_rank = rank_comparer.first_cut_ranks[row][col]
        print(f"New cut-rank is {cut_rank}")
        if cut_rank != partition.cut_rank:
            raise Exception("Unexpected rank after swap")

    report(partition, "Final")


def verify_random_replaces(partition : GraphPartition, rank_calculation_methods: list[str], iterations: int, target : str) -> None:
    """Runs random replaces in the partition's rows or columns and verifies that the cut-rank calculated by different methods match."""

    rank_collectors = [rank_collector_from_name(method_name, partition) for method_name in rank_calculation_methods]
    rank_comparer = CutRankCalculatorComparer(partition)
    cut_rank = partition.cut_rank
    report(partition, "Initial")

    for _ in range(iterations):

        rows, columns = partition.rows_and_columns_copy()
        not_in_partition = list(set(partition.graph.nodes) - set(rows) - set(columns))

        if target == "row":
            columns = not_in_partition
        elif target == "column":
            rows = not_in_partition
        else:
            raise Exception(f"Unknown replace target: '{target}'")
        
        rank_comparer.reset(rows, columns)
        for coll in rank_collectors:
            rank_comparer.calculate_and_compare(coll, log=False)
        if rank_comparer.is_reset():
            raise Exception("No ranks have been calculated")

        row = random.choice(rows)
        column = random.choice(columns)
        if target == "row":
            print(f"Replacing row {row} with unused node {column}")
        elif target == "column":
            print(f"Replacing column {column} with unused node {row}")

        partition.apply_swap(row, column)
        cut_rank = rank_comparer.first_cut_ranks[row][column]
        print(f"New cut-rank is {cut_rank}")
        if cut_rank != partition.cut_rank:
            raise Exception("Unexpected rank after swap")

    report(partition, "Final")


def report(partition : GraphPartition, title: str) -> None:
    print()
    print(f"{title} cut-rank is {partition.cut_rank}")
    print(f"{title} partition:")
    for i, subset in enumerate(partition.subsets):
        print(f"Set {i+1} = {subset}")


if __name__ == "__main__":
    unittest.main()
