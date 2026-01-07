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

                    graph_partition = self.createPartition(partial_partition=True)
                    rank_calculation_methods = ["gauss", method]
                    
                    iterations = 100 if method != "validate" else 10
                    verify_random_replaces(graph_partition, rank_calculation_methods, iterations, target)
    
    
    def test_rank_calculation_for_3_subsets_matches_gauss(self):

        for method in ["single", "row", "all", "validate"]:
            for partial in [False, True]:
                with self.subTest(method=method, partial=partial):

                    graph_partition = self.createPartition(subsets = 3, partial_partition = partial)
                    rank_calculation_methods = ["gauss", method]
                    
                    iterations = 100 if method != "validate" else 10
                    verify_random_swaps(graph_partition, rank_calculation_methods, iterations)
    
    
    def createPartition(self, size: int = 20, subsets: int = 2, partial_partition: bool = False) -> GraphPartition:
        graph_setup = f"r{size}"
        graph_adj_matrix = graph_from_description(graph_setup)
        
        if partial_partition:
            subset_sizes = [1.0 / (subsets + 1)] * subsets
        else:
            subset_sizes = [1.0 / subsets] * (subsets - 1)
            subset_sizes.append(1.0 - sum(subset_sizes)) 

        graph_partition = random_partition(graph_adj_matrix, subset_sizes)
        return graph_partition



def verify_random_swaps(partition : GraphPartition, rank_calculation_methods: list[str], iterations: int) -> None:
    """Runs random swaps on the partition and verifies that the cut-rank calculated by different methods match."""

    rank_collectors = [rank_collector_from_name(method_name, partition) for method_name in rank_calculation_methods]
    rank_comparer = CutRankCalculatorComparer(partition)
    cut_rank = partition.cut_rank
    report(partition, "Initial")

    for _ in range(iterations):

        subset1, subset2 = None, None
        while subset1 == subset2:
            subset1 = random.choice(partition.subsets)[:]
            subset2 = random.choice(partition.subsets)[:]
        
        rank_comparer.reset(subset1, subset2)
        for coll in rank_collectors:
            rank_comparer.calculate_and_compare(coll, log=False)
        if rank_comparer.is_reset():
            raise Exception("No ranks have been calculated")

        node1 = random.choice(subset1)
        node2 = random.choice(subset2)

        print(f"Swapping node {node1} and node {node2}")
        partition.apply_swap(node1, node2)
        cut_rank = rank_comparer.first_cut_ranks[node1][node2]
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

        subset = random.choice(partition.subsets)[:]
        not_in_partition = [index for (index, subset) in enumerate(partition.subset_index) if subset == -1]

        if target == "row":
            rows, columns = subset, not_in_partition
        elif target == "column":
            rows, columns = not_in_partition, subset
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
