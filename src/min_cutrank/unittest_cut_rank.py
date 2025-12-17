import unittest
import sys
import random
from min_cutrank.partition_builder import random_partition
from min_cutrank.graph_partition import GraphPartition
from min_cutrank.test_tools import CutRankCalculatorComparer, graph_from_description, rank_collector_from_name

class TestCutRankCalculation(unittest.TestCase):
    """Tests for verifying that different cut-rank calculation methods give the same results."""

    def test_rank_calculation_matches_gauss(self):

        for method in ["single", "row", "all"]:
            for partial in [False, True]:
                with self.subTest(method=method, partial=partial):
            
                    graph_partition = self.createPartition(partial_partition=partial)
                    rank_calculation_methods = ["gauss", method]
                    
                    verify_random_swaps(graph_partition, rank_calculation_methods, 100)
    
    
    def createPartition(self, size: int = 20, partial_partition: bool = False) -> GraphPartition:
        graph_setup = "r20"
        graph_adj_matrix = graph_from_description(graph_setup)
        graph_partition = random_partition(graph_adj_matrix, 0.3, 0.3) if partial_partition \
            else random_partition(graph_adj_matrix, 0.5)
        return graph_partition



def verify_random_swaps(partition : GraphPartition, rank_calculation_methods: list[str], iterations: int) -> None:
    """Runs random swaps on the partition and verifies that the cut-rank calculated by different methods match."""

    rank_collectors = [rank_collector_from_name(method_name, partition) for method_name in rank_calculation_methods]
    rank_comparer = CutRankCalculatorComparer(partition)
    cut_rank = partition.cut_rank
    print(f"Initial cut-rank is {cut_rank}")
    print("Initial partition:")
    print(f"Set 1 = {partition.rows}")
    print(f"Set 2 = {partition.columns}")

    for _ in range(iterations):

        rank_comparer.reset()
        for coll in rank_collectors:
            rank_comparer.calculate_and_compare(coll, log=False)
        if rank_comparer.is_reset():
            raise Exception("No ranks have been calculated")

        row = random.choice(partition.rows)
        col = random.choice(partition.columns)

        print(f"Swapping row {row} and column {col}")
        partition.apply_swap(row, col)
        cut_rank = rank_comparer.first_cut_ranks[row][col]
        print(f"New cut-rank is {cut_rank}")
        if cut_rank != partition.cut_rank:
            raise Exception("Unexpected rank after swap")

    print()
    print(f"Final rank is {cut_rank}")
    print("Final partition:")
    print(f"Set 1 = {partition.rows}")
    print(f"Set 2 = {partition.columns}")


if __name__ == "__main__":
    unittest.main()
