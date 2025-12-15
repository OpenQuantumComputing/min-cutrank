import unittest
import sys
import random
from min_cutrank.partition_builder import random_partition
from min_cutrank.graph_partition import GraphPartition
from min_cutrank.test_tools import CutRankCalculatorComparer, graph_from_description, rank_collector_from_name

class TestCutRankCalculation(unittest.TestCase):
    """Tests for verifying that different cut-rank calculation methods give the same results."""

    def test_single_rank_calculation_matches_gauss(self):
        graph_setup = "r20"
        graph_adj_matrix = graph_from_description(graph_setup)
        graph_partition = random_partition(graph_adj_matrix, 0.5)
        
        rank_calculation_methods = ["gauss", "single"]
        
        run_random_min_rank(graph_partition, rank_calculation_methods, 100)
    
    
    def test_row_rank_calculation_matches_gauss(self):
        graph_setup = "r20"
        graph_adj_matrix = graph_from_description(graph_setup)
        graph_partition = random_partition(graph_adj_matrix, 0.5)
        
        rank_calculation_methods = ["gauss", "row"]
        
        run_random_min_rank(graph_partition, rank_calculation_methods, 100)
    

    def test_all_rank_calculation_matches_gauss(self):
        graph_setup = "r20"
        graph_adj_matrix = graph_from_description(graph_setup)
        graph_partition = random_partition(graph_adj_matrix, 0.5)
        
        rank_calculation_methods = ["gauss", "all"]
        
        run_random_min_rank(graph_partition, rank_calculation_methods, 100)



def run_random_min_rank(partition : GraphPartition, rank_calculation_methods: list[str], iterations: int) -> None:

    rank_collectors = [rank_collector_from_name(method_name, partition) for method_name in rank_calculation_methods]
    rank_comparer = CutRankCalculatorComparer(partition)
    local_minimum_found = False
    cut_rank = partition.cut_rank
    iteration = 0


    for _ in range(iterations):

        rank_comparer.reset()
        for coll in rank_collectors:
            rank_comparer.calculate_and_compare(coll)
        if rank_comparer.is_reset():
            raise Exception("No ranks have been calculated")

        row = random.choice(partition.rows)
        col = random.choice(partition.columns)

        partition.apply_swap(row, col)
        cut_rank = rank_comparer.first_cut_ranks[row][col]
        if cut_rank != partition.cut_rank:
            raise Exception("Unexpected rank after swap")

    print()
    print(f"Stopped after {iteration} iterations")
    print(f"Final rank is {cut_rank}")
    print("Final partition:")
    print(f"Set 1 = {partition.rows}")
    print(f"Set 2 = {partition.columns}")


if __name__ == "__main__":
    unittest.main()
