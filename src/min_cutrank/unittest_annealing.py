import unittest
import sys
import random
from min_cutrank.compare_grid_annealing import compare_grid_annealing
from min_cutrank.grid_annealing_success import grid_annealing_success
from min_cutrank.sparse_annealing import run_sparse_annealing
from min_cutrank.test_annealing import run_annealing
from min_cutrank.test_cut_rank import run_greedy

class TestAnnealing(unittest.TestCase):
    """Tests for runnning the top level programs. They do not verify any actual results, 
    just that the code runs without errors."""

    def test_run_greedy_validate(self):
        
        run_greedy("-s 12345 -g r20 -p 0.5 -m validate".split(" "))
    
    
    def test_run_sparse_annealing(self):

        run_sparse_annealing("-s 12345 -r 10-15 -c 2.0 -n 5 -p 0.5".split(" "))
    

    def test_run_annealing(self):

        run_annealing("-s 12345 -g g5x6 -m gauss,formula".split(" "))


    def test_run_compare_grid_annealing(self):

        compare_grid_annealing("-s 12345 -r 5-6".split(" "))


    def test_run_grid_annealing_success(self):

        grid_annealing_success("-s 12345 -r 5-6 -n 2 ".split(" "))




if __name__ == "__main__":
    unittest.main()
