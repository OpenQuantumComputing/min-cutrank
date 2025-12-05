# min-cutrank 
This repo contains a scalable metaheuristic for finding graph bipartitions with minimal cut rank, a key metric for distributed measurement-based quantum computing. The implementation combines matrix-based cut-rank evaluations with swap heuristics and simulated annealing to balance optimality and runtime on both grid and sparse topologies. Benchmarks and testing utilities demonstrate how the approach scales across graph families and highlight practical trade-offs when tuning annealing schedules. For more details, see the [paper](./Effective%20Partitioning%20for%20Distributed%20Measurement-Based%20Quantum%20Computing.pdf).

## Graph partitioning algorithm for minizing cut-rank by matrix investigations

- The core object for calculating cut-ranks is GraphPartition from graph_partition.py. The constructor takes two arguments:
  - The adjacency matrix, as a list of lists of int, like [[0, 1, 0, 0], [1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 1, 0]]
  - The initial partition of the nodes, as a list of booleans where elemnet n is True iff node n belongs to partition set 1, like [True, False, False, True]
- Use 'single_swap_cut_rank' from swap_rank_calculator.py to find the cut-rank of one single swap of a specific row and a specific column node. It has time complexity O(n).
- Use 'row_swap_cut_rank' from swap_rank_calculator.py to find the cut-ranks for all swapping combinations of a specific row and any column. It has time complexity O(n^2), but should be faster than doing 'single_swap_cut_rank' for all swaps.
- Use 'all_swap_cut_rank' from swap_rank_calculator.py to find the cut_ranks for all swapping combinations of any row and any column. It has time complexity O(n^2).
- Use the 'apply_swap' method on a GraphPartition object to apply a swap and update all necessary matrices for further swap cut-rank calculations. It should have time complexity O(n^2).

## Annealing algorithm

- Use the method 'cut_rank_annealing_row_formula' from cut_rank_annealing.py to run the annealing algorithm using matrix investigations for the cut-ranks.

## Programs related to cut-rank calculations (see each file for more information)

### Testing

- test_cut_rank.py: Test program for verifying the swap cut-rank formulas and for validating the variables in the GraphPartition object.
- test_annealing.py: Test program for the annealing algorithm.

### Collecting computational results

- compare_grid_annealing.py: Program collecting time measures on the two grid annealing algorithms on NxN grids for a range of N.
- grid_annealing_success.py: Program testing how successful the annealing algorithgm is on NxN grids for a range of N.
- sparse_annealing.py:  Program testing the annealing algorithgm on random sparse graphs of N nodes and c/N probability for each edge for given input constant c

## Results

See results\overview.txt for details
