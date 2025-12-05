import sys
import getopt
import random
from test_tools import parse_int, parse_float, graph_from_description, run_greedy_min_rank
from partition_builder import random_partition


if __name__=="__main__":

    """
    Test program for verifying the swap cut-rank formulas and for validating the variables in the GraphPartition object.

    The program starts with a random partition of a given size over a graph and runs a greedy algorithm that applies swaps as long as there are improvments possible.
    At each iteration, the cut ranks for all possible swaps are calculated by the different cut-rank calculation methods selected by the paramters,
    if two different methods give different rank for the same swap, an exception is raised.
    One of the cut-rank methods will also do a validitation test on all the variables in the GraphPartition object.

    Parameters:
    -s N        The random seed. If omited, no seed is set for the random function. Random numbers can be used to build the graph, and is used to select which swap to apply if several give the same rank improvement.
    -g Graph    The graph setup. See 'graph_from_description' for details
    -p P        The size of the first partition set as a portion of the number of all nodes. Default is 0.5.
    -m Methods  Lists the cut-rank methods to be tested, separated by comma. The alternatives are
                'gauss' calculates each swap cut-rank by Gauss-Jordan elimination on the adjacency matrix
                'single' calculates the swap cut-ranks by #(set 1) x #(set 2) separate calls to 'single_swap_cut_rank'
                'row' calculates the swap cut-ranks by #(set 1) separate calls to 'row_swap_cut_ranks'
                'all' calculates the swap cut-ranks by one single call to 'all_swap_cut_ranks'
                'apply' calculates the swap cut-ranks by actually applying the swaps to the GraphPartition object
                'validate' does the same as 'apply', but also validates all the variables of the GraphPartition object after the swap has been applied
    """

    opt_arguments = sys.argv[1:]

    seed = None
    graph_setup = None
    set_portion = 0.5
    rank_calculation_methods = []

    options = "s:g:p:m:"
    long_options = ["seed=", "graph=", "partition_portion=", "methods="]

    try:
        arguments, values = getopt.getopt(opt_arguments, options, long_options)

        for argument, value in arguments:

            if argument in ("-s", "--seed"):
                seed = parse_int(value, None)
            elif argument in ("-g", "--graph"):
                graph_setup = value
            elif argument in ("-p", "--partition_portion"):
                set_portion = parse_float(value, 0.5)
            elif argument in ("-m", "--methods"):
                rank_calculation_methods = value.split(",")
            elif argument in ("-d", "--directly"):
                ranks_directly = True
            elif argument in ("-a", "--apply"):
                ranks_by_apply = True
                validate = value == "validate"

        if graph_setup == None:
            print("Graph setup missing")

        else:
            if seed != None:
                random.seed(seed)
            graph_adj_matrix = graph_from_description(graph_setup)
            graph_partition = random_partition(graph_adj_matrix, set_portion)

            run_greedy_min_rank(graph_partition, rank_calculation_methods)

    except getopt.error as err:
        print(str(err))
