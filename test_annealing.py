import numpy as np
import sys
import time
import getopt
import random
from cut_rank_annealing import cut_rank_annealing_direct, cut_rank_annealing_row_formula
from test_tools import parse_bool,parse_int, parse_float, graph_from_description, temperatures_from_description, clone_partition
from partition_builder import random_partition


def test_annealing_method(annealing_method, name : str, partition, temperatures, log : bool) -> None:
    print(f"Testing annealing method '{name}'")
    partition_copy = clone_partition(partition)
    start = time.time()
    annealing_method(partition_copy, temperatures, log)
    end = time.time()
    print(f"Annealing method '{name}' completed at cut-rank {partition_copy.cut_rank} in {end - start} sec")


if __name__=="__main__":

    """
    Test program for the annealing algorithm.

    The program starts with a random partition of a given size over a graph and runs each of the specified annealing algorithm implementations on the graph.
    The same graph partition, temperatures and randomness initialization is used on all algorithms, so they are all expected to give the same intermediate and finel outcomes.
    TO-DO: Compare ranks after each row selection for the different algorithms, raise Exception at deviations.

    Parameters:
    -s N        The random seed. If omited, no seed is set for the random function. Random numbers can be used to build the graph, and are used to select which swaps to apply during the annealing algorithm.
    -g Graph    The graph setup. See 'graph_from_description' for details
    -p P        The size of the first partition set as a portion of the number of all nodes. Default is 0.5.
    -t Temp     The temperature setup. See 'temperatures_from_description'. Default is '1e0.1s10', i.e. 10 temperatures on a linear range from 1.0 to 0.1
    -m Methods  Lists the annealing algorithms to be used, separated by comma. The alternatives are
                'gauss' calculates each swap cut-rank by Gauss-Jordan elimination on the adjacency matrix
                'formula' calculates the swap cut-ranks for each selected element in the first partition set by one single call to 'row_swap_cut_ranks'
    -l Bool     Whether the rank at the beginning and after each temperature sweep should be logged to the console.
    """

    opt_arguments = sys.argv[1:]

    seed = None
    graph_setup = None
    set_portion = 0.5

    cut_rank_methods = []
    temperatures = np.linspace(1.0, 0.1, 10)
    log = True

    options = "s:g:p:t:m:l:"
    long_options = ["seed=", "graph=", "partition_portion=", "temperatures=", "cut_rank_methods=", "log="]

    try:
        arguments, values = getopt.getopt(opt_arguments, options, long_options)

        for argument, value in arguments:

            if argument in ("-s", "--seed"):
                seed = parse_int(value, None)
            elif argument in ("-g", "--graph"):
                graph_setup = value
            elif argument in ("-p", "--partition_portion"):
                set_portion = parse_float(value, 0.5)
            elif argument in ("-t", "--temperatures"):
                temperatures = temperatures_from_description(value)
            elif argument in ("-m", "--cut_rank_methods"):
                cut_rank_methods = value.split(",")
            elif argument in ("-l", "--log"):
                log = parse_bool(value, False)

        if graph_setup == None:
            print("Graph setup missing")

        else:
            if seed != None:
                random.seed(seed)
            graph_adj_matrix = graph_from_description(graph_setup)
            graph_partition = random_partition(graph_adj_matrix, set_portion)

            seed_algo = random.randint(0, 65535)
            for cut_rank_m in cut_rank_methods:
                random.seed(seed_algo)

                if cut_rank_m == "gauss":
                    test_annealing_method(cut_rank_annealing_direct, "Gauss-Jordan elimination cut-rank calculation", graph_partition, temperatures, log)

                elif cut_rank_m == "formula":
                    test_annealing_method(cut_rank_annealing_row_formula, "Formula for all cut-ranks on row", graph_partition, temperatures, log)

                else:
                    print(f"Unknown cut-rank annealing method: '{cut_rank_m}'")

    except getopt.error as err:
        print(str(err))
