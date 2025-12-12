import numpy as np
import sys
import getopt
import random
import time
from test_tools import parse_int, parse_float, temperatures_from_description, clone_partition
from partition_builder import random_partition, grid_graph
from cut_rank_annealing import cut_rank_annealing_direct, cut_rank_annealing_row_formula


def run_annealing_method(annealing_method, partition, temperatures) -> None:
    partition_copy = clone_partition(partition)
    start_time = time.time()
    annealing_method(partition_copy, temperatures, False)
    time_total = time.time() - start_time
    print(f"Completed in {time_total} sec")
    return time_total


if __name__=="__main__":

    """
    Program collecting time measures on the two grid annealing algorithms on NxN grids for a range of N.

    The two algorithms to be compared are using the different swap cut-rank calculation methods, either Gauss-Jordan elimination on the adjacency matrix, or by matrix inspection with calls to 'row_swap_cut_ranks'
    The same partition and random initialization is used for both algorithms, so the swap selections made during the annealing algorithm should be the same for the two,
    the only difference is the calculation of the cut-rank.
    The result is stored on an output file if given.

    Parameters:
    -s N        The random seed. If omited, no seed is set for the random function. Random numbers can be used to build the graph, and are used to select which swaps to apply during the annealing algorithm.
    -r Range    The range of N for the NxN grids tested
    -p P        The size of the first partition set as a portion of the number of all nodes. Default is 0.5.
    -t Temp     The temperature setup. See 'temperatures_from_description'. Default is '1e0.1s10', i.e. 10 temperatures on a linear range from 1.0 to 0.1
    -o Outfile  The path to the output file. If absent, not output is written to file.
    """

    opt_arguments = sys.argv[1:]

    seed = None

    range_limits = []
    set_portion = 0.5
    temperatures = np.linspace(1.0, 0.1, 10)
    file_path_out = None

    options = "s:r:p:t:o:"
    long_options = ["seed=", "range=", "partition_portion=", "temperatures=", "output_file="]

    try:
        arguments, values = getopt.getopt(opt_arguments, options, long_options)

        for argument, value in arguments:

            if argument in ("-s", "--seed"):
                seed = parse_int(value, None)
            elif argument in ("-r", "--range"):
                range_limits = [int(s) for s in value.split("-")]
            elif argument in ("-p", "--partition_portion"):
                set_portion = parse_float(value, 0.5)
            elif argument in ("-t", "--temperatures"):
                temperatures = temperatures_from_description(value)
            elif argument in ("-o", "--output_file"):
                file_path_out = value.replace("\\","/")

        if len(range_limits) != 2:
            print("Grid size range is missing")

        else:
            min_size = range_limits[0]
            max_size = range_limits[1]
            above_max_size = max_size + 1
            nmb_sizes = above_max_size - min_size

            if nmb_sizes < 0:
                print("Minimum grid size can not be greater than maximum")

            elif min_size < 3:
                print("Minimum grid must be at least 3")

            else:
                print(f"Running {nmb_sizes} grid sizes from {min_size}x{min_size} to {max_size}x{max_size} for algorithm comparisons")
                if file_path_out == None:
                    print("Output file not given, only writing to console")
                else:
                    print(f"Output will be stored on {file_path_out}")

                if seed != None:
                    random.seed(seed)
                results = {}
                for size in range(min_size, above_max_size):
                    seed_algo = random.randint(0, 65535)
                    adj_mat = grid_graph(size, size)
                    partition = random_partition(adj_mat, set_portion)

                    random.seed(seed_algo)
                    print(f"Running annealing algorithm on grid with size {size}x{size} by matrix investigations")
                    time_formula = run_annealing_method(cut_rank_annealing_row_formula, partition, temperatures)

                    random.seed(seed_algo)
                    print(f"Running annealing algorithm on grid with size {size}x{size} by direct rank calculation")
                    time_direct = run_annealing_method(cut_rank_annealing_direct, partition, temperatures)

                    results[size] = {"formula" : time_formula, "direct" : time_direct}

                if file_path_out != None:
                    with open(file_path_out, "w") as outfile:
                        outfile.write("Size\tRank by matrix inspection\tRank by Gauss-Jordan elimination\n")
                        for size in range(min_size, above_max_size):
                            outfile.write(str(size) + "\t" + str(results[size]["formula"]) + "\t" + str(results[size]["direct"]) + "\n")
                        print(f"Results written to {file_path_out}")

    except getopt.error as err:
        print(str(err))
