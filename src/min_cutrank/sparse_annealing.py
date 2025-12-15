import numpy as np
import sys
import getopt
import random
import time
from test_tools import parse_int, parse_float, temperatures_from_description
from partition_builder import random_partition, random_graph
from cut_rank_annealing import cut_rank_annealing_row_formula


if __name__=="__main__":

    """
    Program testing the annealing algorithgm on random sparse graphs of N nodes and c/N probability for each edge for given input constant c

    Erdös-Rényi 

    For each N, the annealing algorithm is run on a specific number of random Erdös-Rényi graphs G(N,p) of N vertices, and where p = c/N is the probability that each edge in the complete N-graph appears.
    The value of c is the same for all N, the startup partition is a random partition of specific size given as a portion of all nodes, and the operations in the algorithm are single element swaps
    of pairs of elements from the two partition sets. For these algorithm runs, the average final cut-rank and the total time spent on building and running the algorithm for each N are collected.
    The results are stored on an output file if given.

    Parameters:
    -s N        The random seed. If omited, no seed is set for the random function. Random numbers are used to select the initial graph partitions, and to select which swaps to apply during the annealing algorithm.
    -r Range    The range of N, the number of nodes in the graphs to be tested
    -c C-factor N times the probability for the appearance of each edge in the complete N-graph
    -n Samples  The number of graphs to run the algorithm on for each N.
    -p P        The size of the first partition set as a portion of the number of all nodes. Default is 0.5.
    -t Temp     The temperature setup. See 'temperatures_from_description'. Default is '1e0.1s10', i.e. 10 temperatures on a linear range from 1.0 to 0.1
    -o Outfile  The path to the output file. If absent, not output is written to file.
    """

    opt_arguments = sys.argv[1:]

    seed = None

    range_limits = []
    edge_probability_factor = -1.0
    set_portion = 0.5
    temperatures = np.linspace(1.0, 0.1, 10)
    samples = -1
    file_path_out = None

    options = "s:r:c:n:p:t:o:"
    long_options = ["seed=", "range=", "edge_probability_denominator=", "samples=", "partition_portion=", "temperatures=", "output_file="]

    try:
        arguments, values = getopt.getopt(opt_arguments, options, long_options)

        for argument, value in arguments:

            if argument in ("-s", "--seed"):
                seed = parse_int(value, None)
            elif argument in ("-r", "--range"):
                range_limits = [int(s) for s in value.split("-")]
            elif argument in ("-c", "--edge_probability_denominator"):
                edge_probability_factor = parse_float(value, -1.0)
            elif argument in ("-n", "--samples"):
                samples = parse_int(value, None)
            elif argument in ("-p", "--partition_portion"):
                set_portion = parse_float(value, 0.5)
            elif argument in ("-t", "--temperatures"):
                temperatures = temperatures_from_description(value)
            elif argument in ("-o", "--output_file"):
                file_path_out = value.replace("\\","/")

        if len(range_limits) != 2:
            print("Graph size range is missing")

        elif samples <= 0:
            print("Number of samples must be positive")

        elif edge_probability_factor < 0:
            print("c-value in edge probability c/N must be positive")

        else:
            min_size = range_limits[0]
            max_size = range_limits[1]
            above_max_size = max_size + 1
            nmb_sizes = above_max_size - min_size

            if nmb_sizes < 0:
                print("Minimum graph size can not be greater than maximum")

            elif min_size < 2:
                print("Minimum graph size must be at least 2")

            else:
                print(f"Running {nmb_sizes} graph sizes from {min_size} to {max_size} with {samples} samples for each size")
                if file_path_out == None:
                    print("Output file not given, only writing to console")
                else:
                    print(f"Output will be stored on {file_path_out}")

                if seed != None:
                    random.seed(seed)
                results = {}
                for size in range(min_size, above_max_size):

                    edge_prob = edge_probability_factor / size
                    sum_rank = 0
                    print(f"Starting running {samples} of graph size {size} and edge probability {edge_prob}")
                    start_time = time.time()

                    for _ in range(samples):
                        adj_mat = random_graph(size, edge_prob)
                        partition = random_partition(adj_mat, set_portion)
                        cut_rank_annealing_row_formula(partition, temperatures, False)
                        cut_rank = partition.cut_rank

                        sum_rank += cut_rank

                    time_total = time.time() - start_time
                    print(f"Completed in {time_total} sec, average rank was {(sum_rank / samples)}")
                    results[size] = {"sumRanks" : sum_rank, "time" : time_total}

                if file_path_out != None:
                    with open(file_path_out, "w") as outfile:
                        outfile.write("Size\tSamples\tAvg rank\tSeconds\n")
                        for size in range(min_size, above_max_size):
                            outfile.write(str(size) + "\t" + str(samples) + "\t" + str(results[size]["sumRanks"] / samples) + "\t" + str(results[size]["time"]) + "\n")
                        print(f"Results written to {file_path_out}")



    except getopt.error as err:
        print(str(err))
