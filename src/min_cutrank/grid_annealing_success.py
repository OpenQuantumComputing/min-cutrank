import numpy as np
import sys
import getopt
import random
import time
from test_tools import parse_int, parse_float, temperatures_from_description
from partition_builder import random_partition, grid_graph
from cut_rank_annealing import cut_rank_annealing_row_formula


if __name__=="__main__":

    """
    Program testing how successful the annealing algorithgm is on NxN grids for a range of N.

    For each N, the annealing algorithm is run on a specific number of NxN grids with a random partition. For these algorithm runs, the average final cut-rank, the number of times the final cut-rank reaches the known optimal value N,
    and the total time spent on building and running the algorithm for all the NxN grids are collected.
    The result is stored on an output file if given.

    Parameters:
    -s N        The random seed. If omited, no seed is set for the random function. Random numbers are used to select the initial graph partitions, and to select which swaps to apply during the annealing algorithm.
    -r Range    The range of N for the NxN grids tested
    -n Samples  The number of NxN grids to run the algorithm on for each N.
    -p P        The size of the first partition set as a portion of the number of all nodes. Default is 0.5.
    -t Temp     The temperature setup. See 'temperatures_from_description'. Default is '1e0.1s10', i.e. 10 temperatures on a linear range from 1.0 to 0.1
    -o Outfile  The path to the output file. If absent, not output is written to file.
    """

    opt_arguments = sys.argv[1:]

    seed = None

    range_limits = []
    set_portion = 0.5
    temperatures = np.linspace(1.0, 0.1, 10)
    samples = -1
    file_path_out = None

    options = "s:r:n:p:t:o:"
    long_options = ["seed=", "range=", "samples=", "partition_portion=", "temperatures=", "output_file="]

    try:
        arguments, values = getopt.getopt(opt_arguments, options, long_options)

        for argument, value in arguments:

            if argument in ("-s", "--seed"):
                seed = parse_int(value, None)
            elif argument in ("-r", "--range"):
                range_limits = [int(s) for s in value.split("-")]
            elif argument in ("-n", "--samples"):
                samples = parse_int(value, None)
            elif argument in ("-p", "--partition_portion"):
                set_portion = parse_float(value, 0.5)
            elif argument in ("-t", "--temperatures"):
                temperatures = temperatures_from_description(value)
            elif argument in ("-o", "--output_file"):
                file_path_out = value.replace("\\","/")

        if len(range_limits) != 2:
            print("Grid size range is missing")

        elif samples <= 0:
            print("Number of samples must be positive")

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
                print(f"Running {nmb_sizes} grid sizes from {min_size}x{min_size} to {max_size}x{max_size} with {samples} samples for each size")
                if file_path_out == None:
                    print("Output file not given, only writing to console")
                else:
                    print(f"Output will be stored on {file_path_out}")

                if seed != None:
                    random.seed(seed)
                results = {}
                for size in range(min_size, above_max_size):
                    adj_mat = grid_graph(size, size)
                    nmb_exp_rank = 0
                    sum_rank = 0
                    print(f"Starting running {samples} of grid size {size}x{size}")
                    start_time = time.time()

                    for _ in range(samples):
                        partition = random_partition(adj_mat, set_portion)
                        cut_rank_annealing_row_formula(partition, temperatures, False)
                        cut_rank = partition.cut_rank

                        if cut_rank < size:
                            print(f"Got rank {cut_rank}, below expected minimum {size}")
                            raise Exception("Unexpected rank")
                        elif cut_rank == size:
                            nmb_exp_rank += 1
                        
                        sum_rank += cut_rank

                    time_total = time.time() - start_time
                    print(f"Completed in {time_total} sec, {nmb_exp_rank}/{samples} got minimum rank, average rank was {(sum_rank / samples)}")
                    results[size] = {"success" : nmb_exp_rank, "sumRanks" : sum_rank, "time" : time_total}

                if file_path_out != None:
                    with open(file_path_out, "w") as outfile:
                        outfile.write("Size\tSamples\tSuccess\tAvg rank\tSeconds\n")
                        for size in range(min_size, above_max_size):
                            outfile.write(str(size) + "\t" + str(samples) + "\t" + str(results[size]["success"]) + "\t" + str(results[size]["sumRanks"] / samples) + "\t" + str(results[size]["time"]) + "\n")
                        print(f"Results written to {file_path_out}")



    except getopt.error as err:
        print(str(err))
