import numpy as np
import random
from matrix_tools import copy_matrix, rank_matrix_positions
from swap_rank_calculator import row_swap_cut_ranks

from graph_partition import GraphPartition


def cut_rank_annealing_direct(partition : GraphPartition, temperatures, log: bool) -> None:

    rows = partition.rows[:]
    cols = partition.columns[:]
    cut_rank = partition.cut_rank
    nmb_rows = len(rows)
    nmb_cols = len(cols)
    if log:
        print(f"Starting with cut-rank {cut_rank}")

    for temp in temperatures:
        limits = [np.exp(-1.0 / temp), np.exp(-2.0 / temp)]

        for i in range(nmb_rows):
            for j in range(nmb_cols):

                rows[i], cols[j] = cols[j], rows[i]

                copy_matrix(partition.adjacencies, partition.buffer, rows, cols)
                base_rows, _ = rank_matrix_positions(partition.buffer, rows, cols)
                new_cut_rank = len(base_rows)
                delta_rank = new_cut_rank - cut_rank

                if delta_rank <= 0 or random.random() < limits[delta_rank - 1]:
                    cut_rank = new_cut_rank
                else:
                    rows[i], cols[j] = cols[j], rows[i]

        if log:
            print(f"Cut-rank is {cut_rank} after sweep with temperature {temp}")


def cut_rank_annealing_row_formula(partition : GraphPartition, temperatures, log: bool) -> None:

    rows = partition.rows[:]
    cols = partition.columns[:]
    row_ranks = [-1] * partition.nmb_nodes
    cut_rank = partition.cut_rank
    nmb_rows = len(rows)
    nmb_cols = len(cols)
    if log:
        print(f"Starting with cut-rank {cut_rank}")

    for temp in temperatures:
        limits = [np.exp(-1.0 / temp), np.exp(-2.0 / temp)]

        for i in range(nmb_rows):
            row = rows[i]
            for n in partition.nodes:
                row_ranks[n] = -1
            row_swap_cut_ranks(partition, row, row_ranks)
            swap_col = -1
            for j in range(nmb_cols):

                new_cut_rank = row_ranks[cols[j]]
                if new_cut_rank < 0:
                    raise Exception("Cut-rank not calculated")
                delta_rank = new_cut_rank - cut_rank
                if delta_rank <= 0 or random.random() < limits[delta_rank - 1]:
                    swap_col = cols[j]
                    rows[i], cols[j] = cols[j], rows[i]
                    cut_rank = new_cut_rank

            if swap_col >= 0:
                partition.apply_swap(row, swap_col)
            if partition.cut_rank != cut_rank:
                raise Exception("Partition cut-rank does not fit with directly calculated rank")

        if log:
            print(f"Cut-rank is {cut_rank} after sweep with temperature {temp}")
