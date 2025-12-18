from min_cutrank.graph_partition import GraphPartition


def all_swap_cut_ranks(partition : GraphPartition, rows_to_swap: list[int], columns_to_swap: list[int], rows_are_in_partition: bool, columns_are_in_partition: bool, ranks : list[list[int]]) -> None:
    """Finds the cut-ranks for the partitions obtained by swapping any of the given rows with any of the given columns in the given graph partition.
    
    args:
        - partition: 'GraphPartition' The graph partition.
        - rows_to_swap: 'list[int]' The rows to be swapped.
        - columns_to_swap: 'list[int]' The columns to be swapped.
        - rows_are_in_partition: 'bool' If true, rows_to_swap are the partition's rows.
            If false, rows_to_swap are disjunct from the partition's rows and columns.
        - columns_are_in_partition: 'bool' If true, columns_to_swap are the partition's columns.
            If false, columns_to_swap are disjunct from the partition's rows and columns.
        - ranks: 'list[list[int]]' A matrix where position [i][j] will hold the cut-rank after swapping node i and j. 
            Only positions where i is in rows_to_swap and j is in columns_to_swap will be affected.
    """

    old_rank = partition.cut_rank
    nmb_nodes = partition.graph.nmb_nodes
    base_rows = partition.base_rows if rows_are_in_partition else []
    free_rows = partition.free_rows if rows_are_in_partition else rows_to_swap    
    base_columns = partition.base_columns if columns_are_in_partition else []
    free_columns = partition.free_columns if columns_are_in_partition else columns_to_swap
    
    # Preprocessing on rows
    s1_k1 = [-1] * nmb_nodes
    s2 = [False] * nmb_nodes
    q4_952_0 = [False] * nmb_nodes
    q4_952_1 = [False] * nmb_nodes
    for i in partition.base_rows:
        k1 = next((k1 for k1 in partition.free_rows if partition.adj_b_inverse[k1][i] == 1), -1)
        s1_k1[i] = k1
        if columns_are_in_partition:
            if k1 >= 0 and partition.adj_b_inv_adj[k1][i] == 1:
                s2[i] = any(k2 != k1 and partition.adj_b_inv_adj[k2][i] != partition.adj_b_inverse[k2][i] for k2 in partition.free_rows)
            else:
                s2[i] = any(partition.adj_b_inv_adj[k2][i] == 1 for k2 in partition.free_rows)
            q4_952_0[i] = any(partition.adj_b_inv_adj[k][i] == 1 for k in partition.free_rows)
            q4_952_1[i] = any(partition.adj_b_inv_adj[k][i] != partition.adj_b_inverse[k][i] for k in partition.free_rows)
    if columns_are_in_partition:
        for i in free_rows:
            s2[i] = any(k2 != i and partition.adj_b_inv_adj[k2][i] == 1 for k2 in partition.free_rows)

    # Preprocessing on columns
    t1_l1 = [-1] * nmb_nodes
    t2 = [False] * nmb_nodes
    q5_952_0 = [False] * nmb_nodes
    q5_952_1 = [False] * nmb_nodes
    for j in partition.base_columns:
        l1 = next((l1 for l1 in partition.free_columns if partition.b_inverse_adj[j][l1] == 1), -1)
        t1_l1[j] = l1
        if rows_are_in_partition:
            if l1 >= 0 and partition.adj_b_inv_adj[j][l1] == 1:
                t2[j] = any(l2 != l1 and partition.adj_b_inv_adj[j][l2] != partition.b_inverse_adj[j][l2] for l2 in partition.free_columns)
            else:
                t2[j] = any(partition.adj_b_inv_adj[j][l2] == 1 for l2 in partition.free_columns)
            q5_952_0[j] = any(partition.adj_b_inv_adj[j][l] == 1 for l in partition.free_columns)
            q5_952_1[j] = any(partition.adj_b_inv_adj[j][l] != partition.b_inverse_adj[j][l] for l in partition.free_columns)
    if rows_are_in_partition:
        for j in free_columns:
            t2[j] = any(l2 != j and partition.adj_b_inv_adj[j][l2] == 1 for l2 in partition.free_columns)

    # Ranks for i in X^D and j in Y^D
    for i in free_rows:
        for j in free_columns:
            if s2[i]:
                if t2[j]:
                    ranks[i][j] = old_rank + 2
                else:
                    ranks[i][j] = old_rank + 1
            else:
                if t2[j]:
                    ranks[i][j] = old_rank + 1
                else:
                    if rows_are_in_partition and columns_are_in_partition and partition.adj_b_inv_adj[j][i] == 1:
                        ranks[i][j] = old_rank + 1
                    else:
                        ranks[i][j] = old_rank

    # Ranks for i in X^B and j in Y^D
    for i in base_rows:
        k1 = s1_k1[i]
        for j in free_columns:
            if k1 >= 0:
                if s2[i]:
                    if t2[j]:
                        ranks[i][j] = old_rank + 2
                    else:
                        ranks[i][j] = old_rank + 1
                else:
                    if t2[j]:
                        ranks[i][j] = old_rank + 1
                    else:
                        if columns_are_in_partition and partition.adj_b_inv_adj[j][i] != (partition.adj_b_inverse[j][i] & partition.adj_b_inv_adj[k1][i]):
                            ranks[i][j] = old_rank + 1
                        else:
                            ranks[i][j] = old_rank
            else:
                if partition.adj_b_inverse[j][i] == 1:
                    if s2[i]:
                        ranks[i][j] = old_rank + 1
                    else:
                        ranks[i][j] = old_rank
                else:
                    if s2[i]:
                        if t2[j]:
                            ranks[i][j] = old_rank + 1
                        else:
                            ranks[i][j] = old_rank
                    else:
                        if t2[j]:
                            ranks[i][j] = old_rank
                        else:
                            if columns_are_in_partition and partition.adj_b_inv_adj[j][i] == 1:
                                ranks[i][j] = old_rank
                            else:
                                ranks[i][j] = old_rank - 1

    # Ranks for i in X^D and j in Y^B
    for i in free_rows:
        for j in base_columns:
            l1 = t1_l1[j]
            if l1 >= 0:
                if t2[j]:
                    if s2[i]:
                        ranks[i][j] = old_rank + 2
                    else:
                        ranks[i][j] = old_rank + 1
                else:
                    if s2[i]:
                        ranks[i][j] = old_rank + 1
                    else:
                        if rows_are_in_partition and columns_are_in_partition and partition.adj_b_inv_adj[j][i] != (partition.b_inverse_adj[j][i] & partition.adj_b_inv_adj[j][l1]):
                            ranks[i][j] = old_rank + 1
                        else:
                            ranks[i][j] = old_rank
            else:
                if partition.b_inverse_adj[j][i] == 1:
                    if t2[j]:
                        ranks[i][j] = old_rank + 1
                    else:
                        ranks[i][j] = old_rank
                else:
                    if t2[j]:
                        if s2[i]:
                            ranks[i][j] = old_rank + 1
                        else:
                            ranks[i][j] = old_rank
                    else:
                        if s2[i]:
                            ranks[i][j] = old_rank
                        else:
                            if rows_are_in_partition and partition.adj_b_inv_adj[j][i] == 1:
                                ranks[i][j] = old_rank
                            else:
                                ranks[i][j] = old_rank - 1

    # Ranks for i in X^B and j in Y^B
    for i in base_rows:
        k1 = s1_k1[i]
        for j in base_columns:
            l1 = t1_l1[j]

            if partition.base_inverse[j][i] == 1:
                # Case 6

                if k1 >= 0 and l1 >= 0:
                    # Case 6.1
                    if s2[i]:
                        if t2[j]:
                            ranks[i][j] = old_rank + 2
                        else:
                            ranks[i][j] = old_rank + 1
                    else:
                        if t2[j]:
                            ranks[i][j] = old_rank + 1
                        else:
                            if ((partition.adj_b_inv_adj[k1][i] & partition.adj_b_inv_adj[j][l1]) ^ (partition.adj_b_inv_adj[k1][i] & partition.adj_b_inverse[j][i]) ^ (partition.adj_b_inv_adj[j][l1] & partition.b_inverse_adj[j][i])) != partition.adj_b_inv_adj[j][i]:
                                ranks[i][j] = old_rank + 1
                            else:
                                ranks[i][j] = old_rank

                else:
                    # Case 6.2
                    q4 = q4_952_1[i] if partition.b_inverse_adj[j][i] == 1 else q4_952_0[i]
                    q5 = q5_952_1[j] if partition.adj_b_inverse[j][i] == 1 else q5_952_0[j]
                    if q4:
                        if q5:
                            ranks[i][j] = old_rank + 1
                        else:
                            ranks[i][j] = old_rank
                    else:
                        if q5:
                            ranks[i][j] = old_rank
                        else:
                            if partition.adj_b_inv_adj[j][i] != (partition.adj_b_inverse[j][i] & partition.b_inverse_adj[j][i]):
                                ranks[i][j] = old_rank
                            else:
                                ranks[i][j] = old_rank - 1

            else:
                # Case 7

                    if k1 >= 0:

                        if l1 >= 0:
                            # Case 7.1
                            if s2[i]:
                                if t2[j]:
                                    ranks[i][j] = old_rank + 2
                                else:
                                    ranks[i][j] = old_rank + 1
                            else:
                                if t2[j]:
                                    ranks[i][j] = old_rank + 1
                                else:
                                    if ((partition.adj_b_inv_adj[k1][i] & partition.adj_b_inverse[j][i]) ^ (partition.adj_b_inv_adj[j][l1] & partition.b_inverse_adj[j][i])) != partition.adj_b_inv_adj[j][i]:
                                        ranks[i][j] = old_rank + 1
                                    else:
                                        ranks[i][j] = old_rank

                        else:
                            # Case 7.2
                            if t2[j]:
                                if partition.b_inverse_adj[j][i] == 1:
                                    ranks[i][j] = old_rank + 1
                                else:
                                    if s2[i]:
                                        ranks[i][j] = old_rank + 1
                                    else:
                                        ranks[i][j] = old_rank
                            else:
                                if partition.b_inverse_adj[j][i] == 1:
                                    ranks[i][j] = old_rank
                                else:
                                    if s2[i]:
                                        ranks[i][j] = old_rank
                                    else:
                                        if (partition.adj_b_inv_adj[k1][i] & partition.adj_b_inverse[j][i]) != partition.adj_b_inv_adj[j][i]:
                                            ranks[i][j] = old_rank
                                        else:
                                            ranks[i][j] = old_rank - 1

                    else:

                        if l1 >= 0:
                            # Case 7.3
                            if s2[i]:
                                if partition.adj_b_inverse[j][i] == 1:
                                    ranks[i][j] = old_rank + 1
                                else:
                                    if t2[j]:
                                        ranks[i][j] = old_rank + 1
                                    else:
                                        ranks[i][j] = old_rank
                            else:
                                if partition.adj_b_inverse[j][i] == 1:
                                    ranks[i][j] = old_rank
                                else:
                                    if t2[j]:
                                        ranks[i][j] = old_rank
                                    else:
                                        if (partition.adj_b_inv_adj[j][l1] & partition.b_inverse_adj[j][i]) != partition.adj_b_inv_adj[j][i]:
                                            ranks[i][j] = old_rank
                                        else:
                                            ranks[i][j] = old_rank - 1

                        else:
                            # Case 7.4
                            if partition.adj_b_inverse[j][i] == 1:
                                if partition.b_inverse_adj[j][i] == 1:
                                    ranks[i][j] = old_rank
                                else:
                                    if s2[i]:
                                        ranks[i][j] = old_rank
                                    else:
                                        ranks[i][j] = old_rank - 1
                            else:
                                if partition.b_inverse_adj[j][i] == 1:
                                    if t2[j]:
                                        ranks[i][j] = old_rank
                                    else:
                                        ranks[i][j] = old_rank - 1
                                else:
                                    if s2[i]:
                                        if t2[j]:
                                            ranks[i][j] = old_rank
                                        else:
                                            ranks[i][j] = old_rank - 1
                                    else:
                                        if t2[j]:
                                            ranks[i][j] = old_rank - 1
                                        else:
                                            if partition.adj_b_inv_adj[j][i] == 1:
                                                ranks[i][j] = old_rank - 1
                                            else:
                                                ranks[i][j] = old_rank - 2



def row_swap_cut_ranks(partition : GraphPartition, row : int, columns_to_swap : list[int], row_is_in_partition: bool, columns_are_in_partition: bool, ranks : list[int]) -> None:
    """Finds the cut-ranks for the partitions obtained by swapping a specific row with any of the given columns in the given graph partition.
    
    args:
        - partition: 'GraphPartition' The graph partition.
        - row: 'int' The row to be swapped.
        - columns_to_swap: 'list[int]' The columns to be swapped.
        - row_is_in_partition: 'bool' If true, row will be removed from partition.rows and a swapped column will be added instead.
            If false, row is not among the partition's rows or columns. Thus, it is not in the row base, and no swapped column cannot enter the row base.
        - columns_are_in_partition: 'bool' If true, columns_to_swap are the partition's columns.
            If false, columns_to_swap are disjunct from the partition's rows and columns.
        - ranks: 'list[int]' A list where position [j] will hold the cut-rank after swapping 'row' and 'j'. 
        Only positions where j is a column in the current partition will be affected.
    """

    old_rank = partition.cut_rank
    base_columns = partition.base_columns if columns_are_in_partition else []
    free_columns = partition.free_columns if columns_are_in_partition else columns_to_swap

    if (not partition.base_flag[row]):

        s2 = columns_are_in_partition and any(k2 != row and partition.adj_b_inv_adj[k2][row] == 1 for k2 in partition.free_rows)
        for column in free_columns:
            # row in X^D, column in Y^D
            t2 = row_is_in_partition and any(l2 != column and partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
            if s2:
                if t2:
                    ranks[column] = old_rank + 2
                else:
                    ranks[column] = old_rank + 1
            else:
                if t2:
                    ranks[column] = old_rank + 1
                else:
                    if row_is_in_partition and columns_are_in_partition and partition.adj_b_inv_adj[column][row] == 1:
                        ranks[column] = old_rank + 1
                    else:
                        ranks[column] = old_rank

        for column in base_columns:
            # row in X^D, column in Y^B
            l1 = next((l1 for l1 in partition.free_columns if partition.b_inverse_adj[column][l1] == 1), -1)
            if l1 >= 0:
                if partition.adj_b_inv_adj[column][l1] == 1:
                    t2 = row_is_in_partition and any(l2 != l1 and partition.adj_b_inv_adj[column][l2] != partition.b_inverse_adj[column][l2] for l2 in partition.free_columns)
                else:
                    t2 = row_is_in_partition and any(l2 != l1 and partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                if t2:
                    if s2:
                        ranks[column] = old_rank + 2
                    else:
                        ranks[column] = old_rank + 1
                else:
                    if s2:
                        ranks[column] = old_rank + 1
                    else:
                        if row_is_in_partition and partition.adj_b_inv_adj[column][row] != (partition.b_inverse_adj[column][row] & partition.adj_b_inv_adj[column][l1]):
                            ranks[column] = old_rank + 1
                        else:
                            ranks[column] = old_rank
            else:
                t2 = row_is_in_partition and any(partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                if partition.b_inverse_adj[column][row] == 1:
                    if t2:
                        ranks[column] = old_rank + 1
                    else:
                        ranks[column] = old_rank
                else:
                    if t2:
                        if s2:
                            ranks[column] = old_rank + 1
                        else:
                            ranks[column] = old_rank
                    else:
                        if s2:
                            ranks[column] = old_rank
                        else:
                            if row_is_in_partition and partition.adj_b_inv_adj[column][row] == 1:
                                ranks[column] = old_rank
                            else:
                                ranks[column] = old_rank - 1

    else:

        k1 = next((k1 for k1 in partition.free_rows if partition.adj_b_inverse[k1][row] == 1), -1)
        if k1 >= 0 and partition.adj_b_inv_adj[k1][row] == 1:
            s2 = columns_are_in_partition and any(k2 != k1 and partition.adj_b_inv_adj[k2][row] != partition.adj_b_inverse[k2][row] for k2 in partition.free_rows)
        else:
            s2 = columns_are_in_partition and any(k2 != k1 and partition.adj_b_inv_adj[k2][row] == 1 for k2 in partition.free_rows)
        for column in free_columns:
            # row in X^B, column in Y^D
            if k1 >= 0:
                t2 = any(l2 != column and partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                if s2:
                    if t2:
                        ranks[column] = old_rank + 2
                    else:
                        ranks[column] = old_rank + 1
                else:
                    if t2:
                        ranks[column] = old_rank + 1
                    else:
                        if columns_are_in_partition and partition.adj_b_inv_adj[column][row] != (partition.adj_b_inverse[column][row] & partition.adj_b_inv_adj[k1][row]):
                            ranks[column] = old_rank + 1
                        else:
                            ranks[column] = old_rank
            else:
                if partition.adj_b_inverse[column][row] == 1:
                    if s2:
                        ranks[column] = old_rank + 1
                    else:
                        ranks[column] = old_rank
                else:
                    t2 = any(l2 != column and partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                    if s2:
                        if t2:
                            ranks[column] = old_rank + 1
                        else:
                            ranks[column] = old_rank
                    else:
                        if t2:
                            ranks[column] = old_rank
                        else:
                            if columns_are_in_partition and partition.adj_b_inv_adj[column][row] == 1:
                                ranks[column] = old_rank
                            else:
                                ranks[column] = old_rank - 1
        
        q4_0 = any(partition.adj_b_inv_adj[k][row] == 1 for k in partition.free_rows)
        q4_1 = any(partition.adj_b_inv_adj[k][row] != partition.adj_b_inverse[k][row] for k in partition.free_rows)
        for column in base_columns:
            # row in X^B, column in Y^B
            l1 = next((l1 for l1 in partition.free_columns if partition.b_inverse_adj[column][l1] == 1), -1)
            if (partition.base_inverse[column][row] == 1):

                # Full rank matrix with row and column removed is invertible
                if k1 >= 0 and l1 >= 0:
                    if partition.adj_b_inv_adj[column][l1] == 1:
                        t2 = any(l2 != l1 and partition.adj_b_inv_adj[column][l2] != partition.b_inverse_adj[column][l2] for l2 in partition.free_columns)
                    else:
                        t2 = any(l2 != l1 and partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                    if s2:
                        if t2:
                            ranks[column] = old_rank + 2
                        else:
                            ranks[column] = old_rank + 1
                    else:
                        if t2:
                            ranks[column] = old_rank + 1
                        else:
                            if ((partition.adj_b_inv_adj[k1][row] & partition.adj_b_inv_adj[column][l1]) ^ (partition.adj_b_inv_adj[k1][row] & partition.adj_b_inverse[column][row]) ^ (partition.adj_b_inv_adj[column][l1] & partition.b_inverse_adj[column][row])) != partition.adj_b_inv_adj[column][row]:
                                ranks[column] = old_rank + 1
                            else:
                                ranks[column] = old_rank

                else:
                    q4 = q4_1 if partition.b_inverse_adj[column][row] == 1 else q4_0
                    q5 = any(partition.adj_b_inv_adj[column][l] != (partition.adj_b_inverse[column][row] & partition.b_inverse_adj[column][l]) for l in partition.free_columns)
                    if q4:
                        if q5:
                            ranks[column] = old_rank + 1
                        else:
                            ranks[column] = old_rank
                    else:
                        if q5:
                            ranks[column] = old_rank
                        else:
                            if partition.adj_b_inv_adj[column][row] != (partition.adj_b_inverse[column][row] & partition.b_inverse_adj[column][row]):
                                ranks[column] = old_rank
                            else:
                                ranks[column] = old_rank - 1

            else:

                # Full rank matrix with row and column removed is singular
                if k1 >= 0:

                    if l1 >= 0:

                        # Case k1 >= 0 and l1 >= 0
                        if partition.adj_b_inv_adj[column][l1] == 1:
                            t2 = any(l2 != l1 and partition.adj_b_inv_adj[column][l2] != partition.b_inverse_adj[column][l2] for l2 in partition.free_columns)
                        else:
                            t2 = any(l2 != l1 and partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                        if s2:
                            if t2:
                                ranks[column] = old_rank + 2
                            else:
                                ranks[column] = old_rank + 1
                        else:
                            if t2:
                                ranks[column] = old_rank + 1
                            else:
                                if ((partition.adj_b_inv_adj[k1][row] & partition.adj_b_inverse[column][row]) ^ (partition.adj_b_inv_adj[column][l1] & partition.b_inverse_adj[column][row])) != partition.adj_b_inv_adj[column][row]:
                                    ranks[column] = old_rank + 1
                                else:
                                    ranks[column] = old_rank

                    else:

                        # Case k1 >= 0 and l1 < 0
                        t2 = any(partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                        if t2:
                            if partition.b_inverse_adj[column][row] == 1:
                                ranks[column] = old_rank + 1
                            else:
                                if s2:
                                    ranks[column] = old_rank + 1
                                else:
                                    ranks[column] = old_rank
                        else:
                            if partition.b_inverse_adj[column][row] == 1:
                                ranks[column] = old_rank
                            else:
                                if s2:
                                    ranks[column] = old_rank
                                else:
                                    if (partition.adj_b_inv_adj[k1][row] & partition.adj_b_inverse[column][row]) != partition.adj_b_inv_adj[column][row]:
                                        ranks[column] = old_rank
                                    else:
                                        ranks[column] = old_rank - 1

                else:

                    if l1 >= 0:

                        # Case k1 < 0 and l1 >= 0
                        if s2:
                            if partition.adj_b_inverse[column][row] == 1:
                                ranks[column] = old_rank + 1
                            else:
                                if partition.adj_b_inv_adj[column][l1] == 1:
                                    t2 = any(l2 != l1 and partition.adj_b_inv_adj[column][l2] != partition.b_inverse_adj[column][l2] for l2 in partition.free_columns)
                                else:
                                    t2 = any(l2 != l1 and partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                                if t2:
                                    ranks[column] = old_rank + 1
                                else:
                                    ranks[column] = old_rank
                        else:
                            if partition.adj_b_inverse[column][row] == 1:
                                ranks[column] = old_rank
                            else:
                                if partition.adj_b_inv_adj[column][l1] == 1:
                                    t2 = any(l2 != l1 and partition.adj_b_inv_adj[column][l2] != partition.b_inverse_adj[column][l2] for l2 in partition.free_columns)
                                else:
                                    t2 = any(l2 != l1 and partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                                if t2:
                                    ranks[column] = old_rank
                                else:
                                    if (partition.adj_b_inv_adj[column][l1] & partition.b_inverse_adj[column][row]) != partition.adj_b_inv_adj[column][row]:
                                        ranks[column] = old_rank
                                    else:
                                        ranks[column] = old_rank - 1

                    else:

                        # Case k1 < 0 and l1 < 0
                        if partition.adj_b_inverse[column][row] == 1:
                            if partition.b_inverse_adj[column][row] == 1:
                                ranks[column] = old_rank
                            else:
                                if s2:
                                    ranks[column] = old_rank
                                else:
                                    ranks[column] = old_rank - 1
                        else:
                            if partition.b_inverse_adj[column][row] == 1:
                                t2 = any(partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                                if t2:
                                    ranks[column] = old_rank
                                else:
                                    ranks[column] = old_rank - 1
                            else:
                                t2 = any(partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                                if s2:
                                    if t2:
                                        ranks[column] = old_rank
                                    else:
                                        ranks[column] = old_rank - 1
                                else:
                                    if t2:
                                        ranks[column] = old_rank - 1
                                    else:
                                        if partition.adj_b_inv_adj[column][row] == 1:
                                            ranks[column] = old_rank - 1
                                        else:
                                            ranks[column] = old_rank - 2



def single_swap_cut_rank(partition : GraphPartition, row : int, column : int, row_is_in_partition: bool, column_is_in_partition: bool) -> int:
    """Returns the cut-rank for the partition obtained by swapping the given row and column in the given graph partition.
    
    args:
        - partition: 'GraphPartition' The graph partition.
        - row: 'int' The row to be swapped.
        - column: 'int' The column to be swapped.
        - row_is_in_partition: 'bool' If true, row will be removed from partition.rows and column will be added instead.
            If false, row is not among the partition's rows or columns. Thus, it is not in the row base, and column cannot enter the row base.
        - column_is_in_partition: 'bool' If true, column will be removed from partition.columns and row will be added instead.
            If false, column is not among the partition's rows or columns. Thus, it is not in the column base, and row cannot enter the column base.
    """

    old_rank = partition.cut_rank

    if (not partition.base_flag[column]):

        if (not partition.base_flag[row]):

            # row in X^D, column in Y^D
            s2 = column_is_in_partition and any(k2 != row and partition.adj_b_inv_adj[k2][row] == 1 for k2 in partition.free_rows)
            t2 = row_is_in_partition and any(l2 != column and partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
            if s2:
                if t2:
                    return old_rank + 2
                else:
                    return old_rank + 1
            else:
                if t2:
                    return old_rank + 1
                else:
                    if column_is_in_partition and row_is_in_partition and partition.adj_b_inv_adj[column][row] == 1:
                        return old_rank + 1
                    else:
                        return old_rank

        else:

            # row in X^B, column in Y^D
            k1 = next((k1 for k1 in partition.free_rows if partition.adj_b_inverse[k1][row] == 1), -1)
            if k1 >= 0:
                if partition.adj_b_inv_adj[k1][row] == 1:
                    s2 = column_is_in_partition and any(k2 != k1 and partition.adj_b_inv_adj[k2][row] != partition.adj_b_inverse[k2][row] for k2 in partition.free_rows)
                else:
                    s2 = column_is_in_partition and any(k2 != k1 and partition.adj_b_inv_adj[k2][row] == 1 for k2 in partition.free_rows)
                t2 = any(l2 != column and partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                if s2:
                    if t2:
                        return old_rank + 2
                    else:
                        return old_rank + 1
                else:
                    if t2:
                        return old_rank + 1
                    else:
                        if column_is_in_partition and partition.adj_b_inv_adj[column][row] != (partition.adj_b_inverse[column][row] & partition.adj_b_inv_adj[k1][row]):
                            return old_rank + 1
                        else:
                            return old_rank
            else:
                s2 = column_is_in_partition and any(partition.adj_b_inv_adj[k2][row] == 1 for k2 in partition.free_rows)
                if partition.adj_b_inverse[column][row] == 1:
                    if s2:
                        return old_rank + 1
                    else:
                        return old_rank
                else:
                    t2 = any(l2 != column and partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                    if s2:
                        if t2:
                            return old_rank + 1
                        else:
                            return old_rank
                    else:
                        if t2:
                            return old_rank
                        else:
                            if column_is_in_partition and partition.adj_b_inv_adj[column][row] == 1:
                                return old_rank
                            else:
                                return old_rank - 1

    else:

        if (not partition.base_flag[row]):

            # row in X^D, column in Y^B
            l1 = next((l1 for l1 in partition.free_columns if partition.b_inverse_adj[column][l1] == 1), -1)
            if l1 >= 0:
                if partition.adj_b_inv_adj[column][l1] == 1:
                    t2 = row_is_in_partition and any(l2 != l1 and partition.adj_b_inv_adj[column][l2] != partition.b_inverse_adj[column][l2] for l2 in partition.free_columns)
                else:
                    t2 = row_is_in_partition and any(l2 != l1 and partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                s2 = any(k2 != row and partition.adj_b_inv_adj[k2][row] == 1 for k2 in partition.free_rows)
                if t2:
                    if s2:
                        return old_rank + 2
                    else:
                        return old_rank + 1
                else:
                    if s2:
                        return old_rank + 1
                    else:
                        if row_is_in_partition and partition.adj_b_inv_adj[column][row] != (partition.b_inverse_adj[column][row] & partition.adj_b_inv_adj[column][l1]):
                            return old_rank + 1
                        else:
                            return old_rank
            else:
                t2 = row_is_in_partition and any(partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                if partition.b_inverse_adj[column][row] == 1:
                    if t2:
                        return old_rank + 1
                    else:
                        return old_rank
                else:
                    s2 = any(k2 != row and partition.adj_b_inv_adj[k2][row] == 1 for k2 in partition.free_rows)
                    if t2:
                        if s2:
                            return old_rank + 1
                        else:
                            return old_rank
                    else:
                        if s2:
                            return old_rank
                        else:
                            if row_is_in_partition and partition.adj_b_inv_adj[column][row] == 1:
                                return old_rank
                            else:
                                return old_rank - 1

        else:

            # row in X^B, column in Y^B
            k1 = next((k1 for k1 in partition.free_rows if partition.adj_b_inverse[k1][row] == 1), -1)
            l1 = next((l1 for l1 in partition.free_columns if partition.b_inverse_adj[column][l1] == 1), -1)
            if (partition.base_inverse[column][row] == 1):

                # Full rank matrix with row and column removed is invertible
                if k1 >= 0 and l1 >= 0:
                    if partition.adj_b_inv_adj[k1][row] == 1:
                        s2 = any(k2 != k1 and partition.adj_b_inv_adj[k2][row] != partition.adj_b_inverse[k2][row] for k2 in partition.free_rows)
                    else:
                        s2 = any(k2 != k1 and partition.adj_b_inv_adj[k2][row] == 1 for k2 in partition.free_rows)
                    if partition.adj_b_inv_adj[column][l1] == 1:
                        t2 = any(l2 != l1 and partition.adj_b_inv_adj[column][l2] != partition.b_inverse_adj[column][l2] for l2 in partition.free_columns)
                    else:
                        t2 = any(l2 != l1 and partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                    if s2:
                        if t2:
                            return old_rank + 2
                        else:
                            return old_rank + 1
                    else:
                        if t2:
                            return old_rank + 1
                        else:
                            if ((partition.adj_b_inv_adj[k1][row] & partition.adj_b_inv_adj[column][l1]) ^ (partition.adj_b_inv_adj[k1][row] & partition.adj_b_inverse[column][row]) ^ (partition.adj_b_inv_adj[column][l1] & partition.b_inverse_adj[column][row])) != partition.adj_b_inv_adj[column][row]:
                                return old_rank + 1
                            else:
                                return old_rank

                else:
                    q4 = any(partition.adj_b_inv_adj[k][row] != (partition.adj_b_inverse[k][row] & partition.b_inverse_adj[column][row]) for k in partition.free_rows)
                    q5 = any(partition.adj_b_inv_adj[column][l] != (partition.adj_b_inverse[column][row] & partition.b_inverse_adj[column][l]) for l in partition.free_columns)
                    if q4:
                        if q5:
                            return old_rank + 1
                        else:
                            return old_rank
                    else:
                        if q5:
                            return old_rank
                        else:
                            if partition.adj_b_inv_adj[column][row] != (partition.adj_b_inverse[column][row] & partition.b_inverse_adj[column][row]):
                                return old_rank
                            else:
                                return old_rank - 1

            else:

                # Full rank matrix with row and column removed is singular
                if k1 >= 0:

                    if l1 >= 0:

                        # Case k1 >= 0 and l1 >= 0
                        if partition.adj_b_inv_adj[k1][row] == 1:
                            s2 = any(k2 != k1 and partition.adj_b_inv_adj[k2][row] != partition.adj_b_inverse[k2][row] for k2 in partition.free_rows)
                        else:
                            s2 = any(k2 != k1 and partition.adj_b_inv_adj[k2][row] == 1 for k2 in partition.free_rows)
                        if partition.adj_b_inv_adj[column][l1] == 1:
                            t2 = any(l2 != l1 and partition.adj_b_inv_adj[column][l2] != partition.b_inverse_adj[column][l2] for l2 in partition.free_columns)
                        else:
                            t2 = any(l2 != l1 and partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                        if s2:
                            if t2:
                                return old_rank + 2
                            else:
                                return old_rank + 1
                        else:
                            if t2:
                                return old_rank + 1
                            else:
                                if ((partition.adj_b_inv_adj[k1][row] & partition.adj_b_inverse[column][row]) ^ (partition.adj_b_inv_adj[column][l1] & partition.b_inverse_adj[column][row])) != partition.adj_b_inv_adj[column][row]:
                                    return old_rank + 1
                                else:
                                    return old_rank

                    else:

                        # Case k1 >= 0 and l1 < 0
                        t2 = any(partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                        if t2:
                            if partition.b_inverse_adj[column][row] == 1:
                                return old_rank + 1
                            else:
                                if partition.adj_b_inv_adj[k1][row] == 1:
                                    s2 = any(k2 != k1 and partition.adj_b_inv_adj[k2][row] != partition.adj_b_inverse[k2][row] for k2 in partition.free_rows)
                                else:
                                    s2 = any(k2 != k1 and partition.adj_b_inv_adj[k2][row] == 1 for k2 in partition.free_rows)
                                if s2:
                                    return old_rank + 1
                                else:
                                    return old_rank
                        else:
                            if partition.b_inverse_adj[column][row] == 1:
                                return old_rank
                            else:
                                if partition.adj_b_inv_adj[k1][row] == 1:
                                    s2 = any(k2 != k1 and partition.adj_b_inv_adj[k2][row] != partition.adj_b_inverse[k2][row] for k2 in partition.free_rows)
                                else:
                                    s2 = any(k2 != k1 and partition.adj_b_inv_adj[k2][row] == 1 for k2 in partition.free_rows)
                                if s2:
                                    return old_rank
                                else:
                                    if (partition.adj_b_inv_adj[k1][row] & partition.adj_b_inverse[column][row]) != partition.adj_b_inv_adj[column][row]:
                                        return old_rank
                                    else:
                                        return old_rank - 1

                else:

                    if l1 >= 0:

                        # Case k1 < 0 and l1 >= 0
                        s2 = any(partition.adj_b_inv_adj[k2][row] == 1 for k2 in partition.free_rows)
                        if s2:
                            if partition.adj_b_inverse[column][row] == 1:
                                return old_rank + 1
                            else:
                                if partition.adj_b_inv_adj[column][l1] == 1:
                                    t2 = any(l2 != l1 and partition.adj_b_inv_adj[column][l2] != partition.b_inverse_adj[column][l2] for l2 in partition.free_columns)
                                else:
                                    t2 = any(l2 != l1 and partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                                if t2:
                                    return old_rank + 1
                                else:
                                    return old_rank
                        else:
                            if partition.adj_b_inverse[column][row] == 1:
                                return old_rank
                            else:
                                if partition.adj_b_inv_adj[column][l1] == 1:
                                    t2 = any(l2 != l1 and partition.adj_b_inv_adj[column][l2] != partition.b_inverse_adj[column][l2] for l2 in partition.free_columns)
                                else:
                                    t2 = any(l2 != l1 and partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                                if t2:
                                    return old_rank
                                else:
                                    if (partition.adj_b_inv_adj[column][l1] & partition.b_inverse_adj[column][row]) != partition.adj_b_inv_adj[column][row]:
                                        return old_rank
                                    else:
                                        return old_rank - 1

                    else:

                        # Case k1 < 0 and l1 < 0
                        if partition.adj_b_inverse[column][row] == 1:
                            if partition.b_inverse_adj[column][row] == 1:
                                return old_rank
                            else:
                                s2 = any(partition.adj_b_inv_adj[k2][row] == 1 for k2 in partition.free_rows)
                                if s2:
                                    return old_rank
                                else:
                                    return old_rank - 1
                        else:
                            if partition.b_inverse_adj[column][row] == 1:
                                t2 = any(partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                                if t2:
                                    return old_rank
                                else:
                                    return old_rank - 1
                            else:
                                s2 = any(partition.adj_b_inv_adj[k2][row] == 1 for k2 in partition.free_rows)
                                t2 = any(partition.adj_b_inv_adj[column][l2] == 1 for l2 in partition.free_columns)
                                if s2:
                                    if t2:
                                        return old_rank
                                    else:
                                        return old_rank - 1
                                else:
                                    if t2:
                                        return old_rank - 1
                                    else:
                                        if partition.adj_b_inv_adj[column][row] == 1:
                                            return old_rank - 1
                                        else:
                                            return old_rank - 2
