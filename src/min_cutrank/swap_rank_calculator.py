from min_cutrank.graph_partition import GraphPartition, SubMatrix


def all_swap_cut_ranks(partition : GraphPartition, rows_to_swap: list[int], columns_to_swap: list[int], ranks : list[list[int]]) -> None:
    """Finds the cut-ranks for the partitions obtained by swapping any of the given rows with any of the given columns in the given graph partition.
    
    args:
        - partition: 'GraphPartition' The graph partition.
        - rows_to_swap: 'list[int]' The rows to be swapped.
        - columns_to_swap: 'list[int]' The columns to be swapped.
        - ranks: 'list[list[int]]' A matrix where position [i][j] will hold the cut-rank after swapping node i and j. 
            Only positions where i is in rows_to_swap and j is in columns_to_swap will be affected.
    """
    
    row_subset_index = single([partition.subset_index[row] for row in rows_to_swap], "All rows to swap must be in the same subset or not in any subset")
    column_subset_index = single([partition.subset_index[column] for column in columns_to_swap], "All columns to swap must be in the same subset or not in any subset")
    if (row_subset_index == column_subset_index):
        raise Exception("Rows and columns to swap cannot be in the same subset")

    # Initialize ranks
    for row in rows_to_swap:
        for column in columns_to_swap:
            ranks[row][column] = partition.cut_rank

    # Add contributions from each sub-matrix
    for i in range(len(partition.subsets)):
        if i != row_subset_index and row_subset_index != -1:
            add_all_swap_cut_rank_deltas(partition.matrices[row_subset_index][i], rows_to_swap, columns_to_swap, True, i == column_subset_index, ranks)
        if i != column_subset_index and i != row_subset_index and column_subset_index != -1:
            add_all_swap_cut_rank_deltas(partition.matrices[i][column_subset_index], rows_to_swap, columns_to_swap, False, True, ranks)


def add_all_swap_cut_rank_deltas(matrix : SubMatrix, rows_to_swap: list[int], columns_to_swap: list[int], rows_are_in_matrix: bool, columns_are_in_matrix: bool, ranks : list[list[int]]) -> None:
    """Adds the cange in cut-rank for the sub-matrices obtained by swapping any of the given rows with any of the given columns in the associated graph partition.
    
    args:
        - matrix: 'SubMatrix' The sub-matrix.
        - rows_to_swap: 'list[int]' The rows to be swapped.
        - columns_to_swap: 'list[int]' The columns to be swapped.
        - rows_are_in_matrix: 'bool' If true, each of rows_to_swap will be removed from matrix.rows and each of columns_to_swap will be added instead.
            If false, rows_to_swap are not among the matrix's rows or columns. Thus, they are not in the row base, and columns_to_swap cannot enter the row base.
        - columns_are_in_matrix: 'bool' If true, each of columns_to_swap will be removed from matrix.columns and each of rows_to_swap will be added instead.
            If false, columns_to_swap are not among the matrix's rows or columns. Thus, they are not in the column base, and rows_to_swap cannot enter the column base.
        - ranks: 'list[list[int]]' A matrix where position [i][j] represents the cut-rank after swapping node i and j. 
            Only positions where i is in rows_to_swap and j is in columns_to_swap will be affected.
    """

    nmb_nodes = matrix.partition.graph.nmb_nodes
    base_rows = matrix.base_rows if rows_are_in_matrix else []
    free_rows = matrix.free_rows if rows_are_in_matrix else rows_to_swap    
    base_columns = matrix.base_columns if columns_are_in_matrix else []
    free_columns = matrix.free_columns if columns_are_in_matrix else columns_to_swap
    
    # Preprocessing on rows
    s1_k1 = [-1] * nmb_nodes
    s2 = [False] * nmb_nodes
    q4_952_0 = [False] * nmb_nodes
    q4_952_1 = [False] * nmb_nodes
    for i in matrix.base_rows:
        k1 = next((k1 for k1 in matrix.free_rows if matrix.adj_b_inverse[k1][i] == 1), -1)
        s1_k1[i] = k1
        if columns_are_in_matrix:
            if k1 >= 0 and matrix.adj_b_inv_adj[k1][i] == 1:
                s2[i] = any(k2 != k1 and matrix.adj_b_inv_adj[k2][i] != matrix.adj_b_inverse[k2][i] for k2 in matrix.free_rows)
            else:
                s2[i] = any(matrix.adj_b_inv_adj[k2][i] == 1 for k2 in matrix.free_rows)
            q4_952_0[i] = any(matrix.adj_b_inv_adj[k][i] == 1 for k in matrix.free_rows)
            q4_952_1[i] = any(matrix.adj_b_inv_adj[k][i] != matrix.adj_b_inverse[k][i] for k in matrix.free_rows)
    if columns_are_in_matrix:
        for i in free_rows:
            s2[i] = any(k2 != i and matrix.adj_b_inv_adj[k2][i] == 1 for k2 in matrix.free_rows)

    # Preprocessing on columns
    t1_l1 = [-1] * nmb_nodes
    t2 = [False] * nmb_nodes
    q5_952_0 = [False] * nmb_nodes
    q5_952_1 = [False] * nmb_nodes
    for j in matrix.base_columns:
        l1 = next((l1 for l1 in matrix.free_columns if matrix.b_inverse_adj[j][l1] == 1), -1)
        t1_l1[j] = l1
        if rows_are_in_matrix:
            if l1 >= 0 and matrix.adj_b_inv_adj[j][l1] == 1:
                t2[j] = any(l2 != l1 and matrix.adj_b_inv_adj[j][l2] != matrix.b_inverse_adj[j][l2] for l2 in matrix.free_columns)
            else:
                t2[j] = any(matrix.adj_b_inv_adj[j][l2] == 1 for l2 in matrix.free_columns)
            q5_952_0[j] = any(matrix.adj_b_inv_adj[j][l] == 1 for l in matrix.free_columns)
            q5_952_1[j] = any(matrix.adj_b_inv_adj[j][l] != matrix.b_inverse_adj[j][l] for l in matrix.free_columns)
    if rows_are_in_matrix:
        for j in free_columns:
            t2[j] = any(l2 != j and matrix.adj_b_inv_adj[j][l2] == 1 for l2 in matrix.free_columns)

    # Ranks for i in X^D and j in Y^D
    for i in free_rows:
        for j in free_columns:
            if s2[i]:
                if t2[j]:
                    ranks[i][j] += 2
                else:
                    ranks[i][j] += 1
            else:
                if t2[j]:
                    ranks[i][j] += 1
                else:
                    if rows_are_in_matrix and columns_are_in_matrix and matrix.adj_b_inv_adj[j][i] == 1:
                        ranks[i][j] += 1
                    else:
                        ranks[i][j] += 0

    # Ranks for i in X^B and j in Y^D
    for i in base_rows:
        k1 = s1_k1[i]
        for j in free_columns:
            if k1 >= 0:
                if s2[i]:
                    if t2[j]:
                        ranks[i][j] += 2
                    else:
                        ranks[i][j] += 1
                else:
                    if t2[j]:
                        ranks[i][j] += 1
                    else:
                        if columns_are_in_matrix and matrix.adj_b_inv_adj[j][i] != (matrix.adj_b_inverse[j][i] & matrix.adj_b_inv_adj[k1][i]):
                            ranks[i][j] += 1
                        else:
                            ranks[i][j] += 0
            else:
                if matrix.adj_b_inverse[j][i] == 1:
                    if s2[i]:
                        ranks[i][j] += 1
                    else:
                        ranks[i][j] += 0
                else:
                    if s2[i]:
                        if t2[j]:
                            ranks[i][j] += 1
                        else:
                            ranks[i][j] += 0
                    else:
                        if t2[j]:
                            ranks[i][j] += 0
                        else:
                            if columns_are_in_matrix and matrix.adj_b_inv_adj[j][i] == 1:
                                ranks[i][j] += 0
                            else:
                                ranks[i][j] += -1

    # Ranks for i in X^D and j in Y^B
    for i in free_rows:
        for j in base_columns:
            l1 = t1_l1[j]
            if l1 >= 0:
                if t2[j]:
                    if s2[i]:
                        ranks[i][j] += 2
                    else:
                        ranks[i][j] += 1
                else:
                    if s2[i]:
                        ranks[i][j] += 1
                    else:
                        if rows_are_in_matrix and columns_are_in_matrix and matrix.adj_b_inv_adj[j][i] != (matrix.b_inverse_adj[j][i] & matrix.adj_b_inv_adj[j][l1]):
                            ranks[i][j] += 1
                        else:
                            ranks[i][j] += 0
            else:
                if matrix.b_inverse_adj[j][i] == 1:
                    if t2[j]:
                        ranks[i][j] += 1
                    else:
                        ranks[i][j] += 0
                else:
                    if t2[j]:
                        if s2[i]:
                            ranks[i][j] += 1
                        else:
                            ranks[i][j] += 0
                    else:
                        if s2[i]:
                            ranks[i][j] += 0
                        else:
                            if rows_are_in_matrix and matrix.adj_b_inv_adj[j][i] == 1:
                                ranks[i][j] += 0
                            else:
                                ranks[i][j] += -1

    # Ranks for i in X^B and j in Y^B
    for i in base_rows:
        k1 = s1_k1[i]
        for j in base_columns:
            l1 = t1_l1[j]

            if matrix.base_inverse[j][i] == 1:
                # Case 6

                if k1 >= 0 and l1 >= 0:
                    # Case 6.1
                    if s2[i]:
                        if t2[j]:
                            ranks[i][j] += 2
                        else:
                            ranks[i][j] += 1
                    else:
                        if t2[j]:
                            ranks[i][j] += 1
                        else:
                            if ((matrix.adj_b_inv_adj[k1][i] & matrix.adj_b_inv_adj[j][l1]) ^ (matrix.adj_b_inv_adj[k1][i] & matrix.adj_b_inverse[j][i]) ^ (matrix.adj_b_inv_adj[j][l1] & matrix.b_inverse_adj[j][i])) != matrix.adj_b_inv_adj[j][i]:
                                ranks[i][j] += 1
                            else:
                                ranks[i][j] += 0

                else:
                    # Case 6.2
                    q4 = q4_952_1[i] if matrix.b_inverse_adj[j][i] == 1 else q4_952_0[i]
                    q5 = q5_952_1[j] if matrix.adj_b_inverse[j][i] == 1 else q5_952_0[j]
                    if q4:
                        if q5:
                            ranks[i][j] += 1
                        else:
                            ranks[i][j] += 0
                    else:
                        if q5:
                            ranks[i][j] += 0
                        else:
                            if matrix.adj_b_inv_adj[j][i] != (matrix.adj_b_inverse[j][i] & matrix.b_inverse_adj[j][i]):
                                ranks[i][j] += 0
                            else:
                                ranks[i][j] += -1

            else:
                # Case 7

                    if k1 >= 0:

                        if l1 >= 0:
                            # Case 7.1
                            if s2[i]:
                                if t2[j]:
                                    ranks[i][j] += 2
                                else:
                                    ranks[i][j] += 1
                            else:
                                if t2[j]:
                                    ranks[i][j] += 1
                                else:
                                    if ((matrix.adj_b_inv_adj[k1][i] & matrix.adj_b_inverse[j][i]) ^ (matrix.adj_b_inv_adj[j][l1] & matrix.b_inverse_adj[j][i])) != matrix.adj_b_inv_adj[j][i]:
                                        ranks[i][j] += 1
                                    else:
                                        ranks[i][j] += 0

                        else:
                            # Case 7.2
                            if t2[j]:
                                if matrix.b_inverse_adj[j][i] == 1:
                                    ranks[i][j] += 1
                                else:
                                    if s2[i]:
                                        ranks[i][j] += 1
                                    else:
                                        ranks[i][j] += 0
                            else:
                                if matrix.b_inverse_adj[j][i] == 1:
                                    ranks[i][j] += 0
                                else:
                                    if s2[i]:
                                        ranks[i][j] += 0
                                    else:
                                        if (matrix.adj_b_inv_adj[k1][i] & matrix.adj_b_inverse[j][i]) != matrix.adj_b_inv_adj[j][i]:
                                            ranks[i][j] += 0
                                        else:
                                            ranks[i][j] += -1

                    else:

                        if l1 >= 0:
                            # Case 7.3
                            if s2[i]:
                                if matrix.adj_b_inverse[j][i] == 1:
                                    ranks[i][j] += 1
                                else:
                                    if t2[j]:
                                        ranks[i][j] += 1
                                    else:
                                        ranks[i][j] += 0
                            else:
                                if matrix.adj_b_inverse[j][i] == 1:
                                    ranks[i][j] += 0
                                else:
                                    if t2[j]:
                                        ranks[i][j] += 0
                                    else:
                                        if (matrix.adj_b_inv_adj[j][l1] & matrix.b_inverse_adj[j][i]) != matrix.adj_b_inv_adj[j][i]:
                                            ranks[i][j] += 0
                                        else:
                                            ranks[i][j] += -1

                        else:
                            # Case 7.4
                            if matrix.adj_b_inverse[j][i] == 1:
                                if matrix.b_inverse_adj[j][i] == 1:
                                    ranks[i][j] += 0
                                else:
                                    if s2[i]:
                                        ranks[i][j] += 0
                                    else:
                                        ranks[i][j] += -1
                            else:
                                if matrix.b_inverse_adj[j][i] == 1:
                                    if t2[j]:
                                        ranks[i][j] += 0
                                    else:
                                        ranks[i][j] += -1
                                else:
                                    if s2[i]:
                                        if t2[j]:
                                            ranks[i][j] += 0
                                        else:
                                            ranks[i][j] += -1
                                    else:
                                        if t2[j]:
                                            ranks[i][j] += -1
                                        else:
                                            if matrix.adj_b_inv_adj[j][i] == 1:
                                                ranks[i][j] += -1
                                            else:
                                                ranks[i][j] += -2


def row_swap_cut_ranks(partition : GraphPartition, row : int, columns_to_swap : list[int], ranks : list[int]) -> None:
    """Finds the cut-ranks for the partitions obtained by swapping the given row with any of the given columns in the given graph partition.
    
    args:
        - partition: 'GraphPartition' The graph partition.
        - row: 'int' The row to be swapped.
        - columns_to_swap: 'list[int]' The columns to be swapped.
        - ranks: 'list[int]' A list where position [j] represents the cut-rank after swapping 'row' and 'j'. 
            Only positions where j is among columns_to_swap will be affected.
    """

    row_subset_index = partition.subset_index[row]
    column_subset_index = single([partition.subset_index[column] for column in columns_to_swap], "All columns to swap must be in the same subset or not in any subset")
    if (row_subset_index == column_subset_index):
        raise Exception("Row and columns to swap cannot be in the same subset")

    # Initialize ranks
    for column in columns_to_swap:
        ranks[column] = partition.cut_rank

    # Add contributions from each sub-matrix
    for i in range(len(partition.subsets)):
        if i != row_subset_index and row_subset_index != -1:
            add_row_swap_cut_rank_deltas(partition.matrices[row_subset_index][i], row, columns_to_swap, True, i == column_subset_index, ranks)
        if i != column_subset_index and i != row_subset_index and column_subset_index != -1:
            add_row_swap_cut_rank_deltas(partition.matrices[i][column_subset_index], row, columns_to_swap, False, True, ranks)


def add_row_swap_cut_rank_deltas(matrix : SubMatrix, row : int, columns_to_swap : list[int], row_is_in_matrix: bool, columns_are_in_matrix: bool, ranks : list[int]) -> None:
    """Adds the change in rank for the sub-matrices obtained by swapping the given row with any of the given columns in the associated graph partition.
    
    args:
        - matrix: 'SubMatrix' The sub-matrix.
        - row: 'int' The row to be swapped.
        - columns_to_swap: 'list[int]' The columns to be swapped.
        - row_is_in_matrix: 'bool' If true, row will be removed from matrix.rows and each of columns_to_swap will be added instead.
            If false, row is not among the matrix's rows or columns. Thus, it is not in the row base, and columns_to_swap cannot enter the row base.
        - columns_are_in_matrix: 'bool' If true, each of columns_to_swap will be removed from matrix.columns and row will be added instead.
            If false, columns_to_swap are not among the matrix's rows or columns. Thus, they are not in the column base, and row cannot enter the column base.
        - ranks: 'list[int]' A list where position [j] represents the cut-rank after swapping 'row' and 'j'. 
            Only positions where j is among columns_to_swap will be affected.
    """

    base_columns = matrix.base_columns if columns_are_in_matrix else []
    free_columns = matrix.free_columns if columns_are_in_matrix else columns_to_swap

    if (not matrix.base_flag[row]):

        s2 = columns_are_in_matrix and any(k2 != row and matrix.adj_b_inv_adj[k2][row] == 1 for k2 in matrix.free_rows)
        for column in free_columns:
            # row in X^D, column in Y^D
            t2 = row_is_in_matrix and any(l2 != column and matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
            if s2:
                if t2:
                    ranks[column] += 2
                else:
                    ranks[column] += 1
            else:
                if t2:
                    ranks[column] += 1
                else:
                    if row_is_in_matrix and columns_are_in_matrix and matrix.adj_b_inv_adj[column][row] == 1:
                        ranks[column] += 1
                    else:
                        ranks[column] += 0

        for column in base_columns:
            # row in X^D, column in Y^B
            l1 = next((l1 for l1 in matrix.free_columns if matrix.b_inverse_adj[column][l1] == 1), -1)
            if l1 >= 0:
                if matrix.adj_b_inv_adj[column][l1] == 1:
                    t2 = row_is_in_matrix and any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] != matrix.b_inverse_adj[column][l2] for l2 in matrix.free_columns)
                else:
                    t2 = row_is_in_matrix and any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                if t2:
                    if s2:
                        ranks[column] += 2
                    else:
                        ranks[column] += 1
                else:
                    if s2:
                        ranks[column] += 1
                    else:
                        if row_is_in_matrix and matrix.adj_b_inv_adj[column][row] != (matrix.b_inverse_adj[column][row] & matrix.adj_b_inv_adj[column][l1]):
                            ranks[column] += 1
                        else:
                            ranks[column] += 0
            else:
                t2 = row_is_in_matrix and any(matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                if matrix.b_inverse_adj[column][row] == 1:
                    if t2:
                        ranks[column] += 1
                    else:
                        ranks[column] += 0
                else:
                    if t2:
                        if s2:
                            ranks[column] += 1
                        else:
                            ranks[column] += 0
                    else:
                        if s2:
                            ranks[column] += 0
                        else:
                            if row_is_in_matrix and matrix.adj_b_inv_adj[column][row] == 1:
                                ranks[column] += 0
                            else:
                                ranks[column] += -1

    else:

        k1 = next((k1 for k1 in matrix.free_rows if matrix.adj_b_inverse[k1][row] == 1), -1)
        if k1 >= 0 and matrix.adj_b_inv_adj[k1][row] == 1:
            s2 = columns_are_in_matrix and any(k2 != k1 and matrix.adj_b_inv_adj[k2][row] != matrix.adj_b_inverse[k2][row] for k2 in matrix.free_rows)
        else:
            s2 = columns_are_in_matrix and any(k2 != k1 and matrix.adj_b_inv_adj[k2][row] == 1 for k2 in matrix.free_rows)
        for column in free_columns:
            # row in X^B, column in Y^D
            if k1 >= 0:
                t2 = any(l2 != column and matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                if s2:
                    if t2:
                        ranks[column] += 2
                    else:
                        ranks[column] += 1
                else:
                    if t2:
                        ranks[column] += 1
                    else:
                        if columns_are_in_matrix and matrix.adj_b_inv_adj[column][row] != (matrix.adj_b_inverse[column][row] & matrix.adj_b_inv_adj[k1][row]):
                            ranks[column] += 1
                        else:
                            ranks[column] += 0
            else:
                if matrix.adj_b_inverse[column][row] == 1:
                    if s2:
                        ranks[column] += 1
                    else:
                        ranks[column] += 0
                else:
                    t2 = any(l2 != column and matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                    if s2:
                        if t2:
                            ranks[column] += 1
                        else:
                            ranks[column] += 0
                    else:
                        if t2:
                            ranks[column] += 0
                        else:
                            if columns_are_in_matrix and matrix.adj_b_inv_adj[column][row] == 1:
                                ranks[column] += 0
                            else:
                                ranks[column] += -1
        
        q4_0 = any(matrix.adj_b_inv_adj[k][row] == 1 for k in matrix.free_rows)
        q4_1 = any(matrix.adj_b_inv_adj[k][row] != matrix.adj_b_inverse[k][row] for k in matrix.free_rows)
        for column in base_columns:
            # row in X^B, column in Y^B
            l1 = next((l1 for l1 in matrix.free_columns if matrix.b_inverse_adj[column][l1] == 1), -1)
            if (matrix.base_inverse[column][row] == 1):

                # Full rank matrix with row and column removed is invertible
                if k1 >= 0 and l1 >= 0:
                    if matrix.adj_b_inv_adj[column][l1] == 1:
                        t2 = any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] != matrix.b_inverse_adj[column][l2] for l2 in matrix.free_columns)
                    else:
                        t2 = any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                    if s2:
                        if t2:
                            ranks[column] += 2
                        else:
                            ranks[column] += 1
                    else:
                        if t2:
                            ranks[column] += 1
                        else:
                            if ((matrix.adj_b_inv_adj[k1][row] & matrix.adj_b_inv_adj[column][l1]) ^ (matrix.adj_b_inv_adj[k1][row] & matrix.adj_b_inverse[column][row]) ^ (matrix.adj_b_inv_adj[column][l1] & matrix.b_inverse_adj[column][row])) != matrix.adj_b_inv_adj[column][row]:
                                ranks[column] += 1
                            else:
                                ranks[column] += 0

                else:
                    q4 = q4_1 if matrix.b_inverse_adj[column][row] == 1 else q4_0
                    q5 = any(matrix.adj_b_inv_adj[column][l] != (matrix.adj_b_inverse[column][row] & matrix.b_inverse_adj[column][l]) for l in matrix.free_columns)
                    if q4:
                        if q5:
                            ranks[column] += 1
                        else:
                            ranks[column] += 0
                    else:
                        if q5:
                            ranks[column] += 0
                        else:
                            if matrix.adj_b_inv_adj[column][row] != (matrix.adj_b_inverse[column][row] & matrix.b_inverse_adj[column][row]):
                                ranks[column] += 0
                            else:
                                ranks[column] += -1

            else:

                # Full rank matrix with row and column removed is singular
                if k1 >= 0:

                    if l1 >= 0:

                        # Case k1 >= 0 and l1 >= 0
                        if matrix.adj_b_inv_adj[column][l1] == 1:
                            t2 = any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] != matrix.b_inverse_adj[column][l2] for l2 in matrix.free_columns)
                        else:
                            t2 = any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                        if s2:
                            if t2:
                                ranks[column] += 2
                            else:
                                ranks[column] += 1
                        else:
                            if t2:
                                ranks[column] += 1
                            else:
                                if ((matrix.adj_b_inv_adj[k1][row] & matrix.adj_b_inverse[column][row]) ^ (matrix.adj_b_inv_adj[column][l1] & matrix.b_inverse_adj[column][row])) != matrix.adj_b_inv_adj[column][row]:
                                    ranks[column] += 1
                                else:
                                    ranks[column] += 0

                    else:

                        # Case k1 >= 0 and l1 < 0
                        t2 = any(matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                        if t2:
                            if matrix.b_inverse_adj[column][row] == 1:
                                ranks[column] += 1
                            else:
                                if s2:
                                    ranks[column] += 1
                                else:
                                    ranks[column] += 0
                        else:
                            if matrix.b_inverse_adj[column][row] == 1:
                                ranks[column] += 0
                            else:
                                if s2:
                                    ranks[column] += 0
                                else:
                                    if (matrix.adj_b_inv_adj[k1][row] & matrix.adj_b_inverse[column][row]) != matrix.adj_b_inv_adj[column][row]:
                                        ranks[column] += 0
                                    else:
                                        ranks[column] += -1

                else:

                    if l1 >= 0:

                        # Case k1 < 0 and l1 >= 0
                        if s2:
                            if matrix.adj_b_inverse[column][row] == 1:
                                ranks[column] += 1
                            else:
                                if matrix.adj_b_inv_adj[column][l1] == 1:
                                    t2 = any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] != matrix.b_inverse_adj[column][l2] for l2 in matrix.free_columns)
                                else:
                                    t2 = any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                                if t2:
                                    ranks[column] += 1
                                else:
                                    ranks[column] += 0
                        else:
                            if matrix.adj_b_inverse[column][row] == 1:
                                ranks[column] += 0
                            else:
                                if matrix.adj_b_inv_adj[column][l1] == 1:
                                    t2 = any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] != matrix.b_inverse_adj[column][l2] for l2 in matrix.free_columns)
                                else:
                                    t2 = any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                                if t2:
                                    ranks[column] += 0
                                else:
                                    if (matrix.adj_b_inv_adj[column][l1] & matrix.b_inverse_adj[column][row]) != matrix.adj_b_inv_adj[column][row]:
                                        ranks[column] += 0
                                    else:
                                        ranks[column] += -1

                    else:

                        # Case k1 < 0 and l1 < 0
                        if matrix.adj_b_inverse[column][row] == 1:
                            if matrix.b_inverse_adj[column][row] == 1:
                                ranks[column] += 0
                            else:
                                if s2:
                                    ranks[column] += 0
                                else:
                                    ranks[column] += -1
                        else:
                            if matrix.b_inverse_adj[column][row] == 1:
                                t2 = any(matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                                if t2:
                                    ranks[column] += 0
                                else:
                                    ranks[column] += -1
                            else:
                                t2 = any(matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                                if s2:
                                    if t2:
                                        ranks[column] += 0
                                    else:
                                        ranks[column] += -1
                                else:
                                    if t2:
                                        ranks[column] += -1
                                    else:
                                        if matrix.adj_b_inv_adj[column][row] == 1:
                                            ranks[column] += -1
                                        else:
                                            ranks[column] += -2


def single_swap_cut_rank_delta(partition : GraphPartition, row : int, column : int) -> int:
    """Returns the change in cut-rank for the partition obtained by swapping the given row and column in the given graph partition.
    
    args:
        - partition: 'GraphPartition' The graph partition.
        - row: 'int' The row to be swapped.
        - column: 'int' The column to be swapped.
    """

    row_subset_index = partition.subset_index[row]
    column_subset_index = partition.subset_index[column]
    if (row_subset_index == column_subset_index):
        raise Exception("Rows and columns to swap cannot be in the same subset")

    # Initialize rank delta
    rank_delta = 0

    # Add contributions from each sub-matrix
    for i in range(len(partition.subsets)):
        if i != row_subset_index and row_subset_index != -1:
            rank_delta += single_swap_rank_delta(partition.matrices[row_subset_index][i], row, column, True, i == column_subset_index)
        if i != column_subset_index and i != row_subset_index and column_subset_index != -1:
            rank_delta += single_swap_rank_delta(partition.matrices[i][column_subset_index], row, column, False, True)

    return rank_delta


def single_swap_rank_delta(matrix : SubMatrix, row : int, column : int, row_is_in_matrix: bool, column_is_in_matrix: bool) -> int:
    """Returns the change in rank for the sub-matrix obtained by swapping the given row and column in the associated graph partition.
    
    args:
        - matrix: 'SubMatrix' The sub-matrix.
        - row: 'int' The row to be swapped.
        - column: 'int' The column to be swapped.
        - row_is_in_matrix: 'bool' If true, row will be removed from matrix.rows and column will be added instead.
            If false, row is not among the matrix's rows or columns. Thus, it is not in the row base, and column cannot enter the row base.
        - column_is_in_matrix: 'bool' If true, column will be removed from matrix.columns and row will be added instead.
            If false, column is not among the matrix's rows or columns. Thus, it is not in the column base, and row cannot enter the column base.
    """

    if (not matrix.base_flag[column]):

        if (not matrix.base_flag[row]):

            # row in X^D, column in Y^D
            s2 = column_is_in_matrix and any(k2 != row and matrix.adj_b_inv_adj[k2][row] == 1 for k2 in matrix.free_rows)
            t2 = row_is_in_matrix and any(l2 != column and matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
            if s2:
                if t2:
                    return 2
                else:
                    return 1
            else:
                if t2:
                    return 1
                else:
                    if column_is_in_matrix and row_is_in_matrix and matrix.adj_b_inv_adj[column][row] == 1:
                        return 1
                    else:
                        return 0

        else:

            # row in X^B, column in Y^D
            k1 = next((k1 for k1 in matrix.free_rows if matrix.adj_b_inverse[k1][row] == 1), -1)
            if k1 >= 0:
                if matrix.adj_b_inv_adj[k1][row] == 1:
                    s2 = column_is_in_matrix and any(k2 != k1 and matrix.adj_b_inv_adj[k2][row] != matrix.adj_b_inverse[k2][row] for k2 in matrix.free_rows)
                else:
                    s2 = column_is_in_matrix and any(k2 != k1 and matrix.adj_b_inv_adj[k2][row] == 1 for k2 in matrix.free_rows)
                t2 = any(l2 != column and matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                if s2:
                    if t2:
                        return 2
                    else:
                        return 1
                else:
                    if t2:
                        return 1
                    else:
                        if column_is_in_matrix and matrix.adj_b_inv_adj[column][row] != (matrix.adj_b_inverse[column][row] & matrix.adj_b_inv_adj[k1][row]):
                            return 1
                        else:
                            return 0
            else:
                s2 = column_is_in_matrix and any(matrix.adj_b_inv_adj[k2][row] == 1 for k2 in matrix.free_rows)
                if matrix.adj_b_inverse[column][row] == 1:
                    if s2:
                        return 1
                    else:
                        return 0
                else:
                    t2 = any(l2 != column and matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                    if s2:
                        if t2:
                            return 1
                        else:
                            return 0
                    else:
                        if t2:
                            return 0
                        else:
                            if column_is_in_matrix and matrix.adj_b_inv_adj[column][row] == 1:
                                return 0
                            else:
                                return -1

    else:

        if (not matrix.base_flag[row]):

            # row in X^D, column in Y^B
            l1 = next((l1 for l1 in matrix.free_columns if matrix.b_inverse_adj[column][l1] == 1), -1)
            if l1 >= 0:
                if matrix.adj_b_inv_adj[column][l1] == 1:
                    t2 = row_is_in_matrix and any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] != matrix.b_inverse_adj[column][l2] for l2 in matrix.free_columns)
                else:
                    t2 = row_is_in_matrix and any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                s2 = any(k2 != row and matrix.adj_b_inv_adj[k2][row] == 1 for k2 in matrix.free_rows)
                if t2:
                    if s2:
                        return 2
                    else:
                        return 1
                else:
                    if s2:
                        return 1
                    else:
                        if row_is_in_matrix and matrix.adj_b_inv_adj[column][row] != (matrix.b_inverse_adj[column][row] & matrix.adj_b_inv_adj[column][l1]):
                            return 1
                        else:
                            return 0
            else:
                t2 = row_is_in_matrix and any(matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                if matrix.b_inverse_adj[column][row] == 1:
                    if t2:
                        return 1
                    else:
                        return 0
                else:
                    s2 = any(k2 != row and matrix.adj_b_inv_adj[k2][row] == 1 for k2 in matrix.free_rows)
                    if t2:
                        if s2:
                            return 1
                        else:
                            return 0
                    else:
                        if s2:
                            return 0
                        else:
                            if row_is_in_matrix and matrix.adj_b_inv_adj[column][row] == 1:
                                return 0
                            else:
                                return -1

        else:

            # row in X^B, column in Y^B
            k1 = next((k1 for k1 in matrix.free_rows if matrix.adj_b_inverse[k1][row] == 1), -1)
            l1 = next((l1 for l1 in matrix.free_columns if matrix.b_inverse_adj[column][l1] == 1), -1)
            if (matrix.base_inverse[column][row] == 1):

                # Full rank matrix with row and column removed is invertible
                if k1 >= 0 and l1 >= 0:
                    if matrix.adj_b_inv_adj[k1][row] == 1:
                        s2 = any(k2 != k1 and matrix.adj_b_inv_adj[k2][row] != matrix.adj_b_inverse[k2][row] for k2 in matrix.free_rows)
                    else:
                        s2 = any(k2 != k1 and matrix.adj_b_inv_adj[k2][row] == 1 for k2 in matrix.free_rows)
                    if matrix.adj_b_inv_adj[column][l1] == 1:
                        t2 = any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] != matrix.b_inverse_adj[column][l2] for l2 in matrix.free_columns)
                    else:
                        t2 = any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                    if s2:
                        if t2:
                            return 2
                        else:
                            return 1
                    else:
                        if t2:
                            return 1
                        else:
                            if ((matrix.adj_b_inv_adj[k1][row] & matrix.adj_b_inv_adj[column][l1]) ^ (matrix.adj_b_inv_adj[k1][row] & matrix.adj_b_inverse[column][row]) ^ (matrix.adj_b_inv_adj[column][l1] & matrix.b_inverse_adj[column][row])) != matrix.adj_b_inv_adj[column][row]:
                                return 1
                            else:
                                return 0

                else:
                    q4 = any(matrix.adj_b_inv_adj[k][row] != (matrix.adj_b_inverse[k][row] & matrix.b_inverse_adj[column][row]) for k in matrix.free_rows)
                    q5 = any(matrix.adj_b_inv_adj[column][l] != (matrix.adj_b_inverse[column][row] & matrix.b_inverse_adj[column][l]) for l in matrix.free_columns)
                    if q4:
                        if q5:
                            return 1
                        else:
                            return 0
                    else:
                        if q5:
                            return 0
                        else:
                            if matrix.adj_b_inv_adj[column][row] != (matrix.adj_b_inverse[column][row] & matrix.b_inverse_adj[column][row]):
                                return 0
                            else:
                                return -1

            else:

                # Full rank matrix with row and column removed is singular
                if k1 >= 0:

                    if l1 >= 0:

                        # Case k1 >= 0 and l1 >= 0
                        if matrix.adj_b_inv_adj[k1][row] == 1:
                            s2 = any(k2 != k1 and matrix.adj_b_inv_adj[k2][row] != matrix.adj_b_inverse[k2][row] for k2 in matrix.free_rows)
                        else:
                            s2 = any(k2 != k1 and matrix.adj_b_inv_adj[k2][row] == 1 for k2 in matrix.free_rows)
                        if matrix.adj_b_inv_adj[column][l1] == 1:
                            t2 = any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] != matrix.b_inverse_adj[column][l2] for l2 in matrix.free_columns)
                        else:
                            t2 = any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                        if s2:
                            if t2:
                                return 2
                            else:
                                return 1
                        else:
                            if t2:
                                return 1
                            else:
                                if ((matrix.adj_b_inv_adj[k1][row] & matrix.adj_b_inverse[column][row]) ^ (matrix.adj_b_inv_adj[column][l1] & matrix.b_inverse_adj[column][row])) != matrix.adj_b_inv_adj[column][row]:
                                    return 1
                                else:
                                    return 0

                    else:

                        # Case k1 >= 0 and l1 < 0
                        t2 = any(matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                        if t2:
                            if matrix.b_inverse_adj[column][row] == 1:
                                return 1
                            else:
                                if matrix.adj_b_inv_adj[k1][row] == 1:
                                    s2 = any(k2 != k1 and matrix.adj_b_inv_adj[k2][row] != matrix.adj_b_inverse[k2][row] for k2 in matrix.free_rows)
                                else:
                                    s2 = any(k2 != k1 and matrix.adj_b_inv_adj[k2][row] == 1 for k2 in matrix.free_rows)
                                if s2:
                                    return 1
                                else:
                                    return 0
                        else:
                            if matrix.b_inverse_adj[column][row] == 1:
                                return 0
                            else:
                                if matrix.adj_b_inv_adj[k1][row] == 1:
                                    s2 = any(k2 != k1 and matrix.adj_b_inv_adj[k2][row] != matrix.adj_b_inverse[k2][row] for k2 in matrix.free_rows)
                                else:
                                    s2 = any(k2 != k1 and matrix.adj_b_inv_adj[k2][row] == 1 for k2 in matrix.free_rows)
                                if s2:
                                    return 0
                                else:
                                    if (matrix.adj_b_inv_adj[k1][row] & matrix.adj_b_inverse[column][row]) != matrix.adj_b_inv_adj[column][row]:
                                        return 0
                                    else:
                                        return -1

                else:

                    if l1 >= 0:

                        # Case k1 < 0 and l1 >= 0
                        s2 = any(matrix.adj_b_inv_adj[k2][row] == 1 for k2 in matrix.free_rows)
                        if s2:
                            if matrix.adj_b_inverse[column][row] == 1:
                                return 1
                            else:
                                if matrix.adj_b_inv_adj[column][l1] == 1:
                                    t2 = any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] != matrix.b_inverse_adj[column][l2] for l2 in matrix.free_columns)
                                else:
                                    t2 = any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                                if t2:
                                    return 1
                                else:
                                    return 0
                        else:
                            if matrix.adj_b_inverse[column][row] == 1:
                                return 0
                            else:
                                if matrix.adj_b_inv_adj[column][l1] == 1:
                                    t2 = any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] != matrix.b_inverse_adj[column][l2] for l2 in matrix.free_columns)
                                else:
                                    t2 = any(l2 != l1 and matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                                if t2:
                                    return 0
                                else:
                                    if (matrix.adj_b_inv_adj[column][l1] & matrix.b_inverse_adj[column][row]) != matrix.adj_b_inv_adj[column][row]:
                                        return 0
                                    else:
                                        return -1

                    else:

                        # Case k1 < 0 and l1 < 0
                        if matrix.adj_b_inverse[column][row] == 1:
                            if matrix.b_inverse_adj[column][row] == 1:
                                return 0
                            else:
                                s2 = any(matrix.adj_b_inv_adj[k2][row] == 1 for k2 in matrix.free_rows)
                                if s2:
                                    return 0
                                else:
                                    return -1
                        else:
                            if matrix.b_inverse_adj[column][row] == 1:
                                t2 = any(matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                                if t2:
                                    return 0
                                else:
                                    return -1
                            else:
                                s2 = any(matrix.adj_b_inv_adj[k2][row] == 1 for k2 in matrix.free_rows)
                                t2 = any(matrix.adj_b_inv_adj[column][l2] == 1 for l2 in matrix.free_columns)
                                if s2:
                                    if t2:
                                        return 0
                                    else:
                                        return -1
                                else:
                                    if t2:
                                        return -1
                                    else:
                                        if matrix.adj_b_inv_adj[column][row] == 1:
                                            return -1
                                        else:
                                            return -2


def single(items: list[int], error_message: str) -> int:
    """Returns the single integer value in the list. Raises an exception if the list contains different integers."""
    value = items[0]
    for item in items:
        if item != value:
            raise Exception(error_message)
    return value
