from min_cutrank.graph import Graph
from min_cutrank.matrix_tools import insert_zero_matrix, copy_matrix, rank_matrix_positions, matrix_inverse, add_product_matrix


class SubMatrix:
    """A matrix defined by a subset of rows and and a subset of columns in a larger matrix.
    No row and column in the sub-matrix share the same index.
    """

    partition: 'GraphPartition'
    """The graph partition this sub-matrix belongs to."""

    rows : list[int]
    """The row indices of the sub-matrix in the larger matrix."""

    columns : list[int]
    """The column indices of the sub-matrix in the larger matrix."""

    base_flag : list[bool]
    """Flag telling if an index represents a row or column in the selected invertible sub-sub-matrix of this sub-matrix."""

    rank : int
    """The rank of the sub-matrix."""

    base_rows : list[int]
    """The rows that are part of the invertible sub-sub-matrix of this sub-matrix. 
    Same as rows where base_flag[n] is True. Length should equal rank."""

    base_columns : list[int]
    """The columns that are part of the invertible sub-sub-matrix of this sub-matrix.
    Same as columns where base_flag[n] is True. Length should equal rank."""

    free_rows : list[int]
    """The rows in this sub-matric that are *not* part of the invertible sub-sub-matrix of this sub-matrix.
    Same as rows where base_flag[n] is False."""

    free_columns : list[int]
    """The columns in this sub-matrix that are *not* part of the invertible sub-sub-matrix of this sub-matrix.
    Same as columns where base_flag[n] is False."""

    base_inverse : list[list[int]]
    """A square nmb_nodes x nmb_nodes matrix where the base_columns x base_rows submatrix is 'C^(-1)', the inverse of 
    the base_rows x base_columns invertible sub-submatrix of the adjacency matrix."""

    adj_b_inverse: list[list[int]]
    """A square nmb_nodes x nmb_nodes matrix where the nodes x base_rows submatrix represents 'D = A^{base_columns} * C^(-1)' used in the cut-rank calculations."""

    b_inverse_adj: list[list[int]]
    """A square nmb_nodes x nmb_nodes matrix where the base_columns x nodes submatrix represents 'E = C^(-1) * A_{base_rows}' used in the cut-rank calculations."""

    adj_b_inv_adj: list[list[int]]
    """The square nmb_nodes x nmb_nodes matrix 'F = A^{base_columns} * C^(-1) * A_{base_rows} + A' used in the cut-rank calculations."""

    def __init__(self, partition: 'GraphPartition', rows : list[int], columns : list[int]):
        self.partition = partition
        self.rows = rows[:]
        self.columns = columns[:]

        self._build_matrices()


    def _build_matrices(self) -> None:

        partition = self.partition
        nodes = partition.graph.nodes
        adjacencies = partition.graph.adjacencies

        self.base_inverse = partition._empty_matrix()
        copy_matrix(adjacencies, self.base_inverse, self.rows, self.columns)
        (self.base_rows, self.base_columns) = rank_matrix_positions(self.base_inverse, self.rows, self.columns)
        self.rank = len(self.base_rows)
        copy_matrix(adjacencies, self.base_inverse, self.base_rows, self.base_columns)
        matrix_inverse(self.base_inverse, self.base_inverse, self.base_rows, self.base_columns)

        self.adj_b_inverse = partition._empty_matrix()
        add_product_matrix(adjacencies, self.base_inverse, self.adj_b_inverse, nodes, self.base_columns, self.base_rows)
        self.b_inverse_adj = partition._empty_matrix()
        add_product_matrix(self.base_inverse, adjacencies, self.b_inverse_adj, self.base_columns, self.base_rows, nodes)
        self.adj_b_inv_adj = partition._empty_matrix()
        copy_matrix(adjacencies, self.adj_b_inv_adj, nodes, nodes)
        add_product_matrix(self.adj_b_inverse, adjacencies, self.adj_b_inv_adj, nodes, self.base_rows, nodes)

        self.base_flag = [False] * self.partition.graph.nmb_nodes
        for r in self.base_rows:
            self.base_flag[r] = True
        for c in self.base_columns:
            self.base_flag[c] = True
        self._build_free_nodes()


    def _build_free_nodes(self) -> None:

        self.free_rows = [row for row in self.rows if not self.base_flag[row]]
        self.free_columns = [col for col in self.columns if not self.base_flag[col]]


    def _reduce_base(self, removed_rows : list[int], removed_cols : list[int]) -> None:
    
        if len(removed_rows) == 0:
            return

        nodes = self.partition.graph.nodes
        buffer = self.partition.buffer
            
        # Set base nodes
        for row in removed_rows:
            self.base_flag[row] = False
        for col in removed_cols:
            self.base_flag[col] = False
        self.base_rows = [row for row in self.rows if self.base_flag[row]]
        self.base_columns = [col for col in self.columns if self.base_flag[col]]
        self.rank = len(self.base_rows)

        # Get Z
        copy_matrix(self.base_inverse, buffer, removed_cols, removed_rows)
        matrix_inverse(buffer, buffer, removed_cols, removed_rows)

        # Store D^(Delta X) * Z in D^(Delta Y), update D and F
        insert_zero_matrix(self.adj_b_inverse, nodes, removed_cols)
        add_product_matrix(self.adj_b_inverse, buffer, self.adj_b_inverse, nodes, removed_rows, removed_cols)
        add_product_matrix(self.adj_b_inverse, self.b_inverse_adj, self.adj_b_inv_adj, nodes, removed_cols, nodes)
        add_product_matrix(self.adj_b_inverse, self.base_inverse, self.adj_b_inverse, nodes, removed_cols, self.base_rows)

        # Store (C^-1)_YN^(Delta X) * Z in D^(Delta Y), update C^-1 and E
        insert_zero_matrix(self.adj_b_inverse, nodes, removed_cols)
        add_product_matrix(self.base_inverse, buffer, self.adj_b_inverse, self.base_columns, removed_rows, removed_cols)
        add_product_matrix(self.adj_b_inverse, self.b_inverse_adj, self.b_inverse_adj, self.base_columns, removed_cols, nodes)
        add_product_matrix(self.adj_b_inverse, self.base_inverse, self.base_inverse, self.base_columns, removed_cols, self.base_rows)


    def _extend_base(self, added_rows : list[int], added_cols : list[int]) -> None:
    
        if len(added_rows) == 0:
            return

        nodes = self.partition.graph.nodes
        adjacencies = self.partition.graph.adjacencies
        buffer = self.partition.buffer

        # Determine new base
        for row in added_rows:
            self.base_flag[row] = True
        for col in added_cols:
            self.base_flag[col] = True
        new_base_rows = [row for row in self.rows if self.base_flag[row]]
        new_base_columns = [col for col in self.columns if self.base_flag[col]]

        # Store Z in (C^-1)_(Delta Y)^(Delta X)
        copy_matrix(adjacencies, buffer, added_rows, added_cols)  # Stores (C_N)_(Delta X)^(Delta Y) in position for Z-inverse
        insert_zero_matrix(buffer, added_rows, self.base_rows)
        add_product_matrix(adjacencies, self.base_inverse, buffer, added_rows, self.base_columns, self.base_rows)
        add_product_matrix(buffer, adjacencies, buffer, added_rows, self.base_rows, added_cols)  # Gives Z-inverse
        insert_zero_matrix(self.base_inverse, added_cols, new_base_rows)
        insert_zero_matrix(self.base_inverse, self.base_columns, added_rows)
        matrix_inverse(buffer, self.base_inverse, added_rows, added_cols)

        # Get new C^1
        insert_zero_matrix(buffer, self.base_columns, added_cols)
        add_product_matrix(self.base_inverse, adjacencies, buffer, self.base_columns, self.base_rows, added_cols)
        add_product_matrix(self.base_inverse, buffer, self.base_inverse, added_cols, added_rows, self.base_rows)
        add_product_matrix(buffer, self.base_inverse, self.base_inverse, self.base_columns, added_cols, new_base_rows)

        # Get new D
        copy_matrix(adjacencies, self.adj_b_inverse, nodes, added_cols)
        add_product_matrix(self.adj_b_inverse, adjacencies, self.adj_b_inverse, nodes, self.base_rows, added_cols)  # D_O * C_XO^(Delta Y) + A^(Delta Y) stored in (Delta Y)-column of D
        insert_zero_matrix(self.adj_b_inverse, nodes, added_rows)
        add_product_matrix(self.adj_b_inverse, self.base_inverse, self.adj_b_inverse, nodes, added_cols, new_base_rows)

        # Get new E
        copy_matrix(adjacencies, self.b_inverse_adj, added_rows, nodes)
        add_product_matrix(adjacencies, self.b_inverse_adj, self.b_inverse_adj, added_rows, self.base_columns, nodes)  # C_(Delta X)^YO * E_0 + A_(Delta X) stored in (Delta X)-row of E
        insert_zero_matrix(self.b_inverse_adj, added_cols, nodes)
        add_product_matrix(self.base_inverse, self.b_inverse_adj, self.b_inverse_adj, new_base_columns, added_rows, nodes)

        # Get new F
        insert_zero_matrix(buffer, added_cols, nodes)
        add_product_matrix(self.base_inverse, self.b_inverse_adj, buffer, added_cols, added_rows, nodes)
        add_product_matrix(self.adj_b_inverse, buffer, self.adj_b_inv_adj, nodes, added_cols, nodes)

        self.base_rows = new_base_rows
        self.base_columns = new_base_columns
        self.rank = len(self.base_rows)


    def copy(self, other: 'SubMatrix') -> None:
        """Copy all partition state from other into self."""
        if self.partition.graph != other.partition.graph:
            raise Exception("Cannot copy partition from different graph")
        nodes = self.partition.graph.nodes
        
        copy_matrix(other.base_inverse, self.base_inverse, nodes, nodes)
        copy_matrix(other.adj_b_inverse, self.adj_b_inverse, nodes, nodes)
        copy_matrix(other.b_inverse_adj, self.b_inverse_adj, nodes, nodes)
        copy_matrix(other.adj_b_inv_adj, self.adj_b_inv_adj, nodes, nodes)
        self.base_flag[:] = other.base_flag[:]
        self.rows[:] = other.rows[:]
        self.columns[:] = other.columns[:]
        self.rank = other.rank
        self.base_rows[:] = other.base_rows[:]
        self.base_columns[:] = other.base_columns[:]
        self.free_rows[:] = other.free_rows[:]
        self.free_columns[:] = other.free_columns[:]


    def are_in_partition(self, rows, columns):
        rows_are_in_partition = all(row in self.rows for row in rows)
        if not rows_are_in_partition:
            if any(row in self.rows for row in rows):
                raise Exception("Rows are partially in partition")
            if any(row in self.columns for row in rows):
                raise Exception("Rows overlap partition columns")
            
        columns_are_in_partition = all(col in self.columns for col in columns)
        if not columns_are_in_partition:
            if any(col in self.columns for col in columns):
                raise Exception("Columns are partially in partition")
            if any(col in self.rows for col in columns):
                raise Exception("Columns overlap partition rows")

        return rows_are_in_partition, columns_are_in_partition


    def base_changes(self, row, column, row_is_in_partition: bool, column_is_in_partition: bool) -> tuple[list[int], list[int], list[int], list[int]]:
        """Returns changes to the row and column bases when swapping the given row and column for this partition.
        
        args:
            - row: 'int' The row to be swapped.
            - column: 'int' The column to be swapped.
            - row_is_in_partition: 'bool' If true, row will be removed from self.rows and column will be added instead.
                If false, row is not among self.rows or self.columns. Thus, it is not in the row base, and column cannot enter the row base.
            - column_is_in_partition: 'bool' If true, column will be removed from self.columns and row will be added instead.
                If false, column is not among self.rows or self.columns. Thus, it is not in the column base, and row cannot enter the column base.
        """
        remove_rows : list[int]
        remove_columns : list[int]
        add_rows : list[int]
        add_columns : list[int]

        if (not self.base_flag[column]):

            if (not self.base_flag[row]):

                # row in X^D, column in Y^D
                remove_rows = []
                remove_columns = []

                k2 = next((k2 for k2 in self.free_rows if k2 != row and self.adj_b_inv_adj[k2][row] == 1), -1) if column_is_in_partition else -1
                l2 = next((l2 for l2 in self.free_columns if l2 != column and self.adj_b_inv_adj[column][l2] == 1), -1) if row_is_in_partition else -1
                if k2 >= 0:
                    if l2 >= 0:
                        add_rows = [column, k2]
                        add_columns = [row, l2]
                    else:
                        add_rows = [k2]
                        add_columns = [row]
                else:
                    if l2 >= 0:
                        add_rows = [column]
                        add_columns = [l2]
                    else:
                        if column_is_in_partition and row_is_in_partition and self.adj_b_inv_adj[column][row] == 1:
                            add_rows = [column]
                            add_columns = [row]
                        else:
                            add_rows = []
                            add_columns = []

            else:

                # row in X^B, column in Y^D
                alpha = next(a for a in self.base_columns if self.base_inverse[a][row] == 1)
                remove_rows = [row]
                remove_columns = [alpha]

                k1 = next((k1 for k1 in self.free_rows if self.adj_b_inverse[k1][row] == 1), -1)
                if k1 >= 0:
                    if self.adj_b_inv_adj[k1][row] == 1:
                        k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] != self.adj_b_inverse[k2][row]), -1) if column_is_in_partition else -1
                    else:
                        k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] == 1), -1) if column_is_in_partition else -1
                    l2 = next((l2 for l2 in self.free_columns if l2 != column and self.adj_b_inv_adj[column][l2] == 1), -1)
                    if k2 >= 0:
                        if l2 >= 0:
                            add_rows = [column, k1, k2]
                            add_columns = [row, l2, alpha]
                        else:
                            add_rows = [k1, k2]
                            add_columns = [row, alpha]
                    else:
                        if l2 >= 0:
                            add_rows = [column, k1]
                            add_columns = [l2, alpha]
                        else:
                            if column_is_in_partition and self.adj_b_inv_adj[column][row] != (self.adj_b_inverse[column][row] & self.adj_b_inv_adj[k1][row]):
                                add_rows = [column, k1]
                                add_columns = [row, alpha]
                            else:
                                add_rows = [k1]
                                add_columns = [alpha]
                else:
                    k2 = next((k2 for k2 in self.free_rows if self.adj_b_inv_adj[k2][row] == 1), -1) if column_is_in_partition else -1
                    if self.adj_b_inverse[column][row] == 1:
                        if k2 >= 0:
                            add_rows = [column, k2]
                            add_columns = [row, alpha]
                        else:
                            add_rows = [column]
                            add_columns = [alpha]
                    else:
                        l2 = next((l2 for l2 in self.free_columns if l2 != column and self.adj_b_inv_adj[column][l2] == 1), -1)
                        if k2 >= 0:
                            if l2 >= 0:
                                add_rows = [column, k2]
                                add_columns = [row, l2]
                            else:
                                add_rows = [k2]
                                add_columns = [row]
                        else:
                            if l2 >= 0:
                                add_rows = [column]
                                add_columns = [l2]
                            else:
                                if column_is_in_partition and self.adj_b_inv_adj[column][row] == 1:
                                    add_rows = [column]
                                    add_columns = [row]
                                else:
                                    add_rows = []
                                    add_columns = []

        else:

            if (not self.base_flag[row]):

                # row in X^D, column in Y^B
                beta = next(b for b in self.base_rows if self.base_inverse[column][b] == 1)
                remove_rows = [beta]
                remove_columns = [column]

                l1 = next((l1 for l1 in self.free_columns if self.b_inverse_adj[column][l1] == 1), -1)
                if l1 >= 0:
                    if self.adj_b_inv_adj[column][l1] == 1:
                        l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] != self.b_inverse_adj[column][l2]), -1)
                    else:
                        l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] == 1), -1)
                    k2 = next((k2 for k2 in self.free_rows if k2 != row and self.adj_b_inv_adj[k2][row] == 1), -1)
                    if l2 >= 0 and row_is_in_partition:
                        if k2 >= 0:
                            add_rows = [column, k2, beta]
                            add_columns = [row, l1, l2]
                        else:
                            add_rows = [column, beta]
                            add_columns = [l1, l2]
                    else:
                        if k2 >= 0:
                            add_rows = [k2, beta]
                            add_columns = [row, l1]
                        else:
                            if self.adj_b_inv_adj[column][row] != (self.b_inverse_adj[column][row] & self.adj_b_inv_adj[column][l1]) and row_is_in_partition:
                                add_rows = [column, beta]
                                add_columns = [row, l1]
                            else:
                                add_rows = [beta]
                                add_columns = [l1]
                else:
                    l2 = next((l2 for l2 in self.free_columns if self.adj_b_inv_adj[column][l2] == 1), -1) if row_is_in_partition else -1
                    if self.b_inverse_adj[column][row] == 1:
                        if l2 >= 0:
                            add_rows = [column, beta]
                            add_columns = [row, l2]
                        else:
                            add_rows = [beta]
                            add_columns = [row]
                    else:
                        k2 = next((k2 for k2 in self.free_rows if k2 != row and self.adj_b_inv_adj[k2][row] == 1), -1)
                        if l2 >= 0:
                            if k2 >= 0:
                                add_rows = [column, k2]
                                add_columns = [row, l2]
                            else:
                                add_rows = [column]
                                add_columns = [l2]
                        else:
                            if k2 >= 0:
                                add_rows = [k2]
                                add_columns = [row]
                            else:
                                if self.adj_b_inv_adj[column][row] == 1 and row_is_in_partition:
                                    add_rows = [column]
                                    add_columns = [row]
                                else:
                                    add_rows = []
                                    add_columns = []

            else:

                # row in X^B, column in Y^B
                k1 = next((k1 for k1 in self.free_rows if self.adj_b_inverse[k1][row] == 1), -1)
                l1 = next((l1 for l1 in self.free_columns if self.b_inverse_adj[column][l1] == 1), -1)
                if (self.base_inverse[column][row] == 1):

                    # Full rank matrix with row and column removed is invertible
                    remove_rows = [row]
                    remove_columns = [column]

                    if k1 >= 0 and l1 >= 0:
                        if self.adj_b_inv_adj[k1][row] == 1:
                            k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] != self.adj_b_inverse[k2][row]), -1)
                        else:
                            k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] == 1), -1)
                        if self.adj_b_inv_adj[column][l1] == 1:
                            l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] != self.b_inverse_adj[column][l2]), -1)
                        else:
                            l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] == 1), -1)
                        if k2 >= 0:
                            if l2 >= 0:
                                add_rows = [column, k1, k2]
                                add_columns = [row, l1, l2]
                            else:
                                add_rows = [k1, k2]
                                add_columns = [row, l1]
                        else:
                            if l2 >= 0:
                                add_rows = [column, k1]
                                add_columns = [l1, l2]
                            else:
                                if ((self.adj_b_inv_adj[k1][row] & self.adj_b_inv_adj[column][l1]) ^ (self.adj_b_inv_adj[k1][row] & self.adj_b_inverse[column][row]) ^ (self.adj_b_inv_adj[column][l1] & self.b_inverse_adj[column][row])) != self.adj_b_inv_adj[column][row]:
                                    add_rows = [column, k1]
                                    add_columns = [row, l1]
                                else:
                                    add_rows = [k1]
                                    add_columns = [l1]

                    else:
                        k = next((k for k in self.free_rows if self.adj_b_inv_adj[k][row] != (self.adj_b_inverse[k][row] & self.b_inverse_adj[column][row])), -1)
                        l = next((l for l in self.free_columns if self.adj_b_inv_adj[column][l] != (self.adj_b_inverse[column][row] & self.b_inverse_adj[column][l])), -1)
                        if k >= 0:
                            if l >= 0:
                                add_rows = [column, k]
                                add_columns = [row, l]
                            else:
                                add_rows = [k]
                                add_columns = [row]
                        else:
                            if l >= 0:
                                add_rows = [column]
                                add_columns = [l]
                            else:
                                if self.adj_b_inv_adj[column][row] != (self.adj_b_inverse[column][row] & self.b_inverse_adj[column][row]):
                                    add_rows = [column]
                                    add_columns = [row]
                                else:
                                    add_rows = []
                                    add_columns = []

                else:

                    # Full rank matrix with row and column removed is singular
                    alpha = next(a for a in self.base_columns if self.base_inverse[a][row] == 1)
                    beta = next(b for b in self.base_rows if self.base_inverse[column][b] == 1)
                    remove_rows = [row, beta]
                    remove_columns = [column, alpha]

                    if k1 >= 0:

                        if l1 >= 0:

                            # Case k1 >= 0 and l1 >= 0
                            if self.adj_b_inv_adj[k1][row] == 1:
                                k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] != self.adj_b_inverse[k2][row]), -1)
                            else:
                                k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] == 1), -1)
                            if self.adj_b_inv_adj[column][l1] == 1:
                                l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] != self.b_inverse_adj[column][l2]), -1)
                            else:
                                l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] == 1), -1)
                            if k2 >= 0:
                                if l2 >= 0:
                                    add_rows = [column, k1, k2, beta]
                                    add_columns = [row, l1, l2, alpha]
                                else:
                                    add_rows = [k1, k2, beta]
                                    add_columns = [row, l1, alpha]
                            else:
                                if l2 >= 0:
                                    add_rows = [column, k1, beta]
                                    add_columns = [l1, l2, alpha]
                                else:
                                    if ((self.adj_b_inv_adj[k1][row] & self.adj_b_inverse[column][row]) ^ (self.adj_b_inv_adj[column][l1] & self.b_inverse_adj[column][row])) != self.adj_b_inv_adj[column][row]:
                                        add_rows = [column, k1, beta]
                                        add_columns = [row, l1, alpha]
                                    else:
                                        add_rows = [k1, beta]
                                        add_columns = [l1, alpha]

                        else:

                            # Case k1 >= 0 and l1 < 0
                            l2 = next((l2 for l2 in self.free_columns if self.adj_b_inv_adj[column][l2] == 1), -1)
                            if l2 >= 0:
                                if self.b_inverse_adj[column][row] == 1:
                                    add_rows = [column, k1, beta]
                                    add_columns = [row, l2, alpha]
                                else:
                                    if self.adj_b_inv_adj[k1][row] == 1:
                                        k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] != self.adj_b_inverse[k2][row]), -1)
                                    else:
                                        k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] == 1), -1)
                                    if k2 >= 0:
                                        add_rows = [column, k1, k2]
                                        add_columns = [row, l2, alpha]
                                    else:
                                        add_rows = [column, k1]
                                        add_columns = [l2, alpha]
                            else:
                                if self.b_inverse_adj[column][row] == 1:
                                    add_rows = [k1, beta]
                                    add_columns = [row, alpha]
                                else:
                                    if self.adj_b_inv_adj[k1][row] == 1:
                                        k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] != self.adj_b_inverse[k2][row]), -1)
                                    else:
                                        k2 = next((k2 for k2 in self.free_rows if k2 != k1 and self.adj_b_inv_adj[k2][row] == 1), -1)
                                    if k2 >= 0:
                                        add_rows = [k1, k2]
                                        add_columns = [row, alpha]
                                    else:
                                        if (self.adj_b_inv_adj[k1][row] & self.adj_b_inverse[column][row]) != self.adj_b_inv_adj[column][row]:
                                            add_rows = [column, k1]
                                            add_columns = [row, alpha]
                                        else:
                                            add_rows = [k1]
                                            add_columns = [alpha]

                    else:

                        if l1 >= 0:

                            # Case k1 < 0 and l1 >= 0
                            k2 = next((k2 for k2 in self.free_rows if self.adj_b_inv_adj[k2][row] == 1), -1)
                            if k2 >= 0:
                                if self.adj_b_inverse[column][row] == 1:
                                    add_rows = [column, k2, beta]
                                    add_columns = [row, l1, alpha]
                                else:
                                    if self.adj_b_inv_adj[column][l1] == 1:
                                        l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] != self.b_inverse_adj[column][l2]), -1)
                                    else:
                                        l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] == 1), -1)
                                    if l2 >= 0:
                                        add_rows = [column, k2, beta]
                                        add_columns = [row, l1, l2]
                                    else:
                                        add_rows = [k2, beta]
                                        add_columns = [row, l1]
                            else:
                                if self.adj_b_inverse[column][row] == 1:
                                    add_rows = [column, beta]
                                    add_columns = [l1, alpha]
                                else:
                                    if self.adj_b_inv_adj[column][l1] == 1:
                                        l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] != self.b_inverse_adj[column][l2]), -1)
                                    else:
                                        l2 = next((l2 for l2 in self.free_columns if l2 != l1 and self.adj_b_inv_adj[column][l2] == 1), -1)
                                    if l2 >= 0:
                                        add_rows = [column, beta]
                                        add_columns = [l1, l2]
                                    else:
                                        if (self.adj_b_inv_adj[column][l1] & self.b_inverse_adj[column][row]) != self.adj_b_inv_adj[column][row]:
                                            add_rows = [column, beta]
                                            add_columns = [row, l1]
                                        else:
                                            add_rows = [beta]
                                            add_columns = [l1]

                        else:

                            # Case k1 < 0 and l1 < 0
                            if self.adj_b_inverse[column][row] == 1:
                                if self.b_inverse_adj[column][row] == 1:
                                    add_rows = [column, beta]
                                    add_columns = [row, alpha]
                                else:
                                    k2 = next((k2 for k2 in self.free_rows if self.adj_b_inv_adj[k2][row] == 1), -1)
                                    if k2 >= 0:
                                        add_rows = [column, k2]
                                        add_columns = [row, alpha]
                                    else:
                                        add_rows = [column]
                                        add_columns = [alpha]
                            else:
                                if self.b_inverse_adj[column][row] == 1:
                                    l2 = next((l2 for l2 in self.free_columns if self.adj_b_inv_adj[column][l2] == 1), -1)
                                    if l2 >= 0:
                                        add_rows = [column, beta]
                                        add_columns = [row, l2]
                                    else:
                                        add_rows = [beta]
                                        add_columns = [row]
                                else:
                                    k2 = next((k2 for k2 in self.free_rows if self.adj_b_inv_adj[k2][row] == 1), -1)
                                    l2 = next((l2 for l2 in self.free_columns if self.adj_b_inv_adj[column][l2] == 1), -1)
                                    if k2 >= 0:
                                        if l2 >= 0:
                                            add_rows = [column, k2]
                                            add_columns = [row, l2]
                                        else:
                                            add_rows = [k2]
                                            add_columns = [row]
                                    else:
                                        if l2 >= 0:
                                            add_rows = [column]
                                            add_columns = [l2]
                                        else:
                                            if self.adj_b_inv_adj[column][row] == 1:
                                                add_rows = [column]
                                                add_columns = [row]
                                            else:
                                                add_rows = []
                                                add_columns = []
                                                
        return remove_rows, remove_columns, add_rows, add_columns


    def apply_swap(self, row : int, column : int, row_is_in_partition: bool, column_is_in_partition: bool) -> int:

        remove_rows, remove_columns, add_rows, add_columns = self.base_changes(row, column, row_is_in_partition, column_is_in_partition)

        # Set new set of rows and columns
        if row_is_in_partition:
            row_idx = self.rows.index(row)
            self.rows[row_idx] = column
        if column_is_in_partition:
            col_idx = self.columns.index(column)
            self.columns[col_idx] = row

        # Apply reduction and extension
        self._reduce_base(remove_rows, remove_columns)
        self._extend_base(add_rows, add_columns)

        # Update set of fre rows and free columns
        self._build_free_nodes()

        return len(add_rows) - len(remove_rows)
