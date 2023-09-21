#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


typedef struct matrix {
    unsigned int m;
    unsigned int n;
    double* data;
} matrix;


matrix* zmatrix(unsigned int m, unsigned int n);

matrix* imatrix(unsigned int n);

int matrix_set_r(matrix* A, double* data, int length);

matrix* matrix_set(unsigned int m, unsigned int n, double* data, int length);

matrix* matrix_copy(matrix* A);

void print_matrix(matrix* A);

void free_matrix(matrix* A);


int matrix_equals(matrix* A, matrix* B);

matrix* matrix_col_get(matrix* A, unsigned int j);

matrix* matrix_row_get(matrix* A, unsigned int i);

int matrix_col_cut_r(matrix* A, unsigned int j);

int matrix_row_cut_r(matrix* A, unsigned int i);

int matrix_col_put_r(matrix* A, unsigned int j, double* data, int length);

int matrix_row_put_r(matrix* A, unsigned int i, double* data, int length);


int matrix_col_add_mult_r(matrix* A, unsigned int j, double s, matrix* col);

int matrix_row_add_mult_r(matrix* A, unsigned int i, double s, matrix* row);

int matrix_col_swap_r(matrix* A, unsigned int j1, unsigned int j2);

int matrix_row_swap_r(matrix* A, unsigned int i1, unsigned int i2);

matrix* matrix_concat_top_down(matrix** Alist, unsigned int Anum);

matrix* matrix_concat_left_right(matrix** Alist, unsigned int Anum);


void matrix_sclr_mult_r(double s, matrix* A);

matrix* matrix_sclr_mult(double s, matrix* A);

void matrix_T_r(matrix* A);

matrix* matrix_T(matrix* A);

int matrix_add_r(matrix* A, matrix* B);

matrix* matrix_add(matrix* A, matrix* B);

int matrix_minus_r(matrix* A, matrix* B);

matrix* matrix_minus(matrix* A, matrix* B);

matrix* matrix_mtrx_mult(matrix* A, matrix* B);


double matrix_det(matrix* A);

void matrix_ref_r(matrix* A);

matrix* matrix_ref(matrix* A);

void matrix_rref_r(matrix* A);

matrix* matrix_rref(matrix* A);

matrix* matrix_inverse(matrix* A);


matrix* zmatrix(unsigned int m, unsigned int n) {
    matrix* A = (matrix*) malloc(sizeof(matrix));
    A->m = m;
    A->n = n;
    A->data = (double*) calloc(m * n, sizeof(double));
    return A;
}

matrix* imatrix(unsigned int n) {
    matrix* I = zmatrix(n,n);
    for(int i = 0; i < n; i++) {
        I->data[i * (n + 1)] = 1.0;
    }
    return I;
}

matrix* matrix_copy(matrix* A){
    matrix* A2 = zmatrix(A->m, A->n);
    for (int i = 0; i < A->m * A->n; i++) {
        A2->data[i] = A->data[i];
    }
    return A2;
}

int matrix_set_r(matrix* A, double* data, int length) {
    if (A->m * A->n != length) {
        fprintf(stderr, "Wrong number of entries. %s.\n", strerror(22));
        return 0;
    }
    for (int i = 0; i < A->m * A->n; i++) {
        if (fabs(data[i]) < 0.0000005) {
            A->data[i] = 0;
        } else {
            A->data[i] = data[i];
        }
    }
    return 1;
}

matrix* matrix_set(unsigned int m, unsigned int n, double* data, int length) {
    matrix* A = zmatrix(m, n);
    if(matrix_set_r(A, data, length)) {
        return A;
    }
    free_matrix(A);
    return NULL;
}

void print_matrix(matrix* A) {
    for (int i = 0; i < A->m; i++) {
        printf("[");
        for (int j = 0; j < A->n; j++) {
            printf(" %f ", A->data[i * A->n + j]);
        }
        printf("]\n");
    }
    printf("\n");
}

void free_matrix(matrix* A) {
    free(A->data);
    free(A);
}


int matrix_equals(matrix* A, matrix* B) {
    if (A->m != B->m || A->n != B->n) {
        return 0;
    }
    for (int i = 0; i < A->m * A->n; i++) {
        if (A->data[i] != B->data[i]) {
            return 0;
        }
    }
    return 1;
}

matrix* matrix_col_get(matrix* A, unsigned int j) {
    if (j >= A->n) {
        fprintf(stderr, "No such column. %s.\n", strerror(22));
        return NULL;
    }
    matrix* Aj = zmatrix(A->m, 1);
    for (int i = 0; i < A->m; i++) {
        Aj->data[i] = A->data[A->n * i + j];
    }
    return Aj;
}

matrix* matrix_row_get(matrix* A, unsigned int i) {
    if (i >= A->m) {
        fprintf(stderr, "No such row. %s.\n", strerror(22));
        return NULL;
    }
    matrix* ai = zmatrix(1, A->n);
    for (int j = 0; j < A->n; j++) {
        ai->data[j] = A->data[A->n * i + j];
    }
    return ai;
}

int matrix_col_cut_r(matrix* A, unsigned int j) {
    if (j >= A->n) {
        fprintf(stderr, "No such column. %s.\n", strerror(22));
        return 0;
    }

    for (int i = 0; i < A->m * (A->n - 1); i++) {
        A->data[i] = A->data[i + (i + A->n - 1 - j) / (A->n - 1)];
    }

    A->n--;
    A->data = (double*) realloc(A->data, A->m * A->n * sizeof(double));
    return 1;
}

int matrix_row_cut_r(matrix* A, unsigned int i) {
    if (i >= A->m) {
        fprintf(stderr, "No such row. %s.\n", strerror(22));
        return 0;
    }
    for (int j = i * A->n; j < (A->m - 1) * A->n; j++) {
        A->data[j] = A->data[j + A->n];
    }

    A->m--;
    A->data = (double*) realloc(A->data, A->m * A->n * sizeof(double));
    return 1;
}

int matrix_col_put_r(matrix* A, unsigned int j, double* data, int length) {
    if (j >= A->n) {
        fprintf(stderr, "No such column. %s.\n", strerror(22));
        return 0;
    }
    if (length != A->m) {
        fprintf(stderr, "Wrong column length. %s.\n", strerror(22));
        return 0;
    }
    for (int i = 0; i < A->m; i++) {
        if (fabs(data[i]) < 0.0000005) {
            A->data[i * A->n + j] = 0;
        } else {
            A->data[i * A->n + j] = data[i];
        }
    }
    return 1;
}

int matrix_row_put_r(matrix* A, unsigned int i, double* data, int length) {
    if (i >= A->m) {
        fprintf(stderr, "No such row. %s.\n", strerror(22));
        return 0;
    }
    if (length != A->n) {
        fprintf(stderr, "Wrong row length. %s.\n", strerror(22));
        return 0;
    }
    for (int j = 0; j < A->n; j++) {
        if (fabs(data[j]) < 0.0000005) {
            A->data[i * A->n + j] = 0;
        } else {
            A->data[i * A->n + j] = data[j];
        }
    }
    return 1;
}


int matrix_col_add_mult_r(matrix* A, unsigned int j, double s, matrix* col) {
    if (j >= A->n) {
        fprintf(stderr, "No such column. %s.\n", strerror(22));
        return 0;
    }
    if (col->m != A->m) {
        fprintf(stderr, "Wrong column length. %s.\n", strerror(22));
        return 0;
    }
    if (col->n != 1) {
        fprintf(stderr, "Not a column. %s.\n", strerror(22));
        return 0;
    }

    if (fabs(s) < 0.0000005) {
        return 1;
    }
    for (int i = 0; i < A->m; i++) {
        A->data[i * A->n + j] += s * col->data[i];
        if (fabs(A->data[i * A->n + j]) < 0.0000005) {
            A->data[i * A->n + j] = 0;
        }
    }
    return 1;
}

int matrix_row_add_mult_r(matrix* A, unsigned int i, double s, matrix* row) {
    if (i >= A->m) {
        fprintf(stderr, "No such row. %s.\n", strerror(22));
        return 0;
    }
    if (row->n != A->n) {
        fprintf(stderr, "Wrong row length. %s.\n", strerror(22));
        return 0;
    }
    if (row->m != 1) {
        fprintf(stderr, "Not a row. %s.\n", strerror(22));
        return 0;
    }

    if (fabs(s) < 0.0000005) {
        return 1;
    }
    for (int j = 0; j < A->n; j++) {
        A->data[i * A->n + j] += s * row->data[j];
        if (fabs(A->data[i * A->n + j]) < 0.0000005) {
            A->data[i * A->n + j] = 0;
        }

    }
    return 1;
}

int matrix_col_swap_r(matrix* A, unsigned int j1, unsigned int j2) {
    if (j1 >= A->n || j2 >= A->n) {
        fprintf(stderr, "No such column. %s.\n", strerror(22));
        return 0;
    }
    if (j1 == j2) {
        return 1;
    }

    matrix* Aj1 = matrix_col_get(A, j1);
    matrix* Aj2 = matrix_col_get(A, j2);
    matrix_col_put_r(A, j2, Aj1->data, Aj1->m * Aj1->n);
    matrix_col_put_r(A, j1, Aj2->data, Aj2->m * Aj2->n);
    free_matrix(Aj1);
    free_matrix(Aj2);
    return 1;
}

int matrix_row_swap_r(matrix* A, unsigned int i1, unsigned int i2) {
    if (i1 >= A->m || i2 >= A->m) {
        fprintf(stderr, "No such row. %s.\n", strerror(22));
        return 0;
    }
    if (i1 == i2) {
        return 1;
    }

    matrix* ai1 = matrix_row_get(A, i1);
    matrix* ai2 = matrix_row_get(A, i2);
    matrix_row_put_r(A, i2, ai1->data, ai1->m * ai1->n);
    matrix_row_put_r(A, i1, ai2->data, ai2->m * ai2->n);
    free_matrix(ai1);
    free_matrix(ai2);
    return 1;
}

matrix* matrix_concat_top_down(matrix** Alist, unsigned int Anum) {
    if (!Anum) {
        return NULL;
    }
    if (Alist[0] == NULL) {
        fprintf(stderr, "Invalid matrix. %s.\n", strerror(22));
        return NULL;
    }
    if (Anum == 1) {
        return matrix_copy(Alist[0]);
    }
    int m = Alist[0]->m;
    int n = Alist[0]->n;

    for (int i = 1; i < Anum; i++) {
        if (Alist[i] == NULL) {
            fprintf(stderr, "Invalid matrix. %s.\n", strerror(22));
            return NULL;
        }
        if (Alist[i]->n != n) {
            fprintf(stderr, "Cannot concatenate. %s.\n", strerror(22));
            return NULL;
        }
        m += Alist[i]->m;
    }
    matrix* A = zmatrix(m, n);

    int i = 0;
    int depth = 0;
    for (int j = 0; j < A->m; j++) {
        while (depth == Alist[i]->m) {
            i++;
            depth = 0;
        }
        for (int k = 0; k < A->n; k++) {
            A->data[j * A->n + k] = Alist[i]->data[depth * A->n + k];
        }
        depth++;
    }
    return A;
}

matrix* matrix_concat_left_right(matrix** Alist, unsigned int Anum) {
    if (!Anum) {
        return NULL;
    }
    if (Alist[0] == NULL) {
        fprintf(stderr, "Invalid matrix. %s.\n", strerror(22));
        return NULL;
    }
    if (Anum == 1) {
        return matrix_copy(Alist[0]);
    }
    int m = Alist[0]->m;
    int n = Alist[0]->n;

    for (int i = 1; i < Anum; i++) {
        if (Alist[i] == NULL) {
            fprintf(stderr, "Invalid matrix. %s.\n", strerror(22));
            return NULL;
        }
        if (Alist[i]->m != m) {
            fprintf(stderr, "Cannot concatenate. %s.\n", strerror(22));
            return NULL;
        }
        n += Alist[i]->n;
    }
    matrix* A = zmatrix(m, n);

    int i = 0;
    int depth = 0;
    for (int j = 0; j < A->n; j++) {
        while (depth == Alist[i]->n) {
            i++;
            depth = 0;
        }
        for (int k = 0; k < A->m; k++) {
            A->data[k * A->n + j] = Alist[i]->data[k * Alist[i]->n + depth];
        }
        depth++;
    }
    return A;
}


void matrix_sclr_mult_r(double s, matrix* A) {
    if (fabs(s) < 0.0000005) {
        int m = A->m;
        int n = A->n;
        free_matrix(A);
        A = zmatrix(m, n);
        return;
    }
    for (int i = 0; i < A->m * A->n; i++) {
        A->data[i] *= s;
        if (fabs(A->data[i]) < 0.0000005) {
            A->data[i] = 0;
        }
    }
}

matrix* matrix_sclr_mult(double s, matrix* A) {
    matrix* sA = matrix_copy(A);
    matrix_sclr_mult_r(s, sA);
    return sA;
}

void matrix_T_r(matrix* A) {
    if (A->m > 1 && A->n > 1) {
        for(int i = 0; i < A->m * A->n; i++) {
            int j = (i / A->n) + A->n * (i % A->n);
            if (i < j) {
                double temp = A->data[i];
                A->data[i] = A->data[j];
                A->data[j] = temp;
            }
        }
    }

    if (A->m != A->n) {
        int k = A->m;
        A->m = A->n;
        A->n = k;
    }
}

matrix* matrix_T(matrix* A) {
    matrix* AT = matrix_copy(A);
    matrix_T_r(AT);
    return AT;
}

int matrix_add_r(matrix* A, matrix* B) {
    if (A->m != B->m || A->n != B->n) {
        fprintf(stderr, "Size mismatch. %s.\n", strerror(22));
        return 0;
    }
    for (int i = 0; i < A->m * A->n; i++) {
        A->data[i] += B->data[i];
        if (fabs(A->data[i]) < 0.0000005) {
            A->data[i] = 0;
        }
    }
    return 1;
}

matrix* matrix_add(matrix* A, matrix* B) {
    matrix* C = matrix_copy(A);
    if (!matrix_add_r(C, B)) {
        free_matrix(C);
        return NULL;
    }
    return C;
}

int matrix_minus_r(matrix* A, matrix* B) {
    if (A->m != B->m || A->n != B->n) {
        fprintf(stderr, "Size mismatch. %s.\n", strerror(22));
        return 0;
    }
    for (int i = 0; i < A->m * A->n; i++) {
        A->data[i] -= B->data[i];
        if (fabs(A->data[i]) < 0.0000005) {
            A->data[i] = 0;
        }
    }
    return 1;
}

matrix* matrix_minus(matrix* A, matrix* B) {
    matrix* C = matrix_copy(A);
    if (!matrix_minus_r(C, B)) {
        free_matrix(C);
        return NULL;
    }
    return C;
}

matrix* matrix_mtrx_mult(matrix* A, matrix* B) {
    if (A->n != B->m) {
        fprintf(stderr, "Wrong sizes. %s.\n", strerror(22));
        return NULL;
    }

    matrix* C = zmatrix(A->m, B->n);
    for (int i = 0; i < A->m; i++) {
        for (int j = 0; j < B->n; j++) {
            for (int k = 0; k < A->n; k++) {
                C->data[i * B->n + j] += A->data[i * A->n + k] * B->data[k * B->m + j];
            }
            if (fabs(C->data[i * B->n + j]) < 0.0000005) {
                C->data[i * B->n + j] = 0;
            }
        }
    }
    return C;
}


double matrix_det(matrix* A) {
    if (A->m != A->n) {
        fprintf(stderr, "Not square. %s.\n", strerror(22));
        return 0;
    }

    double ans = 0;
    int* permute = (int*) calloc(A->m, sizeof(int));
    int lessFour = 0;
    int breaker = 0;
    while(!breaker) {
        lessFour++;
        if (lessFour == 4) {
            lessFour = 0;
        }

        double product = lessFour < 2 ? 1 : -1;
        int* table = (int*) calloc(A->m, sizeof(int));
        for (int i = 1; i < A->m; i++) {
            table[i] = i;
        }

        for (int i = 0; i < A->m; i++) {
            product *= A->data[i * A->m + table[permute[i]]];
            for (int j = permute[i]; j < A->m - i - 1; j++) {
                table[j] = table[j + 1];
            }
            table = (int*) realloc(table, (A->m - i - 1) * sizeof(int));
        }
        ans += product;

        breaker = 1;
        for (int i = 1; i < A->m; i++) {
            if (permute[A->m - i - 1] < i) {
                permute[A->m - i - 1]++;
                for (int j = 1; j < i; j++) {
                    permute[A->m - j - 1] = 0;
                }
                breaker = 0;
                break;
            }
        }
    }
    free(permute);
    if (fabs(ans) < 0.0000005) {
        return 0;
    }
    return ans;
}

void matrix_ref_r(matrix* A) {
    int echelon_rows = 0;
    for (int i = 0; i < A->n; i++) {
        for (int j = echelon_rows; j < A->m; j++) {
            if (A->data[j * A->n + i] != 0) {
                matrix* row = matrix_row_get(A, j);
                matrix_sclr_mult_r(1 / row->data[j], row);
                row->data[j] = 1;
                matrix_row_swap_r(A, echelon_rows, j);
                matrix_row_put_r(A, echelon_rows, row->data, row->m * row->n);

                echelon_rows++;
                for (int k = echelon_rows; k < A->m; k++) {
                    matrix_row_add_mult_r(A, k, -A->data[k * A->n + i], row);
                }
                free_matrix(row);
                break;
            }
        }
    }

    for (int i = 0; i < A->m * A->n; i++) {
        if (fabs(A->data[i]) < 0.0000005) {
            A->data[i] = 0;
        }
    }
}

matrix* matrix_ref(matrix* A) {
    matrix* Aref = matrix_copy(A);
    matrix_ref_r(Aref);
    return Aref;
}

void matrix_rref_r(matrix* A) {
    matrix_ref_r(A);
    for (int i = 0; i < A->n; i++) {
        for (int j = A->m - 1; j >= 0; j--) {
            if (A->data[j * A->n + i] == 1) {
                matrix* row = matrix_row_get(A, j);
                for (int k = j - 1; k >=0; k--) {
                    matrix_row_add_mult_r(A, k, -A->data[k * A->n + i], row);
                }
                break;
            }
        }
    }

    for (int i = 0; i < A->m * A->n; i++) {
        if (fabs(A->data[i]) < 0.0000005) {
            A->data[i] = 0;
        }
    }

}

matrix* matrix_rref(matrix* A) {
    matrix* Arref = matrix_copy(A);
    matrix_rref_r(Arref);
    return Arref;
}

matrix* matrix_inverse(matrix* A) {
    if (A->m != A->n) {
        fprintf(stderr, "Not square. %s.\n", strerror(22));
        return NULL;
    }
    matrix** Alist = (matrix**) malloc(2 * sizeof(matrix*));
    Alist[0] = A;
    matrix* I = imatrix(A->m);
    Alist[1] = I;
    matrix* AI = matrix_concat_left_right(Alist, 2);
    free_matrix(I);
    free(Alist);

    matrix_rref_r(AI);
    int invertible = 1;
    for (int i = 0; i < A->m; i++) {
        if (AI->data[2 * i * A->n + i] != 1) {
            invertible = 0;
            break;
        }
    }
    if (!invertible) {
        free_matrix(AI);
        return NULL;
    }

    matrix* Ainv = zmatrix(A->m, A->n);
    for (int i = 0; i < A->m; i++) {
        for (int j = 0; j < A->n; j++) {
            Ainv->data[i * A->n + j] = AI->data[(2 * i + 1) * A->n + j];
        }
    }
    free_matrix(AI);
    return Ainv;
}