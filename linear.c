#include "matrix.c"

int main() {
    matrix* A = zmatrix(3, 3);
    double Adata[] = {1, 2, 0.5, -2, 1.5, -1, 0.5, -2, -1};
    matrix_set_r(A, Adata, 9);

    double Bdata[] = {2, -2, 1.5, -1.5, 1, 1, 0.5, 2, 0.5};
    matrix* B = matrix_set(3, 3, Bdata, 9);

    matrix* C = matrix_mtrx_mult(B, A);
    print_matrix(C);
    free_matrix(A);
    free_matrix(B);

    matrix* Cinv = matrix_inverse(C);
    print_matrix(Cinv);

    matrix* I = matrix_mtrx_mult(C, Cinv);
    print_matrix(I);
    free_matrix(C);
    free_matrix(Cinv);
    free_matrix(I);
    return 0;
}