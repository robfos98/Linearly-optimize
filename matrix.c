#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct matrix {
    unsigned int m;
    unsigned int n;
    double *data;
} matrix;

matrix *zmatrix(unsigned int m, unsigned int n);

void free_matrix(matrix *A);


matrix *zmatrix(unsigned int m, unsigned int n) {
    if (m == 0) {
        fprintf(stderr, "Matrices have rows. %s.\n", strerror(22));
        return NULL;
    }
    if (n == 0) {
        fprintf(stderr, "Matrices have columns. %s.\n", strerror(22));
        return NULL;
    }

    matrix* A = (matrix*) malloc(sizeof(matrix));
    A->m = m;
    A->n = n;
    A->data = (double*) calloc(m * n, sizeof(double));
    return A;
}

void free_matrix(matrix *A) {
    free(A->data);
    free(A);
}


int main() {
    matrix* A = zmatrix(2,3);
    for (int i = 0; i < 6; i++) {
        printf("%d: %f\n", i, A->data[i]);
    }

    free_matrix(A);
    return 0;
}