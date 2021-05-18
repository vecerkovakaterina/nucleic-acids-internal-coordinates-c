#include <stdio.h>

/**
 * Multiplies two matrices
 * lm * mn = ln
 * @param l number of rows of the first matrix
 * @param m number of columns of the first matrix and number of rows of the second matrix
 * @param n number of columns of the second matrix
 * @param res the resulting matrix array
 * @param m1 the first matrix
 * @param m2 the second matrix
 */
static void matrix_multiplication(int l, int m, int n, double res[l][n], double m1[l][m], double m2[l][n]) {
    for (int i = 0; i < l; i++) {
        for (int j = 0; j < n; j++) {
            res[i][j] = 0;
        }
    }

    for (int i = 0; i < l; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < m; k++) {
                res[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
}

/**
 * Transpose matrix
 * @param m number of rows of the matrix to be transposed
 * @param n number of columns of the matrix to be transposed
 * @param transposed th resulting transposed matrix
 * @param matrix the matrix to be transposed
 */
static void transpose_matrix(int m, int n, double transposed[n][m], double matrix[m][n]) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            transposed[j][i] = matrix[i][j];
        }
}

/**
 * Subtract two matrices of the same dimensions
 * m1 - m2 = res
 * @param m number of rows of the matrices
 * @param n number of columns of the matrices
 * @param res the difference matrix
 * @param m1 minuend
 * @param m2 subtrahend
 */
static void subtract_matrices_generic(int m, int n, double res[m][n], double m1[m][n], double m2[m][n]) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            res[i][j] = m1[i][j] - m2[i][j];
        }
    }
}

/**
 * Multiply matrix by scalar
 * @param m number of rows of the matrix
 * @param n number of columns of the matrix
 * @param res the product
 * @param scalar the scalar multiplier
 * @param matrix the matrix multiplicand
 */
static void multiply_matrix_by_scalar(int m, int n, double res[m][n], double scalar, double matrix[m][n]) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            res[i][j] = scalar * matrix[i][j];
        }
    }
}

/**
 * Add up values in columns into a vector
 * @param m number of rows in the matrix
 * @param n number of columns in the matrix
 * @param res resulting vector with sums
 * @param matrix the matrix whose columns are to be added
 */
static void get_sum_of_columns_matrix(int m, int n, double res[n], double matrix[m][n]) {
    for (int i = 0; i < n; i++) {
        res[i] = 0.;
        for (int j = 0; j < m; j++) {
            res[i] += matrix[j][i];
        }
    }
}

/**
 * Repeat the vector in m rows to form a m, n matrix
 * @param m number of rows
 * @param n number of columns
 * @param res the resulting m, n matrix
 * @param matrix the m, 1 matrix of ones
 * @param vec the vector to be repeated
 */
static void spread_vector_into_m_rows(int m, int n, double res[m][n], double matrix[][1], const double vec[n]) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            res[i][j] = matrix[i][0] * vec[j];
        }
    }
}

/**
 * Average the values in matrix columns into a vector
 * @param m number of rows in the matrix
 * @param n number of columns in the matrix
 * @param avg the resulting vector
 * @param matrix the matrix whose columns are to be averaged
 */
static void average_matrix_columns(int m, int n, double avg[n], double matrix[m][n]) {
    for (int i = 0; i < n; i++) {
        avg[i] = 0.;
        for (int j = 0; j < m; j++) {
            avg[i] += matrix[j][i];
        }
        avg[i] /= m;
    }
}

/**
 * Multiply a m, n matrix and a n vector to get a n vector
 * @param m number of rows of the matrix
 * @param n number of columns of the matrix and also length of the vector
 * @param res the resulting n vector
 * @param matrix the matrix multiplier
 * @param vec the vector multiplicand
 */
static void matrix_vector_multiplication(int m, int n, double res[m], double matrix[m][n], const double vec[n]) {
    for (int i = 0; i < n; i++) {
        res[i] = 0;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            res[j] += matrix[i][j] * vec[i];
        }
    }
}

/**
 * Multiply a n vector and a m, n matrix to get a n vector
 * @param m number of rows of the matrix
 * @param n number of columns of the matrix and also length of the vector
 * @param res the resulting n vector
 * @param vec the vector multiplier
 * * @param matrix the matrix multiplicand
 */
static void vector_matrix_multiplication(int m, int n, double res[m], const double vec[n], double matrix[m][n]) {
    for (int i = 0; i < n; i++) {
        res[i] = 0;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            res[j] += vec[i] * matrix[i][j];
        }
    }
}

/**
 * Subtract two vectors v1 - v2 of the same length n
 * @param n the length of the vectors
 * @param diff the difference vector
 * @param v1 the minuend vector
 * @param v2 subtrahend vector
 */
static void subtract_two_vectors(int n, double diff[n], const double v1[n], const double v2[n]) {
    for (int i = 0; i < n; i++) {
        diff[i] = v1[i] - v2[i];
    }
}

/**
 * Return the dot product of two vectors of the same length n
 * @param n the length of the two vectors
 * @param vec1 the first vector
 * @param vec2 the second vector
 * @return the dot product
 */
static double dot_product(int n, const double vec1[], const double vec2[]) {
    double dp = 0.;
    for (int i = 0; i < n; i++) {
        dp += vec1[i] * vec2[i];
    }
    return dp;
}

/**
 * Perform the cross product of two 3D vectors
 * @param cross_product the resulting cross product vector
 * @param vec1 the first vector
 * @param vec2 the second vector
 */
static void cross_product(double cross_product[3], const double vec1[3], const double vec2[3]) {
    cross_product[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
    cross_product[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
    cross_product[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
}

/**
 * Return the magnitude of matrix column
 * @param m number of rows of the matrix
 * @param n number of columns of the matrix
 * @param matrix the matrix whose columns magnitudes are to be found
 * @param column index of the column to return the magnitude of
 * @return magnitude of selected column
 */
static double matrix_column_magnitude(int m, int n, double matrix[m][n], int column) {
    double mag = 0.;
    for (int i = 0; i < m; i++) {
        mag += pow(matrix[i][column], 2);
    }
    return sqrt(mag);
}

/**
 * Return magnitude of vector
 * @param n the length of the vector
 * @param vec the vector whose magnitude is to be calculated
 * @return the magnitude of the vector
 */
static double vector_magnitude(int n, double vec[]) {
    double mag = 0.;
    for (int i = 0; i < n; i++) {
        mag += pow(vec[i], 2);
    }
    return sqrt(mag);
}

/**
 * Normalize vector
 * @param n the length of the vector
 * @param vec the vector to be normalized
 */
static void normalize_vector(int n, double vec[]) {
    double magnitude = vector_magnitude(n, vec);
    for (int i = 0; i < n; i++) {
        vec[i] = vec[i] / magnitude;
    }
}

/**
 * Multiply vector by scalar
 * @param n the length of the vector
 * @param scalar the multiplier
 * @param vector the multiplicand
 * @param multiplied the resulting multiplied vector
 */
static void multiply_vector_by_scalar(double n, double scalar, const double vector[], double multiplied[]) {
    for (int i = 0; i < n; i++) {
        multiplied[i] = vector[i] * scalar;
    }
}

/**
 * Return the trace of a square m, m matrix
 * @param m the dimensions of the matrix
 * @param matrix the matrix whose trace is to be returned
 * @return the trace of the square matrix
 */
static double matrix_trace(int m, double matrix[m][m]) {
    double trace = 0.;
    for (int i = 0; i < m; i++) {
        trace += matrix[i][i];
    }
    return trace;
}

/**
 * Convert the value of radians to degrees
 * @param radians
 * @return degrees
 */
static double radians_to_degrees(double radians) {
    return radians * (180.0 / M_PI);
}

/**
 * Convert the value of degrees to radians
 * @param degrees
 * @return radians
 */
static double degrees_to_radians(double degrees) {
    return degrees * (M_PI / 180.0);
}

/**
 * Find the arithmetic mean of two m, n matrices
 * @param m number of rows of the two matrices
 * @param n number of rows of the two matrices
 * @param avg the average matrix
 * @param m1 the first matrix to be averaged
 * @param m2 the second matrix to be averaged
 */
static void average_two_matrices(int m, int n, double avg[m][n], double m1[m][n], double m2[m][n]) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            avg[i][j] = (m1[i][j] + m2[i][j]) / 2;
        }
    }
}

/**
 * Find the arithmetic mean of two vectors
 * @param n the length of the two vectors
 * @param avg the average vector
 * @param vec1 the first vector to be averaged
 * @param vec2 the second vector to be averaged
 */
static void average_two_vectors(int n, double avg[n], const double vec1[n], const double vec2[n]) {
    for (int i = 0; i < n; i++) {
        avg[i] = (vec1[i] + vec2[i]) / 2;
    }
}
