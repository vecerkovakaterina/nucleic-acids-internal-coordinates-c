#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lh3.h"
#include "lingebra.h"
//todo change strand 1 and strand 2
//todo sort snapshot array length
#define SNAPSHOT_ARRAY_LENGTH 100
int snapshot_len = 0;

typedef struct Purine {
    double N9[3];
    double C8[3];
    double N7[3];
    double C5[3];
    double C6[3];
    double N1[3];
    double C2[3];
    double N3[3];
    double C4[3];
} purine;

typedef struct Pyrimidine {
    double N1[3];
    double C2[3];
    double N3[3];
    double C4[3];
    double C5[3];
    double C6[3];
} pyrimidine;

typedef struct Residuum {
    char nucleotide;
    int is_purine;
    int is_pyrimidine;
    purine purine;
    pyrimidine pyrimidine;
} residuum;

typedef struct Eig {
    double eigen_vector[4];
    double eigen_value;
} eig;

void read_purine(residuum *r, const char res_name[], char atom_name[], double x, double y, double z) {
    r->nucleotide = res_name[1];
    r->is_purine = 1;
    r->is_pyrimidine = 0;
    if (strcmp(atom_name, "N9") == 0) {
        r->purine.N9[0] = x;
        r->purine.N9[2] = z;
        r->purine.N9[1] = y;
    } else if (strcmp(atom_name, "C8") == 0) {
        r->purine.C8[0] = x;
        r->purine.C8[1] = y;
        r->purine.C8[2] = z;
    } else if (strcmp(atom_name, "N7") == 0) {
        r->purine.N7[0] = x;
        r->purine.N7[1] = y;
        r->purine.N7[2] = z;
    } else if (strcmp(atom_name, "C5") == 0) {
        r->purine.C5[0] = x;
        r->purine.C5[1] = y;
        r->purine.C5[2] = z;
    } else if (strcmp(atom_name, "C6") == 0) {
        r->purine.C6[0] = x;
        r->purine.C6[1] = y;
        r->purine.C6[2] = z;
    } else if (strcmp(atom_name, "N1") == 0) {
        r->purine.N1[0] = x;
        r->purine.N1[1] = y;
        r->purine.N1[2] = z;
    } else if (strcmp(atom_name, "C2") == 0) {
        r->purine.C2[0] = x;
        r->purine.C2[1] = y;
        r->purine.C2[2] = z;
    } else if (strcmp(atom_name, "N3") == 0) {
        r->purine.N3[0] = x;
        r->purine.N3[1] = y;
        r->purine.N3[2] = z;
    } else if (strcmp(atom_name, "C4") == 0) {
        r->purine.C4[0] = x;
        r->purine.C4[1] = y;
        r->purine.C4[2] = z;
        snapshot_len++;
    }
}

void read_pyrimidine(residuum *r, const char res_name[], char atom_name[], double x, double y, double z) {
    r->nucleotide = res_name[1];
    r->is_purine = 0;
    r->is_pyrimidine = 1;
    if (strcmp(atom_name, "N1") == 0) {
        r->pyrimidine.N1[0] = x;
        r->pyrimidine.N1[1] = y;
        r->pyrimidine.N1[2] = z;
    } else if (strcmp(atom_name, "C2") == 0) {
        r->pyrimidine.C2[0] = x;
        r->pyrimidine.C2[1] = y;
        r->pyrimidine.C2[2] = z;
        snapshot_len++;
    } else if (strcmp(atom_name, "N3") == 0) {
        r->pyrimidine.N3[0] = x;
        r->pyrimidine.N3[1] = y;
        r->pyrimidine.N3[2] = z;
    } else if (strcmp(atom_name, "C4") == 0) {
        r->pyrimidine.C4[0] = x;
        r->pyrimidine.C4[1] = y;
        r->pyrimidine.C4[2] = z;
    } else if (strcmp(atom_name, "C5") == 0) {
        r->pyrimidine.C5[0] = x;
        r->pyrimidine.C5[1] = y;
        r->pyrimidine.C5[2] = z;
    } else if (strcmp(atom_name, "C6") == 0) {
        r->pyrimidine.C6[0] = x;
        r->pyrimidine.C6[1] = y;
        r->pyrimidine.C6[2] = z;
    }
}

void read_file(FILE *fp, residuum snapshot[]) {

    char line[256];
    char atom[5];
    int no_atom;
    char atom_name[5];
    char res_name[5];
    int res_number = 1;
    double x, y, z, a, b;
    char atom_symbol;

    int residuum_number = 0;

    while (!feof(fp)) {
        do {
            fgets(line, sizeof(line), fp);
            sscanf(line, " %s%d%s%s%d%lf%lf%lf%lf%lf%c", atom, &no_atom, atom_name, res_name, &res_number, &x, &y, &z,
                   &a, &b, &atom_symbol);

            if (res_name[1] == 'A' || res_name[1] == 'G') {
                read_purine(&snapshot[snapshot_len], res_name, atom_name, x, y, z);
            } else if (res_name[1] == 'C' || res_name[1] == 'T' || res_name[1] == 'U') {
                read_pyrimidine(&snapshot[snapshot_len], res_name, atom_name, x, y, z);
            }
        } while ((residuum_number == res_number) && !feof(fp));
        residuum_number++;
    }
}

void create_standard_frame_G(purine *g) {
    g->N9[0] = -1.289;
    g->N9[1] = 4.551;
    g->C8[0] = 0.023;
    g->C8[1] = 4.962;
    g->N7[0] = 0.870;
    g->N7[1] = 3.969;
    g->C5[0] = 0.071;
    g->C5[1] = 2.833;
    g->C6[0] = 0.424;
    g->C6[1] = 1.460;
    g->N1[0] = -0.700;
    g->N1[1] = 0.641;
    g->C2[0] = -1.999;
    g->C2[1] = 1.087;
    g->N3[0] = -2.342;
    g->N3[1] = 2.364;
    g->C4[0] = -1.265;
    g->C4[1] = 3.177;
}

void create_standard_frame_A(purine *a) {
    a->N9[0] = -1.291;
    a->N9[1] = 4.498;
    a->C8[0] = 0.024;
    a->C8[1] = 4.897;
    a->N7[0] = 0.877;
    a->N7[1] = 3.902;
    a->C5[0] = 0.071;
    a->C5[1] = 2.771;
    a->C6[0] = 0.369;
    a->C6[1] = 1.398;
    a->N1[0] = -0.668;
    a->N1[1] = 0.532;
    a->C2[0] = -1.912;
    a->C2[1] = 1.023;
    a->N3[0] = -2.320;
    a->N3[1] = 2.290;
    a->C4[0] = -1.267;
    a->C4[1] = 3.124;
}

void create_standard_frame_C(pyrimidine *c) {
    c->N1[0] = -1.285;
    c->N1[1] = 4.542;
    c->C2[0] = -1.472;
    c->C2[1] = 3.158;
    c->N3[0] = -0.391;
    c->N3[1] = 2.344;
    c->C4[0] = 0.837;
    c->C4[1] = 2.868;
    c->C5[0] = 1.056;
    c->C5[1] = 4.275;
    c->C6[0] = -0.023;
    c->C6[1] = 5.068;
}

void create_standard_frame_T(pyrimidine *t) {
    t->N1[0] = -1.284;
    t->N1[1] = 4.500;
    t->C2[0] = -1.462;
    t->C2[1] = 3.135;
    t->N3[0] = -0.298;
    t->N3[1] = 2.407;
    t->C4[0] = 0.994;
    t->C4[1] = 2.897;
    t->C5[0] = 1.106;
    t->C5[1] = 4.338;
    t->C6[0] = -0.024;
    t->C6[1] = 5.057;
}

void create_standard_frame_U(pyrimidine *u) {
    u->N1[0] = -1.284;
    u->N1[1] = 4.500;
    u->C2[0] = -1.462;
    u->C2[1] = 3.131;
    u->N3[0] = -0.302;
    u->N3[1] = 2.397;
    u->C4[0] = 0.989;
    u->C4[1] = 2.884;
    u->C5[0] = 1.089;
    u->C5[1] = 4.311;
    u->C6[0] = -0.024;
    u->C6[1] = 5.053;
}

void create_standard_frames(purine *G, purine *A, pyrimidine *C, pyrimidine *T, pyrimidine *U) {
    create_standard_frame_G(G);
    create_standard_frame_A(A);
    create_standard_frame_C(C);
    create_standard_frame_T(T);
    create_standard_frame_U(U);
}

double *get_atom_from_purine_struct(purine *r, int i) {
    switch (i) {
        case 0:
            return r->N9;
        case 1:
            return r->C8;
        case 2:
            return r->N7;
        case 3:
            return r->C5;
        case 4:
            return r->C6;
        case 5:
            return r->N1;
        case 6:
            return r->C2;
        case 7:
            return r->N3;
        case 8:
            return r->C4;
        default:
            return r->N9;
    }
}

double *get_atom_from_pyrimidine_struct(pyrimidine *r, int i) {
    switch (i) {
        case 0:
            return r->N1;
        case 1:
            return r->C2;
        case 2:
            return r->N3;
        case 3:
            return r->C4;
        case 4:
            return r->C5;
        case 5:
            return r->C6;
        default:
            return r->N1;
    }
}

void matrix_from_purine_residuum(double matrix[9][3], purine residuum) {
    for (int i = 0; i < 9; i++) {
        double *atomic_coords = get_atom_from_purine_struct(&residuum, i);
        for (int j = 0; j < 3; j++) {
            matrix[i][j] = atomic_coords[j];
        }
    }
}

void matrix_from_pyrimidine_residuum(double matrix[6][3], pyrimidine residuum) {
    for (int i = 0; i < 6; i++) {
        double *atomic_coords = get_atom_from_pyrimidine_struct(&residuum, i);
        for (int j = 0; j < 3; j++) {
            matrix[i][j] = atomic_coords[j];
        }
    }
}

void subtract_matrices(int m, int n, double res[m][n], double m1[m][n], double m2[m][n]) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            res[i][j] = m1[i][j] - m2[i][j];
        }
    }
}

// todo separate into functions
void create_covariance_matrix_purine(double cm[3][3], purine experimental, purine standard) {
    double ones[9][1] = {{1},
                         {1},
                         {1},
                         {1},
                         {1},
                         {1},
                         {1},
                         {1},
                         {1}};

    int no_atoms = 9;
    double standard_matrix[9][3];
    double standard_matrix_t[3][9];
    matrix_from_purine_residuum(standard_matrix, standard);
    transpose_matrix(9, 3, standard_matrix_t, standard_matrix);

    double experimental_matrix[9][3];
    matrix_from_purine_residuum(experimental_matrix, experimental);

    double vec1[3];
    get_sum_of_columns_matrix(9, 3, vec1, experimental_matrix);
    double vec2[9][3];
    spread_vector_into_m_rows(9, 3, vec2, ones, vec1);
    double m1[3][3] = {{0, 0, 0},
                       {0, 0, 0},
                       {0, 0, 0}};
    matrix_multiplication(3, 9, 3, m1, standard_matrix_t, vec2);
    multiply_matrix_by_scalar(3, 3, m1, (1. / no_atoms), m1);
    double m2[3][3] = {{0, 0, 0},
                       {0, 0, 0},
                       {0, 0, 0}};
    matrix_multiplication(3, 9, 3, m2, standard_matrix_t, experimental_matrix);
    double m[3][3];
    subtract_matrices(3, 3, m, m2, m1);
    multiply_matrix_by_scalar(3, 3, cm, (1. / (no_atoms - 1)), m);

    //(1/no_atoms-1)(m2 - (1/no_atoms)*m1)
}

// todo separate into functions
void create_covariance_matrix_pyrimidine(double cm[3][3], pyrimidine experimental, pyrimidine standard) {
    double ones[6][1] = {{1},
                         {1},
                         {1},
                         {1},
                         {1},
                         {1}};

    int no_atoms = 6;
    double standard_matrix[6][3];
    double standard_matrix_t[3][6];
    matrix_from_pyrimidine_residuum(standard_matrix, standard);
    transpose_matrix(6, 3, standard_matrix_t, standard_matrix);

    double experimental_matrix[6][3];
    matrix_from_pyrimidine_residuum(experimental_matrix, experimental);

    double vec1[3];
    get_sum_of_columns_matrix(6, 3, vec1, experimental_matrix);
    double vec2[6][3];
    spread_vector_into_m_rows(6, 3, vec2, ones, vec1);
    double m1[3][3] = {{0, 0, 0},
                       {0, 0, 0},
                       {0, 0, 0}};
    matrix_multiplication(3, 6, 3, m1, standard_matrix_t, vec2);
    multiply_matrix_by_scalar(3, 3, m1, (1. / no_atoms), m1);
    double m2[3][3] = {{0, 0, 0},
                       {0, 0, 0},
                       {0, 0, 0}};
    matrix_multiplication(3, 6, 3, m2, standard_matrix_t, experimental_matrix);
    double m[3][3];
    subtract_matrices(3, 3, m, m2, m1);
    multiply_matrix_by_scalar(3, 3, cm, (1. / (no_atoms - 1)), m);


    //(1/no_atoms-1)(m2 - (1/no_atoms)*m1)
}

void get_covariance_matrices(double cm[][3][3], residuum snapshot[], int len, purine std_A, purine std_G,
                             pyrimidine std_C, pyrimidine std_T, pyrimidine std_U) {
    for (int i = 0; i < len; i++) {
        if (snapshot[i].is_purine) {
            if (snapshot[i].nucleotide == 'A') {
                create_covariance_matrix_purine(cm[i], snapshot[i].purine, std_A);
                continue;
            } else if (snapshot[i].nucleotide == 'G') {
                create_covariance_matrix_purine(cm[i], snapshot[i].purine, std_G);
                continue;
            }
        } else if (snapshot[i].is_pyrimidine) {
            if (snapshot[i].nucleotide == 'C') {
                create_covariance_matrix_pyrimidine(cm[i], snapshot[i].pyrimidine, std_C);
                continue;
            } else if (snapshot[i].nucleotide == 'T') {
                create_covariance_matrix_pyrimidine(cm[i], snapshot[i].pyrimidine, std_T);
                continue;
            } else if (snapshot[i].nucleotide == 'U') {
                create_covariance_matrix_pyrimidine(cm[i], snapshot[i].pyrimidine, std_U);
                continue;
            }
        }
    }
}

void get_symmetric_matrices(double sm[][4][4], double cm[][3][3], int len) {
    for (int i = 0; i < len; i++) {
        sm[i][0][0] = cm[i][0][0] + cm[i][1][1] + cm[i][2][2];
        sm[i][0][1] = cm[i][1][2] - cm[i][2][1];
        sm[i][0][2] = cm[i][2][0] - cm[i][0][2];
        sm[i][0][3] = cm[i][0][1] - cm[i][1][0];

        sm[i][1][0] = cm[i][1][2] - cm[i][2][1];
        sm[i][1][1] = cm[i][0][0] - cm[i][1][1] - cm[i][2][2];
        sm[i][1][2] = cm[i][0][1] + cm[i][1][0];
        sm[i][1][3] = cm[i][2][0] + cm[i][0][2];

        sm[i][2][0] = cm[i][2][0] - cm[i][0][2];
        sm[i][2][1] = cm[i][0][1] + cm[i][1][0];
        sm[i][2][2] = -cm[i][0][0] + cm[i][1][1] - cm[i][2][2];
        sm[i][2][3] = cm[i][1][2] + cm[i][2][1];

        sm[i][3][0] = cm[i][0][1] - cm[i][1][0];
        sm[i][3][1] = cm[i][2][0] + cm[i][0][2];
        sm[i][3][2] = cm[i][1][2] + cm[i][2][1];
        sm[i][3][3] = -cm[i][0][0] - cm[i][1][1] + cm[i][2][2];
    }
}

int find_index_largest_eigenvalue(const double eigenvalues[4]) {
    int index = 0;
    double max = eigenvalues[0];
    for (int i = 0; i < 4; i++) {
        if (eigenvalues[i] > max) {
            max = eigenvalues[i];
            index = i;
        }
    }
    return index;
}

void find_largest_eigenvalue(eig *e, double sm[4][4]) {
    double *mem, *eigvalue_r, *eigvalue_i, *eigvec;
    mem = (double *) calloc(4 * 4 + 10, sizeof(double));
    eigvec = mem;
    eigvalue_r = eigvec + 16;
    eigvalue_i = eigvalue_r + 4;

    n_eigeng(sm[0], 4, eigvalue_r, eigvalue_i, eigvec);
    int index = find_index_largest_eigenvalue(eigvalue_r);
    e->eigen_value = eigvalue_r[index];
    e->eigen_vector[0] = -eigvec[0 + index];
    e->eigen_vector[1] = -eigvec[4 + index];
    e->eigen_vector[2] = -eigvec[8 + index];
    e->eigen_vector[3] = -eigvec[12 + index];
}

void get_rotation_matrices(double rm[][3][3], double sm[][4][4], int len) {
    eig eig;
    for (int i = 0; i < len; i++) {
        find_largest_eigenvalue(&eig, sm[i]);
        rm[i][0][0] = eig.eigen_vector[0] * eig.eigen_vector[0] + eig.eigen_vector[1] * eig.eigen_vector[1]
                      - eig.eigen_vector[2] * eig.eigen_vector[2] - eig.eigen_vector[3] * eig.eigen_vector[3];
        rm[i][0][1] = 2 * (eig.eigen_vector[1] * eig.eigen_vector[2] - eig.eigen_vector[0] * eig.eigen_vector[3]);
        rm[i][0][2] = 2 * (eig.eigen_vector[1] * eig.eigen_vector[3] + eig.eigen_vector[0] * eig.eigen_vector[2]);

        rm[i][1][0] = 2 * (eig.eigen_vector[2] * eig.eigen_vector[1] + eig.eigen_vector[0] * eig.eigen_vector[3]);
        rm[i][1][1] = eig.eigen_vector[0] * eig.eigen_vector[0] - eig.eigen_vector[1] * eig.eigen_vector[1]
                      + eig.eigen_vector[2] * eig.eigen_vector[2] - eig.eigen_vector[3] * eig.eigen_vector[3];
        rm[i][1][2] = 2 * (eig.eigen_vector[2] * eig.eigen_vector[3] - eig.eigen_vector[0] * eig.eigen_vector[1]);

        rm[i][2][0] = 2 * (eig.eigen_vector[3] * eig.eigen_vector[1] - eig.eigen_vector[0] * eig.eigen_vector[2]);
        rm[i][2][1] = 2 * (eig.eigen_vector[3] * eig.eigen_vector[2] + eig.eigen_vector[0] * eig.eigen_vector[1]);
        rm[i][2][2] = eig.eigen_vector[0] * eig.eigen_vector[0] - eig.eigen_vector[1] * eig.eigen_vector[1]
                      - eig.eigen_vector[2] * eig.eigen_vector[2] + eig.eigen_vector[3] * eig.eigen_vector[3];
    }
}

//todo separate into functions
void get_rotation_matrices_origins(double origins[][3], residuum snapshot[], double rotation_matrices[][3][3],
                                   int len, purine std_A, purine std_G, pyrimidine std_C, pyrimidine std_T,
                                   pyrimidine std_U) {
    for (int i = 0; i < len; i++) {
        double experimental_avg[3];
        double standard_avg[3];
        double rotation_matrix_t[3][3];
        double s_ave_dot[3];

        if (snapshot[i].is_purine) {
            double experimental_matrix_purine[9][3];
            double standard_matrix_purine[9][3];
            matrix_from_purine_residuum(experimental_matrix_purine, snapshot[i].purine);
            if (snapshot[i].nucleotide == 'G') {
                matrix_from_purine_residuum(standard_matrix_purine, std_G);
            } else if (snapshot[i].nucleotide == 'A') {
                matrix_from_purine_residuum(standard_matrix_purine, std_A);
            }
            average_matrix_columns(9, 3, experimental_avg, experimental_matrix_purine);
            average_matrix_columns(9, 3, standard_avg, standard_matrix_purine);

        } else if (snapshot[i].is_pyrimidine) {
            double experimental_matrix_pyrimidine[9][3];
            double standard_matrix_pyrimidine[9][3];
            matrix_from_pyrimidine_residuum(experimental_matrix_pyrimidine, snapshot[i].pyrimidine);
            if (snapshot[i].nucleotide == 'C') {
                matrix_from_pyrimidine_residuum(standard_matrix_pyrimidine, std_C);
            } else if (snapshot[i].nucleotide == 'T') {
                matrix_from_pyrimidine_residuum(standard_matrix_pyrimidine, std_T);
            } else if (snapshot[i].nucleotide == 'U') {
                matrix_from_pyrimidine_residuum(standard_matrix_pyrimidine, std_U);
            }
            average_matrix_columns(6, 3, experimental_avg, experimental_matrix_pyrimidine);
            average_matrix_columns(6, 3, standard_avg, standard_matrix_pyrimidine);
        }

        transpose_matrix(3, 3, rotation_matrix_t, rotation_matrices[i]);
        matrix_vector_multiplication(3, 3, s_ave_dot, rotation_matrix_t, standard_avg);
        subtract_two_vectors(3, origins[i], experimental_avg, s_ave_dot);
    }
}

void split_snapshot_into_strands(double snapshot[][3][3], double snapshot_origins[][3], double frames_1[][3][3],
                                 double frames_2[][3][3], double origins_1[][3], double origins_2[][3], int len) {
    memcpy(frames_1, snapshot, len * sizeof(*frames_1));
    memcpy(origins_1, snapshot_origins, len * sizeof(*origins_1));

    int len_strand_2 = len;
    for (int i = 0; i < len; i++) {
        memcpy(frames_2[i], snapshot[len_strand_2 + len - 1], sizeof(*frames_2));
        memcpy(origins_2[i], snapshot_origins[len_strand_2 + len - 1], sizeof(*origins_2));
        len_strand_2--;
    }
}

void rotate_strand_2_x_180_deg(double rotated[][3][3], double to_rotate[][3][3], int len) {
    double x_rot_180[3][3] = {{1., 0.,  0.},
                              {0., -1., 0.},
                              {0., 0.,  -1.}};
    for (int i = 0; i < len; i++) {
        matrix_multiplication(3, 3, 3, rotated[i], to_rotate[i], x_rot_180);
    }
}

#ifndef CMAKE_FINAL_THESIS_C_FRAME_FITTING_H
#define CMAKE_FINAL_THESIS_C_FRAME_FITTING_H

#endif //CMAKE_FINAL_THESIS_C_FRAME_FITTING_H
