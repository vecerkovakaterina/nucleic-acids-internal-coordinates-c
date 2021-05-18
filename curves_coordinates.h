#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void curves_get_rotation_matrices(double intra_rotation_matrices[][3][3],
                                  double frames_1[][3][3],
                                  double frames_2[][3][3],
                                  int strand_len) {

    for (int i = 0; i < strand_len; i++) {
        double frame_2_transposed[3][3];
        transpose_matrix(3, 3, frame_2_transposed, frames_2[i]);
        matrix_multiplication(3, 3, 3, intra_rotation_matrices[i], frames_1[i], frame_2_transposed);
    }
}

void gram_schmidt_columns(double frame[][3][3], int len) {

    for (int i = 0; i < len; i++) {
        double transposed[3][3];
        transpose_matrix(3, 3, transposed, frame[i]);

        //normalize x
        normalize_vector(3, transposed[0]);

        //find y orthogonal to x
        double factor1[3];
        multiply_vector_by_scalar(3, dot_product(3, transposed[1], transposed[0]), transposed[0], factor1);
        for (int j = 0; j < 3; j++) {
            transposed[1][j] = transposed[1][j] - factor1[j];
        }

        //normalize y
        normalize_vector(3, transposed[1]);

        //find z orthogonal to x and y
        double factor2[3];
        multiply_vector_by_scalar(3, dot_product(3, transposed[2], transposed[0]), transposed[0], factor2);
        double factor3[3];
        multiply_vector_by_scalar(3, dot_product(3, transposed[2], transposed[1]), transposed[1], factor3);
        for (int j = 0; j < 3; j++) {
            transposed[2][j] = transposed[2][j] - (factor2[j] + factor3[j]);
        }

        //normalize z
        normalize_vector(3, transposed[2]);

        //transpose and copy back to frame
        transpose_matrix(3, 3, frame[i], transposed);
    }

}

void get_rotation_angle_theta(double *theta, double rotation_matrices[][3][3], int len) {
    for (int i = 0; i < len; i++) {
        theta[i] = radians_to_degrees(acos((matrix_trace(3, rotation_matrices[i]) - 1) / 2));
    }
}

void get_unit_rotation_vector(double u_a[][3], double rotation_matrices[][3][3], int len) {
    for (int i = 0; i < len; i++) {
        u_a[i][0] = rotation_matrices[i][1][2] - rotation_matrices[i][2][1];
        u_a[i][1] = rotation_matrices[i][2][0] - rotation_matrices[i][0][2];
        u_a[i][2] = rotation_matrices[i][0][1] - rotation_matrices[i][1][0];

        normalize_vector(3, u_a[i]);
    }
}

void get_intra_middle_frames(double middle_frames[][3][3], double frames_1[][3][3], double frames_2[][3][3], int len) {
    for (int i = 0; i < len; i++) {
        average_two_matrices(3, 3, middle_frames[i], frames_1[i], frames_2[i]);
    }
}

void get_intra_middle_frames_origins(double middle_frames_origins[][3], double origins_1[][3], double origins_2[][3],
                                     int len) {
    for (int i = 0; i < len; i++) {
        average_two_vectors(3, middle_frames_origins[i], origins_1[i], origins_2[i]);
    }
}

void get_translational_coords(double *coord_1, double *coord_2, double *coord_3, double middle_frames[][3][3],
                              double origins_1[][3], double origins_2[][3], int len) {
    for (int i = 0; i < len; i++) {
        double lambda[3];
        subtract_two_vectors(3, lambda, origins_2[i], origins_1[i]);
        double coords[3];
        vector_matrix_multiplication(3, 3, coords, lambda, middle_frames[i]);
        coord_1[i] = coords[0];
        coord_2[i] = coords[1];
        coord_3[i] = coords[2];
    }
}

void get_rotational_coords(double *coord_1, double *coord_2, double *coord_3, double angle[], double vector[][3],
                           double middle_frames[][3][3], int len) {
    for (int i = 0; i < len; i++) {
        double rotation[3];
        multiply_vector_by_scalar(3, angle[i], vector[i], rotation);
        double coords[3];
        vector_matrix_multiplication(3, 3, coords, rotation, middle_frames[i]);
        coord_1[i] = coords[0];
        coord_2[i] = coords[1];
        coord_3[i] = coords[2];
    }
}

void get_inter_rotation_matrices(double rotation_matrices[][3][3], double middle_frames[][3][3], int len) {
    for (int i = 0; i < len; i++) {
        double frame_2_transposed[3][3];
        transpose_matrix(3, 3, frame_2_transposed, middle_frames[i + 1]);
        matrix_multiplication(3, 3, 3, rotation_matrices[i], middle_frames[i], frame_2_transposed);
    }
}

void get_inter_middle_frames(double inter_middle_frames[][3][3], double frames[][3][3], int len) {
    for (int i = 0; i < len; i++) {
        average_two_matrices(3, 3, inter_middle_frames[i], frames[i], frames[i + 1]);
    }
}

void get_inter_middle_frames_origins(double inter_middle_frames_origins[][3], double origins[][3], int len) {
    for (int i = 0; i < len; i++) {
        average_two_vectors(3, inter_middle_frames_origins[i], origins[i], origins[i + 1]);
    }
}

void get_inter_translational_coords(double *shift, double *slide, double *rise, double middle_frames[][3][3],
                                    double origins[][3], int len) {
    for (int i = 0; i < len; i++) {
        double lambda[3];
        subtract_two_vectors(3, lambda, origins[i + 1], origins[i]);
        double coords[3];
        vector_matrix_multiplication(3, 3, coords, lambda, middle_frames[i]);
        shift[i] = coords[0];
        slide[i] = coords[1];
        rise[i] = coords[2];
    }
}

void get_inter_rotational_coords(double *roll, double *tilt, double *twist, const double *angle, double vector[][3],
                                 double middle_frames[][3][3], int len) {
    for (int i = 0; i < len; i++) {
        double rotation[3];
        multiply_vector_by_scalar(3, angle[i], vector[i], rotation);
        double coords[3];
        vector_matrix_multiplication(3, 3, coords, rotation, middle_frames[i]);
        roll[i] = coords[0];
        tilt[i] = coords[1];
        twist[i] = coords[2];
    }
}

void free_curves_arrays(double irm[][3][3], double ta[], double ua[][3], double imf[][3][3], double imfo[][3],
                        double erm[][3][3], double te[], double ue[][3], double emf[][3][3], double emfo[][3]) {
    free(irm);
    free(ta);
    free(ua);
    free(imf);
    free(imfo);
    free(erm);
    free(te);
    free(ue);
    free(emf);
    free(emfo);
}

//todo separate into functinos
void get_curves_coordinates(double frames_1[][3][3], double frames_2[][3][3], double origins_1[][3],
                            double origins_2[][3], double shear[], double stretch[], double stagger[], double buckle[],
                            double propeller[], double opening[], double shift[], double slide[], double rise[],
                            double roll[], double tilt[], double twist[]) {

    int strand_len = (int) (snapshot_len / 2);
    int inter_len = strand_len - 1;

    double (*intra_rotation_matrices)[3][3] = malloc(strand_len * sizeof(*intra_rotation_matrices));
    curves_get_rotation_matrices(intra_rotation_matrices, frames_1, frames_2, strand_len);

    double *theta_a = malloc(strand_len * sizeof(*theta_a));
    get_rotation_angle_theta(theta_a, intra_rotation_matrices, strand_len);

    double (*unit_rotation_vector_a)[3] = malloc(strand_len * sizeof(*unit_rotation_vector_a));
    get_unit_rotation_vector(unit_rotation_vector_a, intra_rotation_matrices, strand_len);

    double (*intra_middle_frames)[3][3] = malloc(strand_len * sizeof(*intra_middle_frames));
    get_intra_middle_frames(intra_middle_frames, frames_1, frames_2, strand_len);
    double (*intra_middle_frames_origins)[3] = malloc(strand_len * sizeof(*intra_middle_frames_origins));
    get_intra_middle_frames_origins(intra_middle_frames_origins, origins_1, origins_2, strand_len);

    //orthonormalize intra-bp middle frames
    gram_schmidt_columns(intra_middle_frames, strand_len);

    get_translational_coords(shear, stagger, stretch, intra_middle_frames, origins_1, origins_2, strand_len);
    get_rotational_coords(buckle, propeller, opening, theta_a, unit_rotation_vector_a, intra_middle_frames, strand_len);

    //inter-bp
    double (*inter_rotation_matrices)[3][3] = malloc(inter_len * sizeof(*inter_rotation_matrices));
    get_inter_rotation_matrices(inter_rotation_matrices, intra_middle_frames, inter_len);
    double (*theta_e) = malloc(inter_len * sizeof(*theta_e));
    get_rotation_angle_theta(theta_e, inter_rotation_matrices, inter_len);
    double (*unit_rotation_vector_e)[3] = malloc(strand_len * sizeof(*unit_rotation_vector_e));
    get_unit_rotation_vector(unit_rotation_vector_e, inter_rotation_matrices, inter_len);
    double (*inter_middle_frames)[3][3] = malloc(inter_len * sizeof(*inter_middle_frames));
    get_inter_middle_frames(inter_middle_frames, intra_middle_frames, inter_len);
    double (*inter_middle_frames_origins)[3] = malloc(inter_len * sizeof(*inter_middle_frames_origins));
    get_inter_middle_frames_origins(inter_middle_frames_origins, intra_middle_frames_origins, inter_len);

    get_inter_translational_coords(shift, slide, rise, inter_middle_frames, intra_middle_frames_origins, inter_len);
    get_inter_rotational_coords(roll, tilt, twist, theta_e, unit_rotation_vector_e, inter_middle_frames, inter_len);

    free_curves_arrays(intra_rotation_matrices, theta_a, unit_rotation_vector_a, intra_middle_frames,
                       intra_middle_frames_origins, inter_rotation_matrices, theta_e, unit_rotation_vector_e,
                       inter_middle_frames, inter_middle_frames_origins);
}

#ifdef CMAKE_FINAL_THESIS_C_CURVES_COORDINATES_H
#endif //CMAKE_FINAL_THESIS_C_CURVES_COORDINATES_H
