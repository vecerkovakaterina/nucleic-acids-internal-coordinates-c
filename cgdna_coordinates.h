#include <stdlib.h>


void get_rotation_angle_theta_cgdna(double *theta, double rotation_matrices[][3][3], int len) {
    for (int i = 0; i < len; i++) {
        theta[i] = radians_to_degrees(2 * tan(acos((matrix_trace(3, rotation_matrices[i]) - 1) / 2.) / 2));
    }
}

void
prepare_intra_rotation_cgdna(double frames_1[][3][3], double frames_2[][3][3], double theta_a[], double u_vec_a[][3],
                             int len) {
    double (*intra_rotation_matrices)[3][3] = malloc(len * sizeof(*intra_rotation_matrices));
    curves_get_rotation_matrices(intra_rotation_matrices, frames_1, frames_2, len);
    get_rotation_angle_theta_cgdna(theta_a, intra_rotation_matrices, len);
    get_unit_rotation_vector(u_vec_a, intra_rotation_matrices, len);

    free(intra_rotation_matrices);
}

void
prepare_inter_rotation_cgdna(double intra_middle_frames[][3][3], double inter_middle_frames[][3][3], double theta_e[],
                             double u_vec_e[][3], int inter_len) {
    double (*inter_rotation_matrices)[3][3] = malloc(inter_len * sizeof(*inter_rotation_matrices));
    get_inter_rotation_matrices(inter_rotation_matrices, intra_middle_frames, inter_len);
    get_rotation_angle_theta_cgdna(theta_e, inter_rotation_matrices, inter_len);
    get_unit_rotation_vector(u_vec_e, inter_rotation_matrices, inter_len);
    get_inter_middle_frames(inter_middle_frames, intra_middle_frames, inter_len);

    free(inter_rotation_matrices);
}

void get_cgdna_coordinates(double frames_1[][3][3], double frames_2[][3][3], double origins_1[][3],
                           double origins_2[][3], double shear[], double stretch[], double stagger[], double buckle[],
                           double propeller[], double opening[], double shift[], double slide[], double rise[],
                           double roll[], double tilt[], double twist[]) {

    int strand_len = (int) (snapshot_len / 2);
    int inter_len = strand_len - 1;

    double *theta_a = malloc(strand_len * sizeof(*theta_a));
    double (*unit_rotation_vector_a)[3] = malloc(strand_len * sizeof(*unit_rotation_vector_a));
    prepare_intra_rotation_cgdna(frames_1, frames_2, theta_a, unit_rotation_vector_a, strand_len);

    double (*intra_middle_frames)[3][3] = malloc(strand_len * sizeof(*intra_middle_frames));
    double (*intra_middle_frames_origins)[3] = malloc(strand_len * sizeof(*intra_middle_frames_origins));
    get_intra_middle_frames_and_origins(intra_middle_frames, intra_middle_frames_origins, frames_1, frames_2,
                                        origins_1, origins_2, strand_len);

    get_translational_coords(shear, stagger, stretch, intra_middle_frames, origins_1, origins_2, strand_len);
    get_rotational_coords(buckle, propeller, opening, theta_a, unit_rotation_vector_a, intra_middle_frames, strand_len);

    //inter-bp
    double (*theta_e) = malloc(inter_len * sizeof(*theta_e));
    double (*unit_rotation_vector_e)[3] = malloc(strand_len * sizeof(*unit_rotation_vector_e));
    double (*inter_middle_frames)[3][3] = malloc(inter_len * sizeof(*inter_middle_frames));
    prepare_inter_rotation_cgdna(intra_middle_frames, inter_middle_frames, theta_e, unit_rotation_vector_e, inter_len);

    get_inter_translational_coords(shift, slide, rise, inter_middle_frames, intra_middle_frames_origins, inter_len);
    get_inter_rotational_coords(roll, tilt, twist, theta_e, unit_rotation_vector_e, inter_middle_frames, inter_len);

    free_curves_arrays(theta_a, unit_rotation_vector_a, intra_middle_frames,
                       intra_middle_frames_origins, theta_e, unit_rotation_vector_e,
                       inter_middle_frames);
}

#ifndef CMAKE_FINAL_THESIS_C_CGDNA_COORDINATES_H
#define CMAKE_FINAL_THESIS_C_CGDNA_COORDINATES_H

#endif //CMAKE_FINAL_THESIS_C_CGDNA_COORDINATES_H
