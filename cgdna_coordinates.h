#include <stdlib.h>

/**
 * Compute the rotation angle according to the cgDNA definition
 * @param theta the array allocated for rotation agles
 * @param rotation_matrices the array of rotation matrices
 * @param len the length of rotation angle array
 */
void get_rotation_angle_theta_cgdna(double *theta, double rotation_matrices[][3][3], int len) {
    for (int i = 0; i < len; i++) {
        theta[i] = radians_to_degrees(2 * tan(acos((matrix_trace(3, rotation_matrices[i]) - 1) / 2.) / 2));
    }
}

/**
 * Compute the rotation angle and vector describing the rotation between neighbouring base frames
 * @param frames_1 the base frames in strand 1
 * @param frames_2 the base frames in strand 1
 * @param theta_a the array allocated for rotation angles
 * @param u_vec_a the array allocated for rotation vectors
 * @param len the length of one strand
 */
void
prepare_intra_rotation_cgdna(double frames_1[][3][3], double frames_2[][3][3], double theta_a[], double u_vec_a[][3],
                             int len) {
    double (*intra_rotation_matrices)[3][3] = malloc(len * sizeof(*intra_rotation_matrices));
    curves_get_rotation_matrices(intra_rotation_matrices, frames_1, frames_2, len);
    get_rotation_angle_theta_cgdna(theta_a, intra_rotation_matrices, len);
    get_unit_rotation_vector(u_vec_a, intra_rotation_matrices, len);

    free(intra_rotation_matrices);
}

/**
 * Compute the rotation angle and vector describing the rotation between consecutive base pair frames
 * @param intra_middle_frames the array of base middle frames
 * @param inter_middle_frames the array of base pair middle frames
 * @param theta_e the array allocated for rotation angles
 * @param u_vec_e the array allocated for rotation vectors
 * @param inter_len the length of inter base pair middle frames array
 */
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

/**
 * Compute both the translational and rotational coordinates according to the cgDNA definition
 * @param frames_1 the array of base frames in strand 1
 * @param frames_2 the array of base frames in strand 2
 * @param origins_1 the array of base origins in strand 1
 * @param origins_2 the array of base origins in strand 2
 * @param shear the array allocated for shear values
 * @param stretch the array allocated for stretch values
 * @param stagger the array allocated for stagger values
 * @param buckle the array allocated for buckle values
 * @param propeller the array allocated for propeller values
 * @param opening the array allocated for opening values
 * @param shift the array allocated for shift values
 * @param slide the array allocated for slide values
 * @param rise the array allocated for rise values
 * @param roll the array allocated for roll values
 * @param tilt the array allocated for tilt values
 * @param twist the array allocated for twist values
 */
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
