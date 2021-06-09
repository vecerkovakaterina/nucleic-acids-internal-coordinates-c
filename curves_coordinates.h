#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/**
 * Compute rotation matrices describing rotation between neighbouring base frames
 * @param intra_rotation_matrices the array allocated for rotation matrices
 * @param frames_1 the array of base frames in strand 1
 * @param frames_2 the array of base frames in strand 2
 * @param strand_len the length of one strand
 */
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

/**
 * Orthonormalize an array of frame matrices
 * @param frame the array of frames to be orthonormalized
 * @param len the length of one strand
 */
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

/**
 * Compute the rotation angle in degrees
 * @param theta the array allocated for theta values
 * @param rotation_matrices the array of rotation matrices
 * @param len the length of one strand
 */
void get_rotation_angle_theta(double *theta, double rotation_matrices[][3][3], int len) {
    for (int i = 0; i < len; i++) {
        theta[i] = radians_to_degrees(acos((matrix_trace(3, rotation_matrices[i]) - 1) / 2.));
    }
}

/**
 * Compute the rotation vector
 * @param u_a the array allocated for rotation vectors
 * @param rotation_matrices the array of rotation matrices
 * @param len the length of one strand
 */
void get_unit_rotation_vector(double u_a[][3], double rotation_matrices[][3][3], int len) {
    for (int i = 0; i < len; i++) {
        u_a[i][0] = rotation_matrices[i][1][2] - rotation_matrices[i][2][1];
        u_a[i][1] = rotation_matrices[i][2][0] - rotation_matrices[i][0][2];
        u_a[i][2] = rotation_matrices[i][0][1] - rotation_matrices[i][1][0];

        normalize_vector(3, u_a[i]);
    }
}

/**
 * Compute the middle frames as averages between neighbouring base frames
 * @param middle_frames the array allocated for middle frames matrices
 * @param frames_1 the base frame matrices of strand 1
 * @param frames_2 the base frame matrices of strand 2
 * @param len the length of one strand
 */
void get_intra_middle_frames(double middle_frames[][3][3], double frames_1[][3][3], double frames_2[][3][3], int len) {
    for (int i = 0; i < len; i++) {
        average_two_matrices(3, 3, middle_frames[i], frames_1[i], frames_2[i]);
    }
}

/**
 * Compute the origins of the middle frames between neighbouring base frames
 * @param middle_frames_origins the array allocated for origin vectors
 * @param origins_1 the array of origins in strand 1
 * @param origins_2 the array of origins in strand 2
 * @param len the length of one strand
 */
void get_intra_middle_frames_origins(double middle_frames_origins[][3], double origins_1[][3], double origins_2[][3],
                                     int len) {
    for (int i = 0; i < len; i++) {
        average_two_vectors(3, middle_frames_origins[i], origins_1[i], origins_2[i]);
    }
}

/**
 * Compute the translational coordinates between neighbouring base frames
 * @param coord_1 the array allocated for shear values
 * @param coord_2 the array allocated for stagger values
 * @param coord_3 the array allocated for stretch values
 * @param middle_frames the array of middle frames between two base frames
 * @param origins_1 the array of origins in strand 1
 * @param origins_2 the array of origins in strand 2
 * @param len the length of one strand
 */
void get_translational_coords(double *shear, double *stagger, double *stretch, double middle_frames[][3][3],
                              double origins_1[][3], double origins_2[][3], int len) {
    for (int i = 0; i < len; i++) {
        double lambda[3];
        subtract_two_vectors(3, lambda, origins_2[i], origins_1[i]);
        double coords[3];
        vector_matrix_multiplication(3, 3, coords, lambda, middle_frames[i]);
        shear[i] = coords[0];
        stagger[i] = coords[1];
        stretch[i] = coords[2];
    }
}

/**
 * Compute the rotational coordinated between neighbouring base frames
 * @param coord_1 the array allocated for buckle values
 * @param coord_2 the array allocated for propeller values
 * @param coord_3 the array allocated for opening values
 * @param angle the array of rotation angles
 * @param vector the array of rotation vectors
 * @param middle_frames the array of middle frames between two base frames
 * @param len the length of one strand
 */
void get_rotational_coords(double *buckle, double *propeller, double *opening, double angle[], double vector[][3],
                           double middle_frames[][3][3], int len) {
    for (int i = 0; i < len; i++) {
        double rotation[3];
        multiply_vector_by_scalar(3, angle[i], vector[i], rotation);
        double coords[3];
        vector_matrix_multiplication(3, 3, coords, rotation, middle_frames[i]);
        buckle[i] = coords[0];
        propeller[i] = coords[1];
        opening[i] = coords[2];
    }
}

/**
 * Compute rotation matrices describing rotation between consecutive base pair frames
 * @param rotation_matrices the array allocated for rotation matrices
 * @param middle_frames the array of base pair frames
 * @param len the length of base pair middle frames array
 */
void get_inter_rotation_matrices(double rotation_matrices[][3][3], double middle_frames[][3][3], int len) {
    for (int i = 0; i < len; i++) {
        double frame_2_transposed[3][3];
        transpose_matrix(3, 3, frame_2_transposed, middle_frames[i + 1]);
        matrix_multiplication(3, 3, 3, rotation_matrices[i], middle_frames[i], frame_2_transposed);
    }
}

/**
 * Compute the base pair frames as averages of consecutive base pair frames
 * @param inter_middle_frames the array allocated for middle frames matrices
 * @param frames the array of base pair frames
 * @param len the length of base pair middle frames array
 */
void get_inter_middle_frames(double inter_middle_frames[][3][3], double frames[][3][3], int len) {
    for (int i = 0; i < len; i++) {
        average_two_matrices(3, 3, inter_middle_frames[i], frames[i], frames[i + 1]);
    }
}

/**
 * Compute the internal translational coordinates between consecutive base pair frames
 * @param shift the array allocated for shift values
 * @param slide the array allocated for slide values
 * @param rise the array allocated for rise values
 * @param middle_frames the array of base pair middle frames matrices
 * @param origins the array of base pair middle frames origins vectors
 * @param len the length of base pair middle frames array
 */
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

/**
 * Compute the internal rotational coordinates between consecutive base pair frames
 * @param roll the array for roll values
 * @param tilt the array for tilt values
 * @param twist the array for twist values
 * @param angle the array of rotation angles
 * @param vector the array of rotation vectors
 * @param middle_frames the base pair middle frames
 * @param len the length of base pair middle frames array
 */
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

/**
 * Free allocated memory
 * @param ta the pointer to rotation angle between two base frames array
 * @param ua the pointer to rotation vector between two base frames array
 * @param imf the pointer to base middle frames matrices array
 * @param imfo the pointer to base middle frames origins vectors array
 * @param te the pointer to rotation angle theta between two base pair frames array
 * @param ue the pointer to rotation vector between two base pair frames array
 * @param emf the pointer to base pair middle frames matrices array
 */
void free_curves_arrays(double ta[], double ua[][3], double imf[][3][3], double imfo[][3], double te[], double ue[][3],
                        double emf[][3][3]) {
    free(ta);
    free(ua);
    free(imf);
    free(imfo);
    free(te);
    free(ue);
    free(emf);
}

/**
 * Compute rotation angle and vector describing rotation between two base frames
 * @param frames_1 the array of base frames in strand 1
 * @param frames_2 the array of base frames in strand 2
 * @param theta_a the array allocated for rotation angles
 * @param u_vec_a the array allocated for rotation vectors
 * @param len the length of one strand
 */
void prepare_intra_rotation(double frames_1[][3][3], double frames_2[][3][3], double theta_a[], double u_vec_a[][3],
                            int len) {
    double (*intra_rotation_matrices)[3][3] = malloc(len * sizeof(*intra_rotation_matrices));
    curves_get_rotation_matrices(intra_rotation_matrices, frames_1, frames_2, len);
    get_rotation_angle_theta(theta_a, intra_rotation_matrices, len);
    get_unit_rotation_vector(u_vec_a, intra_rotation_matrices, len);

    free(intra_rotation_matrices);
}

/**
 * Compute orthonormal base middle frames and their origins
 * @param middle_frames the array allocated for base middle frames matrices
 * @param middle_origins the array allocated for base middle frames origins vectors
 * @param frames_1 the array of base frames in strand 1
 * @param frames_2 the array of base frames in strand 2
 * @param origins_1 the array of base frame origin vectors in strand 1
 * @param origins_2 the array of base frame origin vectors in strand 2
 * @param len the length of one strand
 */
void get_intra_middle_frames_and_origins(double middle_frames[][3][3], double middle_origins[][3],
                                         double frames_1[][3][3], double frames_2[][3][3], double origins_1[][3],
                                         double origins_2[][3], int len) {
    get_intra_middle_frames(middle_frames, frames_1, frames_2, len);
    get_intra_middle_frames_origins(middle_origins, origins_1, origins_2, len);
    gram_schmidt_columns(middle_frames, len);
}

/**
 * Compute the rotation angle and vector describing the rotation between consecutive base pair frames
 * @param intra_middle_frames the array of base middle frames
 * @param inter_middle_frames  the array of base pair middle frames
 * @param theta_e the array allocated for rotation angles
 * @param u_vec_e the array allocated for rotation vectors
 * @param inter_len the length of base pair middle frames array
 */
void prepare_inter_rotation(double intra_middle_frames[][3][3], double inter_middle_frames[][3][3], double theta_e[],
                            double u_vec_e[][3], int inter_len) {
    double (*inter_rotation_matrices)[3][3] = malloc(inter_len * sizeof(*inter_rotation_matrices));
    get_inter_rotation_matrices(inter_rotation_matrices, intra_middle_frames, inter_len);
    get_rotation_angle_theta(theta_e, inter_rotation_matrices, inter_len);
    get_unit_rotation_vector(u_vec_e, inter_rotation_matrices, inter_len);
    get_inter_middle_frames(inter_middle_frames, intra_middle_frames, inter_len);

    free(inter_rotation_matrices);
}

/**
 * Compute both translational and rotational coordinates according to the Curves+ definition
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
void get_curves_coordinates(double frames_1[][3][3], double frames_2[][3][3], double origins_1[][3],
                            double origins_2[][3], double shear[], double stretch[], double stagger[], double buckle[],
                            double propeller[], double opening[], double shift[], double slide[], double rise[],
                            double roll[], double tilt[], double twist[]) {

    int strand_len = (int) (snapshot_len / 2);
    int inter_len = strand_len - 1;

    double *theta_a = malloc(strand_len * sizeof(*theta_a));
    double (*unit_rotation_vector_a)[3] = malloc(strand_len * sizeof(*unit_rotation_vector_a));
    prepare_intra_rotation(frames_1, frames_2, theta_a, unit_rotation_vector_a, strand_len);

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
    prepare_inter_rotation(intra_middle_frames, inter_middle_frames, theta_e, unit_rotation_vector_e, inter_len);

    get_inter_translational_coords(shift, slide, rise, inter_middle_frames, intra_middle_frames_origins, inter_len);
    get_inter_rotational_coords(roll, tilt, twist, theta_e, unit_rotation_vector_e, inter_middle_frames, inter_len);



    free_curves_arrays(theta_a, unit_rotation_vector_a, intra_middle_frames,
                       intra_middle_frames_origins, theta_e, unit_rotation_vector_e,
                       inter_middle_frames);
}

