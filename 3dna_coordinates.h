#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/**
 * Select columns from all base frame matrices
 * @param axis_array the array allocated for the vectors
 * @param frame the base frame matrix
 * @param axis the axis char (x, y, z)
 * @param index the index of the base frame in the array
 */
void select_frames_axis(double *axis_array,
                        double frame[][3][3],
                        char axis,
                        int index) {
    int a = axis - 120;

    for (int i = 0; i < 3; i++) {
        axis_array[i] = frame[index][i][0 + a];
    }
}

/**
 * Select column from base frame matrix
 * @param axis_array the array allocated for the vector
 * @param frame the base frame
 * @param axis the axis char (x, y, z)
 */
void select_frame_axis(double *axis_array,
                       double frame[3][3],
                       char axis) {
    int a = axis - 120;

    for (int i = 0; i < 3; i++) {
        axis_array[i] = frame[i][0 + a];
    }
}

/**
 * Compute the hinge axes between neighbouring base frames
 * @param array the array allocated for the axes vectors
 * @param len the number of base frames in one strand
 * @param strand_I the base frames in strand 1
 * @param strand_II the base frames in strand 2
 */
void get_hinge_axes_intra_array(double array[][3],
                                int len,
                                double strand_I[][3][3],
                                double strand_II[][3][3]) {
    double z_axis_I[3] = {0.};
    double z_axis_II[3] = {0.};

    for (int i = 0; i < len; i++) {
        select_frames_axis(z_axis_I, strand_I, 'z', i);
        select_frames_axis(z_axis_II, strand_II, 'z', i);
        cross_product(array[i], z_axis_I, z_axis_II);
        normalize_vector(3, array[i]);
    }
}

/**
 * Compute the bucklepropeller array
 * @param len the number of base frames in one strand
 * @param strand_1 the base frames in strand 1
 * @param strand_2 the base frames in strand 2
 * @param bp_array the array allocated for the angle values
 */
void create_bucklepropeller_array(int len,
                                  double strand_1[][3][3],
                                  double strand_2[][3][3],
                                  double *bp_array) {
    double z1[3], z2[3];

    for (int i = 0; i < len; i++) {
        select_frames_axis(z1, strand_1, 'z', i);
        select_frames_axis(z2, strand_2, 'z', i);

        double angle = dot_product(3, z1, z2) / (double) (vector_magnitude(3, z1) * vector_magnitude(3, z2));
        angle = radians_to_degrees(acos(angle));
        bp_array[i] = angle;
    }
}

/**
 * Function to create rotation matrix according to the 3DNA definition
 * @param matrix the array allocated for the rotation matrix
 * @param angle the rotation angle
 * @param hinge_axis the hinge axis
 */
void create_rotation_matrix(double matrix[][3],
                            double angle,
                            const double hinge_axis[3]) {
    angle = degrees_to_radians(angle);
    matrix[0][0] = cos(angle) + (1 - cos(angle)) * hinge_axis[0] * hinge_axis[0];
    matrix[0][1] = (1 - cos(angle)) * hinge_axis[0] * hinge_axis[1] - hinge_axis[2] * sin(angle);
    matrix[0][2] = (1 - cos(angle)) * hinge_axis[0] * hinge_axis[2] + hinge_axis[1] * sin(angle);
    matrix[1][0] = (1 - cos(angle)) * hinge_axis[0] * hinge_axis[1] + hinge_axis[2] * sin(angle);
    matrix[1][1] = cos(angle) + (1 - cos(angle)) * (hinge_axis[1] * hinge_axis[1]);
    matrix[1][2] = (1 - cos(angle)) * hinge_axis[1] * hinge_axis[2] - hinge_axis[0] * sin(angle);
    matrix[2][0] = (1 - cos(angle)) * hinge_axis[0] * hinge_axis[2] - hinge_axis[1] * sin(angle);
    matrix[2][1] = (1 - cos(angle)) * hinge_axis[1] * hinge_axis[2] + hinge_axis[0] * sin(angle);
    matrix[2][2] = cos(angle) + (1 - cos(angle)) * (hinge_axis[2] * hinge_axis[2]);
}

/**
 * Function to rotate base frames
 * @param rotated_frames the array allocated for the rotated base frames
 * @param frames the array of frames to be rotated
 * @param hinge_axes the array of hinge axes
 * @param rotation_matrices the array of rotation matrices
 * @param bp_angle the array of rotation angles
 * @param multiple the multiple of the rotation angle to be used
 * @param strand_len the number of base frames in one strand
 */
void rotate_frames(double rotated_frames[][3][3],
                   double frames[][3][3],
                   double hinge_axes[][3],
                   double rotation_matrices[][3][3],
                   const double bp_angle[],
                   double multiple,
                   int strand_len) {
    for (int i = 0; i < strand_len; i++) {
        create_rotation_matrix(rotation_matrices[i], bp_angle[i] * multiple, hinge_axes[i]);
        matrix_multiplication(3, 3, 3, rotated_frames[i], rotation_matrices[i], frames[i]);
    }
}

/**
 * Function to create middle frames of neighbouring base frames
 * @param result the array allocated for the middle frames
 * @param frames_I the array of base frames in strand 1
 * @param frames_II the array of base frames in strand 1
 * @param strand_len the number of base frames in one strand
 */
void get_middle_frames_array(double result[][3][3],
                             double frames_I[][3][3],
                             double frames_II[][3][3],
                             int strand_len) {
    for (int i = 0; i < strand_len; i++) {
        average_two_matrices(3, 3, result[i], frames_I[i], frames_II[i]);
    }
}

/**
 * Function to compute origins of the base middle frames
 * @param result the array allocated for the origin vectors
 * @param origins_I the array of origins in strand 1
 * @param origins_II the array of origins in strand 2
 * @param strand_len the number of base frames in one strand
 */
void get_middle_frames_origins_array(double result[][3],
                                     double origins_I[][3],
                                     double origins_II[][3],
                                     int strand_len) {
    for (int i = 0; i < strand_len; i++) {
        average_two_vectors(3, result[i], origins_I[i], origins_II[i]);
    }
}

/**
 * Compute shear, stretch, stagger coordinates
 * @param shear the array allocated for shear
 * @param stretch the array allocated for stretch
 * @param stagger the array allocated for stagger
 * @param origins_I the array of origins in strand 1
 * @param origins_II the array of origins in strand 2
 * @param middle_frames the array of base middle frames
 * @param strand_len the number of base frames in one strand
 */
void get_shear_stretch_stagger(double shear[],
                               double stretch[],
                               double stagger[],
                               double origins_I[][3],
                               double origins_II[][3],
                               double middle_frames[][3][3],
                               int strand_len) {
    double difference_origins[3];
    double shear_stretch_stagger[3];
    for (int i = 0; i < strand_len; i++) {
        subtract_two_vectors(3, difference_origins, origins_II[i], origins_I[i]);
        matrix_vector_multiplication(3, 3, shear_stretch_stagger, middle_frames[i], difference_origins);
        shear[i] = shear_stretch_stagger[0];
        stretch[i] = shear_stretch_stagger[1];
        stagger[i] = shear_stretch_stagger[2];
    }
}

/**
 * Determine the sign of the opening/twist/phase angle
 * @param axis_1 the first axis determining the sign
 * @param axis_2 the second axis determining the sign
 * @param middle_frames_array the array of middle frames
 * @param index the index of the middle frame in the array
 * @return -1 or +1 to multiply the angle with
 */
int get_angle_sign(double axis_1[3],
                   double axis_2[3],
                   double middle_frames_array[][3][3],
                   int index) {
    double cross[3];
    cross_product(cross, axis_1, axis_2);

    double middle_frame_z[3];
    select_frames_axis(middle_frame_z, middle_frames_array, 'z', index);

    double sign = dot_product(3, cross, middle_frame_z);
    if (sign > 0) {
        return +1;
    }
    return -1;
}

/**
 * Function to compute the opening angle
 * @param opening_array the array allocated for the angles
 * @param rotated_frames_I the array of rotated base frames in strand 1
 * @param rotated_frames_II the array of rotated base frames in strand 2
 * @param strand_len the number of base frames in one strand
 * @param middle_frames_array the array of base middle frames
 */
void get_opening_angle_array(double opening_array[],
                             double rotated_frames_I[][3][3],
                             double rotated_frames_II[][3][3],
                             int strand_len,
                             double middle_frames_array[][3][3]) {
    double axis_rotated_1[3];
    double axis_rotated_2[3];
    for (int i = 0; i < strand_len; i++) {
        select_frames_axis(axis_rotated_1, rotated_frames_I, 'y', i);
        select_frames_axis(axis_rotated_2, rotated_frames_II, 'y', i);
        opening_array[i] = dot_product(3, axis_rotated_1, axis_rotated_2) /
                           (double) (vector_magnitude(3, axis_rotated_1) * vector_magnitude(3, axis_rotated_2));
        opening_array[i] = radians_to_degrees(acos(opening_array[i]));
        opening_array[i] = opening_array[i] * get_angle_sign(axis_rotated_1, axis_rotated_2, middle_frames_array, i);
    }
}

/**
 * Function to compute the phase angle
 * @param phase the array allocated for the angles
 * @param hinge_axes the array of hinge axes
 * @param middle_frames the array of middle frames
 * @param strand_len the number of base frames in one strand
 */
void get_phase_angle_array(double *phase,
                           double hinge_axes[][3],
                           double middle_frames[][3][3],
                           int strand_len) {
    double middle_frame_y[3];
    for (int i = 0; i < strand_len; i++) {
        select_frames_axis(middle_frame_y, middle_frames, 'y', i);
        phase[i] = dot_product(3, hinge_axes[i], middle_frame_y) /
                   (double) (vector_magnitude(3, hinge_axes[i]) * vector_magnitude(3, middle_frame_y));
        phase[i] = radians_to_degrees(acos(phase[i]));
        phase[i] = phase[i] * get_angle_sign(hinge_axes[i], middle_frame_y, middle_frames, i);
    }
}

/**
 * Compute the buckle and propeller or roll and tilt coordinates
 * @param buckle the array allocated for the buckle/roll
 * @param propeller the array allocated for the propeller/tilt
 * @param bp_array the array of bucklepropeller/rolltilt angles
 * @param phase_angle the array of phase angles
 * @param len the length of the arrays
 */
void get_bp_or_rt_angle_array(double *buckle,
                              double *propeller,
                              const double *bp_array,
                              double *phase_angle,
                              int len) {
    for (int i = 0; i < len; i++) {
        buckle[i] = bp_array[i] * sin(degrees_to_radians(phase_angle[i]));
        propeller[i] = bp_array[i] * cos(degrees_to_radians(phase_angle[i]));
    }
}

/**
 * Compute the hinge axes between consecutive base pair frames
 * @param hinge_axes the array allocated for hinge axes
 * @param base_pair_frames the array of base pair frames
 * @param size the length of the array
 */
void get_hinge_axes_inter_array(double hinge_axes[][3],
                                double base_pair_frames[][3][3],
                                int len) {
    double base_pair_z_1[3];
    double base_pair_z_2[3];
    for (int i = 0; i < len; i++) {
        select_frames_axis(base_pair_z_1, base_pair_frames, 'z', i);
        select_frames_axis(base_pair_z_2, base_pair_frames, 'z', i + 1);

        cross_product(hinge_axes[i], base_pair_z_1, base_pair_z_2);
        normalize_vector(3, hinge_axes[i]);
    }
}

/**
 * Function to compute the rolltilt angles
 * @param rolltilt the array allocated for the rolltilt angles
 * @param base_pair_frames the array of base pair frames
 * @param len the length of the arrays
 */
void get_rolltilt_angle_array(double *rolltilt,
                              double base_pair_frames[][3][3],
                              int len) {
    double base_pair_z_1[3];
    double base_pair_z_2[3];
    for (int i = 0; i < len; i++) {
        select_frames_axis(base_pair_z_1, base_pair_frames, 'z', i);
        select_frames_axis(base_pair_z_2, base_pair_frames, 'z', i + 1);

        rolltilt[i] = dot_product(3, base_pair_z_1, base_pair_z_2) /
                      (double) (vector_magnitude(3, base_pair_z_1) * vector_magnitude(3, base_pair_z_2));
        rolltilt[i] = radians_to_degrees(acos(rolltilt[i]));
    }
}

/**
 * Function to compute the base pair middle frames
 * @param middle_frames the array allocated for the base pair middle frames
 * @param base_pair_frames the array of base pair frames
 * @param hinge_axes the array of hinge axes
 * @param rolltilt the array of rolltilt angles
 * @param len the length of the base pair middle frames array
 */
void get_middle_base_pair_frames_array(double middle_frames[][3][3],
                                       double base_pair_frames[][3][3],
                                       double hinge_axes[][3],
                                       const double rolltilt[],
                                       int len) {
    double rotation_matrix_1[3][3];
    double rotation_matrix_2[3][3];
    double rotated_matrix_1[3][3];
    double rotated_matrix_2[3][3];

    for (int i = 0; i < len; i++) {
        create_rotation_matrix(rotation_matrix_1, rolltilt[i] * 0.5, hinge_axes[i]);
        create_rotation_matrix(rotation_matrix_2, rolltilt[i] * (-0.5), hinge_axes[i]);

        matrix_multiplication(3, 3, 3, rotated_matrix_1, rotation_matrix_1, base_pair_frames[i]);
        matrix_multiplication(3, 3, 3, rotated_matrix_2, rotation_matrix_2, base_pair_frames[i + 1]);

        average_two_matrices(3, 3, middle_frames[i], rotated_matrix_1, rotated_matrix_2);
    }
}

/**
 * Compute the base pair middle frames origins
 * @param origins the array allocated for the origins
 * @param base_pair_frames_origins the array of base pair frame origins
 * @param len the length of the array of middle frames origins
 */
void get_middle_base_pair_frames_origins_array(double origins[][3],
                                               double base_pair_frames_origins[][3],
                                               int len) {
    for (int i = 0; i < len; i++) {
        average_two_vectors(3, origins[i], base_pair_frames_origins[i], base_pair_frames_origins[i + 1]);
    }
}

/**
 * Compute the shift, slide, rise coordinates
 * @param shift the array allocated for shift
 * @param slide the array allocated for slide
 * @param rise the array allocated for rise
 * @param base_pair_origins the array of base pair frames origins
 * @param base_pair_middle_frames the array of base pair middle frames
 * @param len the length of the coordinates arrays
 */
void get_shift_slide_rise(double *shift,
                          double *slide,
                          double *rise,
                          double base_pair_origins[][3],
                          double base_pair_middle_frames[][3][3],
                          int len) {
    double origins_difference[3];
    double shift_slide_rise[3];
    for (int i = 0; i < len; i++) {
        subtract_two_vectors(3, origins_difference, base_pair_origins[i + 1], base_pair_origins[i]);
        matrix_vector_multiplication(3, 3, shift_slide_rise, base_pair_middle_frames[i], origins_difference);
        shift[i] = shift_slide_rise[0];
        slide[i] = shift_slide_rise[1];
        rise[i] = shift_slide_rise[2];
    }
}

/**
 * Compute the twist angle values
 * @param twist the array allocated for twist
 * @param rolltilt the array of rolltilt angles
 * @param hinge_axes the array of hinge axes
 * @param base_pair_frames the array of base pair frames
 * @param middle_frames the array of base pair middle frames
 * @param len the length of the twist array
 */
void get_twist_array(double *twist,
                     const double rolltilt[],
                     double hinge_axes[][3],
                     double base_pair_frames[][3][3],
                     double middle_frames[][3][3],
                     int len) {

    double rotation_matrix_1[3][3];
    double rotation_matrix_2[3][3];
    double rotated_matrix_1[3][3];
    double rotated_matrix_2[3][3];
    double rotated_matrix_y_axis_1[3];
    double rotated_matrix_y_axis_2[3];

    for (int i = 0; i < len; i++) {

        create_rotation_matrix(rotation_matrix_1, rolltilt[i] * 0.5, hinge_axes[i]);
        create_rotation_matrix(rotation_matrix_2, rolltilt[i] * (-0.5), hinge_axes[i]);

        matrix_multiplication(3, 3, 3, rotated_matrix_1, rotation_matrix_1, base_pair_frames[i]);
        matrix_multiplication(3, 3, 3, rotated_matrix_2, rotation_matrix_2, base_pair_frames[i + 1]);

        select_frame_axis(rotated_matrix_y_axis_1, rotated_matrix_1, 'y');
        select_frame_axis(rotated_matrix_y_axis_2, rotated_matrix_2, 'y');

        twist[i] = dot_product(3, rotated_matrix_y_axis_1, rotated_matrix_y_axis_2) / (double)
                (vector_magnitude(3, rotated_matrix_y_axis_1) * vector_magnitude(3, rotated_matrix_y_axis_2));

        twist[i] = radians_to_degrees(acos(twist[i]));
        twist[i] *= get_angle_sign(rotated_matrix_y_axis_1, rotated_matrix_y_axis_2, middle_frames, i);
    }
}

/**
 * Function to free middle frames and origins arrays
 * @param mfi the array of base middle frames
 * @param mfio the array of origins
 */
void free_3dna_arrays(double mfi[][3][3], double mfio[][3]) {
    free(mfi);
    free(mfio);
}

/**
 * Function to free arrays allocated in the process of computing intra base pair coordinates
 * @param ha the array of hinge axes
 * @param bp the array of bucklepropeller angles
 * @param rm the array of rotation matrices
 * @param rm1 the array of rotated matrices in strand 1
 * @param rm2 the array of rotated matrices in strand 2
 * @param pa the array of phase angles
 */
void
free_intra_arrays(double ha[][3], double bp[], double rm[][3][3], double rm1[][3][3], double rm2[][3][3], double pa[]) {
    free(ha);
    free(bp);
    free(rm);
    free(rm1);
    free(rm2);
    free(pa);
}

/**
 * Function to compute intra base frame coordinates
 * @param frames_1 the array of base frames in strand 1
 * @param frames_2 the array of base frames in strand 2
 * @param origins_1 the array of origins in strand 1
 * @param origins_2 the array of origins in strand 2
 * @param middle_frames the array of base middle frames
 * @param middle_frames_origins the array of base middle frame origins
 * @param shear the array allocated for shear
 * @param stretch the array allocated for stretch
 * @param stagger the array allocated for stagger
 * @param buckle the array allocated for buckle
 * @param propeller the array allocated for propeller
 * @param opening the array allocated for opening
 * @param strand_len the number of base frames in one strand
 */
void get_intra_coordinates(double frames_1[][3][3], double frames_2[][3][3], double origins_1[][3],
                           double origins_2[][3], double middle_frames[][3][3], double middle_frames_origins[][3],
                           double shear[], double stretch[], double stagger[], double buckle[], double propeller[],
                           double opening[], int strand_len) {
    double (*hinge_axes_array_intra)[3] = malloc(strand_len * sizeof(*hinge_axes_array_intra));
    double *bucklepropeller_array = malloc(strand_len * sizeof(*bucklepropeller_array));
    double (*rotation_matrices_array)[3][3] = malloc(strand_len * sizeof(*rotation_matrices_array));
    double (*rotated_matrices_strand_I)[3][3] = malloc(strand_len * sizeof(*rotated_matrices_strand_I));
    double (*rotated_matrices_strand_II)[3][3] = malloc(strand_len * sizeof(*rotated_matrices_strand_II));
    double *phase_angle_array_base_frames = malloc(strand_len * sizeof(*phase_angle_array_base_frames));

    get_hinge_axes_intra_array(hinge_axes_array_intra, strand_len, frames_1, frames_2);
    create_bucklepropeller_array(strand_len, frames_1, frames_2, bucklepropeller_array);
    rotate_frames(rotated_matrices_strand_I, frames_1, hinge_axes_array_intra, rotation_matrices_array,
                  bucklepropeller_array, 0.5, strand_len);
    rotate_frames(rotated_matrices_strand_II, frames_2, hinge_axes_array_intra, rotation_matrices_array,
                  bucklepropeller_array, -0.5, strand_len);
    get_middle_frames_array(middle_frames, rotated_matrices_strand_I, rotated_matrices_strand_II, strand_len);
    get_middle_frames_origins_array(middle_frames_origins, origins_1, origins_2, strand_len);
    get_shear_stretch_stagger(shear, stretch, stagger, origins_1, origins_2, middle_frames, strand_len);
    get_opening_angle_array(opening, rotated_matrices_strand_I, rotated_matrices_strand_II, strand_len, middle_frames);
    get_phase_angle_array(phase_angle_array_base_frames, hinge_axes_array_intra, middle_frames, strand_len);
    get_bp_or_rt_angle_array(buckle, propeller, bucklepropeller_array, phase_angle_array_base_frames, strand_len);

    free_intra_arrays(hinge_axes_array_intra, bucklepropeller_array, rotation_matrices_array, rotated_matrices_strand_I,
                      rotated_matrices_strand_II, phase_angle_array_base_frames);
}

/**
 * Function to free arrays allocated in the process of computing the inter base pair coordinates
 * @param ha the array of hinge axes
 * @param rt the array of rolltilt angles
 * @param mf the array of base pair middle frames
 * @param mfo the array of base pair middle frame origins
 * @param pa the array of phase angles
 */
void free_inter_arrays(double ha[][3], double rt[], double mf[][3][3], double mfo[][3], double pa[]) {
    free(ha);
    free(rt);
    free(mf);
    free(mfo);
    free(pa);
}

/**
 * Fucntion to compute inter base pair coordinates
 * @param frames the array of base pair frames
 * @param frames_origins the array of base pair frame origins
 * @param shift the array for shift
 * @param slide the array for slide
 * @param rise the array for rise
 * @param roll the array for roll
 * @param tilt the array for tilt
 * @param twist the array for twist
 * @param inter_len the length of the inter base pair coordinate arrays
 */
void get_inter_coordinates(double frames[][3][3], double frames_origins[][3], double shift[], double slide[],
                           double rise[], double roll[], double tilt[], double twist[], int inter_len) {
    double (*hinge_axes_array_inter)[3] = malloc(inter_len * sizeof(*hinge_axes_array_inter));
    double *rolltilt_array = malloc(inter_len * sizeof(*rolltilt_array));
    double (*middle_base_pair_frames_array)[3][3] = malloc(inter_len * sizeof(*middle_base_pair_frames_array));
    double (*middle_base_pair_frames_origins_array)[3] = malloc(
            inter_len * sizeof(*middle_base_pair_frames_origins_array));
    double *phase_angle_array_base_pair_frames = malloc(inter_len * sizeof(*phase_angle_array_base_pair_frames));

    get_hinge_axes_inter_array(hinge_axes_array_inter, frames, inter_len);
    get_rolltilt_angle_array(rolltilt_array, frames, inter_len);
    get_middle_base_pair_frames_array(middle_base_pair_frames_array, frames, hinge_axes_array_inter, rolltilt_array,
                                      inter_len);
    get_middle_base_pair_frames_origins_array(middle_base_pair_frames_origins_array, frames_origins, inter_len);
    get_shift_slide_rise(shift, slide, rise, frames_origins, middle_base_pair_frames_array, inter_len);
    get_twist_array(twist, rolltilt_array, hinge_axes_array_inter, frames, middle_base_pair_frames_array, inter_len);
    get_phase_angle_array(phase_angle_array_base_pair_frames, hinge_axes_array_inter, middle_base_pair_frames_array,
                          inter_len);
    get_bp_or_rt_angle_array(roll, tilt, rolltilt_array, phase_angle_array_base_pair_frames, inter_len);

    free_inter_arrays(hinge_axes_array_inter, rolltilt_array, middle_base_pair_frames_array,
                      middle_base_pair_frames_origins_array, phase_angle_array_base_pair_frames);
}

/**
 * Function to compute all the internal coordinates according to the 3DNA definition
 * @param frames_1 the array of base frames in strand 1
 * @param frames_2 the array of base frames in strand 2
 * @param origins_1 the array of origins in strand 1
 * @param origins_2 the array of origins in strand 2
 * @param shear the array for shear
 * @param stretch the array for stretch
 * @param stagger the array for stagger
 * @param buckle the array for buckle
 * @param propeller the array for propeller
 * @param opening the array for opening
 * @param shift the array for shift
 * @param slide the array for slide
 * @param rise the array for rise
 * @param roll the array for roll
 * @param tilt the array for tilt
 * @param twist the array for twist
 */
void get_3dna_coordinates(double frames_1[][3][3], double frames_2[][3][3], double origins_1[][3],
                          double origins_2[][3], double shear[], double stretch[], double stagger[], double buckle[],
                          double propeller[], double opening[], double shift[], double slide[], double rise[],
                          double roll[], double tilt[], double twist[]) {

    //separate snapshot base frames into strands
    int strand_len = snapshot_len / 2;
    int inter_len = strand_len - 1;

    //base frame coordinates
    double (*middle_base_frames_array)[3][3] = malloc(strand_len * sizeof(*middle_base_frames_array));
    double (*middle_base_frames_origins_array)[3] = malloc(strand_len * sizeof(*middle_base_frames_origins_array));

    get_intra_coordinates(frames_1, frames_2, origins_1, origins_2, middle_base_frames_array,
                          middle_base_frames_origins_array, shear, stagger, stretch, buckle, propeller, opening,
                          strand_len);

    //base pair frames coordinates
    double (*base_pair_frames)[3][3] = middle_base_frames_array;
    double (*base_pair_frames_origins)[3] = middle_base_frames_origins_array;
    get_inter_coordinates(base_pair_frames, base_pair_frames_origins, shift, slide, rise, roll, tilt, twist, inter_len);

    free_3dna_arrays(middle_base_frames_array, middle_base_frames_origins_array);
}
