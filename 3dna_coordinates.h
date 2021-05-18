#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define SNAPSHOT_ARRAY_LENGTH 100

void select_frames_axis(double *axis_array,
                        double frame[][3][3],
                        char axis,
                        int index) {
    int a = axis - 120;

    for (int i = 0; i < 3; i++) {
        axis_array[i] = frame[index][i][0 + a];
    }
}

void select_frame_axis(double *axis_array,
                       double frame[3][3],
                       char axis) {
    int a = axis - 120;

    for (int i = 0; i < 3; i++) {
        axis_array[i] = frame[i][0 + a];
    }
}

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

void create_bucklepropeller_array(int len,
                                  double strand_1[][3][3],
                                  double strand_2[][3][3],
                                  double *bp_array) {
    double z1[3], z2[3];

    for (int i = 0; i < len; i++) {
        select_frames_axis(z1, strand_1, 'z', i);
        select_frames_axis(z2, strand_2, 'z', i);

        double angle = dot_product(3, z1, z2) / (vector_magnitude(3, z1) * vector_magnitude(3, z2));
        angle = radians_to_degrees(acos(angle));
        bp_array[i] = angle;
    }
}

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

void get_middle_frames_array(double result[][3][3],
                             double frames_I[][3][3],
                             double frames_II[][3][3],
                             int strand_len) {
    for (int i = 0; i < strand_len; i++) {
        average_two_matrices(3, 3, result[i], frames_I[i], frames_II[i]);
    }
}

void get_middle_frames_origins_array(double result[][3],
                                     double origins_I[][3],
                                     double origins_II[][3],
                                     int strand_len) {
    for (int i = 0; i < strand_len; i++) {
        average_two_vectors(3, result[i], origins_I[i], origins_II[i]);
    }
}

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
        //todo check 1- 2 or 2 - 1
        subtract_two_vectors(3, difference_origins, origins_II[i], origins_I[i]);
        matrix_vector_multiplication(3, 3, shear_stretch_stagger, middle_frames[i], difference_origins);
        shear[i] = shear_stretch_stagger[0];
        stretch[i] = shear_stretch_stagger[1];
        stagger[i] = shear_stretch_stagger[2];
    }
}

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
                           (vector_magnitude(3, axis_rotated_1) * vector_magnitude(3, axis_rotated_2));
        opening_array[i] = radians_to_degrees(acos(opening_array[i]));
        opening_array[i] = opening_array[i] * get_angle_sign(axis_rotated_1, axis_rotated_2, middle_frames_array, i);
    }
}

void get_phase_angle_array(double *phase,
                           double hinge_axes[][3],
                           double middle_frames[][3][3],
                           int strand_len) {
    double middle_frame_y[3];
    for (int i = 0; i < strand_len; i++) {
        select_frames_axis(middle_frame_y, middle_frames, 'y', i);
        phase[i] = dot_product(3, hinge_axes[i], middle_frame_y) /
                   (vector_magnitude(3, hinge_axes[i]) * vector_magnitude(3, middle_frame_y));
        phase[i] = radians_to_degrees(acos(phase[i]));
        phase[i] = phase[i] * get_angle_sign(hinge_axes[i], middle_frame_y, middle_frames, i);
    }
}

void get_bp_or_rt_angle_array(double *buckle,
                              double *propeller,
                              const double *bp_array,
                              double *phase_angle,
                              int size) {
    for (int i = 0; i < size; i++) {
        buckle[i] = bp_array[i] * sin(degrees_to_radians(phase_angle[i]));
        propeller[i] = bp_array[i] * cos(degrees_to_radians(phase_angle[i]));
    }
}

void get_hinge_axes_inter_array(double hinge_axes[][3],
                                double base_pair_frames[][3][3],
                                int size) {
    double base_pair_z_1[3];
    double base_pair_z_2[3];
    for (int i = 0; i < size; i++) {
        select_frames_axis(base_pair_z_1, base_pair_frames, 'z', i);
        select_frames_axis(base_pair_z_2, base_pair_frames, 'z', i + 1);

        cross_product(hinge_axes[i], base_pair_z_1, base_pair_z_2);
        normalize_vector(3, hinge_axes[i]);
    }
}

void get_rolltilt_angle_array(double *rolltilt,
                              double base_pair_frames[][3][3],
                              int size) {
    double base_pair_z_1[3];
    double base_pair_z_2[3];
    for (int i = 0; i < size; i++) {
        select_frames_axis(base_pair_z_1, base_pair_frames, 'z', i);
        select_frames_axis(base_pair_z_2, base_pair_frames, 'z', i + 1);

        rolltilt[i] = dot_product(3, base_pair_z_1, base_pair_z_2) /
                      (vector_magnitude(3, base_pair_z_1) * vector_magnitude(3, base_pair_z_2));
        rolltilt[i] = radians_to_degrees(acos(rolltilt[i]));
    }
}

void get_middle_base_pair_frames_array(double middle_frames[][3][3],
                                       double base_pair_frames[][3][3],
                                       double hinge_axes[][3],
                                       const double rolltilt[],
                                       int size) {
    double rotation_matrix_1[3][3];
    double rotation_matrix_2[3][3];
    double rotated_matrix_1[3][3];
    double rotated_matrix_2[3][3];

    for (int i = 0; i < size; i++) {
        create_rotation_matrix(rotation_matrix_1, rolltilt[i] * 0.5, hinge_axes[i]);
        create_rotation_matrix(rotation_matrix_2, rolltilt[i] * (-0.5), hinge_axes[i]);

        matrix_multiplication(3, 3, 3, rotated_matrix_1, rotation_matrix_1, base_pair_frames[i]);
        matrix_multiplication(3, 3, 3, rotated_matrix_2, rotation_matrix_2, base_pair_frames[i + 1]);

        average_two_matrices(3, 3, middle_frames[i], rotated_matrix_1, rotated_matrix_2);
    }
}

void get_middle_base_pair_frames_origins_array(double origins[][3],
                                               double base_pair_frames_origins[][3],
                                               int size) {
    for (int i = 0; i < size; i++) {
        average_two_vectors(3, origins[i], base_pair_frames_origins[i], base_pair_frames_origins[i + 1]);
    }
}

void get_shift_slide_rise(double *shift,
                          double *slide,
                          double *rise,
                          double base_pair_origins[][3],
                          double base_pair_middle_frames[][3][3],
                          int size) {
    double origins_difference[3];
    double shift_slide_rise[3];
    for (int i = 0; i < size; i++) {
        subtract_two_vectors(3, origins_difference, base_pair_origins[i + 1], base_pair_origins[i]);
        matrix_vector_multiplication(3, 3, shift_slide_rise, base_pair_middle_frames[i], origins_difference);
        shift[i] = shift_slide_rise[0];
        slide[i] = shift_slide_rise[1];
        rise[i] = shift_slide_rise[2];
    }
}

void get_twist_array(double *twist,
                     const double rolltilt[],
                     double hinge_axes[][3],
                     double base_pair_frames[][3][3],
                     double middle_frames[][3][3],
                     int size) {

    double rotation_matrix_1[3][3];
    double rotation_matrix_2[3][3];
    double rotated_matrix_1[3][3];
    double rotated_matrix_2[3][3];
    double rotated_matrix_y_axis_1[3];
    double rotated_matrix_y_axis_2[3];

    for (int i = 0; i < size; i++) {

        create_rotation_matrix(rotation_matrix_1, rolltilt[i] * 0.5, hinge_axes[i]);
        create_rotation_matrix(rotation_matrix_2, rolltilt[i] * (-0.5), hinge_axes[i]);

        matrix_multiplication(3, 3, 3, rotated_matrix_1, rotation_matrix_1, base_pair_frames[i]);
        matrix_multiplication(3, 3, 3, rotated_matrix_2, rotation_matrix_2, base_pair_frames[i + 1]);

        select_frame_axis(rotated_matrix_y_axis_1, rotated_matrix_1, 'y');
        select_frame_axis(rotated_matrix_y_axis_2, rotated_matrix_2, 'y');

        twist[i] = dot_product(3, rotated_matrix_y_axis_1, rotated_matrix_y_axis_2)
                   / (vector_magnitude(3, rotated_matrix_y_axis_1) * vector_magnitude(3, rotated_matrix_y_axis_2));

        twist[i] = radians_to_degrees(acos(twist[i]));
        twist[i] *= get_angle_sign(rotated_matrix_y_axis_1, rotated_matrix_y_axis_2, middle_frames, i);
    }
}

void free_3dna_arrays(double hai[][3], double bp[], double rm[][3][3], double rm1[][3][3], double rm2[][3][3],
                      double mfi[][3][3], double mfio[][3], double pai[],
                      double hae[][3], double rt[], double mbp[][3][3], double mbpo[][3], double pae[]){
    free(hai); free(bp); free(rm); free(rm1); free(rm2); free(mfi); free(mfio); free(pai);
    free(hae); free(rt); free(mbp); free(mbpo); free(pae);
}

//todo separate into functions
void get_3dna_coordinates(double frames_1[][3][3], double frames_2[][3][3], double origins_1[][3],
                          double origins_2[][3], double shear[], double stretch[], double stagger[], double buckle[],
                          double propeller[], double opening[], double shift[], double slide[], double rise[],
                          double roll[], double tilt[], double twist[]) {

    //separate snapshot base frames into strands
    int strand_len = snapshot_len / 2;
    int inter_len = strand_len - 1;

    //hinge axis
    double (*hinge_axes_array_intra)[3];
    hinge_axes_array_intra = malloc(strand_len * sizeof(*hinge_axes_array_intra));
    get_hinge_axes_intra_array(hinge_axes_array_intra, strand_len, frames_1, frames_2);

    //bucklepropeller_array angle array
    double *bucklepropeller_array = malloc(strand_len * sizeof(*bucklepropeller_array));
    create_bucklepropeller_array(strand_len, frames_1, frames_2, bucklepropeller_array);

    //rotated matrix array
    double (*rotation_matrices_array)[3][3] = malloc(strand_len * sizeof(*rotation_matrices_array));
    double (*rotated_matrices_strand_I)[3][3] = malloc(strand_len * sizeof(*rotated_matrices_strand_I));
    double (*rotated_matrices_strand_II)[3][3] = malloc(strand_len * sizeof(*rotated_matrices_strand_II));

    rotate_frames(rotated_matrices_strand_I, frames_1, hinge_axes_array_intra, rotation_matrices_array,
                  bucklepropeller_array, 0.5, strand_len);
    rotate_frames(rotated_matrices_strand_II, frames_2, hinge_axes_array_intra, rotation_matrices_array,
                  bucklepropeller_array, -0.5, strand_len);

    //middle frames array
    double (*middle_base_frames_array)[3][3] = malloc(strand_len * sizeof(*middle_base_frames_array));
    get_middle_frames_array(middle_base_frames_array, rotated_matrices_strand_I, rotated_matrices_strand_II,
                            strand_len);
    double (*middle_base_frames_origins_array)[3] = malloc(strand_len * sizeof(*middle_base_frames_origins_array));
    get_middle_frames_origins_array(middle_base_frames_origins_array, origins_1,
                                    origins_2, strand_len);

    //translational coordinates shear, stretch, stagger
    get_shear_stretch_stagger(shear, stretch, stagger, origins_1, origins_2,
                              middle_base_frames_array, strand_len);

    //opening angle array
    get_opening_angle_array(opening, rotated_matrices_strand_I, rotated_matrices_strand_II, strand_len,
                            middle_base_frames_array);

    //phase angle array
    double *phase_angle_array_base_frames = malloc(strand_len * sizeof(*phase_angle_array_base_frames));
    get_phase_angle_array(phase_angle_array_base_frames, hinge_axes_array_intra, middle_base_frames_array,
                          strand_len);

    get_bp_or_rt_angle_array(buckle, propeller, bucklepropeller_array, phase_angle_array_base_frames, strand_len);

    //base pair frames coordinates
    double (*base_pair_frames)[3][3] = middle_base_frames_array;
    double (*base_pair_frames_origins)[3] = middle_base_frames_origins_array;

    //hinge axes array
    double (*hinge_axes_array_inter)[3] = malloc(inter_len * sizeof(*hinge_axes_array_inter));
    get_hinge_axes_inter_array(hinge_axes_array_inter, base_pair_frames, inter_len);

    //rolltilt angle array
    double *rolltilt_array = malloc(inter_len * sizeof(*rolltilt_array));
    get_rolltilt_angle_array(rolltilt_array, base_pair_frames, inter_len);

    //middle base-pair frames array
    double (*middle_base_pair_frames_array)[3][3] = malloc(strand_len * sizeof(*middle_base_pair_frames_array));
    get_middle_base_pair_frames_array(middle_base_pair_frames_array, base_pair_frames, hinge_axes_array_inter,
                                      rolltilt_array, inter_len);

    //middle base-pair frames origins array
    double (*middle_base_pair_frames_origins_array)[3] = malloc(
            strand_len * sizeof(*middle_base_pair_frames_origins_array));
    get_middle_base_pair_frames_origins_array(middle_base_pair_frames_origins_array, base_pair_frames_origins,
                                              inter_len);

    //translational base-pair frames coordinates
    get_shift_slide_rise(shift, slide, rise, base_pair_frames_origins, middle_base_pair_frames_array, inter_len);

    //rotational base-pair frames coordinates
    //twist angle array
    get_twist_array(twist, rolltilt_array, hinge_axes_array_inter, base_pair_frames, middle_base_pair_frames_array,
                    inter_len);

    //phase angle array
    double *phase_angle_array_base_pair_frames = malloc(inter_len * sizeof(*phase_angle_array_base_pair_frames));

    get_phase_angle_array(phase_angle_array_base_pair_frames, hinge_axes_array_inter, middle_base_pair_frames_array,
                          inter_len);

    //roll tilt angle array
    get_bp_or_rt_angle_array(tilt, roll, rolltilt_array, phase_angle_array_base_pair_frames, inter_len);
    free_3dna_arrays(hinge_axes_array_intra, bucklepropeller_array, rotation_matrices_array, rotated_matrices_strand_I,
                     rotated_matrices_strand_II, middle_base_frames_array, middle_base_frames_origins_array,
                     phase_angle_array_base_frames, hinge_axes_array_inter,
                     rolltilt_array, middle_base_pair_frames_array, middle_base_pair_frames_origins_array,
                     phase_angle_array_base_pair_frames);
}

#ifndef CMAKE_FINAL_THESIS_C_3DNA_COORDINATES_H
#define CMAKE_FINAL_THESIS_C_3DNA_COORDINATES_H

#endif //CMAKE_FINAL_THESIS_C_3DNA_COORDINATES_H
