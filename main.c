#include <stdio.h>
#include "frame_fitting.h"
#include "curves_coordinates.h"
#include "cgdna_coordinates.h"
#include "3dna_coordinates.h"

#define BUFFER_SIZE 100
#define FILE_PATH_SIZE 100

/**
 * Allocate memory space for frames and origins arrays
 * @param frames_1 pointer to an array of frames in strand I
 * @param frames_2 pointer to an array of frames in strand II
 * @param origins_1 pointer pointer to an array of origins of frames in strand I
 * @param origins_2 pointer pointer to an array of origins of frames in strand II
 */
void allocate_frames_and_origins(double (**frames_1)[3][3], double (**frames_2)[3][3], double (**origins_1)[3],
                                 double (**origins_2)[3]) {
    *frames_1 = malloc(SNAPSHOT_ARRAY_LENGTH * sizeof(**frames_1));
    *frames_2 = malloc(SNAPSHOT_ARRAY_LENGTH * sizeof(**frames_2));
    *origins_1 = malloc(SNAPSHOT_ARRAY_LENGTH * sizeof(**origins_1));
    *origins_2 = malloc(SNAPSHOT_ARRAY_LENGTH * sizeof(**origins_2));
}

/**
 * Allocate memory space for arrays to be filled with relative coordinates values
 * @param shear array pointer
 * @param stretch array pointer
 * @param stagger array pointer
 * @param buckle array pointer
 * @param propeller array pointer
 * @param opening array pointer
 * @param shift array pointer
 * @param slide array pointer
 * @param rise array pointer
 * @param roll array pointer
 * @param tilt array pointer
 * @param twist array pointer
 */
void allocate_coordinates_arrays(double **shear, double **stretch, double **stagger, double **buckle,
                                 double **propeller, double **opening, double **shift, double **slide,
                                 double **rise, double **roll, double **tilt, double **twist) {
    int strand_len = (int) (snapshot_len / 2);
    int inter_len = strand_len - 1;

    *shear = (double *) malloc(strand_len * sizeof(**shear));
    *stretch = (double *) malloc(strand_len * sizeof(**stretch));
    *stagger = (double *) malloc(strand_len * sizeof(**stagger));
    *buckle = (double *) malloc(strand_len * sizeof(**buckle));
    *propeller = (double *) malloc(strand_len * sizeof(**propeller));
    *opening = (double *) malloc(strand_len * sizeof(**opening));

    *shift = (double *) malloc(inter_len * sizeof(**shift));
    *slide = (double *) malloc(inter_len * sizeof(**slide));
    *rise = (double *) malloc(inter_len * sizeof(**rise));
    *roll = (double *) malloc(inter_len * sizeof(**roll));
    *tilt = (double *) malloc(inter_len * sizeof(**tilt));
    *twist = (double *) malloc(inter_len * sizeof(**twist));
}

/**
 * Append the base pair coordinates to the output file
 * @param pdb_number the number of pdb file
 * @param output_path path to the output file
 * @param shear the shear array
 * @param stretch the stretch array
 * @param stagger the stagger array
 * @param buckle the buckle array
 * @param propeller the propeller array
 * @param opening the opening array
 */
void write_to_bp_output_file(int pdb_number, char output_path[], double shear[], double stretch[], double stagger[],
                             double buckle[], double propeller[], double opening[]){
    int intra_len = (int)(snapshot_len / 2);
    for(int i = 1; i <= intra_len; i++){
        char filename_output[FILE_PATH_SIZE] = "\0";
        snprintf(filename_output, FILE_PATH_SIZE, "%s_bp_prm_%d.out", output_path, i);
        char bp_prm_to_append[BUFFER_SIZE] = "\0";
        snprintf(bp_prm_to_append, BUFFER_SIZE, "%d\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n", pdb_number, shear[i], stretch[i], stagger[i],
                 buckle[i], propeller[i], opening[i]);

        FILE *fp = fopen(filename_output, "a");
        if (fp == NULL){
            printf("ERROR: Specified output directory does not exist! %s\n", output_path);
            return;
        }
        fputs(bp_prm_to_append, fp);
        fclose(fp);
    }
}

/**
 * Append the base pair step coordinates to the output file
 * @param pdb_number the number of pdb file
 * @param output_path path to the output file
 * @param shift the shift array
 * @param slide the slide array
 * @param rise the rise array
 * @param roll the roll array
 * @param tilt the tilt array
 * @param twist the twist array
 */
void write_to_bp_step_output_file(int pdb_number, char output_path[], double shift[], double slide[], double rise[], double roll[],
                             double tilt[], double twist[]){
    int inter_len = (int)(snapshot_len / 2) - 1;
    for(int i = 1; i <= inter_len; i++){
        char filename_output[FILE_PATH_SIZE] = "\0";
        snprintf(filename_output, FILE_PATH_SIZE, "%s_bp_step_%d.out", output_path, i);
        char bp_prm_to_append[BUFFER_SIZE] = "\0";
        snprintf(bp_prm_to_append, BUFFER_SIZE, "%d\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n", pdb_number, shift[i], slide[i],
                 rise[i], roll[i], tilt[i], twist[i]);

        FILE *fp = fopen(filename_output, "a");
        if (fp == NULL){
            printf("ERROR: Specified output directory does not exist! %s\n", output_path);
            return;
        }
        fputs(bp_prm_to_append, fp);
        fclose(fp);
    }
}

/**
 * Run the computation of 3DNA definition of internal coordinates
 * @param filename input pdb file with snapshot
 * @return value computation successful (0) or file not found (1)
 */
int run_3dna(char filename_number[], char output_path[], int pdb_number){
    FILE *fp = fopen(filename_number, "r");
    if (fp == NULL) {
        printf("ERROR: Input file %s not found!\n", filename_number);
        return 1;
    }

    double (*frames_strand_1)[3][3], (*frames_strand_2)[3][3], (*origins_strand_1)[3], (*origins_strand_2)[3];
    allocate_frames_and_origins(&frames_strand_1, &frames_strand_2, &origins_strand_1, &origins_strand_2);
    fit_frames(fp, frames_strand_2, frames_strand_1, origins_strand_2, origins_strand_1);


    double *shear, *stretch, *stagger, *buckle, *propeller, *opening, *shift, *slide, *rise, *roll, *tilt, *twist;
    allocate_coordinates_arrays(&shear, &stretch, &stagger, &buckle, &propeller, &opening, &shift, &slide, &rise,
                                &roll, &tilt, &twist);

    get_3dna_coordinates(frames_strand_1, frames_strand_2, origins_strand_1, origins_strand_2, shear, stretch,
                         stagger, buckle, propeller, opening, shift, slide, rise, roll, tilt, twist);

    write_to_bp_output_file(pdb_number, output_path, shear, stretch, stagger, buckle, propeller, opening);
    write_to_bp_step_output_file(pdb_number, output_path, shift, slide, rise, roll, tilt, twist);
    return 0;
}

/**
 * Run the computation of Curves+ definition of internal coordinates
 * @param filename input pdb file with snapshot
 * @return value computation successful (0) or file not found (1)
 */
int run_curves(char filename[], char output_path[], int pdb_number){
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("ERROR: Input file %s not found!\n", filename);
        return 1;
    }

    double (*frames_strand_1)[3][3], (*frames_strand_2)[3][3], (*origins_strand_1)[3], (*origins_strand_2)[3];
    allocate_frames_and_origins(&frames_strand_1, &frames_strand_2, &origins_strand_1, &origins_strand_2);
    fit_frames(fp, frames_strand_2, frames_strand_1, origins_strand_2, origins_strand_1);

    double *shear, *stretch, *stagger, *buckle, *propeller, *opening, *shift, *slide, *rise, *roll, *tilt, *twist;
    allocate_coordinates_arrays(&shear, &stretch, &stagger, &buckle, &propeller, &opening, &shift, &slide, &rise,
                                &roll, &tilt, &twist);

    get_curves_coordinates(frames_strand_1, frames_strand_2, origins_strand_1, origins_strand_2, shear, stretch,
                           stagger, buckle, propeller, opening, shift, slide, rise, roll, tilt, twist);

    write_to_bp_output_file(pdb_number, output_path, shear, stretch, stagger, buckle, propeller, opening);
    write_to_bp_step_output_file(pdb_number, output_path, shift, slide, rise, roll, tilt, twist);
    return 0;
}

/**
 * Run the computation of cgDNA definition of internal coordinates
 * @param filename input pdb file with snapshot
 * @return value computation successful (0) or file not found (1)
 */
int run_cgdna(char filename[], char output_path[], int pdb_number){
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("ERROR: Input file %s not found!\n", filename);
        return 1;
    }

    double (*frames_strand_1)[3][3], (*frames_strand_2)[3][3], (*origins_strand_1)[3], (*origins_strand_2)[3];
    allocate_frames_and_origins(&frames_strand_1, &frames_strand_2, &origins_strand_1, &origins_strand_2);
    fit_frames(fp, frames_strand_2, frames_strand_1, origins_strand_2, origins_strand_1);

    double *shear, *stretch, *stagger, *buckle, *propeller, *opening, *shift, *slide, *rise, *roll, *tilt, *twist;
    allocate_coordinates_arrays(&shear, &stretch, &stagger, &buckle, &propeller, &opening, &shift, &slide, &rise,
                                &roll, &tilt, &twist);

    get_cgdna_coordinates(frames_strand_1, frames_strand_2, origins_strand_1, origins_strand_2, shear, stretch,
                          stagger, buckle, propeller, opening, shift, slide, rise, roll, tilt, twist);

    write_to_bp_output_file(pdb_number, output_path, shear, stretch, stagger, buckle, propeller, opening);
    write_to_bp_step_output_file(pdb_number, output_path, shift, slide, rise, roll, tilt, twist);
    return 0;
}

int main(int argc, char *argv[]) {

    if (argc != 9 || (strcmp(argv[1], "-t") != 0) || (strcmp(argv[3], "-n") != 0) || (strcmp(argv[5], "-i") != 0) ||
            (strcmp(argv[7], "-o") != 0)) {
        printf("USAGE: ./coordinates -t [3dna|curves|cgdna] -n number_of_snapshots -i pdb_files_path -o output_files_path\n");
        return 1;
    }

    char coords_type[10] = "3dna";
    strcpy(coords_type, argv[2]);
    char *ptr;
    long number_snapshots = strtol(argv[4], &ptr, 10);
    char input_path[FILE_PATH_SIZE] = "./";
    strcpy(input_path, argv[6]);
    char output_path[FILE_PATH_SIZE] = "./";
    strcpy(output_path, argv[8]);
    char *filename = strrchr( input_path, '/');
    filename++;

    if (strcmp(coords_type, "3dna") == 0) {
        for(int i = 1; i <= number_snapshots; i++){
            char filename_number[FILE_PATH_SIZE] = "\0";
            snprintf(filename_number, FILE_PATH_SIZE, "%s%d", input_path, i);

            if(run_3dna(filename_number, output_path, i) == 1){
                return 1;
            }
        }
    } else if (strcmp(coords_type, "curves") == 0) {
        for(int i = 1; i <= number_snapshots; i++){
            char filename_number[FILE_PATH_SIZE] = "\0";
            snprintf(filename_number, FILE_PATH_SIZE, "%s%d", input_path, i);

            if(run_curves(filename_number, output_path, i) == 1){
                return 1;
            }
        }
    } else if (strcmp(coords_type, "cgdna") == 0) {
        for(int i = 1; i <= number_snapshots; i++){
            char filename_number[FILE_PATH_SIZE] = "\0";
            snprintf(filename_number, FILE_PATH_SIZE, "%s%d", input_path, i);

            if(run_cgdna(filename_number, output_path, i) == 1){
                return 1;
            }
        }
    }
    return 0;
}
