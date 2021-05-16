#include <stdio.h>
#include "frame_fitting.h"

void allocate_frames_and_origins(double (**frames_1)[3][3], double (**frames_2)[3][3], double (**origins_1)[3],
                                 double (**origins_2)[3]){
    *frames_1 = malloc(SNAPSHOT_ARRAY_LENGTH * sizeof(**frames_1));
    *frames_2 = malloc(SNAPSHOT_ARRAY_LENGTH * sizeof(**frames_2));
    *origins_1 = malloc(SNAPSHOT_ARRAY_LENGTH * sizeof(**origins_1));
    *origins_2 = malloc(SNAPSHOT_ARRAY_LENGTH * sizeof(**origins_2));
}

void allocate_coordinates_arrays(double **shear, double **stretch, double **stagger, double **buckle,
                                 double **propeller, double **opening, double **shift, double **slide,
                                 double **rise, double **roll, double **tilt, double **twist){
    int strand_len = (int)(snapshot_len / 2);
    int inter_len = strand_len - 1;

    *shear = (double*)malloc(strand_len * sizeof(**shear));
    *stretch = (double*)malloc(strand_len * sizeof(**stretch));
    *stagger = (double*)malloc(strand_len * sizeof(**stagger));
    *buckle = (double*)malloc(strand_len * sizeof(**buckle));
    *propeller = (double*)malloc(strand_len * sizeof(**propeller));
    *opening = (double*)malloc(strand_len * sizeof(**opening));

    *shift = (double*)malloc(inter_len * sizeof(**shift));
    *slide = (double*)malloc(inter_len * sizeof(**slide));
    *rise = (double*)malloc(inter_len * sizeof(**rise));
    *roll = (double*)malloc(inter_len * sizeof(**roll));
    *tilt = (double*)malloc(inter_len * sizeof(**tilt));
    *twist = (double*)malloc(inter_len * sizeof(**twist));
}

int main() {
    //todo switch frames
    //todo read arguments from command line
    //todo n pdb files
    FILE *fp = fopen("../teplota.300.pdb.1", "r"); // ../ for running in ide, without for running in terminal
    if(fp == NULL){
        printf("Error: Input file not found!");
        return 1;
    }

    double (*frames_strand_1)[3][3], (*frames_strand_2)[3][3], (*origins_strand_1)[3], (*origins_strand_2)[3];
    allocate_frames_and_origins(&frames_strand_1, &frames_strand_2, &origins_strand_1, &origins_strand_2);
    fit_frames(fp, frames_strand_1, frames_strand_2, origins_strand_1, origins_strand_2);

    double *shear, *stretch, *stagger, *buckle, *propeller, *opening, *shift, *slide, *rise, *roll, *tilt, *twist;
    allocate_coordinates_arrays(&shear, &stretch, &stagger, &buckle, &propeller, &opening, &shift, &slide, &rise,
                                &roll, &tilt, &twist);
    //todo run functions for 3dna, curves, cgdna
    return 0;
}
