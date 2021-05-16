#include <stdio.h>
#include "frame_fitting.h"

int main() {
    //todo strand_len, snapshot_len global variable
    FILE *fp = fopen("../teplota.300.pdb.1", "r"); // ../ for running in ide, without for running in terminal
    if(fp == NULL){
        printf("Error: Input file not found!");
        return 1;
    }

    double (*frames_strand_1)[3][3] = malloc(SNAPSHOT_ARRAY_LENGTH * sizeof(*frames_strand_1));
    double (*frames_strand_2)[3][3] = malloc(SNAPSHOT_ARRAY_LENGTH * sizeof(*frames_strand_2));
    double (*origins_strand_1)[3] = malloc(SNAPSHOT_ARRAY_LENGTH * sizeof(*origins_strand_1));
    double (*origins_strand_2)[3] = malloc(SNAPSHOT_ARRAY_LENGTH * sizeof(*origins_strand_2));
    fit_frames(fp, frames_strand_1, frames_strand_2, origins_strand_1, origins_strand_2);
    //todo switch frames

    return 0;
}
