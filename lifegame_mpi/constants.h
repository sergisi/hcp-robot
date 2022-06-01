//
// Created by sergi on 28/05/22.
//

#ifndef LIFEGAME_MPI_CONSTANTS_H
#define LIFEGAME_MPI_CONSTANTS_H

const int debug = 0 != 0;
const int DNumIterForPartialResults = 25;
const int cores = 1;
const char *DDefaultOutputFilename = "./Life_%04d.txt";

enum {
    tag_prev_grid = 0,
    tag_next_grid = 1,
    tag_end_grid = 2,
    tag_n = 3,
    tag_m = 4,
    tag_it_max = 5,
    tag_initial_grid = 6
};

#include <stdio.h>

#endif //LIFEGAME_MPI_CONSTANTS_H
