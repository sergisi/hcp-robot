//
// Created by sergi on 28/05/22.
//

#ifndef LIFEGAME_MPI_UPDATE_FILE_H
#define LIFEGAME_MPI_UPDATE_FILE_H

#include "read_life.h"
#include "constants.h"
#include "utilities.h"
#include "debug_tools.h"
#include <mpi.h>
#include <omp.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int
mpi_life_update(int it, int it_max, int m, int8_t *grid, int iproc, int size_to_work, int iproc_prev, int iproc_next,
                MPI_Request *requests);


void compute_life_or_dead(int m, int n, const int8_t *grid, int i, int j, int8_t *s);

int
mpi_life_update(int it, int it_max, int m, int8_t *grid, int iproc, int size_to_work, int iproc_prev, int iproc_next,
                MPI_Request *requests) {
    double start, communicating = 0, computing = 0;
    int n = size_to_work;
    int i;
    int j;
    guardia_iproc(iproc, "Initializing s");
    int8_t *s;
    s = malloc((m + 2) * (n + 2) * sizeof(int8_t));
    guardia_iproc(iproc, "Setting to 0");

    for (it = 1; it <= it_max; it++) {
        start = MPI_Wtime();
        MPI_Irecv(grid, m + 2, MPI_INT8_T, iproc_prev, tag_prev_grid, MPI_COMM_WORLD, &requests[0]);
        MPI_Isend(grid + m + 2, m + 2, MPI_INT8_T, iproc_prev, tag_next_grid, MPI_COMM_WORLD, &requests[1]);
        MPI_Isend(grid + size_to_work * (m + 2), m + 2, MPI_INT8_T, iproc_next, tag_prev_grid, MPI_COMM_WORLD,
                  &requests[2]);
        MPI_Irecv(grid + (size_to_work + 1) * (m + 2), m + 2, MPI_INT8_T, iproc_next, tag_next_grid, MPI_COMM_WORLD,
                  &requests[3]);
        MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
        communicating += MPI_Wtime() - start;
        start = MPI_Wtime();
        debug_grid(m, grid, size_to_work, iproc, it);

#pragma omp parallel for schedule(guided) private(i, j) num_threads(cores)
        for (j = 1; j <= n; j++) { // the long one is N
            for (i = 1; i <= m; i++) {
                compute_life_or_dead(m, n, grid, i, j, s);
            }
        }
        guardia_iproc(iproc, "computed s");
/*
  Any dead cell with 3 live neighbors becomes alive.
  Any living cell with less than 2 or more than 3 neighbors dies.
*/

#pragma omp parallel for schedule(guided) private(i, j) num_threads(cores)
        for (j = 1; j <= n; j++) {
            for (i = 1; i <= m; i++) {
                if (grid[i + j * (m + 2)] == 0) {
                    if (s[i - 1 + (j - 1) * m] == 3) {
                        grid[i + j * (m + 2)] = 1;
                    }
                } else if (grid[i + j * (m + 2)] == 1 && (s[i - 1 + (j - 1) * m] < 2 || 3 < s[i - 1 + (j - 1) * m])) {
                    grid[i + j * (m + 2)] = 0;
                }
            }
        }
        guardia_iproc(iproc, "computed grid");
        debug_update(it, iproc);
        computing += MPI_Wtime() - start;
    }
    printf("[%d]: communicating time %f computing time %f\n", iproc, communicating, computing);
    debug_end(iproc);
    free(s);
    return it;
}

void compute_life_or_dead(int m, int n, const int8_t *grid, int i, int j, int8_t *s) {
    int i_prev, i_next, j_prev, j_next;
    i_prev = (1 < i) ? i - 1 : m;
    i_next = (i < m) ? i + 1 : 1;
    j_prev = j - 1;
    j_next = j + 1;
    // adds as ints, as processors is 64 bits
    s[i - 1 + (j - 1) * m] =
            grid[i_prev + (j_prev) * (m + 2)] + grid[i_prev + j * (m + 2)] + grid[i_prev + (j_next) * (m + 2)] +
            grid[i + (j_prev) * (m + 2)] + grid[i + (j_next) * (m + 2)] + grid[i_next + (j_prev) * (m + 2)] +
            grid[i_next + j * (m + 2)] + grid[i_next + (j_next) * (m + 2)];
}

#endif //LIFEGAME_MPI_UPDATE_FILE_H
