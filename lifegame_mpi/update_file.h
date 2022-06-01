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

int mpi_life_update(int it, int it_max, int m, int *grid, int iproc, int size_to_work, int iproc_prev, int iproc_next,
                    MPI_Request *requests);

void life_update(int m, int n, int grid[], int iproc);

void compute_life_or_dead(int m, int n, const int *grid, int i, int j, int *s);

int mpi_life_update(int it, int it_max, int m, int *grid, int iproc, int size_to_work, int iproc_prev, int iproc_next,
                    MPI_Request *requests) {
    for (it = 1; it <= it_max; it++) {
        MPI_Irecv(grid, m + 2,
                  MPI_INT, iproc_prev, tag_prev_grid, MPI_COMM_WORLD, &requests[0]);
        MPI_Isend(grid + m + 2, m + 2,
                  MPI_INT, iproc_prev, tag_next_grid, MPI_COMM_WORLD, &requests[1]);
        MPI_Isend(grid + size_to_work * (m + 2), m + 2,
                  MPI_INT, iproc_next, tag_prev_grid, MPI_COMM_WORLD, &requests[2]);
        MPI_Irecv(grid + (size_to_work + 1) * (m + 2), m + 2,
                  MPI_INT, iproc_next, tag_next_grid, MPI_COMM_WORLD, &requests[3]);
        MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);

        debug_grid(m, grid, size_to_work, iproc, it);

        life_update(m, size_to_work, grid, iproc);
        debug_update(it, iproc);
    }
    debug_end(iproc);
    return it;
}

/******************************************************************************/

void life_update(int m, int n, int grid[], int iproc)

/******************************************************************************/
/*
  Purpose:

    LIFE_UPDATE updates a Life grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, int M, N, the number of rows and columns
    of interior grid cells.

    Input/output, int GRID[(1+M+1)*(1+N+1)], the data.
*/
{
    int i;
    int j;
    guardia_iproc(iproc, "Initializing s");
    int *s;
    s = malloc((m + 2) * (n + 2) * sizeof(int));
    guardia_iproc(iproc, "Setting to 0");
    memset(s, 0, (m + 2) * (n + 2) * sizeof(int));

    guardia_iproc(iproc, "starting update");

    //#pragma omp parallel for shared(s) private(i, j) num_threads(cores) schedule(guided)
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
//#pragma omp parallel for shared(grid) private(i, j) num_threads(cores) schedule(guided)
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
    free(s);
}

void compute_life_or_dead(int m, int n, const int *grid, int i, int j, int *s) {
    int i_prev, i_next, j_prev, j_next;
    i_prev = (1 < i) ? i - 1 : m;
    i_next = (i < m) ? i + 1 : 1;
    j_prev = j - 1;
    j_next = j + 1;
    s[i - 1 + (j - 1) * m] =
            grid[i_prev + (j_prev) * (m + 2)] + grid[i_prev + j * (m + 2)] + grid[i_prev + (j_next) * (m + 2)] +
            grid[i + (j_prev) * (m + 2)] + grid[i + (j_next) * (m + 2)] + grid[i_next + (j_prev) * (m + 2)] +
            grid[i_next + j * (m + 2)] + grid[i_next + (j_next) * (m + 2)];
}
#endif //LIFEGAME_MPI_UPDATE_FILE_H
