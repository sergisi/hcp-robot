//
// Created by sergi on 28/05/22.
//

#ifndef LIFEGAME_MPI_READ_LIFE_H
#define LIFEGAME_MPI_READ_LIFE_H

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
#include "read_life.h"
#include "mpi_utilities.h"

void mpi_read_initial_file(int argc, char *argv[], char *filename, double prob, int *seed, int *it_max, int *m, int *n,
                           int8_t **grid, double *start_time, int *nproc, int *iproc, int *size_to_work, int *iproc_prev,
                           int *iproc_next);

void recieve_all_data(int nproc, int iproc, int *it_max, int *m, int *n, int8_t **grid);

void initialization(int argc, char *const *argv, char *filename, double prob, int *seed, int nproc, int *it_max, int *m,
                    int *n, int8_t **grid, double *start_time);

int8_t *life_init(char *filename, double prob, int m, int n, int *seed);

void life_read(char *input_filename, int m, int n, int8_t grid[]);

void mpi_read_initial_file(int argc, char *argv[], char *filename, double prob, int *seed, int *it_max, int *m, int *n,
                           int8_t **grid, double *start_time, int *nproc, int *iproc, int *size_to_work, int *iproc_prev,
                           int *iproc_next) {
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, iproc);
    MPI_Comm_size(MPI_COMM_WORLD, nproc);
    if ((*iproc) == 0) {
        initialization(argc, argv, filename, prob, seed, (*nproc), it_max, m, n, grid, start_time);
    } else {
        recieve_all_data((*nproc), (*iproc), it_max, m, n, grid);
    }

    (*iproc_prev) = (*iproc) == 0 ? (*nproc) - 1 : (*iproc) - 1;
    (*iproc_next) = (*iproc) == (*nproc) - 1 ? 0 : (*iproc) + 1;

    debug_context((*it_max), (*m), (*n), (*iproc), (*iproc_prev), (*iproc_next));
    (*size_to_work) = get_size_to_work((*n), (*nproc), (*iproc));
}

void
recieve_all_data(int nproc, int iproc, int *it_max, int *m, int *n, int8_t **grid) {// Data to be received: it_max, n, m,
    (*it_max) = get_guarded_int(it_max, iproc, "it_max", tag_it_max);
    (*n) = get_guarded_int(n, iproc, "n", tag_n);
    (*m) = get_guarded_int(m, iproc, "m", tag_m);
    cores = get_guarded_int(&cores, iproc, "cores", tag_cores);
    int size_to_work = get_size_to_work((*n), nproc, iproc);
    // receive grid
    guardian("[1]: Allocating grid");
    (*grid) = (int8_t *) malloc((size_to_work + 2) * (*m + 2) * sizeof(int8_t));
    guardian("[1]: Allocated grid");
    char msg[200];
    sprintf(msg, "[1]: Receiving grid of: %d x %d", size_to_work, (*m));
    guardian(msg);
    MPI_Recv((*grid) + (*m + 2), size_to_work * (*m + 2), MPI_INT8_T, 0, tag_initial_grid, MPI_COMM_WORLD,
             MPI_STATUSES_IGNORE);
    guardian("[1]: received grid");
    debug_grid(*m, *grid, size_to_work, iproc, 0);
}

void initialization(int argc, char *const *argv, char *filename, double prob, int *seed, int nproc, int *it_max, int *m,
                    int *n, int8_t **grid, double *start_time) {
    timestamp();
    printf("\n");
    printf("LIFE GAME SERIAL\n");
    printf("  C version\n");
    printf("  Carry out a few steps of John Conway's\n");
    printf("  Game of Life.\n");
    printf("  Parameters: lifegame [input-file] [x] [y] [iters] .\n");
    printf("\n");

    if (argc > 5)
        (*it_max) = atoi(argv[5]);
    if (argc > 4)
        (*n) = atoi(argv[4]);
    if (argc > 3)
        (*m) = atoi(argv[3]);
    if (argc > 2) {
        filename = argv[2];
    } else {
        filename = NULL;
    }
    if (argc > 1)
        cores = atoi(argv[1]);
    else {
        printf("Not specified number of cores");
        exit(-1);
    }

    (*grid) = life_init(filename, prob, (*m), (*n), seed);
    debug_grid(*m, *grid, *n, 0, 0);
    (*start_time) = MPI_Wtime();
    guardian("[0]: Sending it_max");
    (*it_max) = mpi_send_int(it_max, nproc, tag_it_max);
    guardian("[0]: Sending n");
    (*n) = mpi_send_int(n, nproc, tag_n);
    guardian("[0]: Sending m");
    (*m) = mpi_send_int(m, nproc, tag_m);
    guardian("[0]: Sending grid");
    cores = mpi_send_int(&cores, nproc, tag_cores);
    send_grid((*m), (*n), (*grid), nproc);
    guardian("[0]: End initial");
}

/******************************************************************************/

int8_t *life_init(char *filename, double prob, int m, int n, int *seed)

/******************************************************************************/
/*
  Purpose:

    LIFE_INIT initializes the life grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    08 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, double PROB, the probability that a grid cell
    should be alive.

    Input, int M, N, the number of rows and columns
    of interior grid cells.

    Input/output, int *SEED, a seed for the random
    number generator.

    Output, int LIFE_INIT[(1+M+1)*(1+N+1)], the initial grid.
*/
{
    int8_t *grid;
    int i;
    int j;
    double r;


    grid = (int8_t *) malloc((m + 2) * (n + 2) * sizeof(int8_t ));
    if (grid == NULL)
        perror("Error malloc grid:");

    if (filename != NULL) {
        /* Read input file */
        printf("Reading Input filename %s\n", filename);
        life_read(filename, m, n, grid);
    } else {

#pragma omp parallel for shared(grid) private(i, j) num_threads(cores) schedule(guided)
        for (j = 0; j <= n + 1; j++) {
            for (i = 0; i <= m + 1; i++) {
                grid[i + j * (m + 2)] = 0;
            }
        }

#pragma omp parallel for shared(grid) private(i, j) num_threads(cores) schedule(guided)
        for (j = 1; j <= n; j++) {
            for (i = 1; i <= m; i++) {
                r = r8_uniform_01(seed);
                if (r <= prob) {
                    grid[i + j * (m + 2)] = 1;
                }
            }
        }
    }

    return grid;
}

/******************************************************************************/

void life_read(char *filename, int m, int n, int8_t grid[])

/******************************************************************************/
/*
  Purpose:

    LIFE_READ reads a file to a grid.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Parameters:

    Input, char *OUTPUT_FILENAME, the output file name.

    Input, int M, N, the number of rows and columns
    of interior grid cells.

    Input, int GRID[(1+M+1)*(1+N+1)], the data.
*/
{
    int i;
    int j;
    FILE *input_unit;
/*
  input the file.
*/
    input_unit = fopen(filename, "rt");
    if (input_unit == NULL)
        perror("Reading input file:");
/*
  Read the data.
*/
    for (j = 1; j <= n; j++) {
        for (i = 1; i <= m; i++) {
            fscanf(input_unit, "%hhd", &(grid[i + j * (m + 2)]));
        }
    }
    /* Set the grid borderline to 0's */
#pragma omp parallel for num_threads(cores) schedule(guided)
    for (j = 0; j <= n + 1; j++) {
        grid[0 + j * (m + 2)] = 0;
        grid[(m + 1) + j * (m + 2)] = 0;

    }
#pragma omp parallel for num_threads(cores) schedule(guided)
    for (i = 1; i <= m; i++) {
        grid[i + 0 * (m + 2)] = 0;
        grid[i + (n + 1) * (m + 2)] = 0;
    }
/*
  Close the file.
*/
    fclose(input_unit);
}

#endif //LIFEGAME_MPI_READ_LIFE_H
