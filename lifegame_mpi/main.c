#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <omp.h>
#include <mpi.h>
#include "debug_tools.h"
#include "utilities.h"
#include "constants.h"
#include "read_life.h"
#include "update_file.h"
#include "mpi_utilities.h"


int main(int argc, char *argv[]);

void life_write(char *output_filename, int m, int n, int8_t grid[]);

void
print_result(char *output_filename, int it, int it_max, int m, int n, int8_t *grid, double start_time, double run_time,
             int nproc);

int8_t *save_file(char *output_filename, int it, int it_max, int m, int n, int8_t **grid, double start_time, double run_time,
               int nproc, int iproc, int size_to_work);

int main(int argc, char *argv[]) {
    char *filename = NULL, output_filename[100];
    int it = 0;
    int it_max = 10;
    int m = 10;
    int n = 10;
    int8_t *grid;
    double prob = 0.20;
    int seed = 123456789;
    double start_time = 0, run_time = 0;
    // variables needed for MPI
    int nproc, iproc, size_to_work = 0, iproc_prev, iproc_next;

    MPI_Request requests[4];
    mpi_read_initial_file(argc, argv, filename, prob, &seed, &it_max, &m, &n, &grid, &start_time, &nproc, &iproc,
                          &size_to_work, &iproc_prev, &iproc_next);

    printf("[%d]: %d cores\n", iproc, cores);
    it = mpi_life_update(it, it_max, m, grid, iproc, size_to_work, iproc_prev, iproc_next, requests);

    grid = save_file(output_filename, it, it_max, m, n, &grid, start_time, run_time, nproc, iproc, size_to_work);

    free(grid);
    MPI_Finalize();
    return 0;

}

int8_t *save_file(char *output_filename, int it, int it_max, int m, int n, int8_t **grid, double start_time, double run_time,
               int nproc, int iproc, int size_to_work) {
    if (iproc == 0) {
        recieve_ended_grids(m, n, grid, nproc);
        guardian("[0]: Recieved all grids");
        print_result(output_filename, it, it_max, m, n, (*grid), start_time, run_time, nproc);
    } else {
        debug_grid(m, *grid, size_to_work, iproc, it_max * 2);
        MPI_Send((*grid) + m + 2, size_to_work * (m + 2), MPI_INT8_T, 0, 0, MPI_COMM_WORLD);
    }
    return (*grid);
}

void
print_result(char *output_filename, int it, int it_max, int m, int n, int8_t *grid, double start_time, double run_time,
             int nproc) {
    run_time = MPI_Wtime() - start_time;
    printf("GoL MxN (%d x %d), time=%f iters=%d openmp-cores=%d mpi-cores=%d\n", m, n, run_time, it_max, cores, nproc);
    sprintf(output_filename, DDefaultOutputFilename, it - 1);
    life_write(output_filename, m, n, grid);

/*
  Terminate.
*/
    printf("\n");
    printf("LIFE GAME SERIAL\n");
    printf("  Normal end of execution.\n");
    printf("\n");
    timestamp();
}


/******************************************************************************/

void life_write(char *output_filename, int m, int n, int8_t grid[])

/******************************************************************************/
/*
  Purpose:

    LIFE_WRITE writes a grid to a file.

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    10 September 2013

  Author:

    John Burkardt

  Parameters:

    Input, char *OUTPUT_FILENAME, the output file name.

    Input, int M, N, the number of rows and columns
    of interior grid cells.

    Input, int8_t GRID[(1+M+1)*(1+N+1)], the data.
*/
{
    int i;
    int j;
    FILE *output_unit;
/*
  Open the file.
*/
    output_unit = fopen(output_filename, "wt");
/*
  Write the data.
*/
    for (j = 1; j <= n; j++) {
        for (i = 1; i <= m; i++) {
            fprintf(output_unit, " %d", grid[i + j * (m + 2)]);
        }
        fprintf(output_unit, "\n");
    }
/*
  Close the file.
*/
    fclose(output_unit);
}

