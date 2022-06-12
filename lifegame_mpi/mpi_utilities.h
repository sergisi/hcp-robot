//
// Created by sergi on 28/05/22.
//

#ifndef LIFEGAME_MPI_MPI_UTILITIES_H
#define LIFEGAME_MPI_MPI_UTILITIES_H

#include "update_file.h"
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
#include "mpi_utilities.h"

void recieve_ended_grids(int m, int n, int8_t **grid, int nproc);

int get_start_n(int n, int nproc, int iproc);

int mpi_recv_int(int *it_max, int tag);

int get_size_to_work(int n, int nproc, int iproc);

void recieve_ended_grids(int m, int n, int8_t **grid, int nproc) {
    MPI_Request requests[nproc - 1];
    for (int i = 1; i < nproc; i++) {
        int startn = get_start_n(n, nproc, i);
        int size_to_work = get_size_to_work(n, nproc, i);
        char msg[200] = {0};
        snprintf(msg, 199, "[0] -> [%d]: Receiving from %d with a total of %d elements", i, startn, size_to_work);
        guardian(msg);
        MPI_Irecv(*grid + (startn + 1) * (m + 2), size_to_work * (m + 2), MPI_INT8_T, i, 0, MPI_COMM_WORLD, &requests[i - 1]);
    }
    MPI_Waitall(nproc - 1, requests, MPI_STATUSES_IGNORE);
}

int get_size_to_work(int n, int nproc, int iproc) {
    int size_to_work, module, division;
    module = n % nproc;
    division = n / nproc;
    size_to_work = division + (iproc < module ? 1 : 0);
    return size_to_work;
}

int get_guarded_int(int *it_max, int iproc, const char *var_name, int tag) {
    (*it_max) = mpi_recv_int(it_max, tag);
    if (debug) {
        char msg[200];
        sprintf(msg, "[%d]: : Receiving %s = %d", iproc, var_name, (*it_max));
        guardian(msg);
    }
    return (*it_max);
}

int send_grid(int m, int n, const int8_t *grid, int nproc) {
    int startn, size_to_work = 0;
    MPI_Request requests[nproc - 1];
    for (int i = 1; i < nproc; i++) {
        size_to_work = get_size_to_work(n, nproc, i);
        startn = get_start_n(n, nproc, i);
        char msg[200];
        sprintf(msg, "[0]: Sending grid of: %i x %i", size_to_work, m);
        guardian(msg);
        MPI_Isend(grid + (startn + 1) * (m + 2), size_to_work * (m + 2), MPI_INT8_T, i, tag_initial_grid, MPI_COMM_WORLD, &requests[i - 1]);
    }
    guardian("[0]: Sended all grids");
    MPI_Waitall(nproc - 1, requests, MPI_STATUSES_IGNORE);
    guardian("[0]: All responses gathered");
    return size_to_work;
}

int get_start_n(int n, int nproc, int iproc) {
    int startn;
    int module, division;
    module = n % nproc;
    division = n / nproc;
    startn = division * iproc + (iproc < module ? iproc : module);
    return startn;
}

int mpi_recv_int(int *it_max, int tag) {
    MPI_Recv(it_max, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    return (*it_max);
}

int mpi_send_int(int *it_max, int nproc, int tag) {
    MPI_Request requests[nproc - 1];
    for (int i = 0; i < nproc - 1; i++) {
        MPI_Isend(it_max, 1, MPI_INT, i + 1, tag, MPI_COMM_WORLD, &requests[i]);
    }
    MPI_Waitall(nproc - 1, requests, MPI_STATUSES_IGNORE);
    return (*it_max);
}

#endif //LIFEGAME_MPI_MPI_UTILITIES_H
