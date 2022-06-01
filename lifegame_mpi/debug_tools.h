//
// Created by sergi on 28/05/22.
//

#ifndef LIFEGAME_MPI_DEBUG_TOOLS_H
#define LIFEGAME_MPI_DEBUG_TOOLS_H

#include <stdio.h>
#include "constants.h"

void debug_grid(int m, const int *grid, int size_to_work, int iproc, int it);

void debug_end(int iproc);

void debug_context(int it_max, int m, int n, int iproc, int iproc_prev, int iproc_next);

void debug_update(int it, int iproc);

void guardian(const char *msg);

void guardia_iproc(int iproc, char mmsg[]);

void debug_grid(int m, const int *grid, int size_to_work, int iproc, int it) {
    if (debug && m < 50) {
        for (int i = 0; i < size_to_work + 2; i++) {
            printf("[%d, %d, %d]\t", iproc, it, i);
            for (int j = 1; j <= m; j++) {
                printf("%d ", grid[j + i * (m + 2)]);
            }
            printf("\n");
        }
    }
}

void debug_end(int iproc) {
    char msg[200] = {0};
    sprintf(msg, "[%d]: ended update", iproc);
    guardian(msg);
}

void debug_context(int it_max, int m, int n, int iproc, int iproc_prev, int iproc_next) {
    char msg[200] = {0};
    sprintf(msg,
            "[%d]: starting Context(n=%d, m=%d, it_max=%d, prev=%d, next=%d)",
            iproc, n, m, it_max, iproc_prev, iproc_next);
    guardian(msg);
}

void debug_update(int it, int iproc) {
    char msg[200] = {0};
    snprintf(msg, 199, "[%d]: update %d", iproc, it);
    guardian(msg);
}

void guardian(const char *msg) {
    if (debug) {
        printf("%s\n", msg);
    }
}

void guardia_iproc(int iproc, char mmsg[]) {
    char msg[200] = {0};
    snprintf(msg, 199, "[%d]: %s", iproc, mmsg);
    guardian(msg);
}
#endif //LIFEGAME_MPI_DEBUG_TOOLS_H
