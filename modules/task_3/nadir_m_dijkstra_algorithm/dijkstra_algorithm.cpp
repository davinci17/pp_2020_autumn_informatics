// Copyright 2020 Nadir Mohammed
#include"../../../modules/task_3/nadir_m_dijkstra_algorithm/dijkstra_algorithm.h"
#include<mpi.h>
#include<iostream>
#include<vector>

int maxPoint = 50;
int minDistance(int dBettPoint[], bool sptSet[]) {
    int min = maxPoint, min_index;
    for (int v = 0; v < kCol; v++) {
        if (sptSet[v] == false && dBettPoint[v] <= min)
            min = dBettPoint[v], min_index = v;
    }
    return min_index;
}

void dijkstra(int graph[kCol][kCol], int src, int* dBettPoint) {
    const int kCol = 9;
    dBettPoint = new int[kCol];
    bool sptSet[kCol];
    for (int i = 0; i < kCol; i++) {
        dBettPoint[i] = maxPoint, sptSet[i] = false; }
    dBettPoint[src] = 0;
    for (int count = 0; count < kCol - 1; count++) {
        int u = minDistance(dBettPoint, sptSet);
        sptSet[u] = true;

        for (int v = 0; v < kCol; v++) {
            if (!sptSet[v] && graph[u][v] && dBettPoint[u] != maxPoint && dBettPoint[u] + graph[u][v] < dBettPoint[v]) {
                dBettPoint[v] = dBettPoint[u] + graph[u][v];
            }
        }
    }
    delete[] dBettPoint;
}

void getParallelDijkstras(int graph[kCol * kCol], int src, int* dBettPoint) {
    const int kCol = 9;
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    temp_operat miniPart, genPart;
    dBettPoint = new int[kCol];
    bool sptSet[kCol];
    int* result = new int[kCol / size];
    for (int i = 0; i < kCol; i++) {
        dBettPoint[i] = maxPoint, sptSet[i] = false;
    }
    dBettPoint[src] = 0;
    int min = maxPoint, rest;
    int* secondArr = new int[kCol * (kCol / size)];

    if (rank == 0) {
        rest = 0;
    } else {
        rest = kCol % size + (kCol / size) * rank;
    }
    MPI_Scatter(&graph[(kCol % size) * kCol], (kCol / size) * kCol,
        MPI_INT, &secondArr[0], (kCol / size) * kCol, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        for (int i = 0; i < kCol * (kCol / size); i++) { secondArr[i] = graph[i]; }
    }
    int num_col = kCol / size;

    for (int i = 0; i < kCol - 1; i++) {
        miniPart.attr_val = -1;
        miniPart.pos = -1;
        for (int j = rest; j < num_col + rest; j++) {
            if (!sptSet[j] && (miniPart.pos == -1 || dBettPoint[j] < dBettPoint[miniPart.pos])) {
                miniPart.pos = j;
                miniPart.attr_val = dBettPoint[miniPart.pos];
            }
        }
        if (miniPart.pos == -1) {
            miniPart.attr_val = min;
        }
        MPI_Allreduce(&miniPart, &genPart, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
        if (genPart.pos == -1 || dBettPoint[genPart.pos] == min) {
            break;
        }
        sptSet[genPart.pos] = true;

        for (int cou = 0; cou < num_col; cou++) {
            if (secondArr[genPart.pos + kCol * cou] != 0 &&
                dBettPoint[genPart.pos] + secondArr[genPart.pos + kCol * cou] < dBettPoint[cou + rest]) {
                dBettPoint[cou + rest] = dBettPoint[genPart.pos] + secondArr[genPart.pos + kCol * cou];
            }
        }
        int sizeResult = kCol % size + kCol / size + (kCol / size);
        if (rank == 0) {
            for (int i = 1; i < size; i++) {
                MPI_Recv(&result[0], kCol / size,
                    MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for (int i = kCol % size + kCol / size; i < sizeResult; i++) {
                    dBettPoint[i] = result[i];
                }
            }
        } else {
            for (int i = 1; i < size; i++) {
                if (rank == i) {
                    MPI_Send(&dBettPoint[kCol % size + rank * (kCol / size)], kCol / size,
                        MPI_INT, 0, 0, MPI_COMM_WORLD);
                }
            }
        }
        MPI_Bcast(&dBettPoint[0], sizeResult, MPI_INT, 0, MPI_COMM_WORLD);
    }
    for (int i = 0; i < kCol; i++) {
        std::cout << i << "\t\t" << dBettPoint[i] << std::endl;
    }
    delete[] secondArr;
    delete[] result;
    delete[] dBettPoint;
}
