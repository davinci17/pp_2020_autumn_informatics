// Copyright 2020 Nadir Mohammed
#include"../../../modules/task_3/nadir_m_dijkstra_algorithm/dijkstra_algorithm.h"
#include<mpi.h>
#include<iostream>
#include<vector>

int maxPoint = INT_MAX;
int minDistance(int dist[], bool sptSet[]) {
    int min = maxPoint, min_index;
    for (int v = 0; v < kCol; v++) {
        if (sptSet[v] == false && dist[v] <= min)
            min = dist[v], min_index = v;
    }
    return min_index;
}

void dijkstra(int graph[kCol][kCol], int src, int* dist) {
    const int kCol = 9;
    dist = new int[kCol];
    bool sptSet[kCol];
    for (int i = 0; i < kCol; i++) {
        dist[i] = maxPoint, sptSet[i] = false; }
    dist[src] = 0;
    for (int count = 0; count < kCol - 1; count++) {
        int u = minDistance(dist, sptSet);
        sptSet[u] = true;

        for (int v = 0; v < kCol; v++) {
            if (!sptSet[v] && graph[u][v] && dist[u] != maxPoint && dist[u] + graph[u][v] < dist[v]) {
                dist[v] = dist[u] + graph[u][v];
            }
        }
    }
}

void getParallelDijkstras(int graph[kCol * kCol], int src, int* dist) {
    const int kCol = 9;
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    temp_operat mini_part, gen_part;
    dist = new int[kCol];
    bool sptSet[kCol];
    int* result = new int[kCol / size];
    for (int i = 0; i < kCol; i++) {
        dist[i] = maxPoint, sptSet[i] = false;
    }
    dist[src] = 0;
    int min = maxPoint, rest;
    int* sub_graph = new int[kCol * (kCol / size)];

    if (rank == 0) {
        rest = 0;
    } else {
        rest = kCol % size + (kCol / size) * rank;
    }
    MPI_Scatter(&graph[(kCol % size) * kCol], (kCol / size) * kCol,
        MPI_INT, &sub_graph[0], (kCol / size) * kCol, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        for (int i = 0; i < kCol * (kCol / size); i++) { sub_graph[i] = graph[i]; }
    }
    int num_col = kCol / size;

    for (int i = 0; i < kCol - 1; i++) {
        mini_part.attr_val = -1;
        mini_part.pos = -1;
        for (int j = rest; j < num_col + rest; j++) {
            if (!sptSet[j] && (mini_part.pos == -1 || dist[j] < dist[mini_part.pos])) {
                mini_part.pos = j;
                mini_part.attr_val = dist[mini_part.pos];
            }
        }
        if (mini_part.pos == -1) {
            mini_part.attr_val = min;
        }
        MPI_Allreduce(&mini_part, &gen_part, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
        if (gen_part.pos == -1 || dist[gen_part.pos] == min) {
            break;
        }
        sptSet[gen_part.pos] = true;

        for (int k = 0; k < num_col; k++) {
            if (sub_graph[gen_part.pos + kCol * k] != 0 &&
                dist[gen_part.pos] + sub_graph[gen_part.pos + kCol * k] < dist[k + rest]) {
                dist[k + rest] = dist[gen_part.pos] + sub_graph[gen_part.pos + kCol * k];
            }
        }
        int sizeResult = kCol % size + kCol / size + (kCol / size);
        if (rank == 0) {
            for (int i = 1; i < size; i++) {
                MPI_Recv(&result[0], kCol / size,
                    MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for (int i = kCol % size + kCol / size; i < sizeResult; i++) {
                    dist[i] = result[i];
                }
            }
        } else {
            for (int i = 1; i < size; i++) {
                if (rank == i) {
                    MPI_Send(&dist[kCol % size + rank * (kCol / size)], kCol / size,
                        MPI_INT, 0, 0, MPI_COMM_WORLD);
                }
            }
        }
        MPI_Bcast(&dist[0], sizeResult, MPI_INT, 0, MPI_COMM_WORLD);
    }
    for (int i = 0; i < kCol; i++) {
        std::cout << i << "\t\t" << dist[i] << std::endl;
    }
    delete[] sub_graph;
    delete[] result;
    delete[] dist;
}
