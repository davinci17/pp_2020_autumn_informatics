// Copyright 2020 Nadir Mohammed
#include<iostream>
#include<mpi.h>
#include<vector>
#include"../../../modules/task_3/nadir_m_dijkstra_algorithm/dijkstra_algorithm.h"


int minDistance(int dist[], bool sptSet[]) {
    int min = INT_MAX, min_index;
    for (int v = 0; v < r_col; v++) {
        if (sptSet[v] == false && dist[v] <= min)
            min = dist[v], min_index = v;
    }
    return min_index;
}

void dijkstra(int graph[r_col][r_col], int src, int* dist) {
    dist[r_col];
    bool sptSet[r_col];
    for (int i = 0; i < r_col; i++) { 
        dist[i] = INT_MAX, sptSet[i] = false; }
    dist[src] = 0;
    for (int count = 0; count < r_col - 1; count++) {
        int u = minDistance(dist, sptSet);
        sptSet[u] = true;

        for (int v = 0; v < r_col; v++) {
            if (!sptSet[v] && graph[u][v] && dist[u] != INT_MAX && dist[u] + graph[u][v] < dist[v]) {
                dist[v] = dist[u] + graph[u][v];
            }
        }
    }
    std::cout << "Vertex \t\t Distance from Source\n" << std::endl;
    for (int i = 0; i < r_col; i++) {
        std::cout << i << "\t\t" << dist[i] << std::endl;
    }
    std::cout << "********************************" << std::endl;
    delete[] dist;
}

void getParallelDijkstras(int graph[r_col * r_col], int src, int* dist) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    class temp_operat {
    public: 
        int attr_val, pos;
    };
    dist[r_col];
    bool sptSet[r_col];
    for (int i = 0; i < r_col; i++) { 
        dist[i] = INT_MAX, sptSet[i] = false; 
    }
    dist[src] = 0;
    int min = INT_MAX;
    int* sub_graph = new int[r_col * (r_col / size)];

    temp_operat mini_part;
    temp_operat  gen_part;
    int rest;
    if (rank == 0) {
        rest = 0;
    } else {
        rest = r_col % size + (r_col / size) * rank;
    }
    MPI_Scatter(&graph[(r_col % size) * r_col], (r_col / size) * r_col,
        MPI_INT, &sub_graph[0], (r_col / size) * r_col, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        for (int i = 0; i < r_col * (r_col / size); i++) { sub_graph[i] = graph[i]; }
    }
    int num_col = r_col / size;

    for (int i = 0; i < r_col - 1; i++) {
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
            if (sub_graph[gen_part.pos + r_col * k] != 0 &&
                dist[gen_part.pos] + sub_graph[gen_part.pos + r_col * k] < dist[k + rest]) {
                dist[k + rest] = dist[gen_part.pos] + sub_graph[gen_part.pos + r_col * k];
            }
        }
        int sizeResult = r_col % size + r_col / size + (r_col / size);
        if (rank == 0) {
            int* result = new int[r_col / size];

            for (int i = 1; i < size; i++) {
                MPI_Recv(&result[0], r_col / size,
                    MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                for (int i = r_col % size + r_col / size; i < sizeResult; i++) {
                    dist[i] = result[i];
                }
            }
        } else {
            for (int i = 1; i < size; i++) {
                if (rank == i) {
                    MPI_Send(&dist[r_col % size + rank * (r_col / size)], r_col / size,
                        MPI_INT, 0, 0, MPI_COMM_WORLD);
                }
            }
        }
        MPI_Bcast(&dist[0], sizeResult, MPI_INT, 0, MPI_COMM_WORLD);
    }
    for (int i = 0; i < r_col; i++) {
        std::cout << i << "\t\t" << dist[i] << std::endl;
    }
    delete[] sub_graph;
}
