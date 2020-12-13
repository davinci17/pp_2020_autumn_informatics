// Copyright 2020 Nadir Mohammed
#ifndef MODULES_TASK_3_NADIR_M_DIJKSTRA_ALGORITHM_DIJKSTRA_ALGORITHM_H_
#define MODULES_TASK_3_NADIR_M_DIJKSTRA_ALGORITHM_DIJKSTRA_ALGORITHM_H_

const int r_col = 9;
int minDistance(int dist[], bool sptSet[]);
void dijkstra(int graph[r_col][r_col], int src, int* dist);
void getParallelDijkstras(int graph[r_col * r_col], int src, int* dist);

#endif // MODULES_TASK_3_NADIR_M_DIJKSTRA_ALGORITHM_DIJKSTRA_ALGORITHM_H_
