// Copyright 2020 Nadir Mohammed
#ifndef MODULES_TASK_3_NADIR_M_DIJKSTRA_ALGORITHM_DIJKSTRA_ALGORITHM_H_
#define MODULES_TASK_3_NADIR_M_DIJKSTRA_ALGORITHM_DIJKSTRA_ALGORITHM_H_

const int kCol = 9;
int minDistance(int dist[], bool sptSet[]);
void dijkstra(int graph[kCol][kCol], int src, int* dist);
void getParallelDijkstras(int graph[kCol * kCol], int src, int* dist);
class temp_operat {
    public:
        int attr_val, pos;
};
#endif  // MODULES_TASK_3_NADIR_M_DIJKSTRA_ALGORITHM_DIJKSTRA_ALGORITHM_H_
